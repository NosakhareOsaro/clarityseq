"""
modules.ensemble.ensemble
==========================
Ensemble variant caller merging GATK4 HaplotypeCaller and DeepVariant v1.8.0.

Error profile complementarity
------------------------------
GATK4 uses graph assembly (local de-novo Kmer graph, ~150 bp haplotype windows):
    - Strengths : complex indels, short-tandem-repeat (STR) alleles, joint
                  genotyping across cohorts, well-calibrated VQSR scores.
    - Weaknesses: systematically mis-calls some read-end artefacts; elevated
                  false-positive rate near low-complexity / segmental-dup regions.

DeepVariant v1.8.0 uses 6-channel CNN pileup images (InceptionV3):
    - Strengths : lower SNV false-positive rate, robust in repetitive regions,
                  SPRQ support (v1.8.0) improves tandem-repeat loci.
    - Weaknesses: single-sample only, less accurate for complex indels >50 bp,
                  cannot directly exploit cohort genotype priors.

Because the two callers have **partially non-overlapping error profiles**,
combining them can improve overall precision (INTERSECTION) or sensitivity
(UNION) relative to either caller alone.

Modes
-----
INTERSECTION (default)
    A variant is retained iff it is called PASS by BOTH callers.
    - Higher precision; fewer false positives.
    - Recommended for clinical reporting (reduces false positive P/LP calls).
    - False negatives relative to single-caller output are partially mitigated
      because individual GATK4 and DeepVariant VCFs are also published to the
      pipeline output directory for downstream review.

UNION
    A variant is retained iff it is called PASS by EITHER caller.
    - Higher sensitivity; more false positives.
    - Use only for research or when maximising recall (e.g. research cohorts,
      carrier frequency studies).
    - NOT recommended for primary clinical reporting.

INFO fields added by this module
---------------------------------
    ENSEMBLE_CALLER  : str   "GATK4_ONLY" | "DV_ONLY" | "BOTH"
    ENSEMBLE_MODE    : str   "INTERSECTION" | "UNION"

Clinical rationale for INTERSECTION as default
----------------------------------------------
In a clinical diagnostic context, a false-positive pathogenic/likely-pathogenic
(P/LP) report causes:
    1. Unnecessary patient distress and potential invasive follow-up procedures.
    2. Preventable healthcare resource expenditure (confirmatory Sanger, etc.).
    3. Potential mismanagement (e.g. prophylactic surgery on false cancer risk).

False negatives (missed pathogenic variants) in the ensemble INTERSECTION call
set are partially mitigated by:
    a) The individual GATK4 single-caller VCF also being published.
    b) The individual DeepVariant single-caller VCF also being published.
    c) Laboratory scientists reviewing single-caller outputs for phenotype-
       prioritised genes when ensemble VCF lacks a candidate.

This module implements the ensemble strategy recommended in:
    - ACGS Best Practice Guidelines v1.2 2024 §4.4 (ensemble calling)
    - Poplin et al. 2018 Nat Biotechnol PMID:30247488 (DeepVariant)
    - Zook et al. 2019 Nat Biotechnol PMID:30936564 (GIAB benchmarking)

Usage (CLI)
-----------
    python ensemble.py \\
        --gatk-vcf sample.gatk4.vcf.gz \\
        --dv-vcf   sample.deepvariant.vcf.gz \\
        --output   sample.ensemble.vcf \\
        --mode     INTERSECTION

    python ensemble.py --help
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Iterator

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VALID_MODES = ("INTERSECTION", "UNION")

# VCF field indices (0-based within a split tab-delimited line)
_CHROM = 0
_POS   = 1
_ID    = 2
_REF   = 3
_ALT   = 4
_QUAL  = 5
_FILT  = 6
_INFO  = 7


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _open(path: Path):
    """Open plain or gzip-compressed VCF transparently."""
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _is_pass(filter_field: str) -> bool:
    """Return True if the FILTER column indicates a PASS variant.

    Both literal "PASS" and "." (missing / not applied) are accepted.
    A trailing semi-colon list of non-PASS filters is treated as failing.
    """
    return filter_field in ("PASS", ".")


# ---------------------------------------------------------------------------
# Core parsing
# ---------------------------------------------------------------------------

def parse_vcf(path: Path) -> dict[str, dict]:
    """Parse a VCF file and return a dict keyed by canonical variant key.

    The canonical key is  ``CHROM:POS:REF:ALT``  (multi-allelic sites are
    split into one key per ALT allele so that per-allele status is tracked).

    Parameters
    ----------
    path:
        Path to a VCF file (plain text or ``.gz`` / ``.bgz`` compressed).

    Returns
    -------
    dict mapping variant_key -> record_dict with keys:
        ``chrom``   (str)   chromosome / contig name
        ``pos``     (str)   1-based position string (kept as str for VCF round-trip)
        ``id``      (str)   VCF ID field (dbSNP rsID etc.)
        ``ref``     (str)   reference allele
        ``alt``     (str)   single ALT allele (already split from multi-allelic)
        ``qual``    (str)   QUAL field
        ``filter``  (str)   FILTER field
        ``info``    (str)   INFO field
        ``rest``    (str)   FORMAT + all sample columns joined by tab
        ``pass``    (bool)  True iff FILTER is PASS or "."
        ``line``    (str)   the original (split) VCF line for that alt allele
    """
    variants: dict[str, dict] = {}

    with _open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue  # malformed line — skip silently

            chrom  = parts[_CHROM]
            pos    = parts[_POS]
            vid    = parts[_ID]
            ref    = parts[_REF]
            alts   = parts[_ALT].split(",")  # handle multi-allelic
            qual   = parts[_QUAL]
            filt   = parts[_FILT]
            info   = parts[_INFO]
            rest   = "\t".join(parts[8:]) if len(parts) > 8 else ""

            is_p = _is_pass(filt)

            for alt in alts:
                key = f"{chrom}:{pos}:{ref}:{alt}"
                variants[key] = {
                    "chrom":  chrom,
                    "pos":    pos,
                    "id":     vid,
                    "ref":    ref,
                    "alt":    alt,
                    "qual":   qual,
                    "filter": filt,
                    "info":   info,
                    "rest":   rest,
                    "pass":   is_p,
                    "line":   line,
                }

    return variants


# ---------------------------------------------------------------------------
# VCF header helpers
# ---------------------------------------------------------------------------

def _build_ensemble_header_lines() -> list[str]:
    """Return the INFO meta-lines introduced by this module."""
    return [
        '##INFO=<ID=ENSEMBLE_CALLER,Number=1,Type=String,'
        'Description="Which caller(s) produced this PASS variant: '
        'GATK4_ONLY, DV_ONLY, or BOTH">',
        '##INFO=<ID=ENSEMBLE_MODE,Number=1,Type=String,'
        'Description="Ensemble strategy applied: INTERSECTION or UNION">',
        '##ensemble_merge_command=' + " ".join(sys.argv),
    ]


def _inject_ensemble_info(info_field: str, caller: str, mode: str) -> str:
    """Append ENSEMBLE_CALLER and ENSEMBLE_MODE tags to an INFO string."""
    extra = f"ENSEMBLE_CALLER={caller};ENSEMBLE_MODE={mode}"
    if info_field in (".", ""):
        return extra
    return f"{info_field};{extra}"


# ---------------------------------------------------------------------------
# Merge logic
# ---------------------------------------------------------------------------

def merge_vcfs(
    gatk_vcf: Path,
    dv_vcf:   Path,
    mode:     str = "INTERSECTION",
) -> Iterator[str]:
    """Merge GATK4 and DeepVariant VCFs and yield output VCF lines.

    This generator yields:
        1. All meta-information lines from GATK4 VCF header (##).
        2. The ENSEMBLE INFO meta-lines.
        3. The #CHROM header line from GATK4 VCF.
        4. Data lines for merged variants (sorted by chrom/pos).

    Parameters
    ----------
    gatk_vcf:
        Path to GATK4 HaplotypeCaller VCF (post-VQSR).
    dv_vcf:
        Path to DeepVariant VCF.
    mode:
        ``"INTERSECTION"`` — keep variants PASS in both callers (default).
        ``"UNION"``        — keep variants PASS in either caller.

    Yields
    ------
    str
        One VCF line per yield (without trailing newline).

    Raises
    ------
    ValueError
        If ``mode`` is not one of the valid values.
    """
    if mode not in VALID_MODES:
        raise ValueError(
            f"Invalid ensemble mode {mode!r}. Must be one of {VALID_MODES}."
        )

    # ── 1. Collect GATK4 header lines ────────────────────────────────────────
    gatk_meta_lines: list[str] = []
    gatk_chrom_line: str | None = None

    with _open(gatk_vcf) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                gatk_meta_lines.append(line)
            elif line.startswith("#CHROM"):
                gatk_chrom_line = line
                break   # data lines follow — stop header scan

    # Emit header
    for hline in gatk_meta_lines:
        yield hline
    for eline in _build_ensemble_header_lines():
        yield eline
    if gatk_chrom_line:
        yield gatk_chrom_line

    # ── 2. Parse data records from both callers ───────────────────────────────
    gatk_variants = parse_vcf(gatk_vcf)
    dv_variants   = parse_vcf(dv_vcf)

    # All unique keys across both callers
    all_keys: set[str] = set(gatk_variants) | set(dv_variants)

    # Collect output records before sorting
    output_records: list[tuple[str, int, str]] = []  # (chrom, pos_int, line)

    for key in all_keys:
        in_gatk = key in gatk_variants and gatk_variants[key]["pass"]
        in_dv   = key in dv_variants   and dv_variants[key]["pass"]

        # Apply mode filter
        if mode == "INTERSECTION":
            if not (in_gatk and in_dv):
                continue   # skip variants not in both callers
        else:  # UNION
            if not (in_gatk or in_dv):
                continue   # neither caller produced a PASS call

        # Determine ENSEMBLE_CALLER tag
        if in_gatk and in_dv:
            caller_tag = "BOTH"
        elif in_gatk:
            caller_tag = "GATK4_ONLY"
        else:
            caller_tag = "DV_ONLY"

        # Use GATK4 record as canonical source when available; else DV
        source = gatk_variants[key] if in_gatk else dv_variants[key]

        new_info = _inject_ensemble_info(source["info"], caller_tag, mode)

        # Reconstruct the VCF line
        fields = [
            source["chrom"],
            source["pos"],
            source["id"],
            source["ref"],
            source["alt"],
            source["qual"],
            source["filter"],
            new_info,
        ]
        if source["rest"]:
            fields.append(source["rest"])

        data_line = "\t".join(fields)

        try:
            pos_int = int(source["pos"])
        except ValueError:
            pos_int = 0

        output_records.append((source["chrom"], pos_int, data_line))

    # ── 3. Sort by chrom (natural) then position ──────────────────────────────
    def _sort_key(record: tuple[str, int, str]) -> tuple[list[int | str], int]:
        chrom = record[0]
        pos   = record[1]
        # Natural sort for chromosomes: chr1 < chr2 < ... < chr10 < chrX < chrY
        parts: list[int | str] = []
        for token in chrom.replace("chr", "").split("_"):
            if token.isdigit():
                parts.append(int(token))
            else:
                parts.append(token)
        return (parts, pos)

    output_records.sort(key=_sort_key)

    for _, _, line in output_records:
        yield line


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Command-line interface for ensemble VCF merging.

    Examples
    --------
    # Clinical report (default INTERSECTION mode)
    python ensemble.py \\
        --gatk-vcf  sample.gatk4.vcf.gz \\
        --dv-vcf    sample.deepvariant.vcf.gz \\
        --output    sample.ensemble.vcf

    # Research (UNION — maximise sensitivity)
    python ensemble.py \\
        --gatk-vcf  sample.gatk4.vcf.gz \\
        --dv-vcf    sample.deepvariant.vcf.gz \\
        --output    sample.ensemble.union.vcf \\
        --mode      UNION
    """
    parser = argparse.ArgumentParser(
        prog="ensemble.py",
        description=(
            "Merge GATK4 HaplotypeCaller and DeepVariant VCFs into an "
            "ensemble call set (INTERSECTION or UNION mode)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--gatk-vcf",
        required=True,
        type=Path,
        metavar="GATK_VCF",
        help="Path to GATK4 HaplotypeCaller VCF (post-VQSR). Plain or .gz.",
    )
    parser.add_argument(
        "--dv-vcf",
        required=True,
        type=Path,
        metavar="DV_VCF",
        help="Path to DeepVariant v1.8.0 VCF. Plain or .gz.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="OUTPUT_VCF",
        help=(
            "Output ensemble VCF path. Use '-' to write to stdout. "
            "Plain text output only (pipe through bgzip if needed)."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=VALID_MODES,
        default="INTERSECTION",
        help=(
            "Ensemble strategy. "
            "INTERSECTION (default): keep variants PASS in BOTH callers — "
            "higher precision, recommended for clinical reporting. "
            "UNION: keep variants PASS in EITHER caller — higher sensitivity, "
            "use for research or maximising recall."
        ),
    )

    args = parser.parse_args()

    # Validate inputs exist
    if not args.gatk_vcf.exists():
        sys.exit(f"ERROR: GATK VCF not found: {args.gatk_vcf}")
    if not args.dv_vcf.exists():
        sys.exit(f"ERROR: DeepVariant VCF not found: {args.dv_vcf}")

    # Open output (stdout or file)
    if str(args.output) == "-":
        out_fh = sys.stdout
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        out_fh = open(args.output, "w", encoding="utf-8")

    try:
        n_written = 0
        for line in merge_vcfs(args.gatk_vcf, args.dv_vcf, args.mode):
            out_fh.write(line + "\n")
            if not line.startswith("#"):
                n_written += 1
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    print(
        f"ensemble.py: wrote {n_written} variant records "
        f"(mode={args.mode}) to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
