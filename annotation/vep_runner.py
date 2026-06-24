"""
VEP v111 runner with MANE Select transcript prioritisation.

Transcript Selection Strategy
------------------------------
Ensembl VEP v111 implements MANE Select-first transcript prioritisation
following the MANE project (Morales et al. 2022 PMID:35356062).

Pick order (``--pick_order`` flag):
  1. mane_select       – MANE Select canonical transcript
  2. mane_plus_clinical – MANE Plus Clinical (additional disease transcripts)
  3. canonical         – Ensembl canonical transcript fallback

This order guarantees that reported HGVSc/HGVSp nomenclature uses the
clinically most relevant transcript, aligning with ACGS 2024 §4.

Reference
---------
Morales J, et al. "A joint NCBI and EMBL-EBI transcript set for clinical
genomics and research." Nature. 2022;604:310–315. PMID:35356062.
DOI:10.1038/s41586-022-04558-8

McLaren W, et al. "The Ensembl Variant Effect Predictor." Genome Biology.
2016;17:122. PMID:27268795. DOI:10.1186/s13059-016-0974-4

VEP Plugins Used
----------------
- AlphaMissense   — Cheng et al. 2023 PMID:37703350
- SpliceAI        — Jaganathan et al. 2019 PMID:30661751
- Pangolin        — Zeng et al. 2022 PMID:35190963
- dbNSFP_4.7      — Liu et al. 2020 PMID:33261662
- gnomAD_v4.1     — Karczewski et al. 2020 PMID:32461654 (v4.1 April 2024)

gnomAD v4.1 note: v4.0 contained an allele number (AN) bug for certain
variants in non-EUR populations; v4.1 (released April 2024) corrects this.
Always use v4.1 or later for clinical reporting.
"""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

# Default VEP pick order per ACGS 2024 §4 / Morales 2022
VEP_PICK_ORDER: str = "mane_select,mane_plus_clinical,canonical"

# VEP v111 mandatory flags for clinical-grade annotation
VEP_CORE_FLAGS: list[str] = [
    "--everything",           # Enable all consequence types
    "--pick",                 # One consequence per variant
    f"--pick_order {VEP_PICK_ORDER}",
    "--format vcf",
    "--output_file STDOUT",
    "--json",                 # Machine-readable JSON output
    "--no_stats",             # Skip HTML statistics (faster)
    "--assembly GRCh38",
    "--fasta /data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "--offline",              # Use local cache (no API calls)
    "--cache",
    "--merged",               # Merged Ensembl/RefSeq cache
    "--hgvs",                 # HGVS nomenclature
    "--hgvsg",                # Genomic HGVS
    "--protein",              # Protein consequence
    "--symbol",               # HGNC gene symbol
    "--canonical",            # Mark canonical transcript
    "--mane",                 # Mark MANE Select/Plus Clinical
    "--af_gnomade",           # gnomAD exome AF (v4.1 via cache)
    "--af_gnomadg",           # gnomAD genome AF (v4.1 via cache)
    "--sift b",               # SIFT score + prediction
    "--polyphen b",           # PolyPhen score + prediction
    "--domains",              # Protein domain annotation
    "--regulatory",           # Regulatory feature overlap
    "--numbers",              # Exon/intron numbering
]

# SpliceAI score threshold for reporting (≥0.2 flag; ≥0.5 clinical significance)
SPLICEAI_HIGH_THRESHOLD: float = 0.5
SPLICEAI_LOW_THRESHOLD: float = 0.2


@dataclass
class AnnotatedVariant:
    """Full annotation record for a single variant consequence.

    Attributes:
        chrom: Chromosome (GRCh38, "chr"-prefixed).
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
        gene_symbol: HGNC gene symbol.
        transcript_id: Ensembl transcript ID (ENST…).
        is_mane_select: True if this consequence is on the MANE Select transcript.
        is_mane_plus_clinical: True if MANE Plus Clinical (not Select).
        hgvsc: HGVSc notation, e.g. "NM_007294.4:c.5266dupC".
        hgvsp: HGVSp notation (one-letter), e.g. "p.Gln1756ProfsTer25".
        consequence_terms: List of SO consequence terms.
        impact: VEP impact (HIGH/MODERATE/LOW/MODIFIER).
        revel_score: REVEL score from dbNSFP plugin.
        alphamissense_score: AlphaMissense score from plugin.
        spliceai_ds_max: Maximum SpliceAI delta score across four categories.
        cadd_phred: CADD Phred score.
        gnomad_af: gnomAD v4.1 genome allele frequency.
        gnomad_ac: gnomAD v4.1 allele count.
        gnomad_an: gnomAD v4.1 allele number.
        clinvar_id: ClinVar variation ID if present.
        raw: Full VEP JSON consequence dict for audit trail.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: Optional[str] = None
    transcript_id: Optional[str] = None
    is_mane_select: bool = False
    is_mane_plus_clinical: bool = False
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    consequence_terms: list[str] = field(default_factory=list)
    impact: Optional[str] = None
    revel_score: Optional[float] = None
    alphamissense_score: Optional[float] = None
    spliceai_ds_max: Optional[float] = None
    cadd_phred: Optional[float] = None
    gnomad_af: Optional[float] = None
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None
    clinvar_id: Optional[str] = None
    raw: dict[str, Any] = field(default_factory=dict)


class VEPRunner:
    """Wrapper for running Ensembl VEP v111 as a subprocess.

    Constructs the VEP command line with all required clinical annotation
    plugins, executes it, and parses the JSON output into AnnotatedVariant
    objects.  MANE Select transcript is always preferred (Morales 2022).

    Args:
        vep_binary: Path to the ``vep`` executable (default: "vep", assumes
            it is on PATH).
        cache_dir: Path to the VEP offline cache directory.
        plugin_dir: Path to the VEP plugins directory.
        alphamissense_tsv: Path to the AlphaMissense TSV for the plugin.
        dbnsfp_db: Path to the dbNSFP v4.7 database file.
        spliceai_snv: Path to SpliceAI SNV scores file.
        spliceai_indel: Path to SpliceAI indel scores file.
        extra_flags: Additional VEP flags to append.

    Example:
        >>> runner = VEPRunner(cache_dir="/data/vep_cache")
        >>> variants = runner.run_vep(
        ...     vcf_path=Path("/tmp/sample.vcf"),
        ...     output_path=Path("/tmp/sample.vep.json"),
        ... )
        >>> for v in variants:
        ...     print(v.gene_symbol, v.hgvsc, v.alphamissense_score)
    """

    def __init__(
        self,
        vep_binary: str = "vep",
        cache_dir: Optional[str | Path] = None,
        plugin_dir: Optional[str | Path] = None,
        alphamissense_tsv: Optional[str | Path] = None,
        dbnsfp_db: Optional[str | Path] = None,
        spliceai_snv: Optional[str | Path] = None,
        spliceai_indel: Optional[str | Path] = None,
        extra_flags: Optional[list[str]] = None,
    ) -> None:
        self._vep = vep_binary
        self._cache_dir = Path(cache_dir) if cache_dir else None
        self._plugin_dir = Path(plugin_dir) if plugin_dir else None
        self._am_tsv = Path(alphamissense_tsv) if alphamissense_tsv else None
        self._dbnsfp = Path(dbnsfp_db) if dbnsfp_db else None
        self._spliceai_snv = Path(spliceai_snv) if spliceai_snv else None
        self._spliceai_indel = Path(spliceai_indel) if spliceai_indel else None
        self._extra_flags = extra_flags or []

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_vep(
        self,
        vcf_path: Path,
        output_path: Optional[Path] = None,
    ) -> list[AnnotatedVariant]:
        """Run VEP v111 on a VCF file and return parsed annotations.

        Args:
            vcf_path: Path to input VCF (uncompressed or bgzipped).
            output_path: Optional path for the JSON output file.
                If None, a temporary file is used.

        Returns:
            List of AnnotatedVariant objects, one per variant-consequence pair.
            Where a variant has multiple consequences the MANE Select
            transcript consequence is returned first.

        Raises:
            subprocess.CalledProcessError: If VEP exits with non-zero status.
            FileNotFoundError: If vcf_path does not exist.
        """
        if not vcf_path.exists():
            raise FileNotFoundError(f"Input VCF not found: {vcf_path}")

        # Use a temp file for output if none specified
        tmp: Optional[tempfile.NamedTemporaryFile] = None
        if output_path is None:
            tmp = tempfile.NamedTemporaryFile(suffix=".vep.json", delete=False)
            output_path = Path(tmp.name)
            tmp.close()

        cmd = self._build_command(vcf_path, output_path)
        logger.info("Running VEP: %s", " ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            logger.error("VEP stderr: %s", result.stderr[:2000])
            raise subprocess.CalledProcessError(
                result.returncode, cmd, result.stdout, result.stderr
            )

        logger.debug("VEP completed. Parsing output: %s", output_path)
        variants = self._parse_vep_json(output_path)

        # Clean up temp file
        if tmp is not None:
            output_path.unlink(missing_ok=True)

        return variants

    # ------------------------------------------------------------------
    # Command construction
    # ------------------------------------------------------------------

    def _build_command(self, vcf_path: Path, output_path: Path) -> list[str]:
        """Assemble the VEP command-line arguments.

        Args:
            vcf_path: Input VCF path.
            output_path: Output JSON path.

        Returns:
            List of command tokens ready for subprocess.run.
        """
        cmd = [self._vep]

        # Core flags
        cmd += [
            "--input_file", str(vcf_path),
            "--output_file", str(output_path),
            "--format", "vcf",
            "--json",
            "--no_stats",
            "--assembly", "GRCh38",
            "--everything",
            "--pick",
            "--pick_order", VEP_PICK_ORDER,
            "--hgvs",
            "--hgvsg",
            "--protein",
            "--symbol",
            "--canonical",
            "--mane",
            "--numbers",
            "--domains",
            "--regulatory",
            "--sift", "b",
            "--polyphen", "b",
            "--af_gnomade",
            "--af_gnomadg",
            "--offline",
            "--cache",
            "--merged",
        ]

        if self._cache_dir:
            cmd += ["--dir_cache", str(self._cache_dir)]
        if self._plugin_dir:
            cmd += ["--dir_plugins", str(self._plugin_dir)]

        # Plugin: AlphaMissense (Cheng et al. 2023 PMID:37703350)
        if self._am_tsv and self._am_tsv.exists():
            cmd += ["--plugin", f"AlphaMissense,file={self._am_tsv}"]

        # Plugin: dbNSFP v4.7 (REVEL, BayesDel, CADD, ESM1b)
        if self._dbnsfp and self._dbnsfp.exists():
            dbnsfp_cols = (
                "REVEL_score,BayesDel_noAF_score,BayesDel_addAF_score,"
                "CADD_phred,CADD_raw,AlphaMissense_score,ESM1b_score"
            )
            cmd += ["--plugin", f"dbNSFP,{self._dbnsfp},{dbnsfp_cols}"]

        # Plugin: SpliceAI (Jaganathan et al. 2019 PMID:30661751)
        if self._spliceai_snv and self._spliceai_snv.exists():
            snv_arg = f"snv={self._spliceai_snv}"
            indel_arg = (
                f",indel={self._spliceai_indel}"
                if self._spliceai_indel and self._spliceai_indel.exists()
                else ""
            )
            cmd += ["--plugin", f"SpliceAI,{snv_arg}{indel_arg}"]

        # Plugin: Pangolin (Zeng et al. 2022 PMID:35190963) — splice scoring
        cmd += ["--plugin", "Pangolin"]

        # Plugin: gnomAD v4.1 custom annotation (April 2024)
        cmd += ["--plugin", "gnomADv4.1"]

        cmd.extend(self._extra_flags)
        return cmd

    # ------------------------------------------------------------------
    # JSON parsing
    # ------------------------------------------------------------------

    def _parse_vep_json(self, json_path: Path) -> list[AnnotatedVariant]:
        """Parse a VEP JSON output file into AnnotatedVariant objects.

        VEP JSON has one record per input variant, with a
        ``transcript_consequences`` array. We extract the picked consequence
        (MANE Select preferred) plus additional metadata.

        Args:
            json_path: Path to the VEP JSON output file.

        Returns:
            List of AnnotatedVariant objects.
        """
        variants: list[AnnotatedVariant] = []

        with open(json_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    logger.warning("Skipping malformed VEP JSON line: %s", exc)
                    continue

                variant = self._parse_record(record)
                if variant:
                    variants.append(variant)

        logger.info("Parsed %d annotated variants from VEP output", len(variants))
        return variants

    def _parse_record(self, record: dict[str, Any]) -> Optional[AnnotatedVariant]:
        """Extract one AnnotatedVariant from a VEP JSON record.

        The VEP --pick flag selects one consequence; we look for the
        picked consequence (``pick == 1``) first, then fall back to the
        first MANE Select consequence.

        Args:
            record: Parsed VEP JSON object for one variant.

        Returns:
            AnnotatedVariant or None if no parseable consequence found.
        """
        chrom = record.get("seq_region_name", "")
        pos = int(record.get("start", 0))
        ref = record.get("allele_string", "/").split("/")[0]
        alt = record.get("allele_string", "/").split("/")[-1]

        # Normalise chromosome
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        tcs: list[dict[str, Any]] = record.get("transcript_consequences", [])
        if not tcs:
            return None

        # Prefer the picked consequence, then MANE Select, then first
        picked = next((t for t in tcs if t.get("pick") == 1), None)
        mane = next((t for t in tcs if t.get("mane_select")), None)
        tc = picked or mane or tcs[0]

        # Extract plugin annotations from extra field (VEP appends under 'extras')
        extras: dict[str, Any] = record.get("extras", {})

        # AlphaMissense score from plugin
        am_score = _safe_float(tc.get("alphamissense_score") or extras.get("AlphaMissense"))

        # SpliceAI maximum delta score
        spliceai_max = _extract_spliceai_max(tc, extras)

        # dbNSFP scores
        revel = _safe_float(tc.get("revel_score") or extras.get("REVEL_score"))
        cadd = _safe_float(tc.get("cadd_phred") or extras.get("CADD_phred"))

        # gnomAD v4.1 genome AF (preferred over exome for population freq)
        gnomad_af = _safe_float(tc.get("gnomADg_AF") or record.get("gnomADg_AF"))
        gnomad_ac = _safe_int(tc.get("gnomADg_AC") or record.get("gnomADg_AC"))
        gnomad_an = _safe_int(tc.get("gnomADg_AN") or record.get("gnomADg_AN"))

        # ClinVar variation ID
        clinvar_id = record.get("colocated_variants", [{}])[0].get("var_synonyms", {}).get(
            "ClinVar", [None]
        )
        if isinstance(clinvar_id, list):
            clinvar_id = clinvar_id[0] if clinvar_id else None

        return AnnotatedVariant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gene_symbol=tc.get("gene_symbol"),
            transcript_id=tc.get("transcript_id"),
            is_mane_select=bool(tc.get("mane_select")),
            is_mane_plus_clinical=bool(tc.get("mane_plus_clinical")),
            hgvsc=tc.get("hgvsc"),
            hgvsp=tc.get("hgvsp"),
            consequence_terms=tc.get("consequence_terms", []),
            impact=tc.get("impact"),
            revel_score=revel,
            alphamissense_score=am_score,
            spliceai_ds_max=spliceai_max,
            cadd_phred=cadd,
            gnomad_af=gnomad_af,
            gnomad_ac=gnomad_ac,
            gnomad_an=gnomad_an,
            clinvar_id=clinvar_id,
            raw=record,
        )


# ---------------------------------------------------------------------------
# Private utilities
# ---------------------------------------------------------------------------


def _safe_float(value: Any) -> Optional[float]:
    """Convert value to float, returning None on failure.

    Args:
        value: Any value to convert.

    Returns:
        Float or None.
    """
    if value is None or value == ".":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value: Any) -> Optional[int]:
    """Convert value to int, returning None on failure.

    Args:
        value: Any value to convert.

    Returns:
        Integer or None.
    """
    if value is None or value == ".":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _extract_spliceai_max(
    tc: dict[str, Any],
    extras: dict[str, Any],
) -> Optional[float]:
    """Extract the maximum SpliceAI delta score from VEP annotation.

    SpliceAI reports four delta scores (DS_AG, DS_AL, DS_DG, DS_DL).
    The maximum is used as the summary score for clinical reporting.

    Args:
        tc: Transcript consequence dict from VEP JSON.
        extras: Extras dict from VEP record level.

    Returns:
        Maximum delta score or None.
    """
    ds_keys = ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
               "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]
    scores: list[float] = []

    for key in ds_keys:
        raw = tc.get(key) or extras.get(key)
        val = _safe_float(raw)
        if val is not None:
            scores.append(val)

    return max(scores) if scores else None
