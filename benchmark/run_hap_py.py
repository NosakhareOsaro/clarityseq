#!/usr/bin/env python3
"""Run hap.py benchmarking against GIAB truth sets.

Wraps hap.py v0.3.15 for automated CI benchmarking of GenomeForge variant calls
against GIAB HG001 (NA12878) chr22 truth set v4.2.1 (NIST).

ACGS 2024 §3.1 validation thresholds (CI acceptance criteria):
  SNP sensitivity   >= 99.0%
  Indel sensitivity >= 98.0%

Truth set: GIAB v4.2.1 high-confidence VCF + BED
  ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/
Reference: Krusche et al. 2019 Nature Biotechnology PMID:30988490 (hap.py)
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

GIAB_HG001_TRUTH_URL = (
    "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/"
    "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)
GIAB_HG001_BED_URL = (
    "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/"
    "HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
)
GIAB_CHR22_REGION = "chr22"


def run_hap_py(
    *,
    truth_vcf: Path,
    query_vcf: Path,
    reference_fasta: Path,
    confidence_bed: Path,
    output_prefix: Path,
    region: str = GIAB_CHR22_REGION,
    threads: int = 8,
    engine: str = "vcfeval",
) -> int:
    """Execute hap.py v0.3.15 and return the process exit code.

    Args:
        truth_vcf: GIAB truth VCF (v4.2.1)
        query_vcf: Pipeline output VCF to benchmark
        reference_fasta: GRCh38 FASTA
        confidence_bed: GIAB high-confidence BED
        output_prefix: Output prefix (hap.py creates <prefix>.summary.csv etc.)
        region: Restrict to this region (default: chr22 for CI speed)
        threads: Parallel threads
        engine: Comparison engine (vcfeval recommended; requires rtg-tools)
    """
    cmd = [
        "hap.py",
        str(truth_vcf),
        str(query_vcf),
        "-r", str(reference_fasta),
        "-f", str(confidence_bed),
        "-o", str(output_prefix),
        "--engine", engine,
        "--threads", str(threads),
        "--pass-only",  # Only evaluate PASS variants
    ]
    if region:
        cmd.extend(["-l", region])

    logger.info("Running hap.py: %s", " ".join(cmd))
    result = subprocess.run(cmd, check=False, capture_output=False)
    return result.returncode


def download_giab_truth(output_dir: Path) -> tuple[Path, Path]:
    """Download GIAB HG001 truth VCF + BED if not already cached.

    Returns:
        Tuple of (truth_vcf, confidence_bed) paths
    """
    import urllib.request

    output_dir.mkdir(parents=True, exist_ok=True)
    truth_vcf = output_dir / "HG001_GRCh38_v4.2.1_truth.vcf.gz"
    truth_bed = output_dir / "HG001_GRCh38_v4.2.1_confidence.bed"

    if not truth_vcf.exists():
        logger.info("Downloading GIAB truth VCF...")
        urllib.request.urlretrieve(GIAB_HG001_TRUTH_URL, truth_vcf)
        urllib.request.urlretrieve(GIAB_HG001_TRUTH_URL + ".tbi", str(truth_vcf) + ".tbi")

    if not truth_bed.exists():
        logger.info("Downloading GIAB confidence BED...")
        urllib.request.urlretrieve(GIAB_HG001_BED_URL, truth_bed)

    return truth_vcf, truth_bed


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Run hap.py benchmarking against GIAB HG001 (ACGS 2024 §3.1)"
    )
    parser.add_argument("--query-vcf", required=True, type=Path, help="Pipeline output VCF")
    parser.add_argument(
        "--truth-vcf",
        type=Path,
        default=None,
        help="GIAB truth VCF (downloaded automatically if omitted)",
    )
    parser.add_argument(
        "--confidence-bed",
        type=Path,
        default=None,
        help="GIAB high-confidence BED (downloaded automatically if omitted)",
    )
    parser.add_argument("--reference", required=True, type=Path, help="GRCh38 FASTA")
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix for hap.py files",
    )
    parser.add_argument(
        "--region",
        default=GIAB_CHR22_REGION,
        help=f"Genomic region (default: {GIAB_CHR22_REGION})",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=Path("/tmp/giab"),
        help="Directory for cached GIAB truth files",
    )
    args = parser.parse_args(argv)

    truth_vcf = args.truth_vcf
    confidence_bed = args.confidence_bed

    if truth_vcf is None or confidence_bed is None:
        truth_vcf, confidence_bed = download_giab_truth(args.download_dir)

    if not args.query_vcf.exists():
        logger.error("Query VCF not found: %s", args.query_vcf)
        return 2

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    return run_hap_py(
        truth_vcf=truth_vcf,
        query_vcf=args.query_vcf,
        reference_fasta=args.reference,
        confidence_bed=confidence_bed,
        output_prefix=args.output_prefix,
        region=args.region,
        threads=args.threads,
    )


if __name__ == "__main__":
    sys.exit(main())
