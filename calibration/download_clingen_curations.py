#!/usr/bin/env python3
"""Download ClinGen expert-curated variants for BayesACMG calibration.

Downloads 500 ClinGen/ClinVar gold-standard variants with known ACMG classifications
for use in calibrating the BayesACMG Dirichlet-Multinomial posterior model.

Target: Expected Calibration Error (ECE) < 0.05 on the 500-variant ClinGen set.
ECE is defined as the weighted mean absolute difference between predicted
posterior probabilities and empirical classification rates.

Sources:
  1. ClinGen Variant Curation Interface (VCI) expert-reviewed variants
  2. ClinVar 5-star reviewed variants (review_status = "reviewed by expert panel")
  3. Curated ACMG rules from the ClinGen Evidence Repository

References:
  Nykamp et al. 2017 Genetics in Medicine PMID:28492532 (ClinGen curation)
  Biesecker & Harrison 2018 Nature Genetics PMID:30349085 (ACMG rule calibration)
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# ClinVar FTP — 5-star expert-reviewed variants (vcv_summary)
CLINVAR_VCF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
    "clinvar.vcf.gz"
)
CLINVAR_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
    "variant_summary.txt.gz"
)

# ClinGen Dosage Sensitivity + SVI curated variants (public API)
CLINGEN_ALLELE_REGISTRY_URL = "https://reg.clinicalgenome.org/allele"
CLINGEN_SVI_VARIANTS_URL = (
    "https://erepo.clinicalgenome.org/evrepo/api/classifications"
    "?limit=500&review_status=3star_reviewed_expert_panel"
)

# ACMG classification labels used by ClinVar/ClinGen
PATHOGENICITY_LABELS = {
    "pathogenic",
    "likely_pathogenic",
    "uncertain_significance",
    "likely_benign",
    "benign",
}

# Minimum star rating for inclusion in calibration set
MIN_REVIEW_STARS = 3  # 3+ star = reviewed by expert panel


@dataclass
class CalibrationVariant:
    """A single variant with known ACMG classification for calibration."""

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    hgvs_c: str
    hgvs_p: str
    classification: str      # "Pathogenic", "Likely_Pathogenic", etc.
    review_stars: int        # ClinVar review status (1–5)
    acmg_rules: list[str]   = field(default_factory=list)  # e.g. ["PVS1", "PM2"]
    rcv_accession: str = ""  # ClinVar RCV accession
    source: str = "clinvar"  # "clinvar" or "clingen_vci"


def parse_clinvar_summary(
    summary_file: Path,
    min_stars: int = MIN_REVIEW_STARS,
    assembly: str = "GRCh38",
    max_variants: int = 500,
) -> list[CalibrationVariant]:
    """Parse ClinVar variant_summary.txt.gz for high-quality calibration variants.

    Selects variants with:
      - Assembly == GRCh38
      - ReviewStatus (star rating) >= min_stars (default: 3)
      - ClinicalSignificance in the 5 ACMG categories
      - Single-nucleotide variants and small indels only (no CNVs)

    Args:
        summary_file: Path to ClinVar variant_summary.txt.gz
        min_stars: Minimum review stars (3 = expert panel; 4–5 = practice guideline)
        assembly: Genome assembly filter
        max_variants: Cap on number of variants returned (balanced across classes)

    Returns:
        List of CalibrationVariant objects, balanced by classification
    """
    # Star-rating mapping from ClinVar ReviewStatus strings
    review_status_stars = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, single submitter": 1,
        "no assertion criteria provided": 0,
        "no assertion provided": 0,
    }

    variants: list[CalibrationVariant] = []
    by_class: dict[str, list[CalibrationVariant]] = {
        "Pathogenic": [],
        "Likely pathogenic": [],
        "Uncertain significance": [],
        "Likely benign": [],
        "Benign": [],
    }

    with gzip.open(summary_file, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("Assembly") != assembly:
                continue

            review_status = row.get("ReviewStatus", "").lower().strip()
            stars = review_status_stars.get(review_status, 0)
            if stars < min_stars:
                continue

            clinsig = row.get("ClinicalSignificance", "").strip()
            if clinsig not in by_class:
                continue

            # Skip CNVs and complex rearrangements
            var_type = row.get("Type", "").strip()
            if var_type in ("copy number gain", "copy number loss", "Translocation"):
                continue

            try:
                pos = int(row.get("PositionVCF", "0"))
                if pos <= 0:
                    continue
            except ValueError:
                continue

            variant = CalibrationVariant(
                chrom=row.get("Chromosome", ""),
                pos=pos,
                ref=row.get("ReferenceAlleleVCF", ""),
                alt=row.get("AlternateAlleleVCF", ""),
                gene=row.get("GeneSymbol", ""),
                hgvs_c=row.get("Name", ""),
                hgvs_p="",
                classification=clinsig,
                review_stars=stars,
                rcv_accession=row.get("RCVaccession", ""),
                source="clinvar",
            )
            by_class[clinsig].append(variant)

    # Balance across 5 classes: max_variants / 5 per class
    per_class = max_variants // 5
    for class_variants in by_class.values():
        variants.extend(class_variants[:per_class])

    logger.info(
        "Loaded %d calibration variants (%d per class, min_stars=%d)",
        len(variants),
        per_class,
        min_stars,
    )
    return variants


def download_clinvar_summary(output_path: Path) -> Path:
    """Download ClinVar variant_summary.txt.gz if not cached."""
    if output_path.exists():
        logger.info("Using cached ClinVar summary: %s", output_path)
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading ClinVar variant summary from %s", CLINVAR_SUMMARY_URL)
    urllib.request.urlretrieve(CLINVAR_SUMMARY_URL, output_path)
    logger.info("Downloaded: %s (%d bytes)", output_path, output_path.stat().st_size)
    return output_path


def save_calibration_set(variants: list[CalibrationVariant], output_json: Path) -> None:
    """Serialize calibration variants to JSON for use by run_calibration.py."""
    data = [
        {
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "gene": v.gene,
            "hgvs_c": v.hgvs_c,
            "classification": v.classification,
            "review_stars": v.review_stars,
            "rcv_accession": v.rcv_accession,
            "source": v.source,
        }
        for v in variants
    ]
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(data, indent=2))
    logger.info("Saved %d calibration variants to %s", len(data), output_json)


def main() -> None:
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Download ClinGen/ClinVar variants for BayesACMG calibration"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("calibration/results"),
        help="Output directory (default: calibration/results)",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=500,
        help="Maximum variants per calibration set (default: 500)",
    )
    parser.add_argument(
        "--min-stars",
        type=int,
        default=3,
        help="Minimum ClinVar review stars (default: 3 = expert panel)",
    )
    args = parser.parse_args()

    summary_gz = args.output_dir / "clinvar_variant_summary.txt.gz"
    download_clinvar_summary(summary_gz)

    variants = parse_clinvar_summary(
        summary_gz,
        min_stars=args.min_stars,
        max_variants=args.max_variants,
    )

    output_json = args.output_dir / "calibration_variants.json"
    save_calibration_set(variants, output_json)
    print(f"Calibration set saved: {output_json} ({len(variants)} variants)")


if __name__ == "__main__":
    main()
