"""CYP2D6 star-allele genotyping using Cyrius 1.1.2.

Cyrius is purpose-built for CYP2D6 because GATK4 HaplotypeCaller cannot reliably
call variants in CYP2D6 due to the highly homologous CYP2D7 pseudogene.
Cyrius uses a depth-based graph model that specifically resolves the CYP2D7
interference, achieving >99% concordance vs. long-read sequencing (Twesigomwe 2022).

Star allele → metaboliser phenotype mapping follows CPIC guidelines:
  - NM  (Normal Metaboliser):    activity score 1.25–2.5
  - IM  (Intermediate Metaboliser): activity score 0.5–1.0
  - PM  (Poor Metaboliser):      activity score 0.0
  - UM  (Ultrarapid Metaboliser): activity score > 2.5

References:
  Twesigomwe et al. 2022 npj Genomic Medicine PMID:35513406 (Cyrius validation)
  Gaedigk et al. 2017 CPT PMID:28375858 (CYP2D6 activity score)
  CPIC guideline: https://cpicpgx.org/guidelines/guideline-for-cyp2d6/
"""

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

logger = logging.getLogger(__name__)

CYRIUS_VERSION = "1.1.2"


class MetaboliserPhenotype(str, Enum):
    """CYP2D6 metaboliser phenotype per CPIC activity score bins."""

    POOR = "PM"              # activity score 0.0
    INTERMEDIATE = "IM"      # activity score 0.5–1.0
    NORMAL = "NM"            # activity score 1.25–2.5
    ULTRARAPID = "UM"        # activity score > 2.5
    INDETERMINATE = "Indeterminate"  # incomplete/conflicting data


# CPIC activity score per allele — partial table (most common alleles)
# Full table: CPIC CYP2D6 allele definition table v3.2 (2023)
_ALLELE_ACTIVITY: dict[str, float] = {
    "*1": 1.0,   # Reference/wild-type
    "*2": 1.0,
    "*3": 0.0,   # Loss-of-function (frameshift)
    "*4": 0.0,   # Loss-of-function (1846G>A splice defect)
    "*5": 0.0,   # Gene deletion
    "*6": 0.0,   # Loss-of-function
    "*7": 0.0,
    "*8": 0.0,
    "*9": 0.5,   # Reduced function
    "*10": 0.25, # Reduced function (major allele in East Asian populations)
    "*14": 0.0,
    "*17": 0.5,  # Reduced function (major allele in African populations)
    "*29": 0.5,
    "*35": 1.0,
    "*36": 0.5,
    "*41": 0.5,  # Reduced function (intron 6 splice variant)
    "*xN": 2.0,  # Gene duplication (suffix for duplicated alleles; e.g. *1xN)
}

# CYP2D6 diplotype activity score bins → phenotype (CPIC 2023)
_ACTIVITY_SCORE_TO_PHENOTYPE: list[tuple[tuple[float, float], MetaboliserPhenotype]] = [
    ((0.0, 0.0),   MetaboliserPhenotype.POOR),
    ((0.25, 1.0),  MetaboliserPhenotype.INTERMEDIATE),
    ((1.25, 2.5),  MetaboliserPhenotype.NORMAL),
    ((2.5, 99.0),  MetaboliserPhenotype.ULTRARAPID),
]


@dataclass
class CYP2D6Result:
    """Genotyping result for a single sample."""

    sample_id: str
    diplotype: str            # e.g. "*1/*4", "*2xN/*5"
    activity_score: float
    phenotype: MetaboliserPhenotype
    gene_copy_number: int
    cyrius_filter: str        # PASS / FAIL / N/A
    raw_cyrius_output: dict = field(default_factory=dict)
    affected_drugs: list[str] = field(default_factory=list)


def activity_score_from_diplotype(diplotype: str) -> float:
    """Compute CPIC activity score from a CYP2D6 diplotype string.

    Handles duplications (*1xN → allele score × copy_number) and gene
    deletions (*5 → 0.0).  Unknown alleles default to 0.0 (conservative).

    Args:
        diplotype: e.g. "*1/*4", "*2xN/*5", "*10/*17"
    """
    if diplotype in ("Indeterminate", "N/A", ""):
        return -1.0

    alleles = diplotype.split("/")
    if len(alleles) != 2:
        logger.warning("Unexpected diplotype format: %s", diplotype)
        return -1.0

    total = 0.0
    for allele in alleles:
        allele = allele.strip()
        # Handle duplications: *1xN means N copies of *1
        if "xN" in allele:
            base = allele.split("x")[0]
            base_score = _ALLELE_ACTIVITY.get(base, 0.0)
            # Cyrius outputs xN as ≥2 copies; use 2 as minimum
            total += base_score * 2.0
        else:
            total += _ALLELE_ACTIVITY.get(allele, 0.0)

    return round(total, 2)


def phenotype_from_activity(activity_score: float) -> MetaboliserPhenotype:
    """Map CPIC activity score to metaboliser phenotype.

    Bins per CPIC CYP2D6 guideline (2023 update):
      0.0        → PM (Poor Metaboliser)
      0.25–1.0   → IM (Intermediate Metaboliser)
      1.25–2.5   → NM (Normal Metaboliser)
      > 2.5      → UM (Ultrarapid Metaboliser)
    """
    if activity_score < 0:
        return MetaboliserPhenotype.INDETERMINATE
    for (lo, hi), phenotype in _ACTIVITY_SCORE_TO_PHENOTYPE:
        if lo <= activity_score <= hi:
            return phenotype
    return MetaboliserPhenotype.INDETERMINATE


def run_cyrius(
    bam: Path,
    reference: Path,
    output_dir: Path,
    sample_id: str,
    genome_build: str = "hg38",
) -> Path:
    """Execute Cyrius 1.1.2 on a BAM file and return the output JSON path.

    Args:
        bam: Coordinate-sorted, indexed BAM (DRAGMAP or BWA-MEM2 aligned)
        reference: GRCh38 FASTA (must match BAM)
        output_dir: Directory for Cyrius output files
        sample_id: Sample identifier (used in output filenames)
        genome_build: "hg38" (GRCh38) or "hg19" — use hg38

    Returns:
        Path to Cyrius output JSON file

    Raises:
        subprocess.CalledProcessError: if Cyrius exits non-zero
        RuntimeError: if output JSON not found after successful run
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_json = output_dir / f"{sample_id}.json"

    if output_json.exists():
        logger.info("Cyrius output already exists, skipping: %s", output_json)
        return output_json

    cmd = [
        "python3", "-m", "stargazer",  # Cyrius entry point
        "--manifest", str(bam),
        "--genome", genome_build,
        "--prefix", str(output_dir / sample_id),
        "--reference", str(reference),
    ]
    logger.info("Running Cyrius %s: %s", CYRIUS_VERSION, " ".join(cmd))
    subprocess.run(cmd, check=True, capture_output=False)

    if not output_json.exists():
        raise RuntimeError(
            f"Cyrius completed but output JSON not found: {output_json}"
        )
    return output_json


def parse_cyrius_output(cyrius_json: Path, sample_id: str) -> CYP2D6Result:
    """Parse Cyrius JSON output into a CYP2D6Result.

    Cyrius JSON format:
      {
        "Sample": "SAMPLE1",
        "Genotype": "*1/*4",
        "Filter": "PASS",
        "Copy_Number": 2
      }
    """
    raw = json.loads(cyrius_json.read_text())

    diplotype = raw.get("Genotype", "Indeterminate")
    cyrius_filter = raw.get("Filter", "N/A")
    copy_number = int(raw.get("Copy_Number", 2))

    activity_score = activity_score_from_diplotype(diplotype)
    phenotype = phenotype_from_activity(activity_score)

    return CYP2D6Result(
        sample_id=sample_id,
        diplotype=diplotype,
        activity_score=max(activity_score, 0.0),
        phenotype=phenotype,
        gene_copy_number=copy_number,
        cyrius_filter=cyrius_filter,
        raw_cyrius_output=raw,
    )


def genotype_sample(
    bam: Path,
    reference: Path,
    output_dir: Path,
    sample_id: str,
) -> CYP2D6Result:
    """End-to-end CYP2D6 genotyping: run Cyrius + parse output.

    This is the main entry point for the Nextflow cyrius/main.nf module.
    Results feed into the pgx/cpic_client.py for drug dosing recommendations.

    Args:
        bam: Aligned BAM (DRAGMAP preferred; BWA-MEM2 also supported)
        reference: GRCh38 reference FASTA
        output_dir: Output directory
        sample_id: Sample name

    Returns:
        CYP2D6Result with diplotype, activity score, and metaboliser phenotype
    """
    cyrius_json = run_cyrius(bam, reference, output_dir, sample_id)
    result = parse_cyrius_output(cyrius_json, sample_id)

    logger.info(
        "CYP2D6 result for %s: %s (activity=%.2f, phenotype=%s, filter=%s)",
        sample_id,
        result.diplotype,
        result.activity_score,
        result.phenotype.value,
        result.cyrius_filter,
    )
    return result
