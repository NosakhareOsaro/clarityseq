"""
multi_ancestry.vqsr_selector
==============================
Select gnomAD v4.1 ancestry-stratified VQSR training sets for GATK VQSR.

VQSR (Variant Quality Score Recalibration) uses population-specific training
variants to build a model of variant quality.  Using the wrong training set
(e.g. EUR-only for an AFR sample) inflates false positive rates.

gnomAD v4.1 provides population-stratified high-confidence variant sets
for VQSR training (April 2024 release, 807,162 individuals):
    - All populations: sites_all_v4.1.vcf.gz
    - Population-specific: sites_{AFR,AMR,EAS,EUR,MID,SAS}_v4.1.vcf.gz

Selection strategy:
    For samples with single ancestry (fraction ≥ 0.8):
        Use the population-specific training set.
    For admixed samples:
        Use the global (all-population) training set.
        Log a warning that VQSR calibration may be less optimal.

References:
    gnomAD v4.1: https://gnomad.broadinstitute.org/downloads
    Chen et al. 2024 (gnomAD v4.1 paper).
    GATK VQSR docs: https://gatk.broadinstitute.org/
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from multi_ancestry.ancestry_assigner import AncestryAssignment

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# gnomAD v4.1 resource paths (configurable via environment or caller)
# ---------------------------------------------------------------------------

_GNOMAD_V4_BASE = "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes"

# gnomAD v4.1 VQSR training site VCF paths (GCS URIs or local paths)
_GNOMAD_V4_TRAINING_SETS: dict[str, str] = {
    "AFR": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.AFR.vcf.bgz",
    "AMR": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.AMR.vcf.bgz",
    "EAS": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.EAS.vcf.bgz",
    "EUR": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.EUR.vcf.bgz",
    "MID": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.MID.vcf.bgz",
    "SAS": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.SAS.vcf.bgz",
    "ADMIXED": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.vcf.bgz",  # all-pop
    "ALL": f"{_GNOMAD_V4_BASE}/gnomad.genomes.v4.1.sites.vcf.bgz",
}


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class VQSRTrainingSet:
    """Selected VQSR training set for a sample.

    Attributes:
        sample_id: Sample identifier.
        population_label: Ancestry-assigned population label.
        training_vcf: Path or GCS URI to the selected training VCF.
        fallback: True if global training set was used (admixed sample).
        gnomad_version: gnomAD version string (``"4.1"``).
        min_af_filter: Minimum allele frequency filter applied.
        notes: Any notes about the selection (e.g. admixed fallback reason).
    """

    sample_id: str
    population_label: str
    training_vcf: str
    fallback: bool = False
    gnomad_version: str = "4.1"
    min_af_filter: float = 0.001
    notes: str = ""


# ---------------------------------------------------------------------------
# Selection functions
# ---------------------------------------------------------------------------


def select_vqsr_training(
    assignment: AncestryAssignment,
    gnomad_resource_base: str | None = None,
    min_af_filter: float = 0.001,
) -> VQSRTrainingSet:
    """Select the gnomAD v4.1 VQSR training set for an ancestry-assigned sample.

    Args:
        assignment: AncestryAssignment from ancestry_assigner.assign_ancestry().
        gnomad_resource_base: Optional override for the gnomAD v4.1 base URI.
            Useful for local testing with downloaded VCFs.
        min_af_filter: Minimum allele frequency for high-confidence training
            sites (default 0.001 = 0.1%).

    Returns:
        VQSRTrainingSet with the selected training VCF path.

    References:
        gnomAD v4.1 (April 2024): https://gnomad.broadinstitute.org/downloads
        GATK VQSR: https://gatk.broadinstitute.org/hc/en-us/articles/360035531612
    """
    training_sets = _GNOMAD_V4_TRAINING_SETS.copy()
    if gnomad_resource_base:
        training_sets = {
            pop: f"{gnomad_resource_base.rstrip('/')}/{Path(uri).name}"
            for pop, uri in training_sets.items()
        }

    label = assignment.primary_label
    is_admixed = assignment.is_admixed
    fallback = False
    notes = ""

    if label not in training_sets or is_admixed:
        label_used = "ADMIXED"
        fallback = True
        notes = (
            f"Sample {assignment.sample_id} is admixed "
            f"(components: {assignment.component_fractions}). "
            "Using global gnomAD v4.1 all-population VQSR training set. "
            "VQSR calibration may be suboptimal for admixed samples."
        )
        logger.warning(
            "Admixed sample %s: using global gnomAD v4.1 training set (fallback).",
            assignment.sample_id,
        )
    else:
        label_used = label
        notes = (
            f"gnomAD v4.1 {label} population-specific training set "
            f"(April 2024, ancestry fraction ≥ 0.80)."
        )
        logger.info(
            "Sample %s: selected gnomAD v4.1 %s VQSR training set.",
            assignment.sample_id,
            label,
        )

    training_vcf = training_sets.get(label_used, training_sets["ALL"])

    return VQSRTrainingSet(
        sample_id=assignment.sample_id,
        population_label=label_used,
        training_vcf=training_vcf,
        fallback=fallback,
        gnomad_version="4.1",
        min_af_filter=min_af_filter,
        notes=notes,
    )


def select_vqsr_training_batch(
    assignments: list[AncestryAssignment],
    gnomad_resource_base: str | None = None,
) -> list[VQSRTrainingSet]:
    """Select VQSR training sets for a batch of ancestry-assigned samples.

    Args:
        assignments: List of AncestryAssignment objects.
        gnomad_resource_base: Optional override for gnomAD v4.1 base URI.

    Returns:
        List of VQSRTrainingSet objects, one per sample.
    """
    results = [
        select_vqsr_training(a, gnomad_resource_base=gnomad_resource_base)
        for a in assignments
    ]
    n_fallback = sum(1 for r in results if r.fallback)
    logger.info(
        "VQSR training selection: %d/%d samples use global fallback (admixed).",
        n_fallback,
        len(results),
    )
    return results


def get_population_vcf_uri(
    population: str,
    gnomad_version: str = "4.1",
    gnomad_resource_base: str | None = None,
) -> str:
    """Get the gnomAD VCF URI for a population and version.

    Args:
        population: Population label (``"AFR"``, ``"EUR"``, etc.).
        gnomad_version: gnomAD version string (default ``"4.1"``).
        gnomad_resource_base: Optional override for gnomAD base URI.

    Returns:
        VCF URI string.

    References:
        gnomAD v4.1 downloads: https://gnomad.broadinstitute.org/downloads
    """
    uri = _GNOMAD_V4_TRAINING_SETS.get(population.upper(), _GNOMAD_V4_TRAINING_SETS["ALL"])
    if gnomad_resource_base:
        uri = f"{gnomad_resource_base.rstrip('/')}/{Path(uri).name}"
    return uri
