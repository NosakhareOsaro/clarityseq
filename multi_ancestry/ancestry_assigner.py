"""
multi_ancestry.ancestry_assigner
===================================
Assign population labels to samples from somalier ancestry results.

For genetically admixed individuals, provides a fallback strategy:
    1. Assign primary ancestry if a single population fraction ≥ 0.8.
    2. For admixed samples (no single fraction ≥ 0.8), assign
       ``"ADMIXED"`` label with component fractions listed.
    3. Log admixed samples for VQSR training set selection.

Population label mapping (somalier → gnomAD v4.1):
    AFR → African/African American
    AMR → Latino/Admixed American
    EAS → East Asian
    EUR → European (Finnish + non-Finnish)
    MID → Middle Eastern (gnomAD v4 new label)
    SAS → South Asian
    ADMIXED → admixed/multi-ancestry (use global VQSR training set)

gnomAD v4.1 ancestry groups (April 2024):
    The v4 release introduced MID (Middle Eastern) as a distinct group.
    Reference: Chen et al. 2024 (gnomAD v4 paper).

References:
    Pedersen et al. 2020 Genome Biology PMID:32620139 (somalier).
    Chen et al. 2024 (gnomAD v4.1).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from multi_ancestry.somalier_runner import SomalierAncestryResult

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Population label constants
# ---------------------------------------------------------------------------

_POPULATION_LABELS = frozenset({"AFR", "AMR", "EAS", "EUR", "MID", "SAS", "ADMIXED"})

# Threshold for single-ancestry assignment
_SINGLE_ANCESTRY_THRESHOLD = 0.80  # 80% from one population → assign that label

# Full names for reporting
_POP_FULL_NAMES: dict[str, str] = {
    "AFR": "African/African American",
    "AMR": "Latino/Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "MID": "Middle Eastern",
    "SAS": "South Asian",
    "ADMIXED": "Admixed/Multi-ancestry",
}


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class AncestryAssignment:
    """Population label assignment for a sample.

    Attributes:
        sample_id: Sample identifier.
        primary_label: Assigned population label.
        full_name: Full population name for reporting.
        is_admixed: True if no single ancestry fraction ≥ 0.8.
        component_fractions: Dict of population → fraction for admixed samples.
        vqsr_training_set: Recommended gnomAD v4.1 VQSR training set identifier.
    """

    sample_id: str
    primary_label: str
    full_name: str
    is_admixed: bool
    component_fractions: dict[str, float] = field(default_factory=dict)
    vqsr_training_set: str = ""


# ---------------------------------------------------------------------------
# Assignment function
# ---------------------------------------------------------------------------


def assign_ancestry(
    somalier_result: SomalierAncestryResult,
    threshold: float = _SINGLE_ANCESTRY_THRESHOLD,
) -> AncestryAssignment:
    """Assign a population label from a somalier ancestry result.

    Args:
        somalier_result: SomalierAncestryResult from somalier ancestry command.
        threshold: Minimum fraction for single-ancestry assignment
            (default 0.8; 80% from one population).

    Returns:
        AncestryAssignment with primary label and VQSR training set recommendation.
    """
    fractions = somalier_result.ancestry_fractions

    # Find population with highest fraction
    if fractions:
        top_pop = max(fractions, key=lambda p: fractions[p])
        top_fraction = fractions[top_pop]
    else:
        # Fall back to somalier's predicted ancestry
        top_pop = somalier_result.predicted_ancestry or "ADMIXED"
        top_fraction = somalier_result.predicted_ancestry_p

    if top_fraction >= threshold:
        # Single ancestry assignment
        label = top_pop.upper()
        if label not in _POPULATION_LABELS:
            label = "ADMIXED"
        is_admixed = False
        components = {}
    else:
        # Admixed: no single dominant ancestry
        label = "ADMIXED"
        is_admixed = True
        # Include all fractions above 10% for reporting
        components = {
            pop: frac
            for pop, frac in fractions.items()
            if frac >= 0.10
        }

    # Recommend VQSR training set
    vqsr_training = _recommend_vqsr_training(label)

    assignment = AncestryAssignment(
        sample_id=somalier_result.sample_id,
        primary_label=label,
        full_name=_POP_FULL_NAMES.get(label, label),
        is_admixed=is_admixed,
        component_fractions=components,
        vqsr_training_set=vqsr_training,
    )

    if is_admixed:
        logger.info(
            "Sample %s assigned as ADMIXED (top fraction %.2f < threshold %.2f). "
            "Components: %s",
            somalier_result.sample_id,
            top_fraction,
            threshold,
            components,
        )
    else:
        logger.info(
            "Sample %s assigned ancestry: %s (fraction %.2f)",
            somalier_result.sample_id,
            label,
            top_fraction,
        )

    return assignment


def assign_ancestry_batch(
    somalier_results: list[SomalierAncestryResult],
    threshold: float = _SINGLE_ANCESTRY_THRESHOLD,
) -> list[AncestryAssignment]:
    """Assign population labels to a batch of somalier ancestry results.

    Args:
        somalier_results: List of SomalierAncestryResult objects.
        threshold: Minimum fraction for single-ancestry assignment.

    Returns:
        List of AncestryAssignment objects, one per sample.
    """
    assignments = [assign_ancestry(r, threshold=threshold) for r in somalier_results]
    n_admixed = sum(1 for a in assignments if a.is_admixed)
    logger.info(
        "Ancestry assignment complete: %d/%d samples are admixed.",
        n_admixed,
        len(assignments),
    )
    return assignments


def _recommend_vqsr_training(label: str) -> str:
    """Recommend a gnomAD v4.1 VQSR training set for a population label.

    Args:
        label: Population label (``"AFR"``, ``"EUR"``, etc.).

    Returns:
        gnomAD v4.1 VQSR training set identifier string.

    References:
        gnomAD v4.1 VQSR resources: https://gnomad.broadinstitute.org/downloads
        Chen et al. 2024 (gnomAD v4.1).
    """
    # gnomAD v4.1 provides ancestry-stratified VQSR training sets
    vqsr_sets: dict[str, str] = {
        "AFR": "gnomad_v4.1_AFR_vqsr_training",
        "AMR": "gnomad_v4.1_AMR_vqsr_training",
        "EAS": "gnomad_v4.1_EAS_vqsr_training",
        "EUR": "gnomad_v4.1_EUR_vqsr_training",
        "MID": "gnomad_v4.1_MID_vqsr_training",
        "SAS": "gnomad_v4.1_SAS_vqsr_training",
        "ADMIXED": "gnomad_v4.1_ALL_vqsr_training",  # global training set for admixed
    }
    return vqsr_sets.get(label, "gnomad_v4.1_ALL_vqsr_training")
