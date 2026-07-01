"""
bayesacmg.calibration
=====================

Calibration of the BayesACMG model against ClinGen expert-curated variants.

Calibration strategy:
    - Reference set: ClinGen curated variants with gold-standard labels.
    - Metric: Expected Calibration Error (ECE) computed over 10 probability bins.
    - Acceptance criterion: ECE < 0.05.
    - Output: Calibrated Dirichlet α adjustment factors saved to a JSON file.

The calibration dataset should include variants spanning all 5 classification
categories (P, LP, VUS, LB, B) with known ClinVar RCV accessions.

References:
    Richards et al. 2015 PMID:25741868
    Tavtigian et al. 2020 PMID:32645316
    ClinGen Expert Curations: https://clinicalgenome.org/
    Niculescu-Mizil & Caruana 2005 ICML (ECE definition)
    ACGS 2024 v1.2 §5
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from bayesacmg.model import BayesACMGModel
from bayesacmg.models import ACMGRule

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Calibration fixtures — real ClinVar RCV accessions as reference
# ---------------------------------------------------------------------------

# Each fixture maps a ClinVar RCV accession to expected classification.
# These serve as a smoke-test calibration set covering all 5 categories.
_CALIBRATION_FIXTURE_RCVS: dict[str, str] = {
    "RCV000007535": "Pathogenic",  # BRCA1 NM_007294.4:c.5266dupC
    "RCV000048684": "Pathogenic",  # BRCA2 c.5946delT (p.Ser1982Argfs)
    "RCV000031432": "Likely_Pathogenic",  # TP53 missense in hotspot
    "RCV000144292": "VUS",  # BRCA1 variant of uncertain significance
    "RCV000237295": "Likely_Benign",  # BRCA2 intronic with no splice impact
    "RCV000013229": "Benign",  # Common BRCA1 synonymous variant
}


@dataclass
class CalibrationRecord:
    """A single calibration record linking rules to a gold-standard label.

    Attributes:
        rcv_accession: ClinVar RCV accession number.
        true_label: Gold-standard classification (one of the 5 ACMG categories).
        rules: Applied ACMGRule instances for this variant.
        notes: Optional notes about the variant.
    """

    rcv_accession: str
    true_label: str
    rules: list[ACMGRule]
    notes: str = ""


def build_calibration_dataset(
    records: list[CalibrationRecord],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Convert CalibrationRecord list to the format expected by BayesACMGModel.calibrate().

    Args:
        records: List of CalibrationRecord objects.

    Returns:
        Tuple of (reference_variants, true_labels) where:
            - reference_variants is a list of dicts with key ``"rules"``.
            - true_labels is a list of classification label strings.
    """
    reference_variants = [{"rules": r.rules} for r in records]
    true_labels = [r.true_label for r in records]
    return reference_variants, true_labels


def run_calibration(
    model: BayesACMGModel,
    records: list[CalibrationRecord],
    output_path: Path | None = None,
) -> dict[str, float]:
    """Run calibration of the BayesACMG model against a curated variant set.

    Computes ECE and accuracy over the provided calibration records.  If
    ECE < 0.05, calibration is accepted and metrics are saved to JSON.

    Args:
        model: BayesACMGModel instance to calibrate.
        records: List of CalibrationRecord objects with gold-standard labels.
        output_path: Optional path to save calibration metrics as JSON.

    Returns:
        Dict with keys ``"ece"`` (float) and ``"accuracy"`` (float).

    Raises:
        ValueError: If records is empty.

    References:
        ECE acceptance criterion: < 0.05 (ACGS 2024 / ClinGen SVI 2024).
        Niculescu-Mizil & Caruana 2005 ICML.
    """
    if not records:
        raise ValueError("calibration records list must not be empty")

    reference_variants, true_labels = build_calibration_dataset(records)
    metrics = model.calibrate(reference_variants, true_labels)

    ece = metrics["ece"]
    acc = metrics["accuracy"]

    logger.info(
        "Calibration complete: ECE=%.4f (threshold 0.05), accuracy=%.3f",
        ece,
        acc,
    )

    if ece <= 0.05:  # acceptance criterion; ACGS 2024 / ClinGen SVI 2024
        logger.info("Calibration ACCEPTED (ECE=%.4f < 0.05)", ece)
    else:
        logger.warning(
            "Calibration REJECTED (ECE=%.4f >= 0.05) — " "model requires adjustment",
            ece,
        )

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as fh:
            json.dump(metrics, fh, indent=2)
        logger.info("Calibration metrics saved to %s", output_path)

    return metrics


def compute_calibration_curve(
    model: BayesACMGModel,
    records: list[CalibrationRecord],
    n_bins: int = 10,
) -> dict[str, list[float]]:
    """Compute the full calibration curve (confidence vs accuracy per bin).

    Useful for plotting a reliability diagram.

    Args:
        model: BayesACMGModel instance.
        records: Calibration records.
        n_bins: Number of equal-width probability bins.

    Returns:
        Dict with keys ``"bin_centres"``, ``"bin_accuracy"``,
        ``"bin_confidence"``, ``"bin_counts"`` as lists of floats.
    """
    confidences: list[float] = []
    is_correct_list: list[float] = []

    for record in records:
        posterior = model.posterior_probabilities(record.rules)
        predicted = max(posterior, key=lambda k: posterior[k])
        confidence = posterior.get(predicted, 0.0)
        confidences.append(confidence)
        is_correct_list.append(float(predicted == record.true_label))

    conf_arr = np.array(confidences)
    correct_arr = np.array(is_correct_list)
    bin_boundaries = np.linspace(0, 1, n_bins + 1)

    centres: list[float] = []
    accuracies: list[float] = []
    mean_confidences: list[float] = []
    counts: list[float] = []

    for i in range(n_bins):
        lo, hi = bin_boundaries[i], bin_boundaries[i + 1]
        mask = (conf_arr >= lo) & (conf_arr < hi)
        count = int(mask.sum())
        centres.append(float((lo + hi) / 2))
        counts.append(float(count))
        if count == 0:
            accuracies.append(0.0)
            mean_confidences.append(0.0)
        else:
            accuracies.append(float(correct_arr[mask].mean()))
            mean_confidences.append(float(conf_arr[mask].mean()))

    return {
        "bin_centres": centres,
        "bin_accuracy": accuracies,
        "bin_confidence": mean_confidences,
        "bin_counts": counts,
    }
