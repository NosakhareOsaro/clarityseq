#!/usr/bin/env python3
"""Run BayesACMG calibration against ClinGen/ClinVar gold-standard variants.

Evaluates BayesACMG posterior probability calibration using the 500-variant
ClinGen set prepared by download_clingen_curations.py.

Target: Expected Calibration Error (ECE) < 0.05.
ECE = Σ_b (|b| / N) * |acc(b) - conf(b)|
  where b is a confidence bin, |b| is the number of variants in the bin,
  N is the total variants, acc(b) is the empirical accuracy, and conf(b)
  is the mean predicted probability in that bin.

A well-calibrated model has ECE < 0.05 (Guo et al. 2017).
BayesACMG targets ECE < 0.05 per the PROJECT_GUIDE specification.

References:
  Guo et al. 2017 ICML (calibration of neural networks / ECE definition)
  Nykamp et al. 2017 PMID:28492532 (ClinGen calibration methodology)
  Tavtigian et al. 2018 PMID:29300386 (Bayesian ACMG framework)
"""

from __future__ import annotations

import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

# ECE threshold per PROJECT_GUIDE §7 (BayesACMG calibration target)
ECE_THRESHOLD = 0.05
N_CALIBRATION_BINS = 10  # Decile bins for ECE computation


@dataclass
class CalibrationMetrics:
    """Summary metrics from a calibration run."""

    n_variants: int
    ece: float                      # Expected Calibration Error
    accuracy: float                 # Overall classification accuracy
    pathogenic_sensitivity: float   # Sensitivity for P/LP
    benign_specificity: float       # Specificity for B/LB
    passed: bool                    # ECE < ECE_THRESHOLD


def _classification_to_binary(classification: str) -> int | None:
    """Map 5-tier ACMG classification to binary P/LP=1, B/LB=0, VUS=None."""
    cls = classification.lower().replace(" ", "_")
    if cls in ("pathogenic", "likely_pathogenic"):
        return 1
    if cls in ("benign", "likely_benign"):
        return 0
    return None  # VUS excluded from binary calibration metrics


def compute_ece(
    probabilities: list[float],
    labels: list[int],
    n_bins: int = N_CALIBRATION_BINS,
) -> float:
    """Compute Expected Calibration Error (ECE) with equal-width bins.

    Args:
        probabilities: Predicted P(Pathogenic) for each variant (0.0–1.0)
        labels: True binary labels (1=P/LP, 0=B/LB)
        n_bins: Number of equal-width bins (default: 10)

    Returns:
        ECE value (lower = better calibrated; target < 0.05)
    """
    n = len(probabilities)
    if n == 0:
        return 0.0

    bin_width = 1.0 / n_bins
    ece = 0.0

    for i in range(n_bins):
        lo = i * bin_width
        hi = (i + 1) * bin_width

        # Variants whose predicted probability falls in this bin
        bin_mask = [j for j, p in enumerate(probabilities) if lo <= p < hi]
        if not bin_mask:
            continue

        bin_size = len(bin_mask)
        mean_confidence = sum(probabilities[j] for j in bin_mask) / bin_size
        empirical_accuracy = sum(labels[j] for j in bin_mask) / bin_size

        ece += (bin_size / n) * abs(mean_confidence - empirical_accuracy)

    return round(ece, 6)


def run_bayesacmg_on_variants(
    calibration_json: Path,
) -> tuple[list[float], list[int], list[str]]:
    """Run BayesACMG on all calibration variants and return predictions.

    Imports BayesACMG programmatically — requires `pip install -e bayesacmg/`.

    Returns:
        Tuple of:
          - probabilities: list[float] — P(Pathogenic) posterior for each variant
          - labels: list[int] — true binary labels (1=P/LP, 0=B/LB)
          - classifications: list[str] — predicted 5-tier classifications
    """
    # Import BayesACMG (installed as editable package)
    try:
        from bayesacmg.model import BayesACMGClassifier
        from bayesacmg.models import VariantInput
    except ImportError as e:
        raise ImportError(
            "BayesACMG not installed. Run: pip install -e bayesacmg/"
        ) from e

    variants_data = json.loads(calibration_json.read_text())
    classifier = BayesACMGClassifier()

    probabilities: list[float] = []
    labels: list[int] = []
    classifications: list[str] = []

    for v in variants_data:
        label = _classification_to_binary(v["classification"])
        if label is None:
            continue  # Skip VUS for binary ECE

        # Minimal VariantInput from ClinVar summary data
        variant_input = VariantInput(
            chrom=v["chrom"],
            pos=v["pos"],
            ref=v["ref"],
            alt=v["alt"],
            gene_symbol=v["gene"],
        )

        try:
            result = classifier.classify(variant_input)
            prob = result.posterior_probability
        except Exception as exc:  # noqa: BLE001
            logger.warning("BayesACMG failed for %s: %s", v.get("rcv_accession"), exc)
            continue

        probabilities.append(prob)
        labels.append(label)
        classifications.append(result.classification)

    return probabilities, labels, classifications


def evaluate_calibration(
    calibration_json: Path,
) -> CalibrationMetrics:
    """Full calibration evaluation pipeline.

    Args:
        calibration_json: Output of download_clingen_curations.py

    Returns:
        CalibrationMetrics with ECE and accuracy metrics
    """
    probabilities, labels, classifications = run_bayesacmg_on_variants(
        calibration_json
    )

    n = len(probabilities)
    if n == 0:
        logger.error("No variants classified — check BayesACMG installation")
        raise RuntimeError("No variants were successfully classified")

    ece = compute_ece(probabilities, labels)

    # Binary accuracy (threshold: prob >= 0.9 → Pathogenic, <= 0.1 → Benign)
    correct = sum(
        1 for p, l in zip(probabilities, labels, strict=True)
        if (p >= 0.9 and l == 1) or (p <= 0.1 and l == 0)
    )
    accuracy = correct / n if n > 0 else 0.0

    # Sensitivity = TP / (TP + FN) for P/LP variants
    path_variants = [(p, l) for p, l in zip(probabilities, labels, strict=True) if l == 1]
    tp = sum(1 for p, _ in path_variants if p >= 0.9)
    sensitivity = tp / len(path_variants) if path_variants else 0.0

    # Specificity = TN / (TN + FP) for B/LB variants
    benign_variants = [(p, l) for p, l in zip(probabilities, labels, strict=True) if l == 0]
    tn = sum(1 for p, _ in benign_variants if p <= 0.1)
    specificity = tn / len(benign_variants) if benign_variants else 0.0

    return CalibrationMetrics(
        n_variants=n,
        ece=ece,
        accuracy=round(accuracy, 4),
        pathogenic_sensitivity=round(sensitivity, 4),
        benign_specificity=round(specificity, 4),
        passed=ece < ECE_THRESHOLD,
    )


def main() -> int:
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Calibrate BayesACMG against ClinGen/ClinVar gold standard"
    )
    parser.add_argument(
        "--calibration-json",
        type=Path,
        default=Path("calibration/results/calibration_variants.json"),
        help="Calibration variants JSON from download_clingen_curations.py",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("calibration/results/calibration_metrics.json"),
        help="Output metrics JSON",
    )
    parser.add_argument(
        "--assert-ece",
        type=float,
        default=ECE_THRESHOLD,
        help=f"Maximum allowed ECE (default: {ECE_THRESHOLD})",
    )
    args = parser.parse_args()

    if not args.calibration_json.exists():
        logger.error(
            "Calibration variants JSON not found: %s\n"
            "Run: python calibration/download_clingen_curations.py",
            args.calibration_json,
        )
        return 2

    logger.info("Running calibration against %s", args.calibration_json)
    metrics = evaluate_calibration(args.calibration_json)

    report = {
        "n_variants": metrics.n_variants,
        "ece": metrics.ece,
        "ece_threshold": args.assert_ece,
        "accuracy": metrics.accuracy,
        "pathogenic_sensitivity": metrics.pathogenic_sensitivity,
        "benign_specificity": metrics.benign_specificity,
        "passed": metrics.passed,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2))

    print(f"\n{'='*60}")
    print("BayesACMG Calibration Results")
    print(f"{'='*60}")
    print(f"  Variants evaluated: {metrics.n_variants}")
    print(f"  ECE:                {metrics.ece:.4f}  (threshold < {args.assert_ece})")
    print(f"  Accuracy:           {metrics.accuracy*100:.1f}%")
    print(f"  P/LP sensitivity:   {metrics.pathogenic_sensitivity*100:.1f}%")
    print(f"  B/LB specificity:   {metrics.benign_specificity*100:.1f}%")
    print(f"{'='*60}")
    status = "PASS" if metrics.passed else "FAIL"
    print(f"  RESULT: {status} (ECE {'<' if metrics.passed else '>='} {args.assert_ece})")
    print(f"{'='*60}\n")

    return 0 if metrics.passed else 1


if __name__ == "__main__":
    sys.exit(main())
