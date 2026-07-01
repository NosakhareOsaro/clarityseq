#!/usr/bin/env python3
"""Parse hap.py benchmark results and assert sensitivity thresholds.

CI acceptance criteria (ACGS 2024 §3.1):
  - SNP sensitivity   >= 99.0% vs GIAB HG001 chr22 truth (hap.py v0.3.15)
  - Indel sensitivity >= 98.0% vs GIAB HG001 chr22 truth (hap.py v0.3.15)

Exits 0 if both thresholds are met; exits 1 with diagnostic output otherwise.
Called from .github/workflows/ci.yml integration-test job.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class BenchmarkResult:
    """Parsed sensitivity metrics from a single hap.py summary.csv row."""

    variant_type: str  # "SNP" or "INDEL"
    truth_total: int
    tp: int
    fp: int
    fn: int
    sensitivity: float  # METRIC.Recall
    precision: float    # METRIC.Precision
    f1: float           # METRIC.F1_Score


def parse_happy_summary(summary_csv: Path) -> list[BenchmarkResult]:
    """Parse hap.py summary.csv and return SNP and INDEL results.

    hap.py summary.csv format (v0.3.15):
      Type, Filter, TRUTH.TOTAL, TRUTH.TP, TRUTH.FN, QUERY.TOTAL,
      QUERY.FP, QUERY.UNK, FP.gt, FP.al, METRIC.Recall,
      METRIC.Precision, METRIC.Frac_NA, METRIC.F1_Score, ...

    Only rows with Filter == "PASS" are used.
    """
    results: list[BenchmarkResult] = []

    with summary_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("Filter", "").strip().upper() != "PASS":
                continue
            vtype = row.get("Type", "").strip().upper()
            if vtype not in {"SNP", "INDEL"}:
                continue

            def safe_float(key: str) -> float:
                val = row.get(key, "0").strip()
                return float(val) if val not in ("", ".", "nan") else 0.0

            def safe_int(key: str) -> int:
                val = row.get(key, "0").strip()
                return int(float(val)) if val not in ("", ".", "nan") else 0

            results.append(
                BenchmarkResult(
                    variant_type=vtype,
                    truth_total=safe_int("TRUTH.TOTAL"),
                    tp=safe_int("TRUTH.TP"),
                    fp=safe_int("QUERY.FP"),
                    fn=safe_int("TRUTH.FN"),
                    sensitivity=safe_float("METRIC.Recall") * 100.0,
                    precision=safe_float("METRIC.Precision") * 100.0,
                    f1=safe_float("METRIC.F1_Score") * 100.0,
                )
            )

    return results


def find_summary_csv(benchmark_dir: Path) -> Path:
    """Locate the hap.py summary.csv in the benchmark output directory.

    hap.py outputs: <prefix>.summary.csv — we accept any file matching
    the glob *.summary.csv, preferring files named 'extended.csv' last.
    """
    candidates = sorted(benchmark_dir.glob("*.summary.csv"))
    if not candidates:
        # Fallback: accept plain summary.csv
        candidates = sorted(benchmark_dir.glob("summary.csv"))
    if not candidates:
        raise FileNotFoundError(
            f"No hap.py summary CSV found in {benchmark_dir}. "
            "Expected files matching *.summary.csv"
        )
    return candidates[0]


def emit_json_report(
    results: list[BenchmarkResult],
    output_path: Path | None,
    passed: bool,
) -> None:
    """Write a JSON benchmark report for CI artifact upload."""
    report = {
        "passed": passed,
        "metrics": [
            {
                "type": r.variant_type,
                "truth_total": r.truth_total,
                "TP": r.tp,
                "FP": r.fp,
                "FN": r.fn,
                "sensitivity_pct": round(r.sensitivity, 4),
                "precision_pct": round(r.precision, 4),
                "f1_pct": round(r.f1, 4),
            }
            for r in results
        ],
    }
    if output_path:
        output_path.write_text(json.dumps(report, indent=2))
    else:
        print(json.dumps(report, indent=2))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assert hap.py benchmark meets ACGS 2024 §3.1 thresholds"
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        metavar="DIR",
        help="Directory containing hap.py output (*.summary.csv)",
    )
    parser.add_argument(
        "--assert-snp-sensitivity",
        type=float,
        default=99.0,
        metavar="PCT",
        help="Minimum SNP sensitivity %% (default: 99.0; ACGS 2024 §3.1)",
    )
    parser.add_argument(
        "--assert-indel-sensitivity",
        type=float,
        default=98.0,
        metavar="PCT",
        help="Minimum Indel sensitivity %% (default: 98.0; ACGS 2024 §3.1)",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        metavar="FILE",
        help="Write JSON benchmark report to FILE (default: stdout)",
    )
    args = parser.parse_args(argv)

    benchmark_dir = args.input
    if not benchmark_dir.is_dir():
        print(f"ERROR: {benchmark_dir} is not a directory", file=sys.stderr)
        return 2

    try:
        summary_csv = find_summary_csv(benchmark_dir)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    results = parse_happy_summary(summary_csv)
    if not results:
        print(f"ERROR: No PASS-filtered SNP or INDEL rows in {summary_csv}", file=sys.stderr)
        return 2

    snp_result = next((r for r in results if r.variant_type == "SNP"), None)
    indel_result = next((r for r in results if r.variant_type == "INDEL"), None)

    passed = True
    failures: list[str] = []

    print(f"\n{'='*60}")
    print("ClaritySeq Benchmark Results (hap.py v0.3.15, GIAB HG001)")
    print(f"{'='*60}")
    print(f"Source: {summary_csv}")
    print()

    if snp_result:
        status = "PASS" if snp_result.sensitivity >= args.assert_snp_sensitivity else "FAIL"
        if status == "FAIL":
            passed = False
            failures.append(
                f"SNP sensitivity {snp_result.sensitivity:.2f}% < "
                f"required {args.assert_snp_sensitivity:.1f}%"
            )
        print(
            f"  SNP  sensitivity: {snp_result.sensitivity:6.2f}%  "
            f"(required ≥ {args.assert_snp_sensitivity:.1f}%)  [{status}]"
        )
        print(
            f"       precision:   {snp_result.precision:6.2f}%   "
            f"F1: {snp_result.f1:.2f}%"
        )
        print(
            f"       TP={snp_result.tp:,}  FP={snp_result.fp:,}  FN={snp_result.fn:,}  "
            f"TOTAL={snp_result.truth_total:,}"
        )
    else:
        print("  SNP: NO DATA (FAIL)")
        passed = False
        failures.append("No SNP PASS rows found in summary CSV")

    print()

    if indel_result:
        status = "PASS" if indel_result.sensitivity >= args.assert_indel_sensitivity else "FAIL"
        if status == "FAIL":
            passed = False
            failures.append(
                f"Indel sensitivity {indel_result.sensitivity:.2f}% < "
                f"required {args.assert_indel_sensitivity:.1f}%"
            )
        print(
            f"  INDEL sensitivity: {indel_result.sensitivity:6.2f}%  "
            f"(required ≥ {args.assert_indel_sensitivity:.1f}%)  [{status}]"
        )
        print(
            f"        precision:   {indel_result.precision:6.2f}%   "
            f"F1: {indel_result.f1:.2f}%"
        )
        print(
            f"        TP={indel_result.tp:,}  FP={indel_result.fp:,}  FN={indel_result.fn:,}  "
            f"TOTAL={indel_result.truth_total:,}"
        )
    else:
        print("  INDEL: NO DATA (FAIL)")
        passed = False
        failures.append("No INDEL PASS rows found in summary CSV")

    print(f"\n{'='*60}")
    if passed:
        print("RESULT: ALL THRESHOLDS MET — CI PASSES")
    else:
        print("RESULT: THRESHOLD(S) NOT MET — CI FAILS")
        for f in failures:
            print(f"  ✗ {f}")
    print(f"{'='*60}\n")

    emit_json_report(results, args.output_json, passed)

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
