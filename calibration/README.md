# calibration/

BayesACMG classifier calibration against ClinGen expert-curated variant set.

## Calibration target

**Expected Calibration Error (ECE) < 0.05** on a held-out set of ClinGen-curated variants (target: 500 variants with expert P/LP/VUS/LB/B classifications).

## Process

1. `download_clingen_curations.py` — downloads ClinGen expert panel curations from ClinGen Evidence Repository API
2. `run_calibration.py` — runs BayesACMG on downloaded variants; computes ECE; generates calibration curve plot

## ECE interpretation

ECE measures the difference between predicted probability and actual classification rate. ECE = 0.05 means predicted 80% probability corresponds to ~80% correct classification rate in the real world (±5%).

## Results

Calibration results stored in `results/`. Each run produces:
- `calibration_curve.png` — reliability diagram
- `calibration_metrics.json` — ECE, MCE, Brier score
- `per_rule_coverage.json` — how often each ACMG rule fires

## Running

```bash
# Download ClinGen curations
python calibration/download_clingen_curations.py --output calibration/clingen_curations.json

# Run calibration
python calibration/run_calibration.py --curations calibration/clingen_curations.json
```
