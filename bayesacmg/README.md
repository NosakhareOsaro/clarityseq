# BayesACMG

Bayesian ACMG/AMP variant classifier implementing:

- **ACGS Best Practice Guidelines 2024 v1.2** (Durkie et al., Feb 2024)
- **Richards et al. 2015** PMID:25741868 — original 28-rule framework
- **Tavtigian et al. 2020** PMID:32645316 — Bayesian point-scoring system
- **ClinGen SVI 2024** recommendations (PM2→Supporting; AlphaMissense approved)
- **Walker et al. 2023** PMID:36898414 — splicing subgroup PP3/BP4/BP7
- **ACGS 2024 §6** — mitochondrial variant classification

## Key changes from ACGS 2020

| Rule | ACGS 2020 | ACGS 2024 / ClinGen SVI 2024 |
|------|-----------|-------------------------------|
| PM2 weight | Moderate (2 pts) | **Supporting (1 pt)** — gnomAD v4.1 rationale |
| PP3/BP4 tool | REVEL primary | **AlphaMissense primary** (≥0.564 / ≤0.340) |
| Splicing | Ad hoc | **Walker 2023 framework** (SpliceAI ≥0.5 → PP3 Strong) |
| Mito variants | Generic rules | **ACGS 2024 §6** (haplogroup first; heteroplasmy %) |

## Install

```bash
cd bayesacmg
pip install -e ".[dev]"
```

## Quickstart

```python
from bayesacmg import classify_variant
result = classify_variant(
    chrom="17", pos=41276045, ref="G", alt="T",
    gnomad_af=0.0, alphamissense_score=0.82
)
print(result.classification)  # "Likely Pathogenic"
print(result.posterior_probability)  # 0.91
print(result.credible_interval_95)   # (0.78, 0.97)
```

## Test

```bash
pytest bayesacmg/tests/ -v --cov=bayesacmg --cov-fail-under=90
```
