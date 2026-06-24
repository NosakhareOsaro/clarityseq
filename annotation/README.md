# annotation/

Variant annotation clients and utilities.

## Tool versions and guidelines

| Tool | Version | Evidence supported | Reference |
|------|---------|-------------------|-----------|
| VEP | v111 | All annotations; MANE Select priority | PMID:35356062 |
| AlphaMissense | 2023 scores | PP3 (≥0.564) / BP4 (≤0.340) PRIMARY | PMID:37703350 |
| dbNSFP | v4.7 | REVEL, BayesDel, CADD, ESM1b | Released May 2024 |
| gnomAD | v4.1 | PM2, BA1, BS1 allele frequencies | April 2024 |
| SpliceAI | Latest | PP3/BP4/BP7 for splice variants | PMID:30267119 |
| Pangolin | Latest | Tissue-specific splicing | Supplements SpliceAI |
| ClinVar | Current | PP5/BP6 (use with caution) | — |

## Key requirements

- **gnomAD**: Always use v4.1, NOT v4.0. gnomAD v4.0 had an allele number calculation bug.
- **VEP pick order**: `mane_select,mane_plus_clinical,canonical` (MANE Select priority)
- **AlphaMissense**: PRIMARY PP3/BP4 tool per ClinGen SVI 2024
- **PM2 threshold**: Apply at SUPPORTING (1 pt) — see bayesacmg/src/bayesacmg/rules/pathogenic.py
