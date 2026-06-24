# reporting/

NHS GMS-style clinical reports (HTML + PDF + JSON-LD audit trail).

## ACGS 2024 v1.2 compliance

This module implements all ACGS 2024 reporting requirements:

1. **MANE Select transcripts** (§5; Morales et al. 2022 PMID:35356062) — all HGVSc/HGVSp uses MANE Select notation
2. **VUS review dates** (§9) — each VUS displays "Review by: [date+2yr]"; stored in reclassification daemon schedule
3. **ClinVar submission flag** (Introduction) — novel P/LP variants → pending_clinvar_submissions queue (NHS mandatory)
4. **Classification scheme citation** — "Classified per ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al.) and ACMG/AMP guidelines (Richards et al. 2015 PMID:25741868)"

## Novel contribution

**P(Pathogenic) [95% HDI] column** — Bayesian posterior probability and 95% highest-density interval from BayesACMG. Example: "LP — 78% [61%–91%]". This is the first clinical WGS report template to show ACMG classification uncertainty in a structured, quantified format.

## Templates

| Template | Description |
|----------|-------------|
| `base.html.j2` | NHS GMS-style base with CSS |
| `variants.html.j2` | Main variant table with all ACGS 2024 columns |
| `mito.html.j2` | Mitochondrial section (ACGS 2024 §6: haplogroup first) |
| `expansions.html.j2` | Repeat expansion results |
| `pgx.html.j2` | CYP2D6 pharmacogenomics |

## Audit trail

JSON-LD audit records:
- gnomAD version: "v4.1"
- VEP cache version: "111"
- AlphaMissense score file date
- BayesACMG version + PM2=Supporting note
- ACGS guidelines version: "2024 v1.2"
