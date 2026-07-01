# ClinVar Submission Workflow

ClaritySeq automates NHS-mandated ClinVar submission for novel P/LP variants.

## Regulatory basis

**ACGS 2024 v1.2 Introduction** (Durkie et al., February 2024):
> "Submission of variants to ClinVar by NHS laboratories in England is now a requirement following completion of the information governance review process."

This replaced the previous voluntary guidance in ACGS 2020.

## Prerequisites

1. **NCBI account** with ClinVar submission access
2. **NCBI API key** — set `NCBI_API_KEY` in `.env`
3. **Org ID** — set `CLINVAR_SUBMISSION_ORG_ID` in `.env` (assigned by NCBI)

## What triggers a submission

Automatic submission is triggered for variants meeting **all** of:
- ACMG class: **P (Pathogenic)** or **LP (Likely Pathogenic)**
- BayesACMG posterior probability: **≥ 0.90**
- Not already in ClinVar (no RCV or SCV accession found)
- Variant passed quality filters (VQSR PASS; coverage ≥ 30×)

## What is submitted

Each submission includes:
- **Gene symbol + MANE Select HGVSc/HGVSp** (e.g., `BRCA2:NM_000059.4:c.6406AG>T (p.Glu2136Ter)`)
- **ACMG 5-tier classification** (P/LP/VUS/LB/B)
- **Evidence codes** applied (ACGS 2024 v1.2 criteria)
- **BayesACMG posterior probability** (as submission note)
- **gnomAD v4.1 allele frequency** (or "absent")
- **Phenotype** from Phenopackets v2 (if available; OMIM/HPO terms)

## Submission pipeline

```
Pipeline completion
    → BayesACMG P ≥ 0.90
    → ClinVar lookup: no existing SCV?
    → INSERT INTO clinvar_submission_queue
    → Celery task: reclassification.clinvar_submitter.submit_variant()
    → NCBI ClinVar API submission
    → Store SCV accession + submission date
    → Notify via Slack/email
```

## Manual review

All automatic submissions are flagged for human review before NCBI processes them. The `clinvar_submission_queue` table has a `requires_review` flag (default: True). Clinical scientists review and approve via the audit interface.

## Testing

```bash
# Dry-run submission (does not call NCBI)
python -m reclassification.clinvar_submitter --dry-run --variant-id 12345

# Test with mock NCBI API
pytest reclassification/tests/test_clinvar_submitter.py -v
```

## Monitoring

Check submission status:
```bash
python -m reclassification.clinvar_submitter --check-status --submission-id SUB123456
```
