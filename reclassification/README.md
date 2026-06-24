# reclassification/

ClinVar reclassification monitoring daemon + NHS-mandated ClinVar submission workflow.

## Components

| File | Description |
|------|-------------|
| `daemon.py` | Celery beat scheduler (weekly Monday 08:00 UTC) |
| `clinvar_diff.py` | Weekly ClinVar FTP diff to detect reclassifications |
| `fhir_task.py` | FHIR R4 Task resources (Genomics Reporting IG v3.0.0) |
| `clinvar_submitter.py` | NHS-mandated ClinVar API submission for novel P/LP variants |
| `notifier.py` | Slack/email notifications for reclassification events |
| `models.py` | SQLAlchemy models including ClinVarSubmissionQueue |

## Regulatory context

- **ACGS 2024 v1.2 §9**: "Laboratories must have processes for periodic reclassification. For VUS variants, review within 2 years is recommended."
- **ACGS 2024 Introduction (NHS mandate)**: "Submission of variants to ClinVar by NHS laboratories in England is now a requirement."
- **VUS review dates**: Stored in `VUSReviewSchedule` table; daemon checks weekly.

## FHIR integration

- IG version: **Genomics Reporting IG v3.0.0**
- Task type: `recontact` (IG code system)
- Variant IDs: **GA4GH VRS v2.0** (24-char computed digests)

## Running the daemon

```bash
# Start Redis first
redis-server --daemonize yes

# Start Celery worker
celery -A reclassification.daemon worker -l info

# Start Celery beat (scheduler)
celery -A reclassification.daemon beat -l info
```
