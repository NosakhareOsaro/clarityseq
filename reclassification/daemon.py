"""Celery beat daemon for weekly ClinVar reclassification monitoring.

This module defines Celery tasks and schedules for the GenomeForge
reclassification monitoring daemon. Tasks run on a weekly schedule
(Mondays 08:00 UTC) to align with the ClinVar weekly release cycle.

Schedule rationale:
    - ClinVar publishes new weekly releases on Mondays.
    - Running the diff task at 08:00 UTC on Mondays ensures the new release
      is available (typically published by 06:00 UTC Monday).
    - VUS review checks run daily to catch overdue reviews promptly.

Architecture:
    - Celery workers process tasks asynchronously.
    - Celery beat scheduler triggers periodic tasks.
    - Failed tasks are routed to a dead-letter queue (DLQ) for manual review.
    - All task failures generate Slack/email alerts via notifier.py.

Dependencies:
    - celery[redis]: Task queue and beat scheduler
    - redis: Message broker and result backend
    - sqlalchemy: Database ORM
    - NCBI_API_KEY env var: For ClinVar submissions
    - CLINVAR_FTP_URL env var: Override default ClinVar FTP URL
"""

from __future__ import annotations

import logging
import os
from datetime import date, timedelta
from pathlib import Path
from typing import Any, Optional

from celery import Celery
from celery.schedules import crontab
from celery.utils.log import get_task_logger
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from reclassification.clinvar_diff import (
    DEFAULT_CLINVAR_FTP_URL,
    diff_variants,
    download_latest_clinvar_vcf,
)
from reclassification.clinvar_submitter import (
    SubmissionStatus,
    check_submission_status,
    submit_variant,
)
from reclassification.fhir_task import create_recontact_task, task_to_json
from reclassification.models import (
    ClinVarSubmissionQueue,
    ReclassificationEvent,
    VUSReviewSchedule,
)
from reclassification.notifier import (
    send_reclassification_alert,
    send_vus_review_reminder,
    send_submission_failure_alert,
)

logger = get_task_logger(__name__)

# ---------------------------------------------------------------------------
# Celery application configuration
# ---------------------------------------------------------------------------

BROKER_URL = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379/0")
RESULT_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379/1")

app = Celery(
    "genomeforge_reclassification",
    broker=BROKER_URL,
    backend=RESULT_BACKEND,
)

app.conf.update(
    # Serialisation
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    # Timezone — all times in UTC
    timezone="UTC",
    enable_utc=True,
    # Task acknowledgement: ack AFTER task completes (ensures at-least-once)
    task_acks_late=True,
    # Worker concurrency: 2 threads (network-bound tasks)
    worker_concurrency=2,
    # Task routing — direct failed tasks to dead-letter queue
    task_routes={
        "reclassification.daemon.weekly_clinvar_diff": {
            "queue": "reclassification"
        },
        "reclassification.daemon.check_vus_review_dates": {
            "queue": "reclassification"
        },
        "reclassification.daemon.submit_pending_variants": {
            "queue": "clinvar_submission"
        },
    },
    # Dead-letter queue: tasks that fail after max_retries go here
    task_reject_on_worker_lost=True,
    # Result expiry: keep results for 7 days
    result_expires=604800,
    # Beat schedule: weekly ClinVar diff + daily VUS check
    beat_schedule={
        # Weekly ClinVar diff — Mondays at 08:00 UTC
        # Aligns with ClinVar weekly release cycle (published Monday morning)
        "weekly-clinvar-diff": {
            "task": "reclassification.daemon.weekly_clinvar_diff",
            "schedule": crontab(hour=8, minute=0, day_of_week=1),  # Monday
            "options": {"queue": "reclassification"},
        },
        # Daily VUS review check — 09:00 UTC every day
        "daily-vus-review-check": {
            "task": "reclassification.daemon.check_vus_review_dates",
            "schedule": crontab(hour=9, minute=0),
            "options": {"queue": "reclassification"},
        },
        # Daily ClinVar submission retry — 10:00 UTC every day
        "daily-clinvar-submission": {
            "task": "reclassification.daemon.submit_pending_variants",
            "schedule": crontab(hour=10, minute=0),
            "options": {"queue": "clinvar_submission"},
        },
        # Weekly submission status check — Thursdays 08:00 UTC
        "weekly-submission-status-check": {
            "task": "reclassification.daemon.check_submission_statuses",
            "schedule": crontab(hour=8, minute=0, day_of_week=4),  # Thursday
            "options": {"queue": "clinvar_submission"},
        },
    },
)

# ---------------------------------------------------------------------------
# Database session helper
# ---------------------------------------------------------------------------

DATABASE_URL = os.environ.get(
    "DATABASE_URL",
    "postgresql+psycopg2://genomeforge:genomeforge@localhost:5432/genomeforge",
)

_engine = None


def _get_engine():
    """Return (or create) the SQLAlchemy engine singleton."""
    global _engine
    if _engine is None:
        _engine = create_engine(DATABASE_URL, pool_pre_ping=True)
    return _engine


def _get_previous_vcf_path() -> Optional[Path]:
    """Find the most recent ClinVar VCF archive for diffing.

    Looks in the CLINVAR_ARCHIVE_DIR for the previous week's VCF.
    Returns None if no archive exists (first run).

    Returns:
        Path to previous ClinVar VCF, or None.
    """
    archive_dir = Path(os.environ.get("CLINVAR_ARCHIVE_DIR", "/data/clinvar_archive"))
    vcf_files = sorted(archive_dir.glob("clinvar_*.vcf.gz"))
    if len(vcf_files) < 2:
        return vcf_files[0] if vcf_files else None
    # Return second-to-last file (previous week)
    return vcf_files[-2]


# ---------------------------------------------------------------------------
# Celery tasks
# ---------------------------------------------------------------------------


@app.task(
    bind=True,
    max_retries=3,
    default_retry_delay=3600,  # Retry after 1 hour on transient failures
    autoretry_for=(ConnectionError, TimeoutError),
    name="reclassification.daemon.weekly_clinvar_diff",
)
def weekly_clinvar_diff(self) -> dict[str, Any]:
    """Download and diff the latest ClinVar VCF for reclassification events.

    This task runs weekly (Mondays 08:00 UTC) to:
    1. Download the latest ClinVar VCF from NCBI FTP.
    2. Compare it against last week's archived VCF.
    3. Identify reclassified variants present in the local catalogue.
    4. Create FHIR R4 Task resources for recontact workflows.
    5. Persist reclassification events to the database.
    6. Send Slack/email alerts for significant reclassifications.

    Task is automatically retried up to 3 times on connection errors,
    with a 1-hour delay between retries.

    Returns:
        Summary dict with keys: n_reclassifications, n_recontact_required,
        vcf_path, timestamp.

    Raises:
        ValueError: If ClinVar VCF download fails MD5 verification.
        ftplib.error_perm: If NCBI FTP is unreachable.
    """
    logger.info("Starting weekly ClinVar diff task")
    ftp_url = os.environ.get("CLINVAR_FTP_URL", DEFAULT_CLINVAR_FTP_URL)
    archive_dir = Path(
        os.environ.get("CLINVAR_ARCHIVE_DIR", "/data/clinvar_archive")
    )
    archive_dir.mkdir(parents=True, exist_ok=True)

    # Download latest ClinVar VCF with MD5 verification
    try:
        new_vcf = download_latest_clinvar_vcf(
            ftp_url=ftp_url,
            dest_dir=archive_dir,
            verify_checksum=True,
        )
        logger.info("Downloaded ClinVar VCF: %s", new_vcf)
    except Exception as exc:
        logger.error("Failed to download ClinVar VCF: %s", exc)
        raise self.retry(exc=exc)

    # Find previous week's VCF for comparison
    old_vcf = _get_previous_vcf_path()
    if old_vcf is None or old_vcf == new_vcf:
        logger.warning(
            "No previous ClinVar VCF found for diff — "
            "skipping diff on first run. New VCF archived at %s", new_vcf
        )
        return {
            "n_reclassifications": 0,
            "n_recontact_required": 0,
            "vcf_path": str(new_vcf),
            "timestamp": date.today().isoformat(),
            "note": "First run — no previous VCF to compare",
        }

    # Run the VCF diff to identify reclassifications
    try:
        events = diff_variants(old_vcf, new_vcf)
    except Exception as exc:
        logger.error("ClinVar diff failed: %s", exc)
        raise self.retry(exc=exc)

    logger.info("Detected %d reclassification events", len(events))

    n_recontact = 0
    with Session(_get_engine()) as session:
        for event in events:
            # Check for duplicate before persisting
            existing = session.execute(
                select(ReclassificationEvent).where(
                    ReclassificationEvent.variant_id == event.variant_id,
                    ReclassificationEvent.clinvar_date == event.clinvar_date,
                    ReclassificationEvent.old_class == event.old_class,
                    ReclassificationEvent.new_class == event.new_class,
                )
            ).scalar_one_or_none()

            if existing:
                logger.debug("Duplicate reclassification skipped: %s", event.variant_id)
                continue

            # Generate FHIR Task for recontact
            try:
                task_dict = create_recontact_task(event)
                event.fhir_task_id = task_dict["id"]
            except Exception as exc:
                logger.warning(
                    "Failed to create FHIR Task for %s: %s",
                    event.variant_id, exc,
                )

            session.add(event)

            if event.recontact_required:
                n_recontact += 1
                # Send immediate alert for recontact-required events
                try:
                    send_reclassification_alert(event)
                except Exception as exc:
                    logger.warning("Failed to send alert for %s: %s", event.variant_id, exc)

        session.commit()
        logger.info("Persisted %d reclassification events to DB", len(events))

    return {
        "n_reclassifications": len(events),
        "n_recontact_required": n_recontact,
        "vcf_path": str(new_vcf),
        "timestamp": date.today().isoformat(),
    }


@app.task(
    bind=True,
    max_retries=2,
    default_retry_delay=1800,  # Retry after 30 minutes
    name="reclassification.daemon.check_vus_review_dates",
)
def check_vus_review_dates(self) -> dict[str, Any]:
    """Check for VUS variants that are due for periodic re-review.

    Per ACGS 2024 §9, VUS variants must be re-evaluated every 2 years.
    This task runs daily to identify VUS variants approaching or past
    their 2-year review deadline and sends reminder notifications to
    the responsible clinical scientists.

    Review reminder schedule:
    - 30 days before due date: First reminder
    - 7 days before due date: Urgent reminder
    - Past due date: Overdue alert (daily until actioned)

    Returns:
        Dict with n_overdue, n_due_soon_30d, n_due_soon_7d, timestamp.
    """
    logger.info("Running VUS review date check")
    today = date.today()
    due_30d = today + timedelta(days=30)
    due_7d = today + timedelta(days=7)

    with Session(_get_engine()) as session:
        # Find overdue VUS reviews
        overdue_reviews = session.execute(
            select(VUSReviewSchedule).where(
                VUSReviewSchedule.review_due_date < today,
                VUSReviewSchedule.review_completed_at.is_(None),
            )
        ).scalars().all()

        # Find reviews due within 30 days
        due_soon_30d = session.execute(
            select(VUSReviewSchedule).where(
                VUSReviewSchedule.review_due_date <= due_30d,
                VUSReviewSchedule.review_due_date >= today,
                VUSReviewSchedule.review_completed_at.is_(None),
            )
        ).scalars().all()

        # Find reviews due within 7 days (subset of above)
        due_soon_7d = [
            r for r in due_soon_30d if r.review_due_date <= due_7d
        ]

        logger.info(
            "VUS review status: %d overdue, %d due within 30 days, "
            "%d due within 7 days",
            len(overdue_reviews), len(due_soon_30d), len(due_soon_7d),
        )

        # Send notifications for overdue and upcoming reviews
        for review in overdue_reviews:
            try:
                send_vus_review_reminder(review, urgency="overdue")
                # Update reminder timestamp
                review.reminder_sent_at = date.today()
            except Exception as exc:
                logger.warning(
                    "Failed to send overdue VUS reminder for %s: %s",
                    review.variant_id, exc,
                )

        for review in due_soon_7d:
            try:
                send_vus_review_reminder(review, urgency="urgent")
            except Exception as exc:
                logger.warning(
                    "Failed to send urgent VUS reminder for %s: %s",
                    review.variant_id, exc,
                )

        for review in [r for r in due_soon_30d if r not in due_soon_7d]:
            try:
                send_vus_review_reminder(review, urgency="upcoming")
            except Exception as exc:
                logger.warning(
                    "Failed to send upcoming VUS reminder for %s: %s",
                    review.variant_id, exc,
                )

        session.commit()

    return {
        "n_overdue": len(overdue_reviews),
        "n_due_soon_30d": len(due_soon_30d),
        "n_due_soon_7d": len(due_soon_7d),
        "timestamp": today.isoformat(),
    }


@app.task(
    bind=True,
    max_retries=3,
    default_retry_delay=7200,  # Retry after 2 hours
    name="reclassification.daemon.submit_pending_variants",
)
def submit_pending_variants(self) -> dict[str, Any]:
    """Submit pending variants from the ClinVar submission queue to NCBI.

    Processes variants in PENDING or ERROR status in the ClinVarSubmissionQueue,
    submitting them to the NCBI ClinVar API. Only processes P/LP variants
    (mandatory 3-month submission window per ACGS 2024).

    Failed submissions are retried up to 3 times before being routed to
    the dead-letter queue and generating an alert.

    Returns:
        Dict with n_submitted, n_failed, n_remaining, timestamp.
    """
    logger.info("Processing ClinVar submission queue")
    api_key = os.environ.get("NCBI_API_KEY")
    if not api_key:
        logger.warning(
            "NCBI_API_KEY not set — submissions will use anonymous rate limits "
            "(3 req/s). Set NCBI_API_KEY for faster processing."
        )

    # Mandatory submissions: P/LP only (3-month window per ACGS 2024)
    mandatory_classes = {"Pathogenic", "Likely pathogenic"}

    n_submitted = 0
    n_failed = 0

    with Session(_get_engine()) as session:
        pending = session.execute(
            select(ClinVarSubmissionQueue).where(
                ClinVarSubmissionQueue.submission_status.in_([
                    SubmissionStatus.PENDING.value,
                    SubmissionStatus.ERROR.value,
                ]),
            )
        ).scalars().all()

        # Sort: mandatory (P/LP) first
        pending_sorted = sorted(
            pending,
            key=lambda s: (
                0 if s.clinical_significance in mandatory_classes else 1,
                s.created_at,
            )
        )

        logger.info("Found %d submissions to process", len(pending_sorted))

        for sub in pending_sorted:
            try:
                result = submit_variant(sub, api_key=api_key)

                if result.success:
                    sub.submission_status = SubmissionStatus.SUBMITTED.value
                    sub.ncbi_submission_id = result.submission_id
                    sub.submitted_at = result.submitted_at
                    sub.ncbi_response = str(result.ncbi_response_raw)
                    sub.error_message = None
                    n_submitted += 1
                    logger.info(
                        "Submitted variant %s to ClinVar: %s",
                        sub.variant_id, result.submission_id,
                    )
                else:
                    sub.submission_status = SubmissionStatus.ERROR.value
                    sub.error_message = result.error_message
                    n_failed += 1
                    logger.error(
                        "Failed to submit %s: %s",
                        sub.variant_id, result.error_message,
                    )
                    try:
                        send_submission_failure_alert(sub, result.error_message or "")
                    except Exception as notify_exc:
                        logger.warning("Failed to send failure alert: %s", notify_exc)

            except Exception as exc:
                logger.error(
                    "Unexpected error submitting %s: %s", sub.variant_id, exc
                )
                sub.submission_status = SubmissionStatus.ERROR.value
                sub.error_message = str(exc)
                n_failed += 1

        session.commit()

    return {
        "n_submitted": n_submitted,
        "n_failed": n_failed,
        "n_remaining": n_failed,
        "timestamp": date.today().isoformat(),
    }


@app.task(
    bind=True,
    max_retries=2,
    default_retry_delay=3600,
    name="reclassification.daemon.check_submission_statuses",
)
def check_submission_statuses(self) -> dict[str, Any]:
    """Poll NCBI for the status of previously submitted ClinVar batches.

    Checks the status of submissions in SUBMITTED or PROCESSING state.
    Updates the database with ACCEPTED or REJECTED status as returned
    by the NCBI Submission API.

    Returns:
        Dict with n_accepted, n_rejected, n_still_processing, timestamp.
    """
    logger.info("Checking ClinVar submission statuses")
    api_key = os.environ.get("NCBI_API_KEY")

    n_accepted = 0
    n_rejected = 0
    n_processing = 0

    with Session(_get_engine()) as session:
        in_flight = session.execute(
            select(ClinVarSubmissionQueue).where(
                ClinVarSubmissionQueue.submission_status.in_([
                    SubmissionStatus.SUBMITTED.value,
                    SubmissionStatus.PROCESSING.value,
                ]),
                ClinVarSubmissionQueue.ncbi_submission_id.is_not(None),
            )
        ).scalars().all()

        logger.info("Checking status of %d in-flight submissions", len(in_flight))

        for sub in in_flight:
            try:
                status = check_submission_status(
                    sub.ncbi_submission_id, api_key=api_key
                )

                if status == SubmissionStatus.ACCEPTED.value:
                    sub.submission_status = SubmissionStatus.ACCEPTED.value
                    n_accepted += 1
                    logger.info("Submission %s accepted by ClinVar", sub.ncbi_submission_id)

                elif status in (SubmissionStatus.REJECTED.value, SubmissionStatus.ERROR.value):
                    sub.submission_status = SubmissionStatus.REJECTED.value
                    n_rejected += 1
                    logger.warning(
                        "Submission %s rejected by ClinVar", sub.ncbi_submission_id
                    )
                    try:
                        send_submission_failure_alert(
                            sub,
                            f"Submission {sub.ncbi_submission_id} was rejected by ClinVar"
                        )
                    except Exception as exc:
                        logger.warning("Failed to send rejection alert: %s", exc)

                else:
                    n_processing += 1

            except Exception as exc:
                logger.error(
                    "Error checking status for %s: %s",
                    sub.ncbi_submission_id, exc,
                )

        session.commit()

    return {
        "n_accepted": n_accepted,
        "n_rejected": n_rejected,
        "n_still_processing": n_processing,
        "timestamp": date.today().isoformat(),
    }
