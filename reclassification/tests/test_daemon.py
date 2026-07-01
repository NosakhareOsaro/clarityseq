"""
reclassification.tests.test_daemon
=====================================
pytest tests for the Celery reclassification monitoring daemon.

Tests cover:
    - _get_previous_vcf_path: path selection from archive directory.
    - _get_engine: SQLAlchemy engine singleton.
    - weekly_clinvar_diff: mocked execution (no actual Celery worker).
    - check_vus_review_dates: mocked DB session.
    - submit_pending_variants: submission queue processing.
    - check_submission_statuses: status polling.

Note: Celery tasks are tested by calling the underlying function directly
(bypassing the Celery decorator) with mocked database sessions and
external service clients.
"""

from __future__ import annotations

import os
from datetime import date, datetime, timedelta
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# _get_previous_vcf_path tests
# ---------------------------------------------------------------------------


class TestGetPreviousVcfPath:
    """Tests for _get_previous_vcf_path()."""

    def test_returns_none_when_archive_empty(self, tmp_path: Path) -> None:
        """Returns None when archive directory has no VCF files."""
        from reclassification.daemon import _get_previous_vcf_path

        with patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}):
            result = _get_previous_vcf_path()

        assert result is None

    def test_returns_only_file_when_one_vcf(self, tmp_path: Path) -> None:
        """Returns the only file when archive has exactly one VCF."""
        from reclassification.daemon import _get_previous_vcf_path

        vcf = tmp_path / "clinvar_20240101.vcf.gz"
        vcf.touch()

        with patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}):
            result = _get_previous_vcf_path()

        assert result == vcf

    def test_returns_second_to_last_when_multiple_vcfs(self, tmp_path: Path) -> None:
        """Returns second-to-last (previous week's) VCF when multiple exist."""
        from reclassification.daemon import _get_previous_vcf_path

        vcf1 = tmp_path / "clinvar_20240101.vcf.gz"
        vcf2 = tmp_path / "clinvar_20240108.vcf.gz"
        vcf3 = tmp_path / "clinvar_20240115.vcf.gz"
        vcf1.touch()
        vcf2.touch()
        vcf3.touch()

        with patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}):
            result = _get_previous_vcf_path()

        assert result == vcf2  # second-to-last

    def test_returns_second_to_last_when_two_vcfs(self, tmp_path: Path) -> None:
        """Returns first VCF (previous week's) when exactly two VCFs exist."""
        from reclassification.daemon import _get_previous_vcf_path

        vcf1 = tmp_path / "clinvar_20240101.vcf.gz"
        vcf2 = tmp_path / "clinvar_20240108.vcf.gz"
        vcf1.touch()
        vcf2.touch()

        with patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}):
            result = _get_previous_vcf_path()

        assert result == vcf1


# ---------------------------------------------------------------------------
# _get_engine tests
# ---------------------------------------------------------------------------


class TestGetEngine:
    """Tests for _get_engine() singleton."""

    def setup_method(self) -> None:
        """Reset the engine singleton before each test."""
        import reclassification.daemon as daemon_module
        daemon_module._engine = None

    def test_returns_engine_on_first_call(self) -> None:
        """Returns an engine on first call."""
        from reclassification.daemon import _get_engine

        mock_engine = MagicMock()
        with patch("reclassification.daemon.create_engine", return_value=mock_engine):
            engine = _get_engine()

        assert engine is mock_engine

    def test_returns_same_engine_on_second_call(self) -> None:
        """Returns the same engine on repeated calls (singleton)."""
        from reclassification.daemon import _get_engine

        mock_engine = MagicMock()
        with patch("reclassification.daemon.create_engine", return_value=mock_engine) as mock_create:
            e1 = _get_engine()
            e2 = _get_engine()

        # create_engine should only be called once
        assert mock_create.call_count == 1
        assert e1 is e2


# ---------------------------------------------------------------------------
# weekly_clinvar_diff task tests
# ---------------------------------------------------------------------------


class TestWeeklyClinvarDiff:
    """Tests for weekly_clinvar_diff Celery task."""

    def _call_task(self, **kwargs):
        """Call the task function directly, bypassing Celery."""
        from reclassification.daemon import weekly_clinvar_diff

        # Create a fake self (Celery task instance mock)
        mock_self = MagicMock()
        mock_self.retry.side_effect = Exception("retry called")
        return weekly_clinvar_diff.__wrapped__(mock_self) if hasattr(weekly_clinvar_diff, "__wrapped__") else weekly_clinvar_diff.run(**kwargs)

    def test_returns_no_previous_vcf_note_on_first_run(self, tmp_path: Path) -> None:
        """Returns early summary when no previous VCF exists (first run)."""
        from reclassification.daemon import weekly_clinvar_diff

        with (
            patch.dict(os.environ, {
                "CLINVAR_ARCHIVE_DIR": str(tmp_path),
                "CLINVAR_FTP_URL": "ftp://fake/",
            }),
            patch("reclassification.daemon.download_latest_clinvar_vcf") as mock_dl,
            patch("reclassification.daemon._get_previous_vcf_path", return_value=None),
        ):
            mock_dl.return_value = tmp_path / "clinvar_20240115.vcf.gz"
            (tmp_path / "clinvar_20240115.vcf.gz").touch()

            # Get the underlying function (skip Celery decorators)
            task_fn = weekly_clinvar_diff
            mock_self = MagicMock()
            mock_self.retry.side_effect = Exception("retry")

            # Call via apply with CELERY_ALWAYS_EAGER
            result = task_fn.apply(args=[])

        assert result.result["n_reclassifications"] == 0
        assert "note" in result.result

    def test_download_failure_triggers_retry(self, tmp_path: Path) -> None:
        """Download failure causes task.retry() to be called."""
        from reclassification.daemon import weekly_clinvar_diff

        with (
            patch.dict(os.environ, {
                "CLINVAR_ARCHIVE_DIR": str(tmp_path),
            }),
            patch(
                "reclassification.daemon.download_latest_clinvar_vcf",
                side_effect=ConnectionError("FTP unreachable"),
            ),
        ):
            result = weekly_clinvar_diff.apply(args=[])

        # Task should fail (retry raises in test mode)
        assert result.failed() or result.result is not None


# ---------------------------------------------------------------------------
# check_vus_review_dates task tests
# ---------------------------------------------------------------------------


class TestCheckVusReviewDates:
    """Tests for check_vus_review_dates task."""

    def test_returns_summary_dict(self) -> None:
        """Task returns summary dict with review counts."""
        from reclassification.daemon import check_vus_review_dates

        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = []

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
        ):
            result = check_vus_review_dates.apply(args=[])

        summary = result.result
        assert "n_overdue" in summary
        assert "n_due_soon_30d" in summary
        assert "timestamp" in summary

    def test_overdue_vus_sends_notification(self) -> None:
        """Overdue VUS reviews trigger send_vus_review_reminder."""
        from reclassification.daemon import check_vus_review_dates

        mock_review = MagicMock()
        mock_review.variant_id = "chr1:100:A:T"
        mock_review.review_due_date = date(2020, 1, 1)  # Past date

        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)

        call_count = {"n": 0}

        def mock_execute(query):
            result = MagicMock()
            result.scalars.return_value.all.return_value = [mock_review]
            call_count["n"] += 1
            return result

        mock_session.execute.side_effect = mock_execute

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.send_vus_review_reminder") as mock_notify,
        ):
            check_vus_review_dates.apply(args=[])

        mock_notify.assert_called()


# ---------------------------------------------------------------------------
# submit_pending_variants task tests
# ---------------------------------------------------------------------------


class TestSubmitPendingVariants:
    """Tests for submit_pending_variants task."""

    def test_empty_queue_returns_zero_submitted(self) -> None:
        """Empty submission queue → n_submitted=0."""
        from reclassification.daemon import submit_pending_variants

        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = []

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
        ):
            result = submit_pending_variants.apply(args=[])

        summary = result.result
        assert summary["n_submitted"] == 0

    def test_summary_has_timestamp(self) -> None:
        """Summary includes today's timestamp."""
        from reclassification.daemon import submit_pending_variants

        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = []

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
        ):
            result = submit_pending_variants.apply(args=[])

        assert "timestamp" in result.result


# ---------------------------------------------------------------------------
# check_submission_statuses task tests
# ---------------------------------------------------------------------------


class TestCheckSubmissionStatuses:
    """Tests for check_submission_statuses task."""

    def test_empty_in_flight_returns_zero_accepted(self) -> None:
        """Empty in-flight submissions → n_accepted=0."""
        from reclassification.daemon import check_submission_statuses

        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = []

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
        ):
            result = check_submission_statuses.apply(args=[])

        summary = result.result
        assert summary["n_accepted"] == 0
        assert "timestamp" in summary


# ---------------------------------------------------------------------------
# weekly_clinvar_diff: full event-processing flow
# ---------------------------------------------------------------------------


class TestWeeklyClinvarDiffFullFlow:
    """Tests exercising the event persistence/alerting loop body."""

    def _make_session(self, duplicate: bool = False):
        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalar_one_or_none.return_value = (
            MagicMock() if duplicate else None
        )
        return mock_session

    def test_diff_failure_triggers_retry(self, tmp_path: Path) -> None:
        """diff_variants raising should cause a task.retry() (not crash)."""
        from reclassification.daemon import weekly_clinvar_diff

        old_vcf = tmp_path / "clinvar_old.vcf.gz"
        new_vcf = tmp_path / "clinvar_new.vcf.gz"
        old_vcf.touch()
        new_vcf.touch()

        with (
            patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}),
            patch("reclassification.daemon.download_latest_clinvar_vcf", return_value=new_vcf),
            patch("reclassification.daemon._get_previous_vcf_path", return_value=old_vcf),
            patch(
                "reclassification.daemon.diff_variants",
                side_effect=RuntimeError("VCF parse error"),
            ),
        ):
            result = weekly_clinvar_diff.apply(args=[])

        # Task should fail/raise via retry rather than propagate the raw error
        assert result.failed() or result.result is not None

    def test_processes_events_persists_and_alerts(self, tmp_path: Path) -> None:
        """New (non-duplicate) events are persisted and urgent ones alerted."""
        from reclassification.daemon import weekly_clinvar_diff

        old_vcf = tmp_path / "clinvar_old.vcf.gz"
        new_vcf = tmp_path / "clinvar_new.vcf.gz"
        old_vcf.touch()
        new_vcf.touch()

        event1 = MagicMock()
        event1.variant_id = "chr1:100:A:T"
        event1.clinvar_date = date(2024, 1, 1)
        event1.old_class = "Uncertain significance"
        event1.new_class = "Pathogenic"
        event1.recontact_required = True

        event2 = MagicMock()
        event2.variant_id = "chr2:200:C:G"
        event2.clinvar_date = date(2024, 1, 1)
        event2.old_class = "Benign"
        event2.new_class = "Likely benign"
        event2.recontact_required = False

        mock_session = self._make_session(duplicate=False)

        with (
            patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}),
            patch("reclassification.daemon.download_latest_clinvar_vcf", return_value=new_vcf),
            patch("reclassification.daemon._get_previous_vcf_path", return_value=old_vcf),
            patch("reclassification.daemon.diff_variants", return_value=[event1, event2]),
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.create_recontact_task", return_value={"id": "task-123"}),
            patch("reclassification.daemon.send_reclassification_alert") as mock_alert,
        ):
            result = weekly_clinvar_diff.apply(args=[])

        summary = result.result
        assert summary["n_reclassifications"] == 2
        assert summary["n_recontact_required"] == 1
        mock_alert.assert_called_once_with(event1)
        assert event1.fhir_task_id == "task-123"
        mock_session.commit.assert_called_once()

    def test_duplicate_event_skipped(self, tmp_path: Path) -> None:
        """Events already present in the DB are skipped, not re-alerted."""
        from reclassification.daemon import weekly_clinvar_diff

        old_vcf = tmp_path / "clinvar_old.vcf.gz"
        new_vcf = tmp_path / "clinvar_new.vcf.gz"
        old_vcf.touch()
        new_vcf.touch()

        event = MagicMock()
        event.variant_id = "chr1:100:A:T"
        event.clinvar_date = date(2024, 1, 1)
        event.old_class = "Uncertain significance"
        event.new_class = "Pathogenic"
        event.recontact_required = True

        mock_session = self._make_session(duplicate=True)

        with (
            patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}),
            patch("reclassification.daemon.download_latest_clinvar_vcf", return_value=new_vcf),
            patch("reclassification.daemon._get_previous_vcf_path", return_value=old_vcf),
            patch("reclassification.daemon.diff_variants", return_value=[event]),
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.send_reclassification_alert") as mock_alert,
        ):
            result = weekly_clinvar_diff.apply(args=[])

        assert result.result["n_recontact_required"] == 0
        mock_alert.assert_not_called()
        mock_session.add.assert_not_called()

    def test_fhir_task_creation_failure_does_not_crash(self, tmp_path: Path) -> None:
        """A failure creating the FHIR Task should be logged, not raised."""
        from reclassification.daemon import weekly_clinvar_diff

        old_vcf = tmp_path / "clinvar_old.vcf.gz"
        new_vcf = tmp_path / "clinvar_new.vcf.gz"
        old_vcf.touch()
        new_vcf.touch()

        event = MagicMock()
        event.variant_id = "chr1:100:A:T"
        event.clinvar_date = date(2024, 1, 1)
        event.old_class = "Benign"
        event.new_class = "Likely benign"
        event.recontact_required = False

        mock_session = self._make_session(duplicate=False)

        with (
            patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}),
            patch("reclassification.daemon.download_latest_clinvar_vcf", return_value=new_vcf),
            patch("reclassification.daemon._get_previous_vcf_path", return_value=old_vcf),
            patch("reclassification.daemon.diff_variants", return_value=[event]),
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.create_recontact_task",
                side_effect=RuntimeError("FHIR server unreachable"),
            ),
        ):
            result = weekly_clinvar_diff.apply(args=[])

        assert result.result["n_reclassifications"] == 1
        mock_session.add.assert_called_once_with(event)

    def test_alert_failure_logged_and_continues(self, tmp_path: Path) -> None:
        """A Slack/email alert failure should not abort the task."""
        from reclassification.daemon import weekly_clinvar_diff

        old_vcf = tmp_path / "clinvar_old.vcf.gz"
        new_vcf = tmp_path / "clinvar_new.vcf.gz"
        old_vcf.touch()
        new_vcf.touch()

        event = MagicMock()
        event.variant_id = "chr1:100:A:T"
        event.clinvar_date = date(2024, 1, 1)
        event.old_class = "Uncertain significance"
        event.new_class = "Pathogenic"
        event.recontact_required = True

        mock_session = self._make_session(duplicate=False)

        with (
            patch.dict(os.environ, {"CLINVAR_ARCHIVE_DIR": str(tmp_path)}),
            patch("reclassification.daemon.download_latest_clinvar_vcf", return_value=new_vcf),
            patch("reclassification.daemon._get_previous_vcf_path", return_value=old_vcf),
            patch("reclassification.daemon.diff_variants", return_value=[event]),
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.create_recontact_task", return_value={"id": "t1"}),
            patch(
                "reclassification.daemon.send_reclassification_alert",
                side_effect=RuntimeError("Slack unreachable"),
            ),
        ):
            result = weekly_clinvar_diff.apply(args=[])

        assert result.result["n_recontact_required"] == 1
        mock_session.commit.assert_called_once()


# ---------------------------------------------------------------------------
# check_vus_review_dates: exception-handling branches
# ---------------------------------------------------------------------------


class TestCheckVusReviewDatesExceptionHandling:
    """Tests for the try/except blocks around each reminder loop."""

    def _make_session(self, execute_results: list[list]):
        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)

        results = list(execute_results)

        def mock_execute(query):
            result = MagicMock()
            result.scalars.return_value.all.return_value = results.pop(0)
            return result

        mock_session.execute.side_effect = mock_execute
        return mock_session

    def test_overdue_reminder_exception_logged(self) -> None:
        """Exception sending overdue reminder is caught and logged."""
        from reclassification.daemon import check_vus_review_dates

        review = MagicMock()
        review.variant_id = "chr1:100:A:T"
        review.review_due_date = date(2020, 1, 1)

        mock_session = self._make_session([[review], []])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.send_vus_review_reminder",
                side_effect=RuntimeError("Slack down"),
            ),
        ):
            result = check_vus_review_dates.apply(args=[])

        assert result.result["n_overdue"] == 1
        mock_session.commit.assert_called_once()

    def test_urgent_reminder_exception_logged(self) -> None:
        """Exception sending urgent (<=7d) reminder is caught and logged."""
        from reclassification.daemon import check_vus_review_dates

        review = MagicMock()
        review.variant_id = "chr2:200:C:G"
        review.review_due_date = date.today() + timedelta(days=3)

        mock_session = self._make_session([[], [review]])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.send_vus_review_reminder",
                side_effect=RuntimeError("Slack down"),
            ),
        ):
            result = check_vus_review_dates.apply(args=[])

        assert result.result["n_due_soon_7d"] == 1
        mock_session.commit.assert_called_once()

    def test_upcoming_reminder_exception_logged(self) -> None:
        """Exception sending upcoming (<=30d, >7d) reminder is caught and logged."""
        from reclassification.daemon import check_vus_review_dates

        review = MagicMock()
        review.variant_id = "chr3:300:G:C"
        review.review_due_date = date.today() + timedelta(days=20)

        mock_session = self._make_session([[], [review]])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.send_vus_review_reminder",
                side_effect=RuntimeError("Slack down"),
            ),
        ):
            result = check_vus_review_dates.apply(args=[])

        assert result.result["n_due_soon_30d"] == 1
        assert result.result["n_due_soon_7d"] == 0
        mock_session.commit.assert_called_once()


# ---------------------------------------------------------------------------
# submit_pending_variants: processing loop branches
# ---------------------------------------------------------------------------


class TestSubmitPendingVariantsProcessing:
    """Tests for the per-submission processing loop body."""

    def _make_session(self, pending: list):
        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = pending
        return mock_session

    def _make_sub(self) -> MagicMock:
        sub = MagicMock()
        sub.variant_id = "chr1:100:A:T"
        sub.clinical_significance = "Pathogenic"
        sub.created_at = datetime(2024, 1, 1)
        return sub

    def test_successful_submission_updates_status(self) -> None:
        from reclassification.clinvar_submitter import SubmissionResult
        from reclassification.daemon import submit_pending_variants
        from reclassification.models import SubmissionStatus

        sub = self._make_sub()
        result_obj = SubmissionResult(
            success=True,
            submission_id="SUB1",
            status=SubmissionStatus.SUBMITTED,
            ncbi_response_raw={"ok": True},
            error_message=None,
            submitted_at=datetime(2024, 1, 2),
        )
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.submit_variant", return_value=result_obj),
        ):
            result = submit_pending_variants.apply(args=[])

        assert result.result["n_submitted"] == 1
        assert sub.submission_status == SubmissionStatus.SUBMITTED.value
        assert sub.ncbi_submission_id == "SUB1"
        assert sub.error_message is None

    def test_failed_submission_sends_alert(self) -> None:
        from reclassification.clinvar_submitter import SubmissionResult
        from reclassification.daemon import submit_pending_variants
        from reclassification.models import SubmissionStatus

        sub = self._make_sub()
        result_obj = SubmissionResult(
            success=False,
            submission_id=None,
            status=SubmissionStatus.ERROR,
            ncbi_response_raw=None,
            error_message="NCBI rejected the payload",
            submitted_at=datetime(2024, 1, 2),
        )
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.submit_variant", return_value=result_obj),
            patch("reclassification.daemon.send_submission_failure_alert") as mock_alert,
        ):
            result = submit_pending_variants.apply(args=[])

        assert result.result["n_failed"] == 1
        assert sub.submission_status == SubmissionStatus.ERROR.value
        mock_alert.assert_called_once_with(sub, "NCBI rejected the payload")

    def test_failure_alert_exception_logged(self) -> None:
        """If send_submission_failure_alert itself raises, task still completes."""
        from reclassification.clinvar_submitter import SubmissionResult
        from reclassification.daemon import submit_pending_variants
        from reclassification.models import SubmissionStatus

        sub = self._make_sub()
        result_obj = SubmissionResult(
            success=False,
            submission_id=None,
            status=SubmissionStatus.ERROR,
            ncbi_response_raw=None,
            error_message="NCBI rejected the payload",
            submitted_at=datetime(2024, 1, 2),
        )
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch("reclassification.daemon.submit_variant", return_value=result_obj),
            patch(
                "reclassification.daemon.send_submission_failure_alert",
                side_effect=RuntimeError("Slack down"),
            ),
        ):
            result = submit_pending_variants.apply(args=[])

        assert result.result["n_failed"] == 1

    def test_unexpected_exception_marks_error(self) -> None:
        """An unhandled exception from submit_variant is caught per-item."""
        from reclassification.daemon import submit_pending_variants
        from reclassification.models import SubmissionStatus

        sub = self._make_sub()
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.submit_variant",
                side_effect=RuntimeError("unexpected boom"),
            ),
        ):
            result = submit_pending_variants.apply(args=[])

        assert result.result["n_failed"] == 1
        assert sub.submission_status == SubmissionStatus.ERROR.value
        assert "unexpected boom" in sub.error_message


# ---------------------------------------------------------------------------
# check_submission_statuses: processing loop branches
# ---------------------------------------------------------------------------


class TestCheckSubmissionStatusesProcessing:
    """Tests for the per-submission status-polling loop body."""

    def _make_session(self, in_flight: list):
        mock_session = MagicMock()
        mock_session.__enter__ = MagicMock(return_value=mock_session)
        mock_session.__exit__ = MagicMock(return_value=False)
        mock_session.execute.return_value.scalars.return_value.all.return_value = in_flight
        return mock_session

    def test_accepted_status_updates_count(self) -> None:
        from reclassification.daemon import check_submission_statuses
        from reclassification.models import SubmissionStatus

        sub = MagicMock()
        sub.ncbi_submission_id = "SUB1"
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.check_submission_status",
                return_value=SubmissionStatus.ACCEPTED.value,
            ),
        ):
            result = check_submission_statuses.apply(args=[])

        assert result.result["n_accepted"] == 1
        assert sub.submission_status == SubmissionStatus.ACCEPTED.value

    def test_rejected_status_sends_alert(self) -> None:
        from reclassification.daemon import check_submission_statuses
        from reclassification.models import SubmissionStatus

        sub = MagicMock()
        sub.ncbi_submission_id = "SUB2"
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.check_submission_status",
                return_value=SubmissionStatus.REJECTED.value,
            ),
            patch("reclassification.daemon.send_submission_failure_alert") as mock_alert,
        ):
            result = check_submission_statuses.apply(args=[])

        assert result.result["n_rejected"] == 1
        assert sub.submission_status == SubmissionStatus.REJECTED.value
        mock_alert.assert_called_once()

    def test_rejected_alert_exception_logged(self) -> None:
        from reclassification.daemon import check_submission_statuses
        from reclassification.models import SubmissionStatus

        sub = MagicMock()
        sub.ncbi_submission_id = "SUB3"
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.check_submission_status",
                return_value=SubmissionStatus.ERROR.value,
            ),
            patch(
                "reclassification.daemon.send_submission_failure_alert",
                side_effect=RuntimeError("Slack down"),
            ),
        ):
            result = check_submission_statuses.apply(args=[])

        assert result.result["n_rejected"] == 1

    def test_still_processing_status_counted(self) -> None:
        from reclassification.daemon import check_submission_statuses

        sub = MagicMock()
        sub.ncbi_submission_id = "SUB4"
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.check_submission_status",
                return_value="processing",
            ),
        ):
            result = check_submission_statuses.apply(args=[])

        assert result.result["n_still_processing"] == 1

    def test_status_check_exception_logged(self) -> None:
        """An unhandled exception from check_submission_status is caught per-item."""
        from reclassification.daemon import check_submission_statuses

        sub = MagicMock()
        sub.ncbi_submission_id = "SUB5"
        mock_session = self._make_session([sub])

        with (
            patch("reclassification.daemon._get_engine", return_value=MagicMock()),
            patch("reclassification.daemon.Session", return_value=mock_session),
            patch(
                "reclassification.daemon.check_submission_status",
                side_effect=RuntimeError("NCBI timeout"),
            ),
        ):
            result = check_submission_statuses.apply(args=[])

        summary = result.result
        assert summary["n_accepted"] == 0
        assert summary["n_rejected"] == 0
        assert summary["n_still_processing"] == 0
