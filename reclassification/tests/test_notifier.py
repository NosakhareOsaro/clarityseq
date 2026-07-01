"""
reclassification.tests.test_notifier
=======================================
pytest tests for the Slack and email notification system.

Tests cover:
    - _get_slack_webhook: environment variable reading.
    - _get_smtp_config: SMTP config from env vars.
    - _send_slack_message: webhook posting with mocked requests.
    - _send_email: SMTP sending with mocked smtplib.
    - send_reclassification_alert: Slack and email integration.
    - send_vus_review_reminder: VUS review notification.
    - send_clinvar_submission_failure: submission failure notification.
"""

from __future__ import annotations

import os
from datetime import date, datetime
from unittest.mock import MagicMock, patch

import pytest

from reclassification.notifier import (
    _get_slack_webhook,
    _get_smtp_config,
    _send_email,
    _send_slack_message,
)
from reclassification.models import (
    ClinicalSignificance,
    ClinVarSubmissionQueue,
    ReclassificationEvent,
    VUSReviewSchedule,
)


# ---------------------------------------------------------------------------
# _get_slack_webhook tests
# ---------------------------------------------------------------------------


class TestGetSlackWebhook:
    """Tests for _get_slack_webhook()."""

    def test_returns_none_when_not_configured(self) -> None:
        """Returns None when SLACK_WEBHOOK_URL is not set."""
        with patch.dict(os.environ, {}, clear=True):
            os.environ.pop("SLACK_WEBHOOK_URL", None)
            result = _get_slack_webhook()
        assert result is None

    def test_returns_configured_url(self) -> None:
        """Returns the configured webhook URL."""
        with patch.dict(os.environ, {"SLACK_WEBHOOK_URL": "https://hooks.slack.com/test"}):
            result = _get_slack_webhook()
        assert result == "https://hooks.slack.com/test"


# ---------------------------------------------------------------------------
# _get_smtp_config tests
# ---------------------------------------------------------------------------


class TestGetSmtpConfig:
    """Tests for _get_smtp_config()."""

    def test_default_host_is_localhost(self) -> None:
        """Default SMTP host is localhost."""
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SMTP_HOST", None)
            config = _get_smtp_config()
        assert config["host"] == "localhost"

    def test_default_port_is_587(self) -> None:
        """Default SMTP port is 587."""
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SMTP_PORT", None)
            config = _get_smtp_config()
        assert config["port"] == 587

    def test_custom_host_applied(self) -> None:
        """Custom SMTP_HOST env var is reflected in config."""
        with patch.dict(os.environ, {"SMTP_HOST": "smtp.nhs.uk"}):
            config = _get_smtp_config()
        assert config["host"] == "smtp.nhs.uk"

    def test_to_addrs_split_by_comma(self) -> None:
        """Comma-separated NOTIFICATION_EMAIL_TO is split into a list."""
        with patch.dict(
            os.environ,
            {"NOTIFICATION_EMAIL_TO": "a@nhs.uk, b@nhs.uk , c@nhs.uk"},
        ):
            config = _get_smtp_config()
        assert config["to_addrs"] == ["a@nhs.uk", "b@nhs.uk", "c@nhs.uk"]

    def test_tls_default_is_true(self) -> None:
        """Default SMTP_USE_TLS is True."""
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SMTP_USE_TLS", None)
            config = _get_smtp_config()
        assert config["use_tls"] is True

    def test_tls_false_when_configured(self) -> None:
        """SMTP_USE_TLS=false disables TLS."""
        with patch.dict(os.environ, {"SMTP_USE_TLS": "false"}):
            config = _get_smtp_config()
        assert config["use_tls"] is False


# ---------------------------------------------------------------------------
# _send_slack_message tests
# ---------------------------------------------------------------------------


class TestSendSlackMessage:
    """Tests for _send_slack_message()."""

    def test_returns_false_when_webhook_not_configured(self) -> None:
        """Returns False when SLACK_WEBHOOK_URL is not set."""
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SLACK_WEBHOOK_URL", None)
            result = _send_slack_message("Test message")
        assert result is False

    def test_returns_true_on_successful_post(self) -> None:
        """Returns True on successful Slack webhook post."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()

        with (
            patch.dict(os.environ, {"SLACK_WEBHOOK_URL": "https://hooks.slack.com/test"}),
            patch("requests.post", return_value=mock_response),
        ):
            result = _send_slack_message("Test message")

        assert result is True

    def test_returns_false_on_request_exception(self) -> None:
        """Returns False when requests.post raises an exception."""
        import requests

        with (
            patch.dict(os.environ, {"SLACK_WEBHOOK_URL": "https://hooks.slack.com/test"}),
            patch("requests.post", side_effect=requests.RequestException("timeout")),
        ):
            result = _send_slack_message("Test message")

        assert result is False

    def test_attachments_included_in_payload(self) -> None:
        """Attachments are included in the Slack payload."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()

        attachments = [{"color": "#FF0000", "title": "URGENT", "text": "Recontact required"}]

        with (
            patch.dict(os.environ, {"SLACK_WEBHOOK_URL": "https://hooks.slack.com/test"}),
            patch("requests.post", return_value=mock_response) as mock_post,
        ):
            _send_slack_message("Test", attachments=attachments)

        call_kwargs = mock_post.call_args[1]
        import json
        payload = json.loads(call_kwargs["data"])
        assert "attachments" in payload

    def test_channel_included_in_payload(self) -> None:
        """SLACK_CHANNEL env var is used in the payload."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()

        with (
            patch.dict(os.environ, {
                "SLACK_WEBHOOK_URL": "https://hooks.slack.com/test",
                "SLACK_CHANNEL": "#my-channel",
            }),
            patch("requests.post", return_value=mock_response) as mock_post,
        ):
            _send_slack_message("Test message")

        import json
        payload = json.loads(mock_post.call_args[1]["data"])
        assert payload["channel"] == "#my-channel"


# ---------------------------------------------------------------------------
# _send_email tests
# ---------------------------------------------------------------------------


class TestSendEmail:
    """Tests for _send_email()."""

    def test_returns_false_when_no_recipients(self) -> None:
        """Returns False when no email recipients are configured."""
        with patch.dict(os.environ, {"NOTIFICATION_EMAIL_TO": ""}):
            result = _send_email("Subject", "<p>HTML</p>", "Plain text")
        assert result is False

    def test_returns_true_on_successful_send(self) -> None:
        """Returns True when SMTP sendmail succeeds."""
        import smtplib

        mock_server = MagicMock()
        mock_smtp = MagicMock(return_value=mock_server)

        with (
            patch.dict(os.environ, {
                "NOTIFICATION_EMAIL_TO": "recipient@nhs.uk",
                "SMTP_USE_TLS": "true",
                "SMTP_USERNAME": "",
            }),
            patch("smtplib.SMTP", mock_smtp),
        ):
            result = _send_email("Test Subject", "<p>Body</p>", "Body")

        assert result is True
        mock_server.sendmail.assert_called_once()
        mock_server.quit.assert_called_once()

    def test_returns_false_on_smtp_exception(self) -> None:
        """Returns False when SMTP raises an exception."""
        import smtplib

        with (
            patch.dict(os.environ, {"NOTIFICATION_EMAIL_TO": "recipient@nhs.uk"}),
            patch("smtplib.SMTP", side_effect=smtplib.SMTPException("Connection refused")),
        ):
            result = _send_email("Subject", "<p>HTML</p>", "Text")

        assert result is False

    def test_login_called_when_username_configured(self) -> None:
        """SMTP login is invoked when SMTP_USERNAME is set."""
        mock_server = MagicMock()
        mock_smtp = MagicMock(return_value=mock_server)

        with (
            patch.dict(os.environ, {
                "NOTIFICATION_EMAIL_TO": "recipient@nhs.uk",
                "SMTP_USE_TLS": "true",
                "SMTP_USERNAME": "smtp-user",
                "SMTP_PASSWORD": "smtp-pass",
            }),
            patch("smtplib.SMTP", mock_smtp),
        ):
            result = _send_email("Test Subject", "<p>Body</p>", "Body")

        assert result is True
        mock_server.login.assert_called_once_with("smtp-user", "smtp-pass")

    def test_uses_ssl_when_tls_disabled(self) -> None:
        """Uses SMTP_SSL when SMTP_USE_TLS=false."""
        import smtplib

        mock_server = MagicMock()
        mock_smtp_ssl = MagicMock(return_value=mock_server)

        with (
            patch.dict(os.environ, {
                "NOTIFICATION_EMAIL_TO": "r@nhs.uk",
                "SMTP_USE_TLS": "false",
                "SMTP_USERNAME": "",
            }),
            patch("smtplib.SMTP_SSL", mock_smtp_ssl),
        ):
            result = _send_email("Subject", "<p>HTML</p>", "Text")

        mock_smtp_ssl.assert_called_once()
        assert result is True


# ---------------------------------------------------------------------------
# send_reclassification_alert tests
# ---------------------------------------------------------------------------


class TestSendReclassificationAlert:
    """Tests for send_reclassification_alert()."""

    def _make_event(self, recontact: bool = True) -> ReclassificationEvent:
        """Create a minimal ReclassificationEvent for testing."""
        return ReclassificationEvent(
            variant_id="chr17:43094692:G:A",
            old_class=ClinicalSignificance.VUS,
            new_class=ClinicalSignificance.PATHOGENIC,
            recontact_required=recontact,
            clinvar_date=date(2024, 12, 13),
            detected_at=datetime(2024, 12, 13),
            clinvar_accession="RCV000012345",
        )

    def test_sends_slack_and_email(self) -> None:
        """send_reclassification_alert calls both Slack and email helpers."""
        from reclassification.notifier import send_reclassification_alert

        event = self._make_event(recontact=True)

        with (
            patch("reclassification.notifier._send_slack_message", return_value=True) as mock_slack,
            patch("reclassification.notifier._send_email", return_value=True) as mock_email,
        ):
            send_reclassification_alert(event)

        mock_slack.assert_called_once()
        mock_email.assert_called_once()

    def test_urgent_colour_for_recontact(self) -> None:
        """Recontact required → SLACK_COLOUR_URGENT used in attachment."""
        from reclassification.notifier import (
            SLACK_COLOUR_URGENT,
            send_reclassification_alert,
        )

        event = self._make_event(recontact=True)
        captured_attachments = []

        def capture_slack(text, attachments=None):
            if attachments:
                captured_attachments.extend(attachments)
            return True

        with (
            patch("reclassification.notifier._send_slack_message", side_effect=capture_slack),
            patch("reclassification.notifier._send_email", return_value=True),
        ):
            send_reclassification_alert(event)

        assert any(a.get("color") == SLACK_COLOUR_URGENT for a in captured_attachments)

    def test_fhir_task_id_included_in_message(self) -> None:
        """When set, the event's fhir_task_id should appear in the Slack text."""
        from reclassification.notifier import send_reclassification_alert

        event = self._make_event(recontact=True)
        event.fhir_task_id = "task-abc-123"
        captured_text = {}

        def capture_slack(text, attachments=None):
            captured_text["text"] = text
            if attachments:
                captured_text["attachments"] = attachments
            return True

        with (
            patch("reclassification.notifier._send_slack_message", side_effect=capture_slack),
            patch("reclassification.notifier._send_email", return_value=True),
        ):
            send_reclassification_alert(event)

        message_text = captured_text["attachments"][0]["text"]
        assert "task-abc-123" in message_text

    def test_non_urgent_event_does_not_send_email(self) -> None:
        """Non-recontact events should only trigger Slack, not email."""
        from reclassification.notifier import send_reclassification_alert

        event = self._make_event(recontact=False)

        with (
            patch("reclassification.notifier._send_slack_message", return_value=True) as mock_slack,
            patch("reclassification.notifier._send_email", return_value=True) as mock_email,
        ):
            send_reclassification_alert(event)

        mock_slack.assert_called_once()
        mock_email.assert_not_called()


# ---------------------------------------------------------------------------
# send_vus_review_reminder tests
# ---------------------------------------------------------------------------


class TestSendVusReviewReminder:
    """Tests for send_vus_review_reminder()."""

    def _make_review(self, due_date: date = date(2026, 8, 1)) -> VUSReviewSchedule:
        review = MagicMock(spec=VUSReviewSchedule)
        review.variant_id = "chr17:43094692:G:A"
        review.patient_gms_id = "GMS-0001"
        review.review_due_date = due_date
        review.initial_classification_date = date(2024, 8, 1)
        return review

    def test_overdue_urgency_sends_slack_message(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_URGENT, send_vus_review_reminder

        review = self._make_review(due_date=date(2020, 1, 1))
        captured = {}

        def capture_slack(text, attachments=None):
            captured["text"] = text
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack) as mock_slack:
            send_vus_review_reminder(review, urgency="overdue")

        mock_slack.assert_called_once()
        assert captured["attachments"][0]["color"] == SLACK_COLOUR_URGENT
        assert "OVERDUE" in captured["attachments"][0]["title"]

    def test_urgent_urgency_sends_slack_message(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_WARNING, send_vus_review_reminder

        review = self._make_review(due_date=date.today())
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_vus_review_reminder(review, urgency="urgent")

        assert captured["attachments"][0]["color"] == SLACK_COLOUR_WARNING
        assert "7 DAYS" in captured["attachments"][0]["title"]

    def test_upcoming_urgency_is_default(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_INFO, send_vus_review_reminder

        review = self._make_review(due_date=date.today())
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_vus_review_reminder(review)  # default urgency="upcoming"

        assert captured["attachments"][0]["color"] == SLACK_COLOUR_INFO
        assert "30 DAYS" in captured["attachments"][0]["title"]

    def test_message_includes_variant_and_patient(self) -> None:
        from reclassification.notifier import send_vus_review_reminder

        review = self._make_review()
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_vus_review_reminder(review, urgency="upcoming")

        body = captured["attachments"][0]["text"]
        assert review.variant_id in body
        assert review.patient_gms_id in body


# ---------------------------------------------------------------------------
# send_submission_failure_alert tests
# ---------------------------------------------------------------------------


class TestSendSubmissionFailureAlert:
    """Tests for send_submission_failure_alert()."""

    def _make_submission(self, clinical_significance: str = "Pathogenic") -> ClinVarSubmissionQueue:
        sub = MagicMock(spec=ClinVarSubmissionQueue)
        sub.variant_id = "chr17:43094692:G:A"
        sub.gene_symbol = "BRCA1"
        sub.clinical_significance = clinical_significance
        sub.ncbi_submission_id = None
        return sub

    def test_mandatory_pathogenic_submission_uses_urgent_colour(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_URGENT, send_submission_failure_alert

        sub = self._make_submission(clinical_significance="Pathogenic")
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_submission_failure_alert(sub, "NCBI rejected the submission")

        attachment = captured["attachments"][0]
        assert attachment["color"] == SLACK_COLOUR_URGENT
        assert "MANDATORY" in attachment["text"]
        assert "NCBI rejected the submission" in attachment["text"]

    def test_likely_pathogenic_is_also_mandatory(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_URGENT, send_submission_failure_alert

        sub = self._make_submission(clinical_significance="Likely pathogenic")
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_submission_failure_alert(sub, "timeout")

        assert captured["attachments"][0]["color"] == SLACK_COLOUR_URGENT

    def test_non_mandatory_vus_uses_warning_colour(self) -> None:
        from reclassification.notifier import SLACK_COLOUR_WARNING, send_submission_failure_alert

        sub = self._make_submission(clinical_significance="Uncertain significance")
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_submission_failure_alert(sub, "some error")

        attachment = captured["attachments"][0]
        assert attachment["color"] == SLACK_COLOUR_WARNING
        assert "MANDATORY" not in attachment["text"]

    def test_message_includes_gene_and_ncbi_id(self) -> None:
        from reclassification.notifier import send_submission_failure_alert

        sub = self._make_submission(clinical_significance="Pathogenic")
        sub.ncbi_submission_id = "SUB999"
        captured = {}

        def capture_slack(text, attachments=None):
            captured["attachments"] = attachments
            return True

        with patch("reclassification.notifier._send_slack_message", side_effect=capture_slack):
            send_submission_failure_alert(sub, "rejected")

        body = captured["attachments"][0]["text"]
        assert "BRCA1" in body
        assert "SUB999" in body
