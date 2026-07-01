"""Slack and email notification system for reclassification events.

Sends alerts to clinical scientists and genomics lab teams when:
- A variant reclassification requires patient recontact (ACGS 2024 §9).
- A VUS review date is approaching or overdue.
- A ClinVar submission fails or is rejected by NCBI.

Configuration (environment variables):
    SLACK_WEBHOOK_URL: Slack incoming webhook URL for the genomics-alerts channel.
    SLACK_CHANNEL: Slack channel name (default: #genomics-reclassification).
    NOTIFICATION_EMAIL_FROM: Sender address for email notifications.
    NOTIFICATION_EMAIL_TO: Comma-separated list of recipient addresses.
    SMTP_HOST: SMTP server hostname (default: localhost).
    SMTP_PORT: SMTP port (default: 587).
    SMTP_USERNAME: SMTP authentication username.
    SMTP_PASSWORD: SMTP authentication password.
    SMTP_USE_TLS: 'true' to use STARTTLS (default: true).
"""

from __future__ import annotations

import json
import logging
import os
import smtplib
from datetime import date
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from typing import Literal, Optional

import requests

from reclassification.models import (
    ClinVarSubmissionQueue,
    ReclassificationEvent,
    VUSReviewSchedule,
)

logger = logging.getLogger(__name__)

# Severity colours for Slack message attachments
SLACK_COLOUR_URGENT = "#FF0000"    # Red — recontact required
SLACK_COLOUR_WARNING = "#FFA500"   # Orange — approaching deadline
SLACK_COLOUR_INFO = "#36A64F"      # Green — informational

VUS_URGENCY_TYPE = Literal["overdue", "urgent", "upcoming"]


def _get_slack_webhook() -> Optional[str]:
    """Return the configured Slack webhook URL or None."""
    return os.environ.get("SLACK_WEBHOOK_URL")


def _get_smtp_config() -> dict:
    """Return SMTP configuration from environment variables."""
    return {
        "host": os.environ.get("SMTP_HOST", "localhost"),
        "port": int(os.environ.get("SMTP_PORT", "587")),
        "username": os.environ.get("SMTP_USERNAME", ""),
        "password": os.environ.get("SMTP_PASSWORD", ""),
        "use_tls": os.environ.get("SMTP_USE_TLS", "true").lower() == "true",
        "from_addr": os.environ.get(
            "NOTIFICATION_EMAIL_FROM", "clarityseq@nhs.uk"
        ),
        "to_addrs": [
            addr.strip()
            for addr in os.environ.get(
                "NOTIFICATION_EMAIL_TO", "genomics-lab@nhs.uk"
            ).split(",")
            if addr.strip()
        ],
    }


def _send_slack_message(
    text: str,
    attachments: Optional[list[dict]] = None,
) -> bool:
    """Post a message to the configured Slack webhook.

    Args:
        text: Main message text (supports Slack markdown).
        attachments: Optional list of Slack attachment dicts with
            'color', 'title', 'text', 'fields' keys.

    Returns:
        True if the message was sent successfully, False otherwise.
    """
    webhook_url = _get_slack_webhook()
    if not webhook_url:
        logger.warning(
            "SLACK_WEBHOOK_URL not configured — Slack notification skipped"
        )
        return False

    payload: dict = {
        "channel": os.environ.get("SLACK_CHANNEL", "#genomics-reclassification"),
        "username": "ClaritySeq Reclassification Daemon",
        "icon_emoji": ":dna:",
        "text": text,
    }
    if attachments:
        payload["attachments"] = attachments

    try:
        response = requests.post(
            webhook_url,
            data=json.dumps(payload),
            headers={"Content-Type": "application/json"},
            timeout=10,
        )
        response.raise_for_status()
        logger.debug("Slack notification sent successfully")
        return True
    except requests.RequestException as exc:
        logger.error("Failed to send Slack notification: %s", exc)
        return False


def _send_email(
    subject: str,
    body_html: str,
    body_text: str,
) -> bool:
    """Send an email notification via SMTP.

    Args:
        subject: Email subject line.
        body_html: HTML body content.
        body_text: Plain text body content (for non-HTML clients).

    Returns:
        True if email was sent successfully, False otherwise.
    """
    smtp_conf = _get_smtp_config()
    if not smtp_conf["to_addrs"]:
        logger.warning("No email recipients configured — email notification skipped")
        return False

    msg = MIMEMultipart("alternative")
    msg["Subject"] = subject
    msg["From"] = smtp_conf["from_addr"]
    msg["To"] = ", ".join(smtp_conf["to_addrs"])

    # Attach plain text first, then HTML (MIME spec: last part preferred)
    msg.attach(MIMEText(body_text, "plain"))
    msg.attach(MIMEText(body_html, "html"))

    try:
        if smtp_conf["use_tls"]:
            server = smtplib.SMTP(smtp_conf["host"], smtp_conf["port"])
            server.ehlo()
            server.starttls()
        else:
            server = smtplib.SMTP_SSL(smtp_conf["host"], smtp_conf["port"])

        if smtp_conf["username"]:
            server.login(smtp_conf["username"], smtp_conf["password"])

        server.sendmail(
            smtp_conf["from_addr"],
            smtp_conf["to_addrs"],
            msg.as_string(),
        )
        server.quit()
        logger.info("Email notification sent to %s", smtp_conf["to_addrs"])
        return True
    except smtplib.SMTPException as exc:
        logger.error("Failed to send email notification: %s", exc)
        return False


def send_reclassification_alert(event: ReclassificationEvent) -> None:
    """Send Slack and email alerts for a variant reclassification event.

    Generates a human-readable notification for clinical scientists when a
    variant reclassification is detected. Urgency is determined by whether
    patient recontact is required (ACGS 2024 §9).

    Args:
        event: ReclassificationEvent object with full reclassification details.
    """
    # Determine urgency
    is_urgent = event.recontact_required
    colour = SLACK_COLOUR_URGENT if is_urgent else SLACK_COLOUR_INFO
    urgency_label = "URGENT — Recontact Required" if is_urgent else "Informational"

    title = (
        f"ClinVar Reclassification: {event.variant_id}"
    )
    message_text = (
        f"*{urgency_label}*\n"
        f"Variant `{event.variant_id}` has been reclassified in ClinVar.\n"
        f"*Old classification:* {event.old_class}\n"
        f"*New classification:* {event.new_class}\n"
        f"*ClinVar date:* {event.clinvar_date}\n"
        f"*Detected:* {event.detected_at.date() if event.detected_at else 'today'}"
    )
    if event.clinvar_accession:
        message_text += f"\n*ClinVar accession:* {event.clinvar_accession}"
    if event.fhir_task_id:
        message_text += f"\n*FHIR Task ID:* {event.fhir_task_id}"
    if is_urgent:
        message_text += (
            "\n\n:warning: *Action required:* Patient recontact is required per "
            "ACGS 2024 §9. Please review the associated FHIR Task."
        )

    slack_attachments = [
        {
            "color": colour,
            "title": title,
            "text": message_text,
            "fields": [
                {
                    "title": "Variant ID",
                    "value": event.variant_id,
                    "short": True,
                },
                {
                    "title": "Classification Change",
                    "value": f"{event.old_class} → {event.new_class}",
                    "short": True,
                },
                {
                    "title": "ClinVar Date",
                    "value": str(event.clinvar_date),
                    "short": True,
                },
                {
                    "title": "Recontact Required",
                    "value": "Yes" if is_urgent else "No",
                    "short": True,
                },
            ],
            "footer": "ClaritySeq Reclassification Daemon | ACGS 2024 §9",
            "ts": int(
                event.detected_at.timestamp() if event.detected_at else 0
            ),
        }
    ]

    # Send Slack notification
    _send_slack_message(
        text=f":rotating_light: ClinVar Reclassification Detected: {event.variant_id}",
        attachments=slack_attachments,
    )

    # Send email for urgent (recontact-required) events only
    if is_urgent:
        subject = (
            f"[URGENT] ClinVar Reclassification Requires Patient Recontact: "
            f"{event.variant_id}"
        )
        body_html = f"""
        <html><body>
        <h2 style="color: red;">Variant Reclassification — Recontact Required</h2>
        <p>A variant reclassification has been detected that requires patient recontact
        per ACGS 2024 Best Practice Guidelines §9.</p>
        <table border="1" cellpadding="5">
            <tr><th>Field</th><th>Value</th></tr>
            <tr><td>Variant ID</td><td><code>{event.variant_id}</code></td></tr>
            <tr><td>Old Classification</td><td>{event.old_class}</td></tr>
            <tr><td>New Classification</td><td><strong>{event.new_class}</strong></td></tr>
            <tr><td>ClinVar Date</td><td>{event.clinvar_date}</td></tr>
            <tr><td>ClinVar Accession</td><td>{event.clinvar_accession or 'N/A'}</td></tr>
            <tr><td>FHIR Task ID</td><td>{event.fhir_task_id or 'Pending'}</td></tr>
        </table>
        <p>Please action the FHIR recontact Task within the timeframe required
        by ACGS 2024 §9.</p>
        <hr><p><small>Generated by ClaritySeq Reclassification Daemon</small></p>
        </body></html>
        """
        body_text = (
            f"Variant Reclassification — Recontact Required\n\n"
            f"Variant ID: {event.variant_id}\n"
            f"Old Classification: {event.old_class}\n"
            f"New Classification: {event.new_class}\n"
            f"ClinVar Date: {event.clinvar_date}\n"
            f"ClinVar Accession: {event.clinvar_accession or 'N/A'}\n"
            f"FHIR Task ID: {event.fhir_task_id or 'Pending'}\n\n"
            f"Please action the FHIR recontact Task per ACGS 2024 §9."
        )
        _send_email(subject=subject, body_html=body_html, body_text=body_text)

    logger.info(
        "Reclassification alert sent for variant %s (urgent=%s)",
        event.variant_id, is_urgent,
    )


def send_vus_review_reminder(
    review: VUSReviewSchedule,
    urgency: VUS_URGENCY_TYPE = "upcoming",
) -> None:
    """Send a reminder notification for an upcoming or overdue VUS review.

    Per ACGS 2024 §9, VUS variants require re-evaluation every 2 years.
    This function notifies clinical scientists of approaching and overdue
    review dates.

    Args:
        review: VUSReviewSchedule ORM object with variant and due date info.
        urgency: One of 'overdue', 'urgent' (within 7 days), or 'upcoming'
            (within 30 days). Determines message colour and subject line.
    """
    colour_map: dict[VUS_URGENCY_TYPE, str] = {
        "overdue": SLACK_COLOUR_URGENT,
        "urgent": SLACK_COLOUR_WARNING,
        "upcoming": SLACK_COLOUR_INFO,
    }
    label_map: dict[VUS_URGENCY_TYPE, str] = {
        "overdue": "OVERDUE",
        "urgent": "DUE WITHIN 7 DAYS",
        "upcoming": "DUE WITHIN 30 DAYS",
    }

    colour = colour_map[urgency]
    label = label_map[urgency]
    days_remaining = (review.review_due_date - date.today()).days

    message_text = (
        f"*VUS Review {label}*\n"
        f"Variant `{review.variant_id}` (patient: {review.patient_gms_id}) "
        f"is due for re-evaluation.\n"
        f"*Review due date:* {review.review_due_date}\n"
        f"*Days remaining:* {days_remaining} day(s)\n"
        f"*Initial classification date:* {review.initial_classification_date}\n"
        f"Per ACGS 2024 §9: VUS must be re-evaluated every 2 years."
    )

    _send_slack_message(
        text=f":clock1: VUS Review Reminder: {review.variant_id}",
        attachments=[
            {
                "color": colour,
                "title": f"VUS Review {label}: {review.variant_id}",
                "text": message_text,
                "footer": "ClaritySeq Reclassification Daemon | ACGS 2024 §9",
            }
        ],
    )

    logger.info(
        "VUS review reminder sent for variant %s (urgency=%s, due=%s)",
        review.variant_id, urgency, review.review_due_date,
    )


def send_submission_failure_alert(
    submission: ClinVarSubmissionQueue,
    error_message: str,
) -> None:
    """Send an alert when a ClinVar submission fails or is rejected.

    Failure to submit mandatory P/LP variants within 3 months is a
    compliance issue under the NHS GMS participation agreement. This alert
    ensures clinical scientists are aware of submission failures promptly.

    Args:
        submission: ClinVarSubmissionQueue ORM object for the failed submission.
        error_message: Human-readable error description.
    """
    is_mandatory = submission.clinical_significance in {
        "Pathogenic", "Likely pathogenic"
    }
    colour = SLACK_COLOUR_URGENT if is_mandatory else SLACK_COLOUR_WARNING
    urgency = "MANDATORY P/LP SUBMISSION FAILED" if is_mandatory else "Submission Failed"

    message_text = (
        f"*{urgency}*\n"
        f"ClinVar submission for variant `{submission.variant_id}` has failed.\n"
        f"*Gene:* {submission.gene_symbol}\n"
        f"*Classification:* {submission.clinical_significance}\n"
        f"*Error:* {error_message}\n"
        f"*NCBI Submission ID:* {submission.ncbi_submission_id or 'Not assigned'}"
    )
    if is_mandatory:
        message_text += (
            "\n\n:warning: *Compliance alert:* This is a mandatory submission "
            "under the NHS GMS participation agreement (ACGS 2024 Introduction). "
            "P/LP variants must be submitted within 3 months of report issue."
        )

    _send_slack_message(
        text=f":x: ClinVar Submission Failure: {submission.variant_id}",
        attachments=[
            {
                "color": colour,
                "title": f"ClinVar Submission Failure: {submission.variant_id}",
                "text": message_text,
                "footer": "ClaritySeq ClinVar Submitter",
            }
        ],
    )

    logger.warning(
        "Submission failure alert sent for %s (mandatory=%s): %s",
        submission.variant_id, is_mandatory, error_message,
    )
