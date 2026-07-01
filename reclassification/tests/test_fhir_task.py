"""Tests for FHIR R4 Task generation for variant reclassification recontact.

Validates generated Task resources against:
- HL7 FHIR Genomics Reporting IG v3.0.0 structural requirements.
- GA4GH VRS v2.0 identifier format (ga4gh:VA.<24-char digest>).
- Task workflow lifecycle (status, intent, priority).
- Clinical recontact urgency determination.
"""

from __future__ import annotations

import json
import re
from datetime import date, datetime, timezone
from unittest.mock import MagicMock

import pytest

from reclassification.fhir_task import (
    FHIR_TASK_PROFILE,
    VRS_IDENTIFIER_SYSTEM,
    _build_task_input_reclassification,
    _build_task_input_variant,
    _vrs_identifier_from_variant_id,
    create_recontact_task,
    task_to_json,
    validate_task_structure,
)
from reclassification.models import ClinicalSignificance, ReclassificationEvent

# VRS v2.0 identifier pattern: ga4gh:VA. followed by 24 base64url chars
VRS_PATTERN = re.compile(r"^ga4gh:VA\.[A-Za-z0-9_-]{24}$")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def vus_to_pathogenic_event() -> ReclassificationEvent:
    """Reclassification event: VUS → Pathogenic (requires recontact)."""
    event = MagicMock(spec=ReclassificationEvent)
    event.variant_id = "chr17:43094692:G:A"
    event.old_class = ClinicalSignificance.VUS.value
    event.new_class = ClinicalSignificance.PATHOGENIC.value
    event.clinvar_date = date(2024, 12, 9)
    event.clinvar_accession = "RCV000112345"
    event.detected_at = datetime(2024, 12, 9, 8, 0, 0, tzinfo=timezone.utc)
    event.recontact_required = True
    return event


@pytest.fixture
def pathogenic_to_benign_event() -> ReclassificationEvent:
    """Reclassification event: Pathogenic → Benign (requires recontact)."""
    event = MagicMock(spec=ReclassificationEvent)
    event.variant_id = "chr13:32340300:AT:A"
    event.old_class = ClinicalSignificance.PATHOGENIC.value
    event.new_class = ClinicalSignificance.BENIGN.value
    event.clinvar_date = date(2024, 11, 18)
    event.clinvar_accession = "RCV000567890"
    event.detected_at = datetime(2024, 11, 18, 8, 0, 0, tzinfo=timezone.utc)
    event.recontact_required = True
    return event


@pytest.fixture
def benign_to_likely_benign_event() -> ReclassificationEvent:
    """Reclassification event: Benign → Likely benign (no recontact)."""
    event = MagicMock(spec=ReclassificationEvent)
    event.variant_id = "chr1:100000:A:G"
    event.old_class = ClinicalSignificance.BENIGN.value
    event.new_class = ClinicalSignificance.LIKELY_BENIGN.value
    event.clinvar_date = date(2024, 10, 7)
    event.clinvar_accession = None
    event.detected_at = datetime(2024, 10, 7, 8, 0, 0, tzinfo=timezone.utc)
    event.recontact_required = False
    return event


# ---------------------------------------------------------------------------
# Tests: _vrs_identifier_from_variant_id
# ---------------------------------------------------------------------------


class TestVrsIdentifier:
    """Tests for VRS v2.0 identifier generation."""

    def test_format_matches_vrs_v2_pattern(self):
        vrs_id = _vrs_identifier_from_variant_id("chr17:43094692:G:A")
        assert VRS_PATTERN.match(vrs_id), (
            f"VRS identifier {vrs_id!r} does not match "
            f"expected pattern ga4gh:VA.<24-char-b64url>"
        )

    def test_identifier_is_deterministic(self):
        """Same variant_id always produces the same VRS identifier."""
        vrs1 = _vrs_identifier_from_variant_id("chr17:43094692:G:A")
        vrs2 = _vrs_identifier_from_variant_id("chr17:43094692:G:A")
        assert vrs1 == vrs2

    def test_different_variants_different_identifiers(self):
        """Different variants must produce distinct VRS identifiers."""
        vrs1 = _vrs_identifier_from_variant_id("chr17:43094692:G:A")
        vrs2 = _vrs_identifier_from_variant_id("chr13:32340300:AT:A")
        assert vrs1 != vrs2

    def test_identifier_length(self):
        """VRS identifier after 'ga4gh:VA.' prefix should be 24 chars."""
        vrs_id = _vrs_identifier_from_variant_id("chrX:99999:C:T")
        suffix = vrs_id.replace("ga4gh:VA.", "")
        assert len(suffix) == 24, (
            f"VRS digest should be 24 chars, got {len(suffix)}"
        )

    def test_identifier_uses_base64url_charset(self):
        """VRS digest characters should be in URL-safe base64 charset."""
        vrs_id = _vrs_identifier_from_variant_id("chr22:30000000:G:GATTAC")
        suffix = vrs_id.replace("ga4gh:VA.", "")
        # Base64url: A-Z, a-z, 0-9, _, -
        assert re.match(r"^[A-Za-z0-9_-]+$", suffix)


# ---------------------------------------------------------------------------
# Tests: _build_task_input_variant
# ---------------------------------------------------------------------------


class TestBuildTaskInputVariant:
    """Tests for FHIR Task.input variant component construction."""

    def test_input_has_type_coding(self):
        result = _build_task_input_variant(
            "chr17:43094692:G:A", "ga4gh:VA.abc123def456ghi789jkl012"
        )
        codings = result["type"]["coding"]
        assert len(codings) >= 1
        assert any(c["code"] == "69548-6" for c in codings)

    def test_input_contains_vrs_identifier(self):
        vrs_id = "ga4gh:VA.abc123def456ghi789jkl012"
        result = _build_task_input_variant("chr17:43094692:G:A", vrs_id)
        identifier = result["valueReference"]["identifier"]
        assert identifier["system"] == VRS_IDENTIFIER_SYSTEM
        assert identifier["value"] == vrs_id

    def test_vrs_extension_present(self):
        """FHIR extension for VRS identifier should be present."""
        vrs_id = "ga4gh:VA.abc123def456ghi789jkl012"
        result = _build_task_input_variant("test", vrs_id)
        extensions = result.get("extension", [])
        vrs_ext = next(
            (e for e in extensions if "vrs-identifier" in e.get("url", "")),
            None,
        )
        assert vrs_ext is not None
        assert vrs_ext["valueString"] == vrs_id


# ---------------------------------------------------------------------------
# Tests: create_recontact_task
# ---------------------------------------------------------------------------


class TestCreateRecontactTask:
    """Tests for FHIR R4 Task generation per Genomics Reporting IG v3.0.0."""

    def test_resource_type_is_task(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        assert task["resourceType"] == "Task"

    def test_task_has_uuid_id(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        # UUID format: 8-4-4-4-12 hex chars
        uuid_pattern = re.compile(
            r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
        )
        assert uuid_pattern.match(task["id"]), (
            f"Task.id {task['id']!r} is not a valid UUID"
        )

    def test_urgent_priority_for_upgrade(self, vus_to_pathogenic_event):
        """Upgrade to P/LP should produce urgent priority Task."""
        task = create_recontact_task(vus_to_pathogenic_event)
        assert task["priority"] == "urgent"

    def test_routine_priority_for_downgrade(self, pathogenic_to_benign_event):
        """Downgrade from P to Benign should produce routine priority."""
        task = create_recontact_task(pathogenic_to_benign_event)
        assert task["priority"] == "routine"

    def test_routine_priority_for_non_actionable(self, benign_to_likely_benign_event):
        task = create_recontact_task(benign_to_likely_benign_event)
        assert task["priority"] == "routine"

    def test_task_status_is_requested(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        assert task["status"] == "requested"

    def test_task_intent_is_proposal(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        assert task["intent"] == "proposal"

    def test_genomics_reporting_profile_in_meta(self, vus_to_pathogenic_event):
        """Task.meta.profile must reference the Genomics Reporting IG v3.0.0."""
        task = create_recontact_task(vus_to_pathogenic_event)
        profiles = task["meta"]["profile"]
        assert FHIR_TASK_PROFILE in profiles, (
            f"Expected Genomics Reporting IG profile {FHIR_TASK_PROFILE!r} "
            f"in Task.meta.profile"
        )

    def test_task_has_recontact_code(self, vus_to_pathogenic_event):
        """Task.code must include the recontact action code."""
        task = create_recontact_task(vus_to_pathogenic_event)
        codings = task["code"]["coding"]
        recontact_codes = [c["code"] for c in codings]
        assert "recontact" in recontact_codes

    def test_task_input_contains_vrs_identifier(self, vus_to_pathogenic_event):
        """First Task.input must contain a GA4GH VRS v2.0 identifier."""
        task = create_recontact_task(vus_to_pathogenic_event)
        first_input = task["input"][0]
        vrs_value = first_input["valueReference"]["identifier"]["value"]
        assert VRS_PATTERN.match(vrs_value), (
            f"Task.input VRS identifier {vrs_value!r} does not match "
            f"VRS v2.0 pattern"
        )

    def test_task_has_at_least_two_inputs(self, vus_to_pathogenic_event):
        """Task must have variant input and reclassification reason input."""
        task = create_recontact_task(vus_to_pathogenic_event)
        assert len(task["input"]) >= 2

    def test_task_description_includes_variant(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        assert "chr17:43094692:G:A" in task["description"]

    def test_task_description_includes_clinvar_accession(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        assert "RCV000112345" in task["description"]

    def test_custom_patient_reference(self, vus_to_pathogenic_event):
        task = create_recontact_task(
            vus_to_pathogenic_event,
            patient_reference="Patient/gms-patient-12345",
        )
        assert task["for"]["reference"] == "Patient/gms-patient-12345"

    def test_custom_requester_reference(self, vus_to_pathogenic_event):
        task = create_recontact_task(
            vus_to_pathogenic_event,
            requester_reference="Organization/rvh-genomics-lab",
        )
        assert task["requester"]["reference"] == "Organization/rvh-genomics-lab"

    def test_custom_due_date(self, vus_to_pathogenic_event):
        task = create_recontact_task(
            vus_to_pathogenic_event,
            due_date="2025-03-01",
        )
        assert task["restriction"]["period"]["end"] == "2025-03-01"

    def test_task_note_contains_vrs_id(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        notes_text = " ".join(n["text"] for n in task.get("note", []))
        assert "ga4gh:VA." in notes_text

    def test_task_serialises_to_valid_json(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        json_str = task_to_json(task)
        parsed = json.loads(json_str)
        assert parsed["resourceType"] == "Task"


# ---------------------------------------------------------------------------
# Tests: validate_task_structure
# ---------------------------------------------------------------------------


class TestValidateTaskStructure:
    """Tests for basic Task structural validation."""

    def test_valid_task_has_no_errors(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        errors = validate_task_structure(task)
        assert errors == [], f"Valid Task has unexpected errors: {errors}"

    def test_missing_resource_type_flagged(self):
        task = {"id": "123", "status": "requested", "intent": "proposal",
                "code": {}, "input": [{}]}
        errors = validate_task_structure(task)
        assert any("resourceType" in e for e in errors)

    def test_wrong_resource_type_flagged(self):
        task = {"resourceType": "Observation", "id": "123",
                "status": "requested", "intent": "proposal",
                "code": {}, "input": [{}]}
        errors = validate_task_structure(task)
        assert any("resourceType" in e for e in errors)

    def test_invalid_status_flagged(self):
        task = {"resourceType": "Task", "id": "123",
                "status": "not-a-valid-status", "intent": "proposal",
                "code": {}, "input": [{}]}
        errors = validate_task_structure(task)
        assert any("status" in e for e in errors)

    def test_missing_vrs_identifier_flagged(self, vus_to_pathogenic_event):
        """Task with non-VRS first input should fail validation."""
        task = create_recontact_task(vus_to_pathogenic_event)
        # Corrupt the VRS identifier
        task["input"][0]["valueReference"]["identifier"]["value"] = "not-a-vrs-id"
        errors = validate_task_structure(task)
        assert any("VRS" in e for e in errors)

    def test_invalid_intent_flagged(self):
        task = {"resourceType": "Task", "id": "123", "status": "requested",
                "intent": "not-a-valid-intent",
                "code": {}, "input": [{}]}
        errors = validate_task_structure(task)
        assert any("intent" in e for e in errors)

    def test_empty_input_list_flagged(self, vus_to_pathogenic_event):
        task = create_recontact_task(vus_to_pathogenic_event)
        task["input"] = []
        errors = validate_task_structure(task)
        assert any("input" in e for e in errors)
