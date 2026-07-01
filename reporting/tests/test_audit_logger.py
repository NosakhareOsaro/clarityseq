"""
reporting.tests.test_audit_logger
===================================
pytest tests for the JSON-LD audit trail writer.

Tests cover:
    - write_audit_log: structured JSON-LD output with expected keys/values
      (data sources incl. gnomad_version, classification scheme, variants).
    - _get_hostname: hostname lookup and failure fallback to "unknown".
    - append_audit_event: appending post-report events to an existing
      audit file, and FileNotFoundError when the file does not exist.

References:
    JSON-LD 1.1: https://www.w3.org/TR/json-ld11/
    PROV-O: https://www.w3.org/TR/prov-o/
    ACGS 2024 v1.2 §5 — classification audit trail.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from reporting.audit_logger import (
    _get_hostname,
    append_audit_event,
    write_audit_log,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_audit_data() -> dict:
    """Return a minimal audit_data dict as produced by ReportGenerator."""
    return {
        "report_date": "2024-12-13",
        "patient_id": "TEST-PATIENT-001",
        "sample_id": "LAB-2024-12345",
        "pipeline_version": "1.0.0",
        "assembly": "GRCh38",
        "tools": {
            "vep": "111",
            "gnomad": "4.1",
            "alphamissense": "Cheng et al. 2023 PMID:37703350",
        },
        "classification_scheme": {
            "framework": "ACGS 2024 v1.2",
            "point_system": "Tavtigian et al. 2020 PMID:32645316",
            "pm2_weight": "Supporting (1 pt) — ClinGen SVI 2024",
            "pp3_bp4_primary": "AlphaMissense (Cheng 2023): ≥0.564→PP3, ≤0.340→BP4",
            "mane_select": "Morales et al. 2022 PMID:35356062",
        },
        "variants": [
            {
                "gene": "BRCA1",
                "transcript": "NM_007294.4",
                "hgvsc": "NM_007294.4:c.5266dupC",
                "hgvsp": "NP_009225.1:p.Gln1756ProfsTer25",
                "acmg_class": "Pathogenic",
                "rules_applied": ["PVS1", "PM2"],
                "gnomad_af": None,
                "alphamissense_score": None,
                "posterior_p": 0.99,
                "vus_review_date": None,
                "pending_clinvar_submission": False,
            }
        ],
    }


# ---------------------------------------------------------------------------
# write_audit_log tests
# ---------------------------------------------------------------------------


class TestWriteAuditLog:
    """Tests for write_audit_log() JSON-LD output structure."""

    def test_writes_file(self, tmp_path: Path, sample_audit_data: dict) -> None:
        """write_audit_log creates the output file."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        assert out.exists()

    def test_json_ld_context_present(self, tmp_path: Path, sample_audit_data: dict) -> None:
        """Output contains the @context and @type JSON-LD keys."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert "@context" in data
        assert data["@type"] == "prov:Activity"
        assert data["@context"]["@vocab"] == "https://schema.org/"

    def test_id_uses_sample_id_and_report_date(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """@id field embeds sample_id and report_date."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert data["@id"] == "urn:genomeforge:audit:LAB-2024-12345:2024-12-13"

    def test_data_sources_contains_gnomad_version(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """data_sources.gnomad_version reflects the tools.gnomad value."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert data["data_sources"]["gnomad_version"] == "4.1"
        assert data["data_sources"]["vep_version"] == "111"
        assert data["data_sources"]["genome_assembly"] == "GRCh38"
        assert data["data_sources"]["gnomad_individuals"] == 807162

    def test_data_sources_defaults_when_tools_missing(self, tmp_path: Path) -> None:
        """Missing 'tools' key falls back to default vep/gnomad versions."""
        out = tmp_path / "audit.jsonld"
        write_audit_log({"sample_id": "S1", "report_date": "2024-01-01"}, out)
        data = json.loads(out.read_text())
        assert data["data_sources"]["gnomad_version"] == "4.1"
        assert data["data_sources"]["vep_version"] == "111"

    def test_sample_block_contains_patient_and_sample_id(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """sample block records patient_id and sample_id."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert data["sample"]["patient_id"] == "TEST-PATIENT-001"
        assert data["sample"]["sample_id"] == "LAB-2024-12345"

    def test_variants_list_mapped_correctly(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """Each variant dict is mapped into the JSON-LD variants array."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert len(data["variants"]) == 1
        v = data["variants"][0]
        assert v["gene"] == "BRCA1"
        assert v["acmg_class"] == "Pathogenic"
        assert v["rules_applied"] == ["PVS1", "PM2"]
        assert v["@type"] == "genomeforge:VariantClassification"

    def test_classification_scheme_default_when_missing(self, tmp_path: Path) -> None:
        """Missing classification_scheme falls back to the default dict."""
        out = tmp_path / "audit.jsonld"
        write_audit_log({"sample_id": "S1", "report_date": "2024-01-01"}, out)
        data = json.loads(out.read_text())
        assert "Supporting" in data["classification_scheme"]["pm2_weight"]

    def test_system_info_contains_hostname_and_user(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """system_info block records hostname and user."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        data = json.loads(out.read_text())
        assert "hostname" in data["system_info"]
        assert "user" in data["system_info"]

    def test_pretty_false_produces_compact_json(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """pretty=False writes JSON without indentation (single line)."""
        out = tmp_path / "audit_compact.jsonld"
        write_audit_log(sample_audit_data, out, pretty=False)
        text = out.read_text()
        assert "\n" not in text.strip("\n")

    def test_creates_parent_directories(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """write_audit_log creates missing parent directories."""
        out = tmp_path / "nested" / "dir" / "audit.jsonld"
        write_audit_log(sample_audit_data, out)
        assert out.exists()


# ---------------------------------------------------------------------------
# _get_hostname tests
# ---------------------------------------------------------------------------


class TestGetHostname:
    """Tests for _get_hostname() including the failure fallback."""

    def test_returns_hostname_string(self) -> None:
        """Returns a non-empty string under normal conditions."""
        hostname = _get_hostname()
        assert isinstance(hostname, str)
        assert hostname != ""

    def test_returns_unknown_on_socket_error(self) -> None:
        """When socket.gethostname() raises, returns 'unknown'."""
        with patch("reporting.audit_logger.socket.gethostname", side_effect=OSError("boom")):
            assert _get_hostname() == "unknown"


# ---------------------------------------------------------------------------
# append_audit_event tests
# ---------------------------------------------------------------------------


class TestAppendAuditEvent:
    """Tests for append_audit_event() post-report event logging."""

    def test_raises_file_not_found_when_missing(self, tmp_path: Path) -> None:
        """FileNotFoundError is raised when audit_path does not exist."""
        missing = tmp_path / "does_not_exist.jsonld"
        with pytest.raises(FileNotFoundError, match="Audit file not found"):
            append_audit_event(missing, "clinvar_submission", {"accession": "SCV000123"})

    def test_appends_event_to_existing_file(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """Appending an event adds it to the 'events' list in the audit file."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)

        append_audit_event(
            out,
            "clinvar_submission",
            {"accession": "SCV000123", "status": "submitted"},
        )

        data = json.loads(out.read_text())
        assert "events" in data
        assert len(data["events"]) == 1
        event = data["events"][0]
        assert event["@type"] == "genomeforge:clinvar_submission"
        assert event["data"]["accession"] == "SCV000123"
        assert "timestamp" in event

    def test_multiple_events_accumulate(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """Multiple calls to append_audit_event accumulate events in order."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)

        append_audit_event(out, "vus_reclassification", {"new_class": "Likely_Pathogenic"})
        append_audit_event(out, "report_amendment", {"reason": "typo fix"})

        data = json.loads(out.read_text())
        assert len(data["events"]) == 2
        assert data["events"][0]["@type"] == "genomeforge:vus_reclassification"
        assert data["events"][1]["@type"] == "genomeforge:report_amendment"

    def test_preserves_existing_audit_fields(
        self, tmp_path: Path, sample_audit_data: dict
    ) -> None:
        """Appending an event does not clobber existing audit fields."""
        out = tmp_path / "audit.jsonld"
        write_audit_log(sample_audit_data, out)

        append_audit_event(out, "clinvar_submission", {"accession": "SCV000999"})

        data = json.loads(out.read_text())
        assert data["sample"]["sample_id"] == "LAB-2024-12345"
        assert data["@type"] == "prov:Activity"
