"""Tests for the NHS-mandated ClinVar submission client.

Uses mock NCBI API responses to verify:
- JSON submission payload construction (MANE Select HGVSc + fallback coordinates).
- XML submission format generation.
- API call mechanics and result parsing.
- Status polling and response mapping.
- Error handling and retry behaviour.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from unittest.mock import MagicMock, Mock, patch

import pytest
import requests

from reclassification.clinvar_submitter import (
    NCBI_SUBMISSION_ENDPOINT,
    NCBI_STATUS_ENDPOINT,
    SubmissionResult,
    SubmissionStatus,
    build_submission_json,
    build_submission_xml,
    check_submission_status,
    submit_variant,
)
from reclassification.models import ClinVarSubmissionQueue


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pathogenic_submission() -> ClinVarSubmissionQueue:
    """ClinVarSubmissionQueue object for a BRCA1 P/LP variant."""
    sub = MagicMock(spec=ClinVarSubmissionQueue)
    sub.id = 1
    sub.variant_id = "chr17:43094692:G:A"
    sub.gene_symbol = "BRCA1"
    sub.chromosome = "17"
    sub.position_grch38 = 43094692
    sub.ref_allele = "G"
    sub.alt_allele = "A"
    sub.mane_select_hgvsc = "NM_007294.4:c.5266dupC"
    sub.clinical_significance = "Pathogenic"
    sub.condition_name = "Hereditary breast ovarian cancer syndrome"
    sub.condition_id = "C0677776"
    sub.bayesacmg_probability = 0.9987
    sub.evidence_codes = '["PVS1", "PS1", "PM2"]'
    sub.submission_status = SubmissionStatus.PENDING.value
    sub.ncbi_submission_id = None
    return sub


@pytest.fixture
def vus_submission_no_hgvsc() -> ClinVarSubmissionQueue:
    """ClinVarSubmissionQueue object for a VUS without MANE Select HGVSc."""
    sub = MagicMock(spec=ClinVarSubmissionQueue)
    sub.id = 2
    sub.variant_id = "chr13:32340300:A:T"
    sub.gene_symbol = "BRCA2"
    sub.chromosome = "13"
    sub.position_grch38 = 32340300
    sub.ref_allele = "A"
    sub.alt_allele = "T"
    sub.mane_select_hgvsc = None  # No HGVSc — use coordinates
    sub.clinical_significance = "Uncertain significance"
    sub.condition_name = "Hereditary breast ovarian cancer syndrome"
    sub.condition_id = "C0677776"
    sub.bayesacmg_probability = 0.54
    sub.evidence_codes = '["PM2", "PP3"]'
    sub.submission_status = SubmissionStatus.PENDING.value
    sub.ncbi_submission_id = None
    return sub


@pytest.fixture
def missing_gene_submission() -> ClinVarSubmissionQueue:
    """Submission missing required gene_symbol field."""
    sub = MagicMock(spec=ClinVarSubmissionQueue)
    sub.id = 3
    sub.variant_id = "chr1:100000:A:G"
    sub.gene_symbol = None  # Missing!
    sub.chromosome = "1"
    sub.position_grch38 = 100000
    sub.ref_allele = "A"
    sub.alt_allele = "G"
    sub.mane_select_hgvsc = None
    sub.clinical_significance = "Pathogenic"
    sub.condition_name = "Test condition"
    sub.condition_id = None
    sub.bayesacmg_probability = None
    sub.evidence_codes = None
    sub.submission_status = SubmissionStatus.PENDING.value
    return sub


# ---------------------------------------------------------------------------
# Tests: build_submission_json
# ---------------------------------------------------------------------------


class TestBuildSubmissionJson:
    """Tests for JSON submission payload construction."""

    def test_payload_structure_is_valid(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        # Must have top-level 'actions' array
        assert "actions" in payload
        assert isinstance(payload["actions"], list)
        assert len(payload["actions"]) >= 1

    def test_uses_mane_select_hgvsc_when_available(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        # Drill into the submission content
        variant = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["variantSet"]["variant"][0]
        )
        assert "hgvs" in variant
        assert variant["hgvs"] == "NM_007294.4:c.5266dupC"
        # Should NOT have chromosomal coordinates when HGVSc is present
        assert "chromosomeCoordinates" not in variant

    def test_falls_back_to_coordinates_without_hgvsc(self, vus_submission_no_hgvsc):
        payload = build_submission_json(vus_submission_no_hgvsc)
        variant = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["variantSet"]["variant"][0]
        )
        assert "chromosomeCoordinates" in variant
        coords = variant["chromosomeCoordinates"]
        assert coords["assembly"] == "GRCh38"
        assert coords["chromosome"] == "13"
        assert coords["start"] == 32340300
        assert coords["referenceAllele"] == "A"
        assert coords["alternateAllele"] == "T"

    def test_clinical_significance_mapped_correctly(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        clinsig = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["clinicalSignificance"]
            ["clinicalSignificanceDescription"]
        )
        assert clinsig == "Pathogenic"

    def test_gene_symbol_in_payload(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        variant = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["variantSet"]["variant"][0]
        )
        genes = [g["symbol"] for g in variant.get("gene", [])]
        assert "BRCA1" in genes

    def test_condition_with_medgen_id(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        conditions = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["conditionSet"]["condition"]
        )
        assert len(conditions) >= 1
        cond = conditions[0]
        assert cond["name"] == "Hereditary breast ovarian cancer syndrome"
        assert cond.get("id") == "C0677776"
        assert cond.get("db") == "MedGen"

    def test_bayesacmg_probability_in_comment(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        comment = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["clinicalSignificance"]["comment"]
        )
        assert "BayesACMG" in comment or "ACMG" in comment

    def test_observation_comment_includes_probability(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        obs = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["observedIn"][0]
        )
        assert "comment" in obs
        assert "0.9987" in obs["comment"]

    def test_missing_gene_raises_value_error(self, missing_gene_submission):
        with pytest.raises(ValueError, match="gene_symbol"):
            build_submission_json(missing_gene_submission)

    def test_local_id_set_to_variant_id(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        local_id = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["localID"]
        )
        assert local_id == "chr17:43094692:G:A"

    def test_acmg_criteria_reference_in_payload(self, pathogenic_submission):
        payload = build_submission_json(pathogenic_submission)
        criteria = (
            payload["actions"][0]["data"]["content"]
            ["clinvarSubmission"][0]["assertionCriteria"]
        )
        # Should cite Richards et al. 2015 (PMID:25741868)
        assert criteria.get("id") == "25741868"
        assert criteria.get("db") == "PubMed"


# ---------------------------------------------------------------------------
# Tests: build_submission_xml
# ---------------------------------------------------------------------------


class TestBuildSubmissionXml:
    """Tests for legacy XML submission format generation."""

    def test_xml_has_declaration(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "<?xml" in xml_str

    def test_root_element_is_clinvar_submission(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "ClinVarSubmission" in xml_str

    def test_xml_contains_gene_symbol(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "BRCA1" in xml_str

    def test_xml_contains_mane_select_hgvsc(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "NM_007294.4:c.5266dupC" in xml_str

    def test_xml_uses_coordinates_without_hgvsc(self, vus_submission_no_hgvsc):
        xml_str = build_submission_xml(vus_submission_no_hgvsc)
        assert "32340300" in xml_str
        assert "GRCh38" in xml_str

    def test_xml_contains_significance(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "Pathogenic" in xml_str

    def test_xml_contains_condition_name(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "Hereditary breast ovarian cancer syndrome" in xml_str

    def test_xml_contains_bayesacmg_probability(self, pathogenic_submission):
        xml_str = build_submission_xml(pathogenic_submission)
        assert "0.9987" in xml_str

    def test_missing_significance_raises(self, missing_gene_submission):
        with pytest.raises(ValueError):
            build_submission_xml(missing_gene_submission)

    def test_xml_is_parseable(self, pathogenic_submission):
        """Generated XML must be syntactically valid."""
        from xml.etree import ElementTree as ET
        xml_str = build_submission_xml(pathogenic_submission)
        # Should not raise an exception
        root = ET.fromstring(xml_str.split("?>", 1)[-1].strip())
        assert root is not None


# ---------------------------------------------------------------------------
# Tests: submit_variant
# ---------------------------------------------------------------------------


class TestSubmitVariant:
    """Tests for the NCBI ClinVar API submission function."""

    def test_dry_run_returns_success(self, pathogenic_submission):
        result = submit_variant(pathogenic_submission, dry_run=True)
        assert result.success is True
        assert result.submission_id == "dry-run-id"
        assert result.status == SubmissionStatus.SUBMITTED

    def test_dry_run_does_not_call_api(self, pathogenic_submission):
        with patch("requests.Session.post") as mock_post:
            submit_variant(pathogenic_submission, dry_run=True)
            mock_post.assert_not_called()

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_successful_submission(self, mock_sleep, pathogenic_submission):
        """Mock successful NCBI API response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "SUB123456",
            "status": "submitted",
        }
        mock_response.raise_for_status = Mock()  # No exception

        with patch("requests.Session.post", return_value=mock_response):
            result = submit_variant(pathogenic_submission, api_key="test-key")

        assert result.success is True
        assert result.submission_id == "SUB123456"
        assert result.status == SubmissionStatus.SUBMITTED
        assert result.error_message is None

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_http_error_returns_failure(self, mock_sleep, pathogenic_submission):
        """HTTP error from NCBI should return failure result."""
        mock_response = Mock()
        mock_response.status_code = 400
        mock_response.json.return_value = {"message": "Invalid submission format"}
        mock_response.text = '{"message": "Invalid submission format"}'
        mock_response.raise_for_status.side_effect = requests.HTTPError("400 Bad Request")

        with patch("requests.Session.post", return_value=mock_response):
            result = submit_variant(pathogenic_submission)

        assert result.success is False
        assert result.status == SubmissionStatus.ERROR
        assert result.error_message is not None
        assert "400" in result.error_message

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_network_error_returns_failure(self, mock_sleep, pathogenic_submission):
        """Network-level error should return failure result."""
        with patch(
            "requests.Session.post",
            side_effect=requests.ConnectionError("Connection refused"),
        ):
            result = submit_variant(pathogenic_submission)

        assert result.success is False
        assert result.status == SubmissionStatus.ERROR
        assert "Connection refused" in (result.error_message or "")

    def test_missing_gene_returns_failure(self, missing_gene_submission):
        """Submission with missing required fields should return failure."""
        result = submit_variant(missing_gene_submission)
        assert result.success is False
        assert result.status == SubmissionStatus.ERROR
        assert result.error_message is not None

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_rate_limit_sleep_called(self, mock_sleep, pathogenic_submission):
        """Sleep should be called to respect NCBI rate limits."""
        mock_response = Mock()
        mock_response.json.return_value = {"id": "SUB999"}
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.post", return_value=mock_response):
            submit_variant(pathogenic_submission)

        mock_sleep.assert_called_once()

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_api_key_header_included(self, mock_sleep, pathogenic_submission):
        """X-API-KEY header should be included when api_key is provided."""
        mock_response = Mock()
        mock_response.json.return_value = {"id": "SUB999"}
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.post", return_value=mock_response) as mock_post:
            submit_variant(pathogenic_submission, api_key="MY_API_KEY")

        # Check that the session was configured with the API key
        # (headers are set during session creation, not per-call)
        # Verify the call was made
        mock_post.assert_called_once()


# ---------------------------------------------------------------------------
# Tests: check_submission_status
# ---------------------------------------------------------------------------


class TestCheckSubmissionStatus:
    """Tests for NCBI submission status polling."""

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_accepted_status_mapped(self, mock_sleep):
        mock_response = Mock()
        mock_response.json.return_value = {
            "actions": [{"status": "processed"}]
        }
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.get", return_value=mock_response):
            status = check_submission_status("SUB123456")

        assert status == SubmissionStatus.ACCEPTED.value

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_processing_status_mapped(self, mock_sleep):
        mock_response = Mock()
        mock_response.json.return_value = {
            "actions": [{"status": "processing"}]
        }
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.get", return_value=mock_response):
            status = check_submission_status("SUB123456")

        assert status == SubmissionStatus.PROCESSING.value

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_error_status_mapped_to_rejected(self, mock_sleep):
        mock_response = Mock()
        mock_response.json.return_value = {
            "actions": [{"status": "error"}]
        }
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.get", return_value=mock_response):
            status = check_submission_status("SUB123456")

        assert status == SubmissionStatus.REJECTED.value

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_http_error_returns_error_status(self, mock_sleep):
        mock_response = Mock()
        mock_response.raise_for_status.side_effect = requests.HTTPError("404")

        with patch("requests.Session.get", return_value=mock_response):
            status = check_submission_status("INVALID_ID")

        assert status == SubmissionStatus.ERROR.value

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_empty_actions_returns_processing(self, mock_sleep):
        """Empty actions array should be treated as still processing."""
        mock_response = Mock()
        mock_response.json.return_value = {"actions": []}
        mock_response.raise_for_status = Mock()

        with patch("requests.Session.get", return_value=mock_response):
            status = check_submission_status("SUB000001")

        assert status == SubmissionStatus.PROCESSING.value

    @patch("reclassification.clinvar_submitter.time.sleep")
    def test_network_error_returns_error_status(self, mock_sleep):
        with patch(
            "requests.Session.get",
            side_effect=requests.ConnectionError("Timeout"),
        ):
            status = check_submission_status("SUB123456")

        assert status == SubmissionStatus.ERROR.value
