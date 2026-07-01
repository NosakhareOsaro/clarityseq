"""
annotation.tests.test_clinvar_client
=======================================
pytest tests for the ClinVar clinical significance lookup client.

Tests cover:
    - ClinVarClient.get_clinvar_data: SPDI-based lookup, HTTP error handling.
    - ClinVarClient.get_by_variation_id: VCV-prefixed and bare ID lookup.
    - ClinVarClient._lookup_by_spdi: real HTTP body incl. 404 handling.
    - ClinVarClient._parse_variation: response-shape variations, submissions,
      condition parsing, RCV fallback.
    - _chrom_to_refseq: known and unknown chromosome mapping.
    - _review_status_to_stars: star-rating derivation from review status text.
    - _suggest_acmg_evidence: PP5/BP6 suggestion logic.

References:
    Landrum et al. 2016 PMID:26582918 (ClinVar).
"""

from __future__ import annotations

from unittest.mock import AsyncMock, MagicMock, patch

import httpx
import pytest

from annotation.clinvar_client import (
    ClinVarClient,
    ClinVarData,
    ClinVarSubmission,
    _chrom_to_refseq,
    _review_status_to_stars,
    _suggest_acmg_evidence,
)


# ---------------------------------------------------------------------------
# Dataclass defaults
# ---------------------------------------------------------------------------


class TestClinVarDataclasses:
    """Sanity checks for the ClinVarSubmission / ClinVarData dataclasses."""

    def test_submission_requires_submitter_and_classification(self) -> None:
        sub = ClinVarSubmission(submitter="LabCorp", classification="Pathogenic")
        assert sub.submitter == "LabCorp"
        assert sub.condition is None

    def test_clinvar_data_defaults(self) -> None:
        data = ClinVarData()
        assert data.star_rating == 0
        assert data.condition_names == []
        assert data.submissions == []
        assert data.classification is None


# ---------------------------------------------------------------------------
# _chrom_to_refseq
# ---------------------------------------------------------------------------


class TestChromToRefseq:
    """Tests for _chrom_to_refseq()."""

    def test_chr17_maps_to_correct_accession(self) -> None:
        assert _chrom_to_refseq("17") == "000017.11"

    def test_chrx_maps_correctly(self) -> None:
        assert _chrom_to_refseq("X") == "000023.11"

    def test_mt_maps_correctly(self) -> None:
        assert _chrom_to_refseq("MT") == "012920.1"

    def test_lowercase_x_maps_correctly(self) -> None:
        """Lower-case chromosome letters should be upper-cased before lookup."""
        assert _chrom_to_refseq("x") == "000023.11"

    def test_unknown_chrom_falls_back_to_chr1(self) -> None:
        """Unmapped chromosome names fall back to the chr1 accession stub."""
        assert _chrom_to_refseq("UNKNOWN") == "000001.11"


# ---------------------------------------------------------------------------
# _review_status_to_stars
# ---------------------------------------------------------------------------


class TestReviewStatusToStars:
    """Tests for _review_status_to_stars()."""

    def test_practice_guideline_is_four_stars(self) -> None:
        assert _review_status_to_stars("practice guideline") == 4

    def test_expert_panel_is_three_stars(self) -> None:
        assert _review_status_to_stars("reviewed by expert panel") == 3

    def test_multiple_submitters_no_conflicts_is_two_stars(self) -> None:
        assert (
            _review_status_to_stars(
                "criteria provided, multiple submitters, no conflicts"
            )
            == 2
        )

    def test_single_submitter_is_one_star(self) -> None:
        assert (
            _review_status_to_stars("criteria provided, single submitter") == 1
        )

    def test_criteria_provided_alone_is_one_star(self) -> None:
        assert _review_status_to_stars("criteria provided") == 1

    def test_no_assertion_criteria_is_zero_stars(self) -> None:
        """The 0-star status text contains the substring "criteria provided"
        and must not be misclassified as 1 star."""
        assert _review_status_to_stars("no assertion criteria provided") == 0

    def test_no_assertion_provided_is_zero_stars(self) -> None:
        assert _review_status_to_stars("no assertion provided") == 0

    def test_empty_string_is_zero_stars(self) -> None:
        assert _review_status_to_stars("") == 0

    def test_case_insensitive(self) -> None:
        assert _review_status_to_stars("PRACTICE GUIDELINE") == 4


# ---------------------------------------------------------------------------
# _suggest_acmg_evidence
# ---------------------------------------------------------------------------


class TestSuggestAcmgEvidence:
    """Tests for _suggest_acmg_evidence()."""

    def test_pathogenic_with_stars_returns_pp5(self) -> None:
        assert _suggest_acmg_evidence("Pathogenic", 2) == "PP5"

    def test_likely_pathogenic_with_stars_returns_pp5(self) -> None:
        assert _suggest_acmg_evidence("Likely pathogenic", 1) == "PP5"

    def test_pathogenic_likely_pathogenic_combo_returns_pp5(self) -> None:
        assert (
            _suggest_acmg_evidence("Pathogenic/Likely pathogenic", 3) == "PP5"
        )

    def test_benign_with_stars_returns_bp6(self) -> None:
        assert _suggest_acmg_evidence("Benign", 2) == "BP6"

    def test_likely_benign_with_stars_returns_bp6(self) -> None:
        assert _suggest_acmg_evidence("Likely benign", 1) == "BP6"

    def test_zero_stars_returns_none_even_if_pathogenic(self) -> None:
        """No assertion criteria (0 stars) should never suggest PP5/BP6."""
        assert _suggest_acmg_evidence("Pathogenic", 0) is None

    def test_none_classification_returns_none(self) -> None:
        assert _suggest_acmg_evidence(None, 3) is None

    def test_uncertain_significance_returns_none(self) -> None:
        assert _suggest_acmg_evidence("Uncertain significance", 3) is None

    def test_conflicting_returns_none(self) -> None:
        assert _suggest_acmg_evidence("Conflicting", 2) is None


# ---------------------------------------------------------------------------
# ClinVarClient._parse_variation
# ---------------------------------------------------------------------------


MOCK_VARIATION_DIRECT = {
    "variation_id": 12345,
    # Presence of a top-level "classification" key (any value) is what makes
    # _parse_variation treat this dict directly as the "record" rather than
    # looking for a nested "result" key.
    "classification": None,
    "germline_classification": {
        "description": "Pathogenic",
        "review_status": "reviewed by expert panel",
    },
    "accession_list": ["RCV000031282"],
    "trait_set": [{"name": "Hereditary breast and ovarian cancer syndrome"}],
    "submissions": [
        {
            "submitter_name": "Invitae",
            "classification": "Pathogenic",
            "last_evaluated": "2023-01-15",
            "review_status": "criteria provided, single submitter",
            "condition_name": "HBOC",
            "collection_method": "clinical testing",
        }
    ],
    "date_last_updated": "2024-03-01",
}


class TestParseVariation:
    """Tests for ClinVarClient._parse_variation()."""

    def setup_method(self) -> None:
        self.client = ClinVarClient()

    def test_empty_data_returns_none(self) -> None:
        assert self.client._parse_variation({}) is None

    def test_direct_form_parsed(self) -> None:
        """Response with 'classification' at top level should parse directly."""
        result = self.client._parse_variation(MOCK_VARIATION_DIRECT)
        assert result is not None
        assert result.classification == "Pathogenic"
        assert result.star_rating == 3
        assert result.acmg_evidence == "PP5"
        assert result.rcv_accession == "RCV000031282"
        assert result.variation_id == "12345"

    def test_conditions_extracted(self) -> None:
        result = self.client._parse_variation(MOCK_VARIATION_DIRECT)
        assert result is not None
        assert "Hereditary breast and ovarian cancer syndrome" in result.condition_names

    def test_submissions_parsed(self) -> None:
        result = self.client._parse_variation(MOCK_VARIATION_DIRECT)
        assert result is not None
        assert len(result.submissions) == 1
        sub = result.submissions[0]
        assert sub.submitter == "Invitae"
        assert sub.condition == "HBOC"
        assert sub.method == "clinical testing"

    def test_result_wrapped_form_parsed(self) -> None:
        """Response nested under a 'result' key should also parse correctly."""
        wrapped = {"result": MOCK_VARIATION_DIRECT}
        result = self.client._parse_variation(wrapped)
        assert result is not None
        assert result.classification == "Pathogenic"

    def test_empty_result_key_returns_none(self) -> None:
        assert self.client._parse_variation({"result": {}}) is None

    def test_classification_as_dict_form(self) -> None:
        """agg_class can itself be a dict with a 'description' key."""
        data = {
            "classification": {"description": "Benign"},
            "germline_classification": {},
            "review_status": "criteria provided, multiple submitters, no conflicts",
        }
        result = self.client._parse_variation(data)
        assert result is not None
        assert result.classification == "Benign"
        assert result.acmg_evidence == "BP6"

    def test_missing_rcv_falls_back_to_accession(self) -> None:
        data = {
            "classification": "Pathogenic",
            "germline_classification": {},
            "accession": "RCV999999999",
            "review_status": "reviewed by expert panel",
        }
        result = self.client._parse_variation(data)
        assert result is not None
        assert result.rcv_accession == "RCV999999999"

    def test_unrecognised_classification_string_kept_as_is(self) -> None:
        """A classification string not in the map is passed through unmodified."""
        data = {
            "classification": "Some Novel Classification",
            "germline_classification": {},
        }
        result = self.client._parse_variation(data)
        assert result is not None
        assert result.classification == "Some Novel Classification"

    def test_default_review_status_used_when_missing(self) -> None:
        """When no review_status is present, falls back by star rating."""
        data = {"classification": "Uncertain significance", "germline_classification": {}}
        result = self.client._parse_variation(data)
        assert result is not None
        assert result.star_rating == 0
        assert result.review_status == "no assertion criteria provided"


# ---------------------------------------------------------------------------
# ClinVarClient.get_clinvar_data
# ---------------------------------------------------------------------------


class TestGetClinVarData:
    """Tests for ClinVarClient.get_clinvar_data()."""

    @pytest.mark.asyncio
    async def test_returns_data_on_success(self) -> None:
        client = ClinVarClient()
        with patch.object(
            client, "_lookup_by_spdi", new_callable=AsyncMock
        ) as mock_lookup:
            mock_lookup.return_value = MOCK_VARIATION_DIRECT
            result = await client.get_clinvar_data("chr17", 43094692, "G", "A")

        assert result is not None
        assert result.classification == "Pathogenic"

    @pytest.mark.asyncio
    async def test_returns_none_when_lookup_returns_none(self) -> None:
        client = ClinVarClient()
        with patch.object(
            client, "_lookup_by_spdi", new_callable=AsyncMock
        ) as mock_lookup:
            mock_lookup.return_value = None
            result = await client.get_clinvar_data("chr17", 43094692, "G", "A")

        assert result is None

    @pytest.mark.asyncio
    async def test_http_error_returns_none(self) -> None:
        client = ClinVarClient()
        with patch.object(
            client, "_lookup_by_spdi", new_callable=AsyncMock
        ) as mock_lookup:
            mock_lookup.side_effect = httpx.HTTPError("connection reset")
            result = await client.get_clinvar_data("chr17", 43094692, "G", "A")

        assert result is None

    @pytest.mark.asyncio
    async def test_spdi_built_from_bare_chrom(self) -> None:
        """SPDI passed to _lookup_by_spdi should use 0-based pos & refseq stub."""
        client = ClinVarClient()
        with patch.object(
            client, "_lookup_by_spdi", new_callable=AsyncMock
        ) as mock_lookup:
            mock_lookup.return_value = None
            await client.get_clinvar_data("chr17", 43094692, "G", "A")

        called_spdi = mock_lookup.call_args[0][0]
        assert called_spdi == "NC_000017.11:43094691:G:A"


# ---------------------------------------------------------------------------
# ClinVarClient.get_by_variation_id
# ---------------------------------------------------------------------------


class TestGetByVariationId:
    """Tests for ClinVarClient.get_by_variation_id()."""

    @pytest.mark.asyncio
    async def test_vcv_prefix_stripped(self) -> None:
        client = ClinVarClient()
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_VARIATION_DIRECT
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            result = await client.get_by_variation_id("VCV000012345")

        assert result is not None
        call_url = mock_get.call_args[0][0]
        # Leading zeros stripped after "VCV" removal
        assert "/variation/12345/" in call_url

    @pytest.mark.asyncio
    async def test_bare_numeric_id_used_directly(self) -> None:
        client = ClinVarClient()
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_VARIATION_DIRECT
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            result = await client.get_by_variation_id("12345")

        assert result is not None
        call_url = mock_get.call_args[0][0]
        assert "/variation/12345/" in call_url

    @pytest.mark.asyncio
    async def test_http_error_returns_none(self) -> None:
        client = ClinVarClient()
        with patch(
            "httpx.AsyncClient.get", side_effect=httpx.HTTPError("timeout")
        ):
            result = await client.get_by_variation_id("VCV000012345")

        assert result is None

    @pytest.mark.asyncio
    async def test_api_key_included_in_params(self) -> None:
        client = ClinVarClient(api_key="secret-key")
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_VARIATION_DIRECT
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            await client.get_by_variation_id("12345")

        call_kwargs = mock_get.call_args
        params = call_kwargs.kwargs.get("params") or call_kwargs[1].get("params")
        assert params["api_key"] == "secret-key"


# ---------------------------------------------------------------------------
# ClinVarClient._lookup_by_spdi
# ---------------------------------------------------------------------------


class TestLookupBySpdi:
    """Tests for ClinVarClient._lookup_by_spdi() real HTTP body."""

    @pytest.mark.asyncio
    async def test_returns_none_on_404(self) -> None:
        """A 404 status (variant absent from ClinVar) returns None, not an error."""
        client = ClinVarClient()
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response

            result = await client._lookup_by_spdi("NC_000017.11:43094691:G:A")

        assert result is None

    @pytest.mark.asyncio
    async def test_returns_json_on_success(self) -> None:
        client = ClinVarClient()
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.raise_for_status = MagicMock()
            mock_response.json.return_value = MOCK_VARIATION_DIRECT
            mock_get.return_value = mock_response

            result = await client._lookup_by_spdi("NC_000017.11:43094691:G:A")

        assert result == MOCK_VARIATION_DIRECT

    @pytest.mark.asyncio
    async def test_api_key_added_to_params_when_set(self) -> None:
        client = ClinVarClient(api_key="my-api-key")
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.raise_for_status = MagicMock()
            mock_response.json.return_value = {}
            mock_get.return_value = mock_response

            await client._lookup_by_spdi("NC_000017.11:43094691:G:A")

        call = mock_get.call_args
        params = call.kwargs.get("params") or call[1].get("params")
        assert params["api_key"] == "my-api-key"
