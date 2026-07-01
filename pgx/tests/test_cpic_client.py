"""
pgx.tests.test_cpic_client
============================
pytest tests for the CPIC REST API client.

Tests cover:
    - _parse_recommendation: parsing raw CPIC API dicts.
    - get_recommendations: caching, HTTP success, HTTP fallback.
    - _get_builtin_cyp2d6_recommendations: built-in fallback logic.
    - DrugRecommendation dataclass construction.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from pgx.cpic_client import (
    DrugRecommendation,
    _fetch_recommendations_from_api,
    _get_builtin_cyp2d6_recommendations,
    _parse_recommendation,
    get_recommendations,
)


# ---------------------------------------------------------------------------
# DrugRecommendation dataclass
# ---------------------------------------------------------------------------


class TestDrugRecommendation:
    """Tests for DrugRecommendation dataclass."""

    def test_construction(self) -> None:
        """DrugRecommendation can be constructed with all fields."""
        rec = DrugRecommendation(
            drug_name="codeine",
            gene="CYP2D6",
            diplotype="*1/*4",
            phenotype="IM",
            classification="Use label recommended dosage",
            implications="Reduced CYP2D6 activity",
            recommendation="Use label recommended codeine dosage.",
            cpic_level="A",
        )
        assert rec.drug_name == "codeine"
        assert rec.gene == "CYP2D6"
        assert rec.diplotype == "*1/*4"
        assert rec.phenotype == "IM"
        assert rec.cpic_level == "A"

    def test_default_citations_empty(self) -> None:
        """Default citations list is empty."""
        rec = DrugRecommendation(
            drug_name="codeine",
            gene="CYP2D6",
            diplotype="*1/*1",
            phenotype="NM",
            classification="Use label",
            implications="Normal",
            recommendation="Use label",
        )
        assert rec.citations == []

    def test_default_cpic_level_is_a(self) -> None:
        """Default CPIC level is A."""
        rec = DrugRecommendation(
            drug_name="codeine",
            gene="CYP2D6",
            diplotype="*1/*1",
            phenotype="NM",
            classification="Use label",
            implications="Normal",
            recommendation="Use label",
        )
        assert rec.cpic_level == "A"


# ---------------------------------------------------------------------------
# _parse_recommendation tests
# ---------------------------------------------------------------------------


class TestParseRecommendation:
    """Tests for _parse_recommendation()."""

    def test_parses_drug_name(self) -> None:
        """Drug name is extracted from raw dict."""
        raw = {
            "drugName": "codeine",
            "phenotype": "NM",
            "classification": "Use label",
            "implications": "Normal CYP2D6",
            "recommendation": "Use label recommended dosage.",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*1/*1")
        assert rec.drug_name == "codeine"

    def test_parses_nested_drug_name(self) -> None:
        """Drug name is extracted from nested drug.name when drugName absent."""
        raw = {
            "drug": {"name": "tramadol"},
            "phenotype": "PM",
            "classification": "Avoid",
            "implications": "Poor metaboliser",
            "recommendation": "Avoid tramadol",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*4/*4")
        assert rec.drug_name == "tramadol"

    def test_parses_phenotype(self) -> None:
        """Phenotype is extracted from raw dict."""
        raw = {
            "drugName": "codeine",
            "phenotype": "PM",
            "classification": "Avoid",
            "implications": "No function",
            "recommendation": "Avoid codeine",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*4/*4")
        assert rec.phenotype == "PM"

    def test_nested_phenotype(self) -> None:
        """Phenotype falls back to phenotypes[gene] dict."""
        raw = {
            "drugName": "codeine",
            "phenotypes": {"CYP2D6": "UM"},
            "classification": "Avoid",
            "implications": "Ultra-rapid",
            "recommendation": "Avoid",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*1xN/*1")
        assert rec.phenotype == "UM"

    def test_cpic_level_default_a(self) -> None:
        """cpicLevel defaults to A when absent."""
        raw = {
            "drugName": "codeine",
            "phenotype": "NM",
            "classification": "Use label",
            "implications": "Normal",
            "recommendation": "Normal dose",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*1/*1")
        assert rec.cpic_level == "A"

    def test_cpic_level_b_parsed(self) -> None:
        """cpicLevel B is parsed correctly."""
        raw = {
            "drugName": "ondansetron",
            "phenotype": "PM",
            "classification": "Alternative",
            "implications": "Increased exposure",
            "recommendation": "Reduce dose by 25%",
            "cpicLevel": "B",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*4/*4")
        assert rec.cpic_level == "B"

    def test_gene_and_diplotype_preserved(self) -> None:
        """Gene and diplotype are preserved from arguments."""
        raw = {
            "drugName": "codeine",
            "phenotype": "NM",
            "classification": "Use label",
            "implications": "Normal",
            "recommendation": "Use label",
        }
        rec = _parse_recommendation(raw, "CYP2D6", "*2/*17")
        assert rec.gene == "CYP2D6"
        assert rec.diplotype == "*2/*17"


# ---------------------------------------------------------------------------
# _get_builtin_cyp2d6_recommendations tests
# ---------------------------------------------------------------------------


class TestGetBuiltinCyp2d6Recommendations:
    """Tests for _get_builtin_cyp2d6_recommendations()."""

    def test_returns_list_of_recommendations(self) -> None:
        """Returns a list of DrugRecommendation objects."""
        recs = _get_builtin_cyp2d6_recommendations("*1/*1")
        assert isinstance(recs, list)
        assert len(recs) >= 1
        assert all(isinstance(r, DrugRecommendation) for r in recs)

    def test_pm_diplotype_avoids_codeine(self) -> None:
        """PM diplotype (*4/*4) generates Avoid codeine recommendation."""
        recs = _get_builtin_cyp2d6_recommendations("*4/*4")
        codeine_recs = [r for r in recs if r.drug_name == "codeine"]
        assert codeine_recs
        assert "Avoid" in codeine_recs[0].recommendation

    def test_nm_diplotype_normal_codeine(self) -> None:
        """NM diplotype (*1/*1) generates normal codeine dosage recommendation."""
        recs = _get_builtin_cyp2d6_recommendations("*1/*1")
        codeine_recs = [r for r in recs if r.drug_name == "codeine"]
        assert codeine_recs
        assert "label" in codeine_recs[0].recommendation.lower()

    def test_im_diplotype_includes_codeine(self) -> None:
        """IM diplotype (*1/*10) generates codeine recommendation."""
        recs = _get_builtin_cyp2d6_recommendations("*1/*10")
        codeine_recs = [r for r in recs if r.drug_name == "codeine"]
        assert codeine_recs
        assert len(codeine_recs[0].recommendation) > 0

    def test_includes_tramadol(self) -> None:
        """Tramadol recommendation is always included."""
        recs = _get_builtin_cyp2d6_recommendations("*1/*1")
        drug_names = [r.drug_name for r in recs]
        assert "tramadol" in drug_names

    def test_all_recommendations_have_cpic_level_a(self) -> None:
        """All built-in recommendations have CPIC level A."""
        recs = _get_builtin_cyp2d6_recommendations("*4/*4")
        assert all(r.cpic_level == "A" for r in recs)

    def test_all_recommendations_have_gene_cyp2d6(self) -> None:
        """All built-in recommendations are for CYP2D6."""
        recs = _get_builtin_cyp2d6_recommendations("*1/*4")
        assert all(r.gene == "CYP2D6" for r in recs)


# ---------------------------------------------------------------------------
# get_recommendations tests
# ---------------------------------------------------------------------------


class TestGetRecommendations:
    """Tests for get_recommendations() with caching."""

    def setup_method(self) -> None:
        """Clear recommendation cache before each test."""
        from pgx import cpic_client
        cpic_client._recommendation_cache.clear()

    def test_returns_drug_recommendations_on_api_success(self) -> None:
        """Returns list of DrugRecommendation on successful API call."""
        mock_api_data = [
            {
                "drugName": "codeine",
                "phenotype": "NM",
                "classification": "Use label recommended dosage",
                "implications": "Normal CYP2D6 activity",
                "recommendation": "Use label recommended codeine dosage.",
                "cpicLevel": "A",
                "guidelineName": "CYP2D6 Codeine Guideline",
                "url": "https://cpicpgx.org/guidelines/guideline-for-codeine/",
            }
        ]

        with patch("pgx.cpic_client._fetch_recommendations_from_api", return_value=mock_api_data):
            recs = get_recommendations("CYP2D6", "*1/*1")

        assert len(recs) == 1
        assert recs[0].drug_name == "codeine"
        assert recs[0].phenotype == "NM"

    def test_falls_back_to_builtin_on_http_error(self) -> None:
        """Falls back to built-in CYP2D6 recommendations on HTTP error."""
        import httpx

        with patch(
            "pgx.cpic_client._fetch_recommendations_from_api",
            side_effect=httpx.RequestError("Connection refused"),
        ):
            recs = get_recommendations("CYP2D6", "*4/*4")

        assert len(recs) >= 1
        assert all(isinstance(r, DrugRecommendation) for r in recs)

    def test_caches_successful_result(self) -> None:
        """Successful API call result is cached; second call uses cache."""
        mock_api_data = [
            {
                "drugName": "codeine",
                "phenotype": "NM",
                "classification": "Use label",
                "implications": "Normal",
                "recommendation": "Use label",
            }
        ]

        with patch(
            "pgx.cpic_client._fetch_recommendations_from_api",
            return_value=mock_api_data,
        ) as mock_fetch:
            get_recommendations("CYP2D6", "*1/*1")
            get_recommendations("CYP2D6", "*1/*1")  # second call

        # Should only call API once (cache hit on second call)
        assert mock_fetch.call_count == 1

    def test_gene_symbol_uppercased(self) -> None:
        """Gene symbol is uppercased in the cache key."""
        mock_api_data: list = []

        with patch(
            "pgx.cpic_client._fetch_recommendations_from_api",
            return_value=mock_api_data,
        ):
            get_recommendations("cyp2d6", "*1/*1")

        from pgx import cpic_client
        assert ("CYP2D6", "*1/*1") in cpic_client._recommendation_cache


# ---------------------------------------------------------------------------
# _fetch_recommendations_from_api tests (mocked httpx)
# ---------------------------------------------------------------------------


class TestFetchRecommendationsFromApi:
    """Tests for _fetch_recommendations_from_api() with mocked httpx.get."""

    def test_returns_list_when_api_returns_list(self) -> None:
        """A raw list response from the CPIC API is returned as-is."""
        mock_response = MagicMock()
        mock_response.json.return_value = [{"drugName": "codeine"}]

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_recommendations_from_api("CYP2D6", "*1/*1")

        assert result == [{"drugName": "codeine"}]

    def test_extracts_data_key_when_dict_response(self) -> None:
        """A dict response extracts the 'data' key as the recommendation list."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"data": [{"drugName": "tramadol"}]}

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_recommendations_from_api("CYP2D6", "*4/*4")

        assert result == [{"drugName": "tramadol"}]

    def test_returns_empty_list_when_data_key_absent(self) -> None:
        """A dict response without a 'data' key returns an empty list."""
        mock_response = MagicMock()
        mock_response.json.return_value = {}

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_recommendations_from_api("CYP2D6", "*1/*4")

        assert result == []

    def test_raise_for_status_called(self) -> None:
        """raise_for_status() is called on the response to surface HTTP errors."""
        mock_response = MagicMock()
        mock_response.json.return_value = []

        with patch("httpx.get", return_value=mock_response):
            _fetch_recommendations_from_api("CYP2D6", "*1/*1")

        mock_response.raise_for_status.assert_called_once()

    def test_passes_gene_and_diplotype_as_params(self) -> None:
        """gene and diplotype are passed as query params to httpx.get."""
        mock_response = MagicMock()
        mock_response.json.return_value = []

        with patch("httpx.get", return_value=mock_response) as mock_get:
            _fetch_recommendations_from_api("CYP2D6", "*1/*4")

        _, kwargs = mock_get.call_args
        assert kwargs["params"] == {"gene": "CYP2D6", "diplotype": "*1/*4"}
