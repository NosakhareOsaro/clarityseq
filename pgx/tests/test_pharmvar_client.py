"""
pgx.tests.test_pharmvar_client
================================
pytest tests for the PharmVar REST API client.

Tests cover:
    - _parse_star_allele: raw dict parsing and activity mapping.
    - StarAllele properties: is_no_function, is_normal_function.
    - get_star_allele_activity: lookup with cache, default fallback.
    - compute_diplotype_activity_score: sum of two allele activities.
    - get_star_alleles_for_gene: caching and error propagation.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from pgx.pharmvar_client import (
    StarAllele,
    _fetch_gene_alleles_from_api,
    _parse_star_allele,
    compute_diplotype_activity_score,
    get_star_allele_activity,
    get_star_alleles_for_gene,
)


# ---------------------------------------------------------------------------
# StarAllele properties
# ---------------------------------------------------------------------------


class TestStarAlleleProperties:
    """Tests for StarAllele.is_no_function and is_normal_function."""

    def test_is_no_function_true(self) -> None:
        """'No Function' status returns True for is_no_function."""
        allele = StarAllele(
            allele_id="4",
            name="*4",
            gene="CYP2D6",
            function="No Function",
            activity_value=0.0,
        )
        assert allele.is_no_function is True

    def test_is_no_function_false_for_normal(self) -> None:
        """'Normal Function' status returns False for is_no_function."""
        allele = StarAllele(
            allele_id="1",
            name="*1",
            gene="CYP2D6",
            function="Normal Function",
            activity_value=1.0,
        )
        assert allele.is_no_function is False

    def test_is_normal_function_true(self) -> None:
        """'Normal Function' status returns True for is_normal_function."""
        allele = StarAllele(
            allele_id="1",
            name="*1",
            gene="CYP2D6",
            function="Normal Function",
            activity_value=1.0,
        )
        assert allele.is_normal_function is True

    def test_is_normal_function_false_for_no_function(self) -> None:
        """'No Function' status returns False for is_normal_function."""
        allele = StarAllele(
            allele_id="4",
            name="*4",
            gene="CYP2D6",
            function="No Function",
            activity_value=0.0,
        )
        assert allele.is_normal_function is False

    def test_is_normal_function_false_for_decreased(self) -> None:
        """'Decreased Function' returns False for is_normal_function."""
        allele = StarAllele(
            allele_id="10",
            name="*10",
            gene="CYP2D6",
            function="Decreased Function",
            activity_value=0.5,
        )
        assert allele.is_normal_function is False


# ---------------------------------------------------------------------------
# _parse_star_allele tests
# ---------------------------------------------------------------------------


class TestParseStarAllele:
    """Tests for _parse_star_allele()."""

    def test_parses_allele_name(self) -> None:
        """Allele name is extracted from alleleName field."""
        raw = {
            "alleleName": "*4",
            "id": "42",
            "functionStatus": "No Function",
        }
        allele = _parse_star_allele(raw)
        assert allele.name == "*4"

    def test_parses_name_fallback(self) -> None:
        """Allele name falls back to name field when alleleName absent."""
        raw = {
            "name": "*1",
            "id": "1",
            "functionStatus": "Normal Function",
        }
        allele = _parse_star_allele(raw)
        assert allele.name == "*1"

    def test_no_function_activity_is_zero(self) -> None:
        """No Function → activity_value=0.0."""
        raw = {"alleleName": "*4", "id": "4", "functionStatus": "No Function"}
        allele = _parse_star_allele(raw)
        assert allele.activity_value == 0.0

    def test_normal_function_activity_is_one(self) -> None:
        """Normal Function → activity_value=1.0."""
        raw = {"alleleName": "*1", "id": "1", "functionStatus": "Normal Function"}
        allele = _parse_star_allele(raw)
        assert allele.activity_value == 1.0

    def test_decreased_function_activity_is_half(self) -> None:
        """Decreased Function → activity_value=0.5."""
        raw = {"alleleName": "*10", "id": "10", "functionStatus": "Decreased Function"}
        allele = _parse_star_allele(raw)
        assert allele.activity_value == 0.5

    def test_increased_function_activity_is_two(self) -> None:
        """Increased Function → activity_value=2.0."""
        raw = {"alleleName": "*1xN", "id": "xn", "functionStatus": "Increased Function"}
        allele = _parse_star_allele(raw)
        assert allele.activity_value == 2.0

    def test_unknown_function_defaults_to_half(self) -> None:
        """Unknown Function defaults to activity_value=0.5 (conservative)."""
        raw = {"alleleName": "*99", "id": "99", "functionStatus": "Unknown Function"}
        allele = _parse_star_allele(raw)
        assert allele.activity_value == 0.5

    def test_haplotype_name_populated(self) -> None:
        """haplotype_name is gene+allele_name."""
        raw = {"alleleName": "*4", "id": "4", "functionStatus": "No Function"}
        allele = _parse_star_allele(raw, gene="CYP2D6")
        assert allele.haplotype_name == "CYP2D6*4"

    def test_defining_variants_empty_by_default(self) -> None:
        """defining_variants is empty when variants key absent."""
        raw = {"alleleName": "*1", "id": "1", "functionStatus": "Normal Function"}
        allele = _parse_star_allele(raw)
        assert allele.defining_variants == []

    def test_defining_variants_parsed(self) -> None:
        """defining_variants list is parsed from variants key."""
        raw = {
            "alleleName": "*4",
            "id": "4",
            "functionStatus": "No Function",
            "variants": [{"pos": 42126611, "ref": "C", "alt": "T"}],
        }
        allele = _parse_star_allele(raw)
        assert len(allele.defining_variants) == 1


# ---------------------------------------------------------------------------
# get_star_allele_activity tests
# ---------------------------------------------------------------------------


class TestGetStarAlleleActivity:
    """Tests for get_star_allele_activity()."""

    def setup_method(self) -> None:
        """Clear PharmVar gene allele cache."""
        from pgx import pharmvar_client
        pharmvar_client._gene_allele_cache.clear()

    def test_returns_activity_from_cached_alleles(self) -> None:
        """Returns activity value from successfully fetched allele list."""
        mock_alleles = [
            StarAllele("1", "*1", "CYP2D6", "Normal Function", 1.0),
            StarAllele("4", "*4", "CYP2D6", "No Function", 0.0),
        ]
        with patch("pgx.pharmvar_client.get_star_alleles_for_gene", return_value=mock_alleles):
            assert get_star_allele_activity("*1") == pytest.approx(1.0)
            assert get_star_allele_activity("*4") == pytest.approx(0.0)

    def test_adds_asterisk_if_missing(self) -> None:
        """Allele name without leading asterisk is normalised to *4."""
        mock_alleles = [
            StarAllele("4", "*4", "CYP2D6", "No Function", 0.0),
        ]
        with patch("pgx.pharmvar_client.get_star_alleles_for_gene", return_value=mock_alleles):
            result = get_star_allele_activity("4")
        assert result == pytest.approx(0.0)

    def test_falls_back_to_default_on_api_error(self) -> None:
        """Falls back to default activity table on API error."""
        with patch(
            "pgx.pharmvar_client.get_star_alleles_for_gene",
            side_effect=Exception("Network error"),
        ):
            # *4 is in the default table with activity 0.0
            result = get_star_allele_activity("*4")
        assert result == pytest.approx(0.0)

    def test_unknown_allele_returns_none_when_api_fails(self) -> None:
        """Returns None for unknown allele when API fails."""
        with patch(
            "pgx.pharmvar_client.get_star_alleles_for_gene",
            side_effect=Exception("Network error"),
        ):
            result = get_star_allele_activity("*9999")
        assert result is None

    def test_default_star1_activity(self) -> None:
        """Default activity for *1 is 1.0."""
        with patch(
            "pgx.pharmvar_client.get_star_alleles_for_gene",
            side_effect=Exception("Unavailable"),
        ):
            result = get_star_allele_activity("*1")
        assert result == pytest.approx(1.0)

    def test_default_star5_activity_is_zero(self) -> None:
        """Default activity for *5 (gene deletion) is 0.0."""
        with patch(
            "pgx.pharmvar_client.get_star_alleles_for_gene",
            side_effect=Exception("Unavailable"),
        ):
            result = get_star_allele_activity("*5")
        assert result == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# compute_diplotype_activity_score tests
# ---------------------------------------------------------------------------


class TestComputeDiplotypeActivityScore:
    """Tests for compute_diplotype_activity_score()."""

    def test_star1_star1_gives_two(self) -> None:
        """*1/*1 (NM/NM) → activity score 2.0."""
        mock_alleles = [StarAllele("1", "*1", "CYP2D6", "Normal Function", 1.0)]
        with patch("pgx.pharmvar_client.get_star_alleles_for_gene", return_value=mock_alleles):
            score = compute_diplotype_activity_score("*1", "*1")
        assert score == pytest.approx(2.0)

    def test_star1_star9_gives_one_and_half(self) -> None:
        """*1/*9 (NM/Decreased) → activity score 1.5."""
        mock_alleles = [
            StarAllele("1", "*1", "CYP2D6", "Normal Function", 1.0),
            StarAllele("9", "*9", "CYP2D6", "Decreased Function", 0.5),
        ]
        with patch("pgx.pharmvar_client.get_star_alleles_for_gene", return_value=mock_alleles):
            score = compute_diplotype_activity_score("*1", "*9")
        assert score == pytest.approx(1.5)

    def test_star9_star9_gives_one(self) -> None:
        """*9/*9 (Decreased/Decreased) → activity score 1.0."""
        mock_alleles = [StarAllele("9", "*9", "CYP2D6", "Decreased Function", 0.5)]
        with patch("pgx.pharmvar_client.get_star_alleles_for_gene", return_value=mock_alleles):
            score = compute_diplotype_activity_score("*9", "*9")
        assert score == pytest.approx(1.0)

    def test_unknown_allele_defaults_to_one(self) -> None:
        """Unknown allele defaults to activity 1.0 (conservative assumption)."""
        with patch(
            "pgx.pharmvar_client.get_star_alleles_for_gene",
            side_effect=Exception("Unavailable"),
        ):
            # *9999 is unknown → fallback to 1.0 each
            score = compute_diplotype_activity_score("*9999", "*9999")
        assert score == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# get_star_alleles_for_gene tests
# ---------------------------------------------------------------------------


class TestGetStarAllelesForGene:
    """Tests for get_star_alleles_for_gene()."""

    def setup_method(self) -> None:
        """Clear cache before each test."""
        from pgx import pharmvar_client
        pharmvar_client._gene_allele_cache.clear()

    def test_returns_star_alleles_on_success(self) -> None:
        """Returns list of StarAllele on successful API call."""
        mock_data = [
            {"alleleName": "*1", "id": "1", "functionStatus": "Normal Function"},
            {"alleleName": "*4", "id": "4", "functionStatus": "No Function"},
        ]
        with patch("pgx.pharmvar_client._fetch_gene_alleles_from_api", return_value=mock_data):
            alleles = get_star_alleles_for_gene("CYP2D6")
        assert len(alleles) == 2
        assert all(isinstance(a, StarAllele) for a in alleles)

    def test_alleles_sorted_by_name(self) -> None:
        """Alleles are sorted alphabetically by name."""
        mock_data = [
            {"alleleName": "*9", "id": "9", "functionStatus": "Decreased Function"},
            {"alleleName": "*1", "id": "1", "functionStatus": "Normal Function"},
        ]
        with patch("pgx.pharmvar_client._fetch_gene_alleles_from_api", return_value=mock_data):
            alleles = get_star_alleles_for_gene("CYP2D6")
        assert alleles[0].name < alleles[1].name

    def test_result_cached(self) -> None:
        """Second call uses cache — API not called twice."""
        mock_data = [{"alleleName": "*1", "id": "1", "functionStatus": "Normal Function"}]
        with patch(
            "pgx.pharmvar_client._fetch_gene_alleles_from_api",
            return_value=mock_data,
        ) as mock_fetch:
            get_star_alleles_for_gene("CYP2D6")
            get_star_alleles_for_gene("CYP2D6")
        assert mock_fetch.call_count == 1

    def test_propagates_http_error(self) -> None:
        """HTTP error from API is propagated (no silent fallback)."""
        import httpx

        with (
            patch(
                "pgx.pharmvar_client._fetch_gene_alleles_from_api",
                side_effect=httpx.RequestError("Connection refused"),
            ),
            pytest.raises(httpx.RequestError),
        ):
            get_star_alleles_for_gene("CYP2D6")


# ---------------------------------------------------------------------------
# _fetch_gene_alleles_from_api tests (mocked httpx)
# ---------------------------------------------------------------------------


class TestFetchGeneAllelesFromApi:
    """Tests for _fetch_gene_alleles_from_api() with mocked httpx.get."""

    def test_returns_list_when_api_returns_list(self) -> None:
        """A raw list response from the PharmVar API is returned as-is."""
        mock_response = MagicMock()
        mock_response.json.return_value = [{"alleleName": "*1"}]

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_gene_alleles_from_api("CYP2D6")

        assert result == [{"alleleName": "*1"}]

    def test_extracts_data_key(self) -> None:
        """A dict response extracts the 'data' key as the allele list."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"data": [{"alleleName": "*4"}]}

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_gene_alleles_from_api("CYP2D6")

        assert result == [{"alleleName": "*4"}]

    def test_extracts_alleles_key_fallback(self) -> None:
        """A dict response falls back to the 'alleles' key when 'data' is absent."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"alleles": [{"alleleName": "*9"}]}

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_gene_alleles_from_api("CYP2D6")

        assert result == [{"alleleName": "*9"}]

    def test_returns_empty_list_when_no_known_key(self) -> None:
        """A dict response with neither 'data' nor 'alleles' returns an empty list."""
        mock_response = MagicMock()
        mock_response.json.return_value = {}

        with patch("httpx.get", return_value=mock_response):
            result = _fetch_gene_alleles_from_api("CYP2D6")

        assert result == []

    def test_passes_gene_param(self) -> None:
        """gene is passed as a query param to httpx.get."""
        mock_response = MagicMock()
        mock_response.json.return_value = []

        with patch("httpx.get", return_value=mock_response) as mock_get:
            _fetch_gene_alleles_from_api("CYP2D6")

        _, kwargs = mock_get.call_args
        assert kwargs["params"] == {"gene": "CYP2D6"}

    def test_raise_for_status_called(self) -> None:
        """raise_for_status() is called on the response to surface HTTP errors."""
        mock_response = MagicMock()
        mock_response.json.return_value = []

        with patch("httpx.get", return_value=mock_response):
            _fetch_gene_alleles_from_api("CYP2D6")

        mock_response.raise_for_status.assert_called_once()
