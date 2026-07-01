"""
annotation.tests.test_clingen_client
=======================================
pytest tests for the ClinGen gene-disease validity client.

Tests cover:
    - ClinGenClient.get_gene_validity: mocked HTTP responses.
    - ClinGenClient.get_gene_validity_by_hgnc: HGNC-based lookup.
    - ClinGenClient.is_valid_disease_gene: threshold comparison.
    - _parse_response: classification parsing.
    - _highest_classification / _validity_index: ordering logic.
    - HTTP error handling.

References:
    Strande et al. 2017 PMID:28552198 (ClinGen gene validity framework).
"""

from __future__ import annotations

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from annotation.clingen_client import (
    ClinGenClient,
    GeneValidity,
    GeneValidityClassification,
    _highest_classification,
    _validity_index,
)


# ---------------------------------------------------------------------------
# _validity_index tests
# ---------------------------------------------------------------------------


class TestValidityIndex:
    """Tests for _validity_index() classification ordering."""

    def test_definitive_has_lowest_index(self) -> None:
        """Definitive is the highest evidence → lowest index."""
        assert _validity_index("Definitive") < _validity_index("Strong")

    def test_strong_before_moderate(self) -> None:
        """Strong has lower index than Moderate."""
        assert _validity_index("Strong") < _validity_index("Moderate")

    def test_moderate_before_limited(self) -> None:
        """Moderate has lower index than Limited."""
        assert _validity_index("Moderate") < _validity_index("Limited")

    def test_limited_before_disputed(self) -> None:
        """Limited has lower index than Disputed."""
        assert _validity_index("Limited") < _validity_index("Disputed")

    def test_unknown_returns_max(self) -> None:
        """Unknown classification returns highest index (weakest)."""
        max_known = _validity_index("Animal Model Only")
        assert _validity_index("TOTALLY_UNKNOWN") > max_known


# ---------------------------------------------------------------------------
# _highest_classification tests
# ---------------------------------------------------------------------------


class TestHighestClassification:
    """Tests for _highest_classification() aggregation function."""

    def test_single_definitive(self) -> None:
        """Single 'Definitive' classification returns Definitive."""
        assert _highest_classification(["Definitive"]) == "Definitive"

    def test_definitive_beats_strong(self) -> None:
        """Definitive + Strong → Definitive (highest evidence)."""
        result = _highest_classification(["Strong", "Definitive", "Limited"])
        assert result == "Definitive"

    def test_empty_returns_none(self) -> None:
        """Empty list returns None."""
        assert _highest_classification([]) is None

    def test_all_limited(self) -> None:
        """All Limited → Limited."""
        assert _highest_classification(["Limited", "Limited"]) == "Limited"


# ---------------------------------------------------------------------------
# ClinGenClient._parse_response tests
# ---------------------------------------------------------------------------


MOCK_BRCA1_RESPONSE = {
    "results": [
        {
            "gene": {"symbol": "BRCA1", "hgnc_id": "HGNC:1100"},
            "disease": {
                "label": "Hereditary Breast and Ovarian Cancer Syndrome",
                "curie": "OMIM:604370",
            },
            "classification": {"label": "Definitive"},
            "mode_of_inheritance": {"label": "Autosomal dominant"},
            "date": "2024-01-15",
            "sopVersion": "7",
            "affiliation": {"id": "10002", "name": "Hereditary Cancer GCEP"},
            "report_url": "https://clinicalgenome.org/curation-activities/gene-disease-validity/affliation/10002/gene/HGNC:1100/",
        }
    ]
}

MOCK_EMPTY_RESPONSE = {"results": []}


class TestParseResponse:
    """Tests for ClinGenClient._parse_response()."""

    def setup_method(self) -> None:
        """Create a client for each test."""
        self.client = ClinGenClient()

    def test_definitive_parsed_correctly(self) -> None:
        """Definitive classification is parsed from response."""
        result = self.client._parse_response("BRCA1", MOCK_BRCA1_RESPONSE)
        assert result is not None
        assert result.highest_classification == "Definitive"
        assert result.gene_symbol == "BRCA1"

    def test_hgnc_id_extracted(self) -> None:
        """HGNC ID is extracted from the gene entry."""
        result = self.client._parse_response("BRCA1", MOCK_BRCA1_RESPONSE)
        assert result is not None
        assert result.hgnc_id == "HGNC:1100"

    def test_definitive_diseases_populated(self) -> None:
        """Definitive disease labels are collected."""
        result = self.client._parse_response("BRCA1", MOCK_BRCA1_RESPONSE)
        assert result is not None
        assert len(result.definitive_diseases) == 1
        assert "Hereditary" in result.definitive_diseases[0]

    def test_is_diagnostic_true_for_definitive(self) -> None:
        """Classification with Definitive label has is_diagnostic=True."""
        result = self.client._parse_response("BRCA1", MOCK_BRCA1_RESPONSE)
        assert result is not None
        assert result.classifications[0].is_diagnostic is True

    def test_empty_results_returns_none(self) -> None:
        """Empty results in response returns None."""
        result = self.client._parse_response("UNKNOWN", MOCK_EMPTY_RESPONSE)
        assert result is None

    def test_lof_mechanism_inferred_from_ad_moi(self) -> None:
        """Autosomal dominant MOI → has_lof_disease_mechanism=True."""
        result = self.client._parse_response("BRCA1", MOCK_BRCA1_RESPONSE)
        assert result is not None
        assert result.has_lof_disease_mechanism is True

    def test_strong_diseases_populated(self) -> None:
        """Strong diseases are extracted separately."""
        response = {
            "results": [
                {
                    "gene": {"symbol": "BRCA2", "hgnc_id": "HGNC:1101"},
                    "disease": {"label": "Fanconi Anemia", "curie": "OMIM:605724"},
                    "classification": {"label": "Strong"},
                    "mode_of_inheritance": {"label": "Autosomal recessive"},
                    "date": "2023-06-01",
                    "sopVersion": "7",
                    "affiliation": {},
                }
            ]
        }
        result = self.client._parse_response("BRCA2", response)
        assert result is not None
        assert len(result.strong_diseases) == 1
        assert result.highest_classification == "Strong"


# ---------------------------------------------------------------------------
# ClinGenClient.get_gene_validity tests
# ---------------------------------------------------------------------------


class TestGetGeneValidity:
    """Tests for ClinGenClient.get_gene_validity()."""

    @pytest.mark.asyncio
    async def test_returns_gene_validity_on_success(self) -> None:
        """Returns GeneValidity object for a known gene."""
        client = ClinGenClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get_method:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_BRCA1_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get_method.return_value = mock_response

            result = await client.get_gene_validity("BRCA1")

        assert result is not None
        assert result.highest_classification == "Definitive"

    @pytest.mark.asyncio
    async def test_returns_none_on_http_error(self) -> None:
        """Returns None when HTTP request fails."""
        client = ClinGenClient()

        import httpx
        with patch("httpx.AsyncClient.get", side_effect=httpx.HTTPError("connection failed")):
            result = await client.get_gene_validity("BRCA1")

        assert result is None

    @pytest.mark.asyncio
    async def test_returns_none_for_unknown_gene(self) -> None:
        """Returns None when no ClinGen curations exist."""
        client = ClinGenClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get_method:
            mock_response = MagicMock()
            mock_response.json.return_value = {"results": []}
            mock_response.raise_for_status = MagicMock()
            mock_get_method.return_value = mock_response

            result = await client.get_gene_validity("UNKNOWN_GENE_XYZ")

        assert result is None


# ---------------------------------------------------------------------------
# ClinGenClient.is_valid_disease_gene tests
# ---------------------------------------------------------------------------


class TestIsValidDiseaseGene:
    """Tests for ClinGenClient.is_valid_disease_gene()."""

    def setup_method(self) -> None:
        """Create a client for each test."""
        self.client = ClinGenClient()

    def _make_validity(self, classification: str) -> GeneValidity:
        """Create a GeneValidity with a given highest_classification."""
        return GeneValidity(
            hgnc_id="HGNC:1100",
            gene_symbol="BRCA1",
            highest_classification=classification,
            classifications=[],
        )

    def test_definitive_is_valid_for_strong_minimum(self) -> None:
        """Definitive classification passes Strong minimum threshold."""
        validity = self._make_validity("Definitive")
        assert self.client.is_valid_disease_gene(validity, min_level="Strong") is True

    def test_strong_is_valid_for_strong_minimum(self) -> None:
        """Strong classification meets Strong minimum threshold."""
        validity = self._make_validity("Strong")
        assert self.client.is_valid_disease_gene(validity, min_level="Strong") is True

    def test_moderate_fails_strong_minimum(self) -> None:
        """Moderate does not meet Strong minimum threshold."""
        validity = self._make_validity("Moderate")
        assert self.client.is_valid_disease_gene(validity, min_level="Strong") is False

    def test_limited_fails_strong_minimum(self) -> None:
        """Limited does not meet Strong minimum threshold."""
        validity = self._make_validity("Limited")
        assert self.client.is_valid_disease_gene(validity, min_level="Strong") is False

    def test_none_validity_returns_false(self) -> None:
        """None validity (gene not in ClinGen) returns False."""
        assert self.client.is_valid_disease_gene(None, min_level="Strong") is False

    def test_none_classification_returns_false(self) -> None:
        """None highest_classification returns False."""
        validity = GeneValidity(
            hgnc_id=None,
            gene_symbol="GENE",
            highest_classification=None,
        )
        assert self.client.is_valid_disease_gene(validity) is False

    def test_moderate_is_valid_for_moderate_minimum(self) -> None:
        """Moderate meets Moderate minimum threshold."""
        validity = self._make_validity("Moderate")
        assert self.client.is_valid_disease_gene(validity, min_level="Moderate") is True

    def test_limited_fails_moderate_minimum(self) -> None:
        """Limited does not meet Moderate minimum threshold."""
        validity = self._make_validity("Limited")
        assert self.client.is_valid_disease_gene(validity, min_level="Moderate") is False


# ---------------------------------------------------------------------------
# GeneValidityClassification dataclass
# ---------------------------------------------------------------------------


class TestGeneValidityClassification:
    """Tests for GeneValidityClassification dataclass."""

    def test_is_diagnostic_true_for_definitive(self) -> None:
        """Definitive classification has is_diagnostic=True."""
        clsf = GeneValidityClassification(
            disease_label="HBOC",
            classification="Definitive",
            is_diagnostic=True,
        )
        assert clsf.is_diagnostic is True

    def test_is_diagnostic_false_for_limited(self) -> None:
        """Limited classification has is_diagnostic=False."""
        clsf = GeneValidityClassification(
            disease_label="Unknown",
            classification="Limited",
            is_diagnostic=False,
        )
        assert clsf.is_diagnostic is False


# ---------------------------------------------------------------------------
# ClinGenClient.get_gene_validity_by_hgnc tests
# ---------------------------------------------------------------------------


class TestGetGeneValidityByHgnc:
    """Tests for ClinGenClient.get_gene_validity_by_hgnc()."""

    @pytest.mark.asyncio
    async def test_returns_gene_validity_on_success(self) -> None:
        """Returns a GeneValidity built from the response's gene symbol."""
        client = ClinGenClient()

        with patch("httpx.AsyncClient.get") as mock_get_method:
            mock_response = MagicMock()
            payload = dict(MOCK_BRCA1_RESPONSE)
            payload["gene"] = {"symbol": "BRCA1"}
            mock_response.json.return_value = payload
            mock_response.raise_for_status = MagicMock()
            mock_get_method.return_value = mock_response

            result = await client.get_gene_validity_by_hgnc("HGNC:1100")

        assert result is not None
        assert result.gene_symbol == "BRCA1"
        assert result.highest_classification == "Definitive"

    @pytest.mark.asyncio
    async def test_falls_back_to_hgnc_id_when_gene_symbol_missing(self) -> None:
        """When the response has no 'gene.symbol', hgnc_id is used as the
        gene_symbol fallback."""
        client = ClinGenClient()

        with patch("httpx.AsyncClient.get") as mock_get_method:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_BRCA1_RESPONSE  # no "gene" key
            mock_response.raise_for_status = MagicMock()
            mock_get_method.return_value = mock_response

            result = await client.get_gene_validity_by_hgnc("HGNC:1100")

        assert result is not None
        assert result.gene_symbol == "HGNC:1100"

    @pytest.mark.asyncio
    async def test_returns_none_on_http_error(self) -> None:
        """Returns None when the HTTP request fails."""
        client = ClinGenClient()

        import httpx

        with patch(
            "httpx.AsyncClient.get", side_effect=httpx.HTTPError("connection failed")
        ):
            result = await client.get_gene_validity_by_hgnc("HGNC:1100")

        assert result is None

    @pytest.mark.asyncio
    async def test_hgnc_id_passed_as_query_param(self) -> None:
        """hgnc_id should be included in the request query params."""
        client = ClinGenClient()

        with patch("httpx.AsyncClient.get") as mock_get_method:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_EMPTY_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get_method.return_value = mock_response

            await client.get_gene_validity_by_hgnc("HGNC:1100")

        call = mock_get_method.call_args
        params = call.kwargs.get("params") or call[1].get("params")
        assert params["hgnc_id"] == "HGNC:1100"
