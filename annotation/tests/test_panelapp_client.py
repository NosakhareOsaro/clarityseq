"""
annotation.tests.test_panelapp_client
========================================
pytest tests for the Genomics England PanelApp client.

Tests cover:
    - PanelAppClient.get_gene_panels: mocked HTTP responses.
    - PanelAppClient.get_panel_by_id: panel metadata lookup.
    - PanelAppClient.get_gene_entry: full gene entry with green_panels.
    - PanelAppClient.search_panels: text-based panel search.
    - _parse_gene_panels: confidence level filtering and sorting.
    - Panel and GeneEntry dataclasses.

References:
    Martin et al. 2019 Nature Genetics PMID:31676867 (PanelApp).
    Gene confidence levels: 3=Green (diagnostic), 2=Amber, 1=Red.
"""

from __future__ import annotations

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from annotation.panelapp_client import (
    GeneEntry,
    Panel,
    PanelAppClient,
    _safe_int,
)


# ---------------------------------------------------------------------------
# Mock API responses
# ---------------------------------------------------------------------------

MOCK_BRCA1_PANELS_RESPONSE = {
    "results": [
        {
            "panel": {
                "id": 510,
                "name": "Hereditary Breast and Ovarian Cancer",
                "disease_group": "Tumour Syndromes",
                "disease_sub_group": "Inherited Breast Cancer",
                "status": "public",
                "version": "2.18",
            },
            "confidence_level": "3",  # Green — diagnostic
            "mode_of_inheritance": "MONOALLELIC",
            "phenotypes": ["Hereditary Breast and Ovarian Cancer", "OMIM:604370"],
        },
        {
            "panel": {
                "id": 285,
                "name": "Rare Disease Tier 1",
                "disease_group": "Rare Disease",
                "disease_sub_group": "Multisystem",
                "status": "public",
                "version": "5.7",
            },
            "confidence_level": "2",  # Amber — not diagnostic
            "mode_of_inheritance": "MONOALLELIC",
            "phenotypes": ["Breast cancer"],
        },
        {
            "panel": {
                "id": 999,
                "name": "Low Confidence Panel",
                "disease_group": "Cancer",
                "disease_sub_group": "",
                "status": "public",
                "version": "1.0",
            },
            "confidence_level": "1",  # Red — below min_confidence=2
            "mode_of_inheritance": "",
            "phenotypes": [],
        },
    ]
}

MOCK_EMPTY_RESPONSE = {"results": []}


# ---------------------------------------------------------------------------
# PanelAppClient._parse_gene_panels tests
# ---------------------------------------------------------------------------


class TestParsePanelGenes:
    """Tests for PanelAppClient._parse_gene_panels()."""

    def setup_method(self) -> None:
        """Create a client with default min_confidence=2."""
        self.client = PanelAppClient(min_confidence=2)

    def test_green_panel_included(self) -> None:
        """Panel with confidence 3 (Green) is included."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        panel_ids = [p.panel_id for p in panels]
        assert 510 in panel_ids  # Green panel

    def test_amber_panel_included(self) -> None:
        """Panel with confidence 2 (Amber) is included at min_confidence=2."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        panel_ids = [p.panel_id for p in panels]
        assert 285 in panel_ids

    def test_red_panel_excluded_at_default_confidence(self) -> None:
        """Panel with confidence 1 (Red) is excluded at min_confidence=2."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        panel_ids = [p.panel_id for p in panels]
        assert 999 not in panel_ids

    def test_green_panels_sorted_first(self) -> None:
        """Green panels (confidence 3) appear before Amber (2) in sorted output."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        if len(panels) >= 2:
            assert panels[0].gene_confidence >= panels[-1].gene_confidence

    def test_is_diagnostic_true_for_green(self) -> None:
        """Panel with confidence 3 has is_diagnostic=True."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        green = [p for p in panels if p.gene_confidence == 3]
        assert all(p.is_diagnostic for p in green)

    def test_is_diagnostic_false_for_amber(self) -> None:
        """Panel with confidence 2 has is_diagnostic=False."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        amber = [p for p in panels if p.gene_confidence == 2]
        assert all(not p.is_diagnostic for p in amber)

    def test_panel_url_populated(self) -> None:
        """Panel URL is set from panel ID."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        green = [p for p in panels if p.panel_id == 510]
        assert green
        assert "510" in (green[0].panel_url or "")

    def test_phenotypes_list_populated(self) -> None:
        """Non-empty phenotype strings are included in phenotypes list."""
        panels = self.client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        green = [p for p in panels if p.panel_id == 510]
        assert green
        assert len(green[0].phenotypes) >= 1

    def test_empty_results_returns_empty_list(self) -> None:
        """Empty results dict returns empty list."""
        panels = self.client._parse_gene_panels(MOCK_EMPTY_RESPONSE)
        assert panels == []

    def test_min_confidence_3_excludes_amber(self) -> None:
        """min_confidence=3 excludes Amber panels."""
        client = PanelAppClient(min_confidence=3)
        panels = client._parse_gene_panels(MOCK_BRCA1_PANELS_RESPONSE)
        panel_ids = [p.panel_id for p in panels]
        assert 285 not in panel_ids  # Amber excluded
        assert 510 in panel_ids      # Green still included


# ---------------------------------------------------------------------------
# PanelAppClient.get_gene_panels tests
# ---------------------------------------------------------------------------


class TestGetGenePanels:
    """Tests for PanelAppClient.get_gene_panels()."""

    @pytest.mark.asyncio
    async def test_returns_panels_on_success(self) -> None:
        """Returns list of Panel objects for a known gene."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_BRCA1_PANELS_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            panels = await client.get_gene_panels("BRCA1")

        assert len(panels) >= 1
        assert all(isinstance(p, Panel) for p in panels)

    @pytest.mark.asyncio
    async def test_returns_empty_list_on_http_error(self) -> None:
        """Returns empty list on HTTP error (no raise)."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get", side_effect=httpx.HTTPError("503")):
            panels = await client.get_gene_panels("BRCA1")

        assert panels == []

    @pytest.mark.asyncio
    async def test_returns_empty_for_unknown_gene(self) -> None:
        """Returns empty list for a gene with no PanelApp entries."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_EMPTY_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            panels = await client.get_gene_panels("UNKNOWN_GENE")

        assert panels == []

    @pytest.mark.asyncio
    async def test_gene_symbol_uppercased_in_url(self) -> None:
        """Gene symbol is uppercased in the API URL."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_EMPTY_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            await client.get_gene_panels("brca1")

        call_url = mock_get.call_args[0][0]
        assert "BRCA1" in call_url


# ---------------------------------------------------------------------------
# PanelAppClient.get_panel_by_id tests
# ---------------------------------------------------------------------------


class TestGetPanelById:
    """Tests for PanelAppClient.get_panel_by_id()."""

    @pytest.mark.asyncio
    async def test_returns_dict_on_success(self) -> None:
        """Returns panel dict on successful response."""
        client = PanelAppClient()
        mock_data = {"id": 510, "name": "Hereditary Breast Cancer"}

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = mock_data
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            result = await client.get_panel_by_id(510)

        assert result is not None
        assert result["id"] == 510

    @pytest.mark.asyncio
    async def test_returns_none_on_http_error(self) -> None:
        """Returns None on HTTP error."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get", side_effect=httpx.HTTPError("404")):
            result = await client.get_panel_by_id(9999)

        assert result is None


# ---------------------------------------------------------------------------
# PanelAppClient.get_gene_entry tests
# ---------------------------------------------------------------------------


class TestGetGeneEntry:
    """Tests for PanelAppClient.get_gene_entry()."""

    @pytest.mark.asyncio
    async def test_returns_gene_entry_with_green_panels(self) -> None:
        """GeneEntry includes green_panels for Green confidence levels."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_BRCA1_PANELS_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            entry = await client.get_gene_entry("BRCA1")

        assert entry is not None
        assert isinstance(entry, GeneEntry)
        assert len(entry.green_panels) >= 1
        assert entry.gene_symbol == "BRCA1"

    @pytest.mark.asyncio
    async def test_returns_none_when_no_panels(self) -> None:
        """Returns None when gene has no PanelApp entries."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_EMPTY_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            entry = await client.get_gene_entry("UNKNOWN_GENE")

        assert entry is None

    @pytest.mark.asyncio
    async def test_highest_confidence_populated(self) -> None:
        """GeneEntry.highest_confidence is the maximum confidence level."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = MOCK_BRCA1_PANELS_RESPONSE
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            entry = await client.get_gene_entry("BRCA1")

        assert entry is not None
        assert entry.highest_confidence == 3  # Green


# ---------------------------------------------------------------------------
# PanelAppClient.search_panels tests
# ---------------------------------------------------------------------------


class TestSearchPanels:
    """Tests for PanelAppClient.search_panels()."""

    @pytest.mark.asyncio
    async def test_returns_results_list(self) -> None:
        """Returns list of panel dicts matching search."""
        client = PanelAppClient()
        mock_data = {
            "results": [
                {"id": 510, "name": "Hereditary Breast and Ovarian Cancer"},
            ]
        }

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = mock_data
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            results = await client.search_panels(name="breast")

        assert len(results) == 1
        assert results[0]["id"] == 510

    @pytest.mark.asyncio
    async def test_returns_empty_on_http_error(self) -> None:
        """Returns empty list on HTTP error."""
        client = PanelAppClient()

        import httpx
        with patch("httpx.AsyncClient.get", side_effect=httpx.HTTPError("503")):
            results = await client.search_panels(name="breast")

        assert results == []

    @pytest.mark.asyncio
    async def test_disease_group_param_included(self) -> None:
        """disease_group filter should be included in the request params."""
        client = PanelAppClient()
        mock_data = {"results": [{"id": 1, "name": "Cancer Panel"}]}

        import httpx
        with patch("httpx.AsyncClient.get") as mock_get:
            mock_response = MagicMock()
            mock_response.json.return_value = mock_data
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            results = await client.search_panels(disease_group="Cancer")

        assert results == mock_data["results"]
        call = mock_get.call_args
        params = call.kwargs.get("params") or call[1].get("params")
        assert params["disease_group"] == "Cancer"

    @pytest.mark.asyncio
    async def test_result_none_returns_empty_list(self) -> None:
        """When _get returns a falsy value, an empty list is returned."""
        client = PanelAppClient()
        with patch.object(client, "_get", new_callable=AsyncMock) as mock_inner:
            mock_inner.return_value = None
            results = await client.search_panels(name="none-such")

        assert results == []


# ---------------------------------------------------------------------------
# PanelAppClient.get_gene_panels — raw falsy response branch
# ---------------------------------------------------------------------------


class TestGetGenePanelsFalsyRaw:
    """Covers the `if not raw: return []` branch (raw is not an HTTPError)."""

    @pytest.mark.asyncio
    async def test_none_response_returns_empty_list(self) -> None:
        client = PanelAppClient()
        with patch.object(client, "_get", new_callable=AsyncMock) as mock_inner:
            mock_inner.return_value = None
            panels = await client.get_gene_panels("BRCA1")

        assert panels == []

    @pytest.mark.asyncio
    async def test_empty_dict_response_returns_empty_list(self) -> None:
        client = PanelAppClient()
        with patch.object(client, "_get", new_callable=AsyncMock) as mock_inner:
            mock_inner.return_value = {}
            panels = await client.get_gene_panels("BRCA1")

        assert panels == []


# ---------------------------------------------------------------------------
# PanelAppClient._parse_gene_panels — missing panel_info branch
# ---------------------------------------------------------------------------


class TestParseGenePanelsMissingPanelInfo:
    """Covers the `if not panel_info: continue` branch."""

    def test_item_without_panel_key_is_skipped(self) -> None:
        client = PanelAppClient(min_confidence=2)
        data = {
            "results": [
                {
                    # No "panel" key at all
                    "confidence_level": "3",
                    "mode_of_inheritance": "MONOALLELIC",
                    "phenotypes": [],
                },
                {
                    "panel": {
                        "id": 510,
                        "name": "Hereditary Breast and Ovarian Cancer",
                        "status": "public",
                        "version": "2.18",
                    },
                    "confidence_level": "3",
                    "mode_of_inheritance": "MONOALLELIC",
                    "phenotypes": [],
                },
            ]
        }

        panels = client._parse_gene_panels(data)

        assert len(panels) == 1
        assert panels[0].panel_id == 510

    def test_item_with_empty_panel_dict_is_skipped(self) -> None:
        client = PanelAppClient(min_confidence=2)
        data = {"results": [{"panel": {}, "confidence_level": "3"}]}

        panels = client._parse_gene_panels(data)

        assert panels == []


# ---------------------------------------------------------------------------
# _safe_int tests
# ---------------------------------------------------------------------------


class TestSafeInt:
    """Tests for the module-level _safe_int() helper."""

    def test_valid_string_converts(self) -> None:
        assert _safe_int("3") == 3

    def test_none_returns_none(self) -> None:
        assert _safe_int(None) is None

    def test_invalid_string_returns_none(self) -> None:
        """A ValueError during int() conversion is caught."""
        assert _safe_int("not_a_number") is None

    def test_uncoercible_type_returns_none(self) -> None:
        """A TypeError during int() conversion (e.g. a dict) is caught."""
        assert _safe_int({"a": 1}) is None
