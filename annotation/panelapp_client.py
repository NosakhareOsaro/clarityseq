"""
Genomics England PanelApp client for gene panel membership lookup.

PanelApp is a crowd-sourced knowledgebase of virtual gene panels for rare
disease and cancer genetics.  It is the authoritative source for NHS GMS
(Genomic Medicine Service) gene panels.

References
----------
Martin AR, et al. "PanelApp crowdsources expert knowledge to establish
consensus diagnostic gene panels." Nature Genetics. 2019;51:1560–1565.
PMID:31676867. DOI:10.1038/s41588-019-0528-2

API Documentation
-----------------
https://panelapp.genomicsengland.co.uk/api/docs/
Base URL: https://panelapp.genomicsengland.co.uk/api/v1/

Gene Confidence Levels
-----------------------
Each gene in a panel has a confidence level (evidence-based rating):
  3 (Green)  — High confidence; included in diagnostic panel
  2 (Amber)  — Moderate confidence; evidence insufficient for diagnostic use
  1 (Red)    — Low confidence / disputed
  0 (None)   — No evidence

For clinical reporting, only Green (3) genes should be considered diagnostic.
Amber (2) genes may be flagged as potentially relevant.

NHS GMS Virtual Panels
-----------------------
The NHS GMS uses specific PanelApp panel IDs for different conditions.
Key panels:
  R59  — Rare and undiagnosed diseases (ID: 285)
  R73  — Intellectual disability (ID: 285)
  R208 — Hereditary breast and ovarian cancer (ID: 510)
  R25  — Familial hypercholesterolaemia (ID: 14)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

# PanelApp API base URL
PANELAPP_API_BASE: str = "https://panelapp.genomicsengland.co.uk/api/v1"

# Confidence level colour mapping
_CONFIDENCE_COLOURS: dict[int, str] = {
    3: "Green",
    2: "Amber",
    1: "Red",
    0: "None",
}

# Minimum confidence level for diagnostic reporting
MIN_DIAGNOSTIC_CONFIDENCE: int = 3  # Green genes only


@dataclass
class Panel:
    """A PanelApp virtual gene panel.

    Attributes:
        panel_id: Numeric PanelApp panel identifier.
        name: Panel display name (e.g. "Rare Disease Tier 1").
        disease_group: Disease group classification.
        disease_sub_group: Disease sub-group.
        status: Panel status ("public", "retired", etc.).
        version: Panel version string.
        gene_confidence: Confidence level for this gene in this panel (0–3).
        gene_confidence_colour: Colour label ("Green", "Amber", "Red").
        mode_of_inheritance: Inheritance mode for gene in this context.
        phenotypes: List of associated phenotype descriptions.
        is_diagnostic: True if gene is Green (level 3) in this panel.
        panel_url: PanelApp web URL for this panel.
    """

    panel_id: int
    name: str
    disease_group: Optional[str] = None
    disease_sub_group: Optional[str] = None
    status: str = "public"
    version: Optional[str] = None
    gene_confidence: int = 0
    gene_confidence_colour: str = "None"
    mode_of_inheritance: Optional[str] = None
    phenotypes: list[str] = field(default_factory=list)
    is_diagnostic: bool = False
    panel_url: Optional[str] = None


@dataclass
class GeneEntry:
    """A gene's entry across all PanelApp panels.

    Attributes:
        hgnc_id: HGNC identifier (e.g. "HGNC:1100").
        gene_symbol: HGNC gene symbol.
        ensembl_gene_id: Ensembl gene ID (ENSG…).
        panels: List of Panel objects the gene appears in.
        green_panels: Subset of panels where confidence == Green (3).
        highest_confidence: Maximum confidence level across all panels.
    """

    hgnc_id: Optional[str]
    gene_symbol: str
    ensembl_gene_id: Optional[str] = None
    panels: list[Panel] = field(default_factory=list)
    green_panels: list[Panel] = field(default_factory=list)
    highest_confidence: int = 0


class PanelAppClient:
    """Client for the Genomics England PanelApp REST API.

    Supports gene panel lookup by HGNC symbol, Ensembl gene ID, or
    PanelApp panel ID.

    Args:
        api_base: PanelApp API base URL.
        http_timeout: HTTP request timeout in seconds.
        min_confidence: Minimum confidence level to include in results.
            Default is 2 (Amber and above); set to 3 for Green-only.

    Example:
        >>> client = PanelAppClient()
        >>> panels = await client.get_gene_panels("BRCA1")
        >>> green = [p for p in panels if p.is_diagnostic]
        >>> print(f"BRCA1 in {len(green)} diagnostic panels")
    """

    def __init__(
        self,
        api_base: str = PANELAPP_API_BASE,
        http_timeout: float = 10.0,
        min_confidence: int = 2,
    ) -> None:
        self._api_base = api_base
        self._http_timeout = http_timeout
        self._min_confidence = min_confidence

    async def get_gene_panels(self, gene_symbol: str) -> list[Panel]:
        """Retrieve all PanelApp panels containing a gene.

        Args:
            gene_symbol: HGNC gene symbol (case-insensitive).

        Returns:
            List of Panel objects ordered by confidence (Green first).
            Only panels with confidence >= min_confidence are returned.
        """
        url = f"{self._api_base}/genes/{gene_symbol.upper()}/"
        try:
            raw = await self._get(url)
        except httpx.HTTPError as exc:
            logger.warning("PanelApp lookup failed for %s: %s", gene_symbol, exc)
            return []

        if not raw:
            return []

        return self._parse_gene_panels(raw)

    async def get_panel_by_id(self, panel_id: int) -> Optional[dict]:
        """Retrieve panel metadata by PanelApp panel ID.

        Args:
            panel_id: Numeric PanelApp panel identifier.

        Returns:
            Raw panel metadata dict or None.
        """
        url = f"{self._api_base}/panels/{panel_id}/"
        try:
            return await self._get(url)
        except httpx.HTTPError as exc:
            logger.warning("PanelApp panel %d lookup failed: %s", panel_id, exc)
            return None

    async def get_gene_entry(self, gene_symbol: str) -> Optional[GeneEntry]:
        """Return a full GeneEntry for a gene across all PanelApp panels.

        Args:
            gene_symbol: HGNC gene symbol.

        Returns:
            GeneEntry with panels and green_panels populated, or None.
        """
        panels = await self.get_gene_panels(gene_symbol)
        if not panels:
            return None

        green = [p for p in panels if p.is_diagnostic]
        max_conf = max((p.gene_confidence for p in panels), default=0)

        return GeneEntry(
            hgnc_id=None,  # Populated from API when available
            gene_symbol=gene_symbol,
            panels=panels,
            green_panels=green,
            highest_confidence=max_conf,
        )

    async def search_panels(
        self,
        name: Optional[str] = None,
        disease_group: Optional[str] = None,
    ) -> list[dict]:
        """Search PanelApp panels by name or disease group.

        Args:
            name: Partial panel name for text search.
            disease_group: Disease group filter.

        Returns:
            List of panel summary dicts.
        """
        url = f"{self._api_base}/panels/"
        params: dict[str, str] = {}
        if name:
            params["name"] = name
        if disease_group:
            params["disease_group"] = disease_group

        try:
            result = await self._get(url, params=params)
            return result.get("results", []) if result else []
        except httpx.HTTPError as exc:
            logger.warning("PanelApp panel search failed: %s", exc)
            return []

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    async def _get(self, url: str, params: Optional[dict] = None) -> dict:
        """Execute a GET request against the PanelApp API.

        Args:
            url: Full API URL.
            params: Optional query parameters.

        Returns:
            Parsed JSON response dict.

        Raises:
            httpx.HTTPError: On HTTP errors.
        """
        async with httpx.AsyncClient(timeout=self._http_timeout) as client:
            resp = await client.get(url, params=params or {})
            resp.raise_for_status()
            return resp.json()  # type: ignore[return-value]

    def _parse_gene_panels(self, data: dict) -> list[Panel]:
        """Parse the PanelApp gene API response into Panel objects.

        Args:
            data: Raw JSON from the genes endpoint.

        Returns:
            List of Panel objects filtered by min_confidence.
        """
        panels: list[Panel] = []
        results = data.get("results", [data])  # single-result endpoints return the gene directly

        # Handle the paginated format and the single-gene format
        gene_data = data if "panel" in data else {}
        panel_items = data.get("results", [])

        # Direct gene endpoint returns a single object with panels nested
        # The actual structure: GET /genes/BRCA1/ returns the gene entry
        # with a 'panel' field for each result in a paginated list
        for item in panel_items:
            panel_info = item.get("panel", {})
            if not panel_info:
                continue

            confidence_level = _safe_int(item.get("confidence_level")) or 0
            if confidence_level < self._min_confidence:
                continue

            colour = _CONFIDENCE_COLOURS.get(confidence_level, "None")
            phenotypes = [p for p in item.get("phenotypes", []) if p]

            panel = Panel(
                panel_id=_safe_int(panel_info.get("id")) or 0,
                name=panel_info.get("name", ""),
                disease_group=panel_info.get("disease_group"),
                disease_sub_group=panel_info.get("disease_sub_group"),
                status=panel_info.get("status", "public"),
                version=str(panel_info.get("version", "")),
                gene_confidence=confidence_level,
                gene_confidence_colour=colour,
                mode_of_inheritance=item.get("mode_of_inheritance"),
                phenotypes=phenotypes,
                is_diagnostic=(confidence_level >= MIN_DIAGNOSTIC_CONFIDENCE),
                panel_url=(
                    f"https://panelapp.genomicsengland.co.uk/panels/{panel_info.get('id')}/"
                ),
            )
            panels.append(panel)

        # Sort: Green (3) first, then Amber (2), then Red (1)
        panels.sort(key=lambda p: p.gene_confidence, reverse=True)
        return panels


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _safe_int(value: object) -> Optional[int]:
    """Convert value to int safely.

    Args:
        value: Input value.

    Returns:
        Integer or None.
    """
    if value is None:
        return None
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
