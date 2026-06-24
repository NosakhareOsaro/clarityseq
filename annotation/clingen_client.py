"""
ClinGen gene-disease validity client.

ClinGen (Clinical Genome Resource) curates gene-disease relationships using
a standardised framework with defined evidence levels.

References
----------
Strande NT, et al. "Evaluating the Clinical Validity of Gene-Disease
Associations: An Evidence-Based Framework Developed by the Clinical Genome
Resource." American Journal of Human Genetics. 2017;100(6):895–906.
PMID:28552198. DOI:10.1016/j.ajhg.2017.04.015

ClinGen Gene Validity API
--------------------------
Base URL: https://search.clinicalgenome.org/kb/gene-validity
REST API: https://clinicalgenome.org/api/gene-validity/

Evidence Classifications
------------------------
ClinGen uses a semi-quantitative scoring system to classify gene-disease
relationships:

  Definitive      — Overwhelming evidence; gene causative for disease
  Strong          — Strong evidence from multiple independent sources
  Moderate        — Moderate evidence; likely disease-causing
  Limited         — Some evidence; insufficient for diagnostic use
  Disputed        — Conflicting evidence in literature
  Refuted         — Evidence does not support causation
  Animal Model    — Evidence only from animal models

For clinical variant interpretation, only Definitive and Strong are
generally considered sufficient for full ACMG criteria application.
Moderate may support variant interpretation with caveats.

ACMG/AMP Integration
--------------------
Gene-disease validity affects several ACMG rules:
- PP2: Missense variant in a gene with low missense tolerance AND where
  missense variants cause disease (requires Definitive/Strong validity)
- BP1: Missense variant in a gene where only truncating variants cause
  disease (requires Definitive/Strong loss-of-function disease mechanism)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

# ClinGen Gene Validity Classification API
_CLINGEN_API_BASE = "https://clinicalgenome.org/api"
_CLINGEN_GV_ENDPOINT = f"{_CLINGEN_API_BASE}/gene-validity"

# ClinGen SVI evidence classification order (highest to lowest)
_VALIDITY_ORDER: list[str] = [
    "Definitive",
    "Strong",
    "Moderate",
    "Limited",
    "Disputed",
    "Refuted",
    "No Known Disease Relationship",
    "Animal Model Only",
]

_DIAGNOSTIC_VALIDITY_LEVELS: frozenset[str] = frozenset({"Definitive", "Strong"})
_SUPPORTIVE_VALIDITY_LEVELS: frozenset[str] = frozenset({"Moderate"})


@dataclass
class GeneValidityClassification:
    """A single ClinGen gene-disease validity classification.

    Attributes:
        disease_label: MeSH/OMIM disease name.
        disease_id: Disease ontology identifier (e.g. "OMIM:113705").
        moi: Mode of inheritance (e.g. "Autosomal dominant").
        classification: ClinGen evidence classification level.
        classification_date: ISO date of classification.
        sop_version: ClinGen SOP version used for curation.
        gcep_id: Gene Curation Expert Panel identifier.
        gcep_name: Gene Curation Expert Panel name.
        report_url: URL to the full ClinGen curation report.
        is_diagnostic: True for Definitive or Strong classifications.
    """

    disease_label: str
    disease_id: Optional[str] = None
    moi: Optional[str] = None
    classification: str = "Limited"
    classification_date: Optional[str] = None
    sop_version: Optional[str] = None
    gcep_id: Optional[str] = None
    gcep_name: Optional[str] = None
    report_url: Optional[str] = None
    is_diagnostic: bool = False


@dataclass
class GeneValidity:
    """Aggregated ClinGen gene-disease validity for a gene.

    Attributes:
        hgnc_id: HGNC gene identifier (e.g. "HGNC:1100").
        gene_symbol: HGNC gene symbol (e.g. "BRCA1").
        classifications: All gene-disease validity curations for this gene.
        highest_classification: The highest confidence classification found.
        definitive_diseases: Diseases with Definitive classification.
        strong_diseases: Diseases with Strong classification.
        has_lof_disease_mechanism: True if any disease association is primarily
            driven by loss-of-function mechanism.
        has_missense_disease_mechanism: True if any disease is caused by
            gain-of-function or missense-specific mechanism.
    """

    hgnc_id: Optional[str]
    gene_symbol: str
    classifications: list[GeneValidityClassification] = field(default_factory=list)
    highest_classification: Optional[str] = None
    definitive_diseases: list[str] = field(default_factory=list)
    strong_diseases: list[str] = field(default_factory=list)
    has_lof_disease_mechanism: bool = False
    has_missense_disease_mechanism: bool = False


class ClinGenClient:
    """Client for querying ClinGen gene-disease validity curations.

    Args:
        api_base: ClinGen API base URL.
        http_timeout: HTTP timeout in seconds.

    Example:
        >>> client = ClinGenClient()
        >>> validity = await client.get_gene_validity("BRCA1")
        >>> if validity:
        ...     print(validity.highest_classification)
        ...     print(validity.definitive_diseases)
    """

    def __init__(
        self,
        api_base: str = _CLINGEN_API_BASE,
        http_timeout: float = 10.0,
    ) -> None:
        self._api_base = api_base
        self._http_timeout = http_timeout

    async def get_gene_validity(self, gene_symbol: str) -> Optional[GeneValidity]:
        """Retrieve ClinGen gene-disease validity curations for a gene.

        Args:
            gene_symbol: HGNC gene symbol (case-insensitive).

        Returns:
            GeneValidity object with all curations, or None if no curations
            exist or the API request fails.
        """
        url = f"{_CLINGEN_GV_ENDPOINT}"
        params = {"gene": gene_symbol.upper(), "format": "json"}

        try:
            async with httpx.AsyncClient(timeout=self._http_timeout) as client:
                resp = await client.get(url, params=params)
                resp.raise_for_status()
                data = resp.json()
        except httpx.HTTPError as exc:
            logger.warning("ClinGen validity lookup failed for %s: %s", gene_symbol, exc)
            return None

        return self._parse_response(gene_symbol, data)

    async def get_gene_validity_by_hgnc(self, hgnc_id: str) -> Optional[GeneValidity]:
        """Retrieve gene-disease validity by HGNC ID.

        Args:
            hgnc_id: HGNC identifier, e.g. "HGNC:1100".

        Returns:
            GeneValidity or None.
        """
        url = f"{_CLINGEN_GV_ENDPOINT}"
        params = {"hgnc_id": hgnc_id, "format": "json"}

        try:
            async with httpx.AsyncClient(timeout=self._http_timeout) as client:
                resp = await client.get(url, params=params)
                resp.raise_for_status()
                data = resp.json()
        except httpx.HTTPError as exc:
            logger.warning("ClinGen HGNC lookup failed for %s: %s", hgnc_id, exc)
            return None

        gene_symbol = data.get("gene", {}).get("symbol", hgnc_id)
        return self._parse_response(gene_symbol, data)

    def is_valid_disease_gene(
        self,
        validity: Optional[GeneValidity],
        min_level: str = "Strong",
    ) -> bool:
        """Check if a gene meets a minimum ClinGen validity threshold.

        Args:
            validity: GeneValidity object from get_gene_validity().
            min_level: Minimum ClinGen classification to accept.
                One of: "Definitive", "Strong", "Moderate", "Limited".

        Returns:
            True if the gene's highest classification is at or above
            the specified minimum level.
        """
        if validity is None or validity.highest_classification is None:
            return False

        min_idx = _validity_index(min_level)
        gene_idx = _validity_index(validity.highest_classification)

        # Lower index = higher evidence in _VALIDITY_ORDER
        return gene_idx <= min_idx

    # ------------------------------------------------------------------
    # Internal parsing
    # ------------------------------------------------------------------

    def _parse_response(self, gene_symbol: str, data: dict) -> Optional[GeneValidity]:
        """Parse ClinGen API response into a GeneValidity object.

        Args:
            gene_symbol: Gene symbol for the result.
            data: Raw JSON response dict.

        Returns:
            GeneValidity or None if no classifications found.
        """
        raw_results = data.get("results", data.get("entities", []))
        if not raw_results:
            return None

        classifications: list[GeneValidityClassification] = []
        hgnc_id: Optional[str] = None

        for item in raw_results:
            gene_info = item.get("gene", {})
            if not hgnc_id:
                hgnc_id = gene_info.get("hgnc_id")

            disease = item.get("disease", {})
            clsf = item.get("classification", {})
            clsf_label = clsf.get("label", "Limited")

            classification = GeneValidityClassification(
                disease_label=disease.get("label", "Unknown"),
                disease_id=disease.get("curie"),
                moi=item.get("mode_of_inheritance", {}).get("label"),
                classification=clsf_label,
                classification_date=item.get("date"),
                sop_version=item.get("sopVersion"),
                gcep_id=item.get("affiliation", {}).get("id"),
                gcep_name=item.get("affiliation", {}).get("name"),
                report_url=item.get("report_url"),
                is_diagnostic=clsf_label in _DIAGNOSTIC_VALIDITY_LEVELS,
            )
            classifications.append(classification)

        if not classifications:
            return None

        # Determine highest classification level
        valid_labels = [c.classification for c in classifications]
        highest = _highest_classification(valid_labels)

        definitive = [
            c.disease_label for c in classifications
            if c.classification == "Definitive"
        ]
        strong = [
            c.disease_label for c in classifications
            if c.classification == "Strong"
        ]

        # Infer disease mechanism from MOI and disease names (heuristic)
        has_lof = any(
            c.moi in ("Autosomal dominant", "X-linked") for c in classifications
        )
        has_missense = any(
            "gain" in (c.disease_label or "").lower()
            or "activating" in (c.disease_label or "").lower()
            for c in classifications
        )

        return GeneValidity(
            hgnc_id=hgnc_id,
            gene_symbol=gene_symbol,
            classifications=classifications,
            highest_classification=highest,
            definitive_diseases=definitive,
            strong_diseases=strong,
            has_lof_disease_mechanism=has_lof,
            has_missense_disease_mechanism=has_missense,
        )


# ---------------------------------------------------------------------------
# Private utilities
# ---------------------------------------------------------------------------


def _validity_index(level: str) -> int:
    """Return the index of a classification level in _VALIDITY_ORDER.

    Args:
        level: Classification level string.

    Returns:
        Index (lower = stronger evidence), or len(_VALIDITY_ORDER) if unknown.
    """
    try:
        return _VALIDITY_ORDER.index(level)
    except ValueError:
        return len(_VALIDITY_ORDER)


def _highest_classification(levels: list[str]) -> Optional[str]:
    """Return the highest-evidence ClinGen classification from a list.

    Args:
        levels: List of classification level strings.

    Returns:
        The highest-evidence classification, or None if the list is empty.
    """
    if not levels:
        return None
    return min(levels, key=_validity_index)
