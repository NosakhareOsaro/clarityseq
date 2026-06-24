"""
pgx.cpic_client
================
CPIC REST API client for drug dosing recommendations.

CPIC (Clinical Pharmacogenomics Implementation Consortium) provides
evidence-based dosing guidelines for drug-gene pairs.
https://cpicpgx.org/

CPIC API v1.0:
    Base URL: https://api.cpicpgx.org/v1/
    Endpoints used:
        GET /recommendation?gene={gene}&diplotype={diplotype}
            — Get dosing recommendations for a gene/diplotype pair.
        GET /guideline?genesymbol={gene}
            — List all guidelines covering a gene.
        GET /pair?genesymbol={gene}
            — Drug-gene pair definitions.

CYP2D6 CPIC guidelines cover (as of 2024):
    - Codeine (PMID:30447227)
    - Tramadol (PMID:30447227)
    - Oxycodone (PMID:30447227)
    - Tamoxifen (PMID:29385237)
    - Tricyclic antidepressants (PMID:25974703)
    - SSRIs: fluvoxamine, paroxetine (PMID:25974703)
    - Atomoxetine (PMID:30289176)
    - Ondansetron, tropisetron (PMID:32148910)

References:
    Relling et al. 2020 CPT PMID:32779747 (CPIC overview).
    CPIC level A guidelines: mandatory clinical use.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Any

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential

logger = logging.getLogger(__name__)

_CPIC_API_BASE = "https://api.cpicpgx.org/v1"
_CACHE_TTL_SECONDS = 86_400  # 24-hour TTL


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class DrugRecommendation:
    """CPIC drug dosing recommendation for a pharmacogene diplotype.

    Attributes:
        drug_name: Drug name (e.g. ``"codeine"``).
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        diplotype: Diplotype string (e.g. ``"*1/*4"``).
        phenotype: Metaboliser phenotype (NM/IM/PM/UM).
        classification: CPIC recommendation classification
            (``"Use label recommended dosage"`` etc.).
        implications: Clinical implications text.
        recommendation: Specific dosing recommendation text.
        cpic_level: CPIC evidence level (A, B, C, D).
        guideline_name: CPIC guideline name.
        guideline_url: URL to full CPIC guideline.
        citations: Literature citations for this recommendation.
    """

    drug_name: str
    gene: str
    diplotype: str
    phenotype: str
    classification: str
    implications: str
    recommendation: str
    cpic_level: str = "A"
    guideline_name: str = ""
    guideline_url: str = ""
    citations: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

# Cache: {(gene, diplotype): (fetched_at, [DrugRecommendation])}
_recommendation_cache: dict[tuple[str, str], tuple[float, list[DrugRecommendation]]] = {}


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    reraise=True,
)
def _fetch_recommendations_from_api(gene: str, diplotype: str) -> list[dict[str, Any]]:
    """Fetch dosing recommendations from the CPIC REST API.

    Args:
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        diplotype: Diplotype string (e.g. ``"*1/*4"``).

    Returns:
        List of raw recommendation dicts from the CPIC API.

    Raises:
        httpx.HTTPStatusError: On non-2xx response.
        httpx.RequestError: On network error.
    """
    url = f"{_CPIC_API_BASE}/recommendation"
    response = httpx.get(
        url,
        params={"gene": gene, "diplotype": diplotype},
        timeout=15.0,
    )
    response.raise_for_status()
    data = response.json()
    if isinstance(data, list):
        return data
    return data.get("data", [])


def _parse_recommendation(raw: dict[str, Any], gene: str, diplotype: str) -> DrugRecommendation:
    """Parse a CPIC API recommendation response dict.

    Args:
        raw: Raw CPIC API recommendation dict.
        gene: Gene symbol.
        diplotype: Diplotype string.

    Returns:
        DrugRecommendation populated from CPIC API response.
    """
    return DrugRecommendation(
        drug_name=raw.get("drugName", raw.get("drug", {}).get("name", "")),
        gene=gene,
        diplotype=diplotype,
        phenotype=raw.get("phenotype", raw.get("phenotypes", {}).get(gene, "")),
        classification=raw.get("classification", ""),
        implications=raw.get("implications", raw.get("implications_text", "")),
        recommendation=raw.get("recommendation", raw.get("recommendation_text", "")),
        cpic_level=raw.get("cpicLevel", "A"),
        guideline_name=raw.get("guidelineName", ""),
        guideline_url=raw.get("url", f"https://cpicpgx.org/guidelines/"),
        citations=raw.get("citations", []),
    )


def get_recommendations(gene: str, diplotype: str) -> list[DrugRecommendation]:
    """Retrieve all CPIC dosing recommendations for a gene/diplotype pair.

    Uses 24-hour in-memory cache.  Falls back to built-in CYP2D6 recommendations
    if the CPIC API is unavailable.

    Args:
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        diplotype: Diplotype string (e.g. ``"*1/*4"``).

    Returns:
        List of DrugRecommendation objects for all drugs with CPIC guidelines
        for this gene/diplotype combination.

    References:
        CPIC API: https://api.cpicpgx.org/
        CYP2D6 guideline: https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/
    """
    cache_key = (gene.upper(), diplotype)
    now = time.time()

    if cache_key in _recommendation_cache:
        fetched_at, cached = _recommendation_cache[cache_key]
        if (now - fetched_at) < _CACHE_TTL_SECONDS:
            return cached

    try:
        raw_recs = _fetch_recommendations_from_api(gene, diplotype)
        recs = [_parse_recommendation(r, gene, diplotype) for r in raw_recs]
        _recommendation_cache[cache_key] = (now, recs)
        logger.info("CPIC: loaded %d recommendations for %s %s", len(recs), gene, diplotype)
        return recs
    except (httpx.HTTPStatusError, httpx.RequestError) as exc:
        logger.warning("CPIC API unavailable (%s); using built-in recommendations.", exc)
        return _get_builtin_cyp2d6_recommendations(diplotype)


def _get_builtin_cyp2d6_recommendations(diplotype: str) -> list[DrugRecommendation]:
    """Return built-in CYP2D6 CPIC recommendations when API is unavailable.

    Covers the most clinically important CYP2D6 drug interactions.
    Based on CPIC CYP2D6 guideline (cpicpgx.org, accessed 2024).

    Args:
        diplotype: CYP2D6 diplotype string (e.g. ``"*1/*4"``).

    Returns:
        List of DrugRecommendation objects for common CYP2D6 drugs.
    """
    # Map diplotype to phenotype for built-in recommendations
    from pgx.cyrius_runner import classify_phenotype

    # Compute rough activity score from diplotype for phenotype classification
    _NO_FUNCTION = {"*3", "*4", "*5", "*6"}
    _DECREASED = {"*9", "*10", "*17", "*41"}

    alleles = diplotype.replace("*", "").split("/") if "/" in diplotype else [diplotype.replace("*", "")]
    score = 0.0
    for a in alleles[:2]:
        star = f"*{a}"
        if star in _NO_FUNCTION:
            score += 0.0
        elif star in _DECREASED:
            score += 0.5
        else:
            score += 1.0

    phenotype = classify_phenotype(score)

    _CODEINE_RECS: dict[str, str] = {
        "NM": "Use label recommended codeine dosage.",
        "IM": "Use label recommended codeine dosage. Monitor for reduced efficacy.",
        "PM": "Avoid codeine. Alternative non-opioid analgesics recommended (CPIC Level A).",
        "UM": "Avoid codeine. Risk of serious toxicity due to rapid morphine accumulation (CPIC Level A).",
    }

    _TRAMADOL_RECS: dict[str, str] = {
        "NM": "Use label recommended tramadol dosage.",
        "IM": "Use lower tramadol doses. Monitor for reduced analgesia.",
        "PM": "Avoid tramadol. Ineffective; no conversion to active O-desmethyltramadol.",
        "UM": "Avoid tramadol. Risk of serious CNS toxicity.",
    }

    return [
        DrugRecommendation(
            drug_name="codeine",
            gene="CYP2D6",
            diplotype=diplotype,
            phenotype=phenotype,
            classification=_CODEINE_RECS.get(phenotype, "Consult CPIC guidelines."),
            implications=(
                f"CYP2D6 {phenotype}: codeine is a prodrug requiring CYP2D6 for "
                "conversion to morphine."
            ),
            recommendation=_CODEINE_RECS.get(phenotype, "Consult CPIC guidelines."),
            cpic_level="A",
            guideline_name="CPIC Guideline for Codeine and CYP2D6",
            guideline_url="https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
            citations=["Crews et al. 2014 PMID:24458010", "CPIC CYP2D6 guideline (2020 update)"],
        ),
        DrugRecommendation(
            drug_name="tramadol",
            gene="CYP2D6",
            diplotype=diplotype,
            phenotype=phenotype,
            classification=_TRAMADOL_RECS.get(phenotype, "Consult CPIC guidelines."),
            implications=(
                f"CYP2D6 {phenotype}: tramadol is a prodrug requiring CYP2D6 for "
                "conversion to active O-desmethyltramadol."
            ),
            recommendation=_TRAMADOL_RECS.get(phenotype, "Consult CPIC guidelines."),
            cpic_level="A",
            guideline_name="CPIC Guideline for Codeine and CYP2D6",
            guideline_url="https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
            citations=["Crews et al. 2014 PMID:24458010"],
        ),
    ]
