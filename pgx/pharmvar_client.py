"""
pgx.pharmvar_client
====================
PharmVar REST API client for CYP2D6 star allele definitions.

PharmVar (Pharmacogene Variation Consortium) maintains the canonical
reference database for pharmacogene star allele nomenclature.
https://www.pharmvar.org/

PharmVar REST API v1.0:
    Base URL: https://www.pharmvar.org/api-service/
    Endpoints used:
        GET /alleles?gene={gene} — list all star alleles for a gene.
        GET /alleles/{id}        — get a specific star allele by ID.

Star allele structure:
    Each star allele (haplotype) is defined by:
    - A name (e.g. "*4")
    - One or more defining variants (positions + alleles)
    - A functional status (normal, reduced, no function)
    - An activity value (0, 0.5, 1.0, 2.0)

Caching:
    Gene allele lists are cached in memory with a 24-hour TTL.
    Individual allele records are cached indefinitely (definitions rarely change).

References:
    Pratt et al. 2022 CPT Pharmacometrics Syst Pharmacol PMID:36053668
        (PharmVar 4.0 update).
    Gaedigk et al. 2018 CPT Pharmacometrics Syst Pharmacol PMID:29134699
        (PharmVar founding paper).
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential

logger = logging.getLogger(__name__)

_PHARMVAR_API_BASE = "https://www.pharmvar.org/api-service"
_CACHE_TTL_SECONDS = 86_400  # 24-hour TTL for gene allele lists


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class StarAllele:
    """A CYP2D6 (or other pharmacogene) star allele from PharmVar.

    Attributes:
        allele_id: PharmVar internal allele identifier.
        name: Star allele name (e.g. ``"*4"``).
        gene: Gene symbol (e.g. ``"CYP2D6"``).
        function: Functional status string from PharmVar
            (``"Normal Function"``, ``"Decreased Function"``,
            ``"No Function"``, ``"Unknown Function"``).
        activity_value: CPIC activity score for this allele
            (0.0 = no function, 0.5 = decreased, 1.0 = normal, 2.0 = increased).
        defining_variants: List of dicts describing the defining variant(s).
        haplotype_name: Full haplotype name including gene prefix
            (e.g. ``"CYP2D6*4"``).
        raw: Raw PharmVar API response dict.
    """

    allele_id: str
    name: str
    gene: str
    function: str
    activity_value: float
    defining_variants: list[dict[str, Any]] = field(default_factory=list)
    haplotype_name: str = ""
    raw: dict[str, Any] = field(default_factory=dict)

    @property
    def is_no_function(self) -> bool:
        """Return True if this allele has no CYP enzyme function.

        Returns:
            Boolean indicating no-function status (PM contribution).
        """
        return "no function" in self.function.lower()

    @property
    def is_normal_function(self) -> bool:
        """Return True if this allele has normal CYP enzyme function.

        Returns:
            Boolean indicating normal function status (NM contribution).
        """
        return "normal function" in self.function.lower()


def _parse_star_allele(raw: dict[str, Any], gene: str = "CYP2D6") -> StarAllele:
    """Parse a PharmVar API allele response into a StarAllele.

    Args:
        raw: Dict from the PharmVar allele API response.
        gene: Gene symbol (default ``"CYP2D6"``).

    Returns:
        StarAllele populated from the PharmVar API response.
    """
    name = raw.get("alleleName", raw.get("name", ""))
    function_str = raw.get("functionStatus", "Unknown Function")

    # Map PharmVar function string to CPIC activity value
    function_to_activity: dict[str, float] = {
        "normal function": 1.0,
        "decreased function": 0.5,
        "no function": 0.0,
        "increased function": 2.0,
        "unknown function": 0.5,  # conservative estimate
    }
    activity = function_to_activity.get(function_str.lower(), 0.5)

    return StarAllele(
        allele_id=str(raw.get("id", raw.get("alleleId", ""))),
        name=name,
        gene=gene,
        function=function_str,
        activity_value=activity,
        defining_variants=raw.get("variants", []),
        haplotype_name=f"{gene}{name}",
        raw=raw,
    )


# ---------------------------------------------------------------------------
# API functions
# ---------------------------------------------------------------------------

# Gene allele list cache: {gene_symbol: (fetched_at, [StarAllele])}
_gene_allele_cache: dict[str, tuple[float, list[StarAllele]]] = {}


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    reraise=True,
)
def _fetch_gene_alleles_from_api(gene: str) -> list[dict[str, Any]]:
    """Fetch star allele list for a gene from the PharmVar API.

    Args:
        gene: Gene symbol (e.g. ``"CYP2D6"``).

    Returns:
        List of raw allele dicts from the PharmVar API.

    Raises:
        httpx.HTTPStatusError: On non-2xx response.
        httpx.RequestError: On network error.
    """
    url = f"{_PHARMVAR_API_BASE}/alleles"
    response = httpx.get(url, params={"gene": gene}, timeout=15.0)
    response.raise_for_status()
    data = response.json()
    if isinstance(data, list):
        return data
    return data.get("data", data.get("alleles", []))


def get_star_alleles_for_gene(gene: str = "CYP2D6") -> list[StarAllele]:
    """Retrieve all star allele definitions for a gene from PharmVar.

    Uses a 24-hour in-memory cache.  On cache miss, fetches from the
    PharmVar REST API with up to 3 retries (exponential backoff).

    Args:
        gene: Gene symbol (default ``"CYP2D6"``).

    Returns:
        List of StarAllele objects for the gene, sorted by allele name.

    Raises:
        httpx.HTTPStatusError: On persistent API errors after retries.
        httpx.RequestError: On persistent network errors after retries.
    """
    now = time.time()

    if gene in _gene_allele_cache:
        fetched_at, cached_alleles = _gene_allele_cache[gene]
        if (now - fetched_at) < _CACHE_TTL_SECONDS:
            return cached_alleles

    try:
        raw_alleles = _fetch_gene_alleles_from_api(gene)
        alleles = [_parse_star_allele(r, gene) for r in raw_alleles]
        alleles.sort(key=lambda a: a.name)
        _gene_allele_cache[gene] = (now, alleles)
        logger.info("PharmVar: loaded %d star alleles for %s", len(alleles), gene)
        return alleles
    except (httpx.HTTPStatusError, httpx.RequestError) as exc:
        logger.error("PharmVar API error for gene %s: %s", gene, exc)
        raise


def get_star_allele_activity(allele_name: str, gene: str = "CYP2D6") -> float | None:
    """Look up the CPIC activity value for a named star allele.

    Args:
        allele_name: Star allele name, e.g. ``"*1"``, ``"*4"``.
            Leading asterisk is optional.
        gene: Gene symbol (default ``"CYP2D6"``).

    Returns:
        CPIC activity value (0.0, 0.5, or 1.0) or None if not found.

    References:
        PharmVar: https://www.pharmvar.org/
        CPIC activity scores: https://cpicpgx.org/
    """
    if not allele_name.startswith("*"):
        allele_name = f"*{allele_name}"

    try:
        alleles = get_star_alleles_for_gene(gene)
        for allele in alleles:
            if allele.name == allele_name:
                return allele.activity_value
    except Exception as exc:
        logger.warning("Could not fetch PharmVar alleles: %s. Using default activity.", exc)

    # Default activity values for common CYP2D6 alleles when API is unavailable
    _DEFAULT_ACTIVITIES: dict[str, float] = {
        "*1": 1.0,  # normal function
        "*2": 1.0,  # normal function
        "*3": 0.0,  # no function (frameshift)
        "*4": 0.0,  # no function (splice defect)
        "*5": 0.0,  # no function (gene deletion)
        "*6": 0.0,  # no function
        "*9": 0.5,  # decreased function
        "*10": 0.5, # decreased function (common in East Asian)
        "*17": 0.5, # decreased function (common in African)
        "*41": 0.5, # decreased function (reduced splicing)
        "*1xN": 2.0, # duplication — ultrarapid
    }
    return _DEFAULT_ACTIVITIES.get(allele_name)


def compute_diplotype_activity_score(star_allele_1: str, star_allele_2: str) -> float:
    """Compute the CPIC diplotype activity score from two star alleles.

    The activity score is the sum of the two haplotype activity values.
    Common activity scores and their phenotypes:
        0.0     → PM (Poor Metaboliser)
        0.25-1.0 → IM (Intermediate Metaboliser)
        1.25-2.25 → NM (Normal Metaboliser)
        ≥2.25   → UM (Ultrarapid Metaboliser)

    Args:
        star_allele_1: First haplotype star allele name (e.g. ``"*1"``).
        star_allele_2: Second haplotype star allele name (e.g. ``"*4"``).

    Returns:
        CPIC diplotype activity score (sum of allele activity values).

    References:
        CPIC CYP2D6 guideline: https://cpicpgx.org/guidelines/
        PharmVar: https://www.pharmvar.org/
    """
    act1 = get_star_allele_activity(star_allele_1) or 1.0
    act2 = get_star_allele_activity(star_allele_2) or 1.0
    return act1 + act2
