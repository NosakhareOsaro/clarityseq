"""
bayesacmg.vcep_client
=====================

ClinGen CSpec registry client for gene-specific VCEP rule specifications.

Why VCEP specifications are essential
--------------------------------------
The ACMG/AMP 28-criterion framework (Richards 2015 PMID:25741868) provides
a general foundation, but gene-specific expert panels (VCEPs — Variant
Curation Expert Panels) issue refined specifications that may OVERRIDE the
general framework.  Examples:

- BRCA1/2 (ENIGMA VCEP): specific AF thresholds, evidence weightings.
- RASopathy genes (KCNQ1, etc.): modified PP3/BP4 thresholds.
- TP53: specific functional assay interpretations.

The ClinGen CSpec registry (https://cspec.genome.network/) provides a
machine-readable API to retrieve these specifications.  Callers MUST query
this registry before applying any rule that VCEPs commonly override (PM2,
PP3/BP4 thresholds, PS3/BS3 criteria).

CSpec registry API
-------------------
Base URL: https://cspec.genome.network/cspec/api/svi/
Per-gene:  https://cspec.genome.network/cspec/api/svi/?gene=<HGNC_SYMBOL>
Response:  JSON array of VCEP specification objects.

Caching strategy
-----------------
Responses are cached in memory with a 24-hour TTL per gene symbol.
On cache miss the API is called with up to 3 retries using exponential
backoff (tenacity).  The cache is not persisted to disk; restart clears it.

PM2 override
-------------
Some VCEPs allow PM2 at Moderate weight for very specific gene-disease pairs
(e.g. PTEN, RYR2).  If a VCEP specification sets ``pm2_weight`` to
``"moderate"``, the caller should override the default Supporting weight
returned by rule_pm2() in pathogenic.py.

References:
    ClinGen CSpec registry: https://cspec.genome.network/cspec/
    Richards et al. 2015 PMID:25741868
    ClinGen SVI Working Group 2024
    ACGS 2024 v1.2 §5.
"""

from __future__ import annotations

import asyncio
import time
from dataclasses import dataclass, field
from typing import Any

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CSPEC_API_BASE = "https://cspec.genome.network/cspec/api/svi/"
_CACHE_TTL_SECONDS = 86_400  # 24-hour TTL per ClinGen CSpec caching guidance
_MAX_RETRIES = 3
_RETRY_WAIT_MIN = 1    # seconds
_RETRY_WAIT_MAX = 10   # seconds


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

@dataclass
class VCEPSpec:
    """Gene-specific VCEP specification retrieved from ClinGen CSpec registry.

    Attributes:
        gene_symbol: HGNC gene symbol.
        vcep_name: Name of the Variant Curation Expert Panel.
        vcep_id: Unique VCEP identifier from the CSpec registry.
        pm2_weight: Override for PM2 weight; ``"supporting"`` (default) or
            ``"moderate"`` (gene-specific VCEP override).
        pp3_threshold_alphamissense: Override for AlphaMissense PP3 threshold.
            None uses the default 0.564 from ClinGen SVI 2024.
        bp4_threshold_alphamissense: Override for AlphaMissense BP4 threshold.
            None uses the default 0.340 from ClinGen SVI 2024.
        custom_thresholds: Dict of rule_id → custom threshold for other rules.
        raw_spec: Raw JSON spec from the CSpec registry API.
        retrieved_at: Unix timestamp when this spec was retrieved.
    """

    gene_symbol: str
    vcep_name: str = ""
    vcep_id: str = ""
    pm2_weight: str = "supporting"              # default ClinGen SVI 2024
    pp3_threshold_alphamissense: float | None = None   # None → use 0.564
    bp4_threshold_alphamissense: float | None = None   # None → use 0.340
    custom_thresholds: dict[str, Any] = field(default_factory=dict)
    raw_spec: dict[str, Any] = field(default_factory=dict)
    retrieved_at: float = 0.0

    @property
    def is_expired(self) -> bool:
        """Return True if this cached spec has exceeded the 24-hour TTL.

        Returns:
            Boolean indicating cache staleness.
        """
        return (time.time() - self.retrieved_at) > _CACHE_TTL_SECONDS


# ---------------------------------------------------------------------------
# In-memory cache
# ---------------------------------------------------------------------------

_spec_cache: dict[str, VCEPSpec] = {}


def _parse_spec(gene_symbol: str, raw: list[dict[str, Any]]) -> VCEPSpec:
    """Parse CSpec registry API response into a VCEPSpec.

    Args:
        gene_symbol: The queried HGNC gene symbol.
        raw: Parsed JSON list of specification objects from the CSpec API.

    Returns:
        VCEPSpec populated from the first matching result, or a default
        spec (pm2_weight="supporting") if no VCEP specification exists.
    """
    if not raw:
        return VCEPSpec(
            gene_symbol=gene_symbol,
            vcep_name="",
            retrieved_at=time.time(),
        )

    spec_data = raw[0]  # take first result; typically one VCEP per gene
    pm2_weight = "supporting"

    # Parse PM2 weight override if present in spec
    for criterion in spec_data.get("criteria", []):
        crit_id = criterion.get("id", "").upper()
        if crit_id == "PM2":
            strength = criterion.get("strength", "").lower()
            if "moderate" in strength:
                pm2_weight = "moderate"   # VCEP override: PM2 at Moderate
            break

    return VCEPSpec(
        gene_symbol=gene_symbol,
        vcep_name=spec_data.get("vcep_name", ""),
        vcep_id=str(spec_data.get("vcep_id", "")),
        pm2_weight=pm2_weight,
        raw_spec=spec_data,
        retrieved_at=time.time(),
    )


@retry(
    stop=stop_after_attempt(_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=_RETRY_WAIT_MIN, max=_RETRY_WAIT_MAX),
    reraise=True,
)
async def _fetch_spec_from_api(gene_symbol: str, client: httpx.AsyncClient) -> list[dict[str, Any]]:
    """Fetch VCEP spec from ClinGen CSpec registry API.

    Applies exponential backoff retry (up to 3 attempts) via tenacity.

    Args:
        gene_symbol: HGNC gene symbol to query.
        client: Active httpx.AsyncClient.

    Returns:
        Parsed JSON list from the CSpec API response.

    Raises:
        httpx.HTTPStatusError: On 4xx/5xx responses after retries exhausted.
        httpx.RequestError: On network-level errors after retries exhausted.
    """
    url = f"{_CSPEC_API_BASE}?gene={gene_symbol}"
    response = await client.get(url, timeout=10.0)
    response.raise_for_status()
    data = response.json()
    # CSpec API returns either a list or a dict with a "results" key
    if isinstance(data, list):
        return data
    return data.get("results", [])


async def get_vcep_spec(
    gene_symbol: str,
    client: httpx.AsyncClient | None = None,
) -> VCEPSpec:
    """Retrieve VCEP specification for a gene from the ClinGen CSpec registry.

    Uses an in-memory cache with 24-hour TTL.  On cache miss (or expiry),
    queries https://cspec.genome.network/cspec/api/svi/?gene=<symbol>.
    On API error, returns a default VCEPSpec with standard thresholds.

    Args:
        gene_symbol: HGNC gene symbol to query (e.g. ``"BRCA1"``).
        client: Optional httpx.AsyncClient.  If None, a transient client is
            created internally (slightly less efficient; prefer passing one
            from the caller's context).

    Returns:
        VCEPSpec for the gene.  If no VCEP specification exists in the CSpec
        registry, returns a default VCEPSpec with pm2_weight="supporting".

    References:
        ClinGen CSpec registry: https://cspec.genome.network/cspec/
        ClinGen SVI Working Group 2024.

    Examples:
        >>> import asyncio, httpx
        >>> async def example():
        ...     async with httpx.AsyncClient() as c:
        ...         spec = await get_vcep_spec("BRCA1", client=c)
        ...     return spec.pm2_weight
        >>> asyncio.run(example())
        'supporting'
    """
    symbol = gene_symbol.upper()

    # Cache hit
    if symbol in _spec_cache and not _spec_cache[symbol].is_expired:
        return _spec_cache[symbol]

    # Cache miss — call the API
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient()

    try:
        raw = await _fetch_spec_from_api(symbol, client)
        spec = _parse_spec(symbol, raw)
    except (httpx.HTTPStatusError, httpx.RequestError, Exception):
        # On any error, return default spec (fail-safe: don't block classification)
        spec = VCEPSpec(
            gene_symbol=symbol,
            vcep_name="",
            retrieved_at=time.time(),
        )
    finally:
        if own_client:
            await client.aclose()

    _spec_cache[symbol] = spec
    return spec


def get_vcep_spec_sync(gene_symbol: str) -> VCEPSpec:
    """Synchronous wrapper around get_vcep_spec for use in non-async contexts.

    Uses asyncio.run() — do not call from within an already-running event loop.

    Args:
        gene_symbol: HGNC gene symbol.

    Returns:
        VCEPSpec for the gene.

    References:
        ClinGen CSpec registry: https://cspec.genome.network/cspec/
    """
    return asyncio.run(get_vcep_spec(gene_symbol))


def clear_cache() -> None:
    """Clear the in-memory VCEP spec cache.

    Useful in tests to prevent cache contamination between test cases.

    Returns:
        None.
    """
    _spec_cache.clear()
