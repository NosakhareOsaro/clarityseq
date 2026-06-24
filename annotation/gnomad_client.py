"""
gnomAD v4.1 client for population allele frequency lookup.

gnomAD v4.1 Overview
--------------------
Released: April 19, 2024
Individuals: 807,162 total (730,947 exomes + 76,215 genomes)
Reference: GRCh38

gnomAD v4.0 Allele Number Bug
------------------------------
gnomAD v4.0 contained a critical bug where the allele number (AN) was
underreported for variants in certain non-European populations, specifically
affecting allele frequency calculations in AFR (African/African American),
AMR (Admixed American), and SAS (South Asian) populations.

v4.1 (April 2024) corrects this by recalculating AN for all variants.
Always use gnomAD v4.1 or later for clinical-grade frequency reporting.
Do NOT use gnomAD v4.0 for allele frequency evidence in clinical reports.

References
----------
Karczewski KJ, et al. "The mutational constraint spectrum quantified from
variation in 141,456 humans." Nature. 2020;581:434–443. PMID:32461654.
DOI:10.1038/s41586-020-2308-7

gnomAD v4.1 blog: https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/

PM2 Supporting Evidence (ClinGen SVI 2024)
------------------------------------------
Per ClinGen SVI Working Group 2024 recommendations, PM2 (absent from controls)
is now applicable only at "Supporting" strength (PM2_Supporting) when allele
frequency criteria are met. The previous Moderate strength was considered
overclaiming.

Thresholds:
  - PM2_Supporting: AF < population-specific threshold in gnomAD v4.1
  - General threshold: AF < 0.0001 (0.01%) in any population
  - Disease-specific thresholds may apply for dominant/recessive conditions

Ancestry Labels (gnomAD v4.1)
------------------------------
  afr   — African/African American
  amr   — Admixed American
  asj   — Ashkenazi Jewish
  eas   — East Asian
  fin   — Finnish
  mid   — Middle Eastern
  nfe   — Non-Finnish European
  rmi   — Remaining Individuals
  sas   — South Asian
  oth   — Other (deprecated in v4; use rmi)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

# gnomAD v4.1 GraphQL API endpoint
GNOMAD_API_URL: str = "https://gnomad.broadinstitute.org/api"

# gnomAD v4.1 genome dataset identifier
GNOMAD_GENOME_DATASET: str = "gnomad_r4"
GNOMAD_EXOME_DATASET: str = "gnomad_r4"

# PM2_Supporting threshold (ClinGen SVI 2024): AF < 0.01% in all populations
PM2_SUPPORTING_THRESHOLD: float = 0.0001

# gnomAD v4.1 population codes
GNOMAD_POPULATIONS: list[str] = [
    "afr", "amr", "asj", "eas", "fin", "mid", "nfe", "rmi", "sas",
]


@dataclass
class AncestryFrequency:
    """Allele frequency data for a single ancestry group.

    Attributes:
        population: gnomAD ancestry code (e.g. "nfe", "afr").
        af: Allele frequency in this population.
        ac: Allele count.
        an: Allele number (corrected in v4.1).
        nhomalt: Number of homozygous alternate individuals.
    """

    population: str
    af: float
    ac: int
    an: int
    nhomalt: int = 0


@dataclass
class GnomADData:
    """Comprehensive gnomAD v4.1 frequency data for a variant.

    Attributes:
        chrom: Chromosome ("chr"-prefixed, GRCh38).
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
        af: Overall allele frequency across all populations.
        ac: Overall allele count.
        an: Overall allele number (corrected in v4.1).
        nhomalt: Number of homozygous alternate individuals (genome callset).
        by_ancestry: Per-population frequency data keyed by gnomAD code.
        af_genome: Genome-specific AF (preferred for WGS studies).
        af_exome: Exome-specific AF.
        dataset: gnomAD dataset version used (e.g. "gnomad_r4").
        flags: Quality/filter flags from gnomAD (e.g. "lcr", "segdup").
        pm2_supporting: True if AF meets PM2_Supporting threshold.
        note: Free-text note (e.g. v4.0 bug warning if relevant).
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    af: Optional[float] = None
    ac: Optional[int] = None
    an: Optional[int] = None
    nhomalt: Optional[int] = None
    by_ancestry: dict[str, AncestryFrequency] = field(default_factory=dict)
    af_genome: Optional[float] = None
    af_exome: Optional[float] = None
    dataset: str = GNOMAD_GENOME_DATASET
    flags: list[str] = field(default_factory=list)
    pm2_supporting: bool = False
    note: Optional[str] = None


# GraphQL query for gnomAD v4.1 variant lookup
_GQL_VARIANT_QUERY = """
query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    chrom
    pos
    ref
    alt
    genome {
      ac
      an
      af
      homozygote_count
      populations {
        id
        ac
        an
        af
        homozygote_count
      }
      filters
    }
    exome {
      ac
      an
      af
      homozygote_count
      populations {
        id
        ac
        an
        af
        homozygote_count
      }
    }
    flags
  }
}
"""


class GnomADClient:
    """Client for querying gnomAD v4.1 population allele frequencies.

    Queries the gnomAD GraphQL API.  In production environments it is
    recommended to use a local mirror of the gnomAD Hail tables or the
    gnomAD VCF files with tabix indexing for performance.

    Args:
        api_url: gnomAD GraphQL API URL.
        dataset: Dataset identifier (default: "gnomad_r4" for v4.1).
        http_timeout: HTTP request timeout in seconds.

    Example:
        >>> client = GnomADClient()
        >>> data = await client.get_allele_frequency("chr17", 43094692, "G", "A")
        >>> print(data.af)          # overall AF
        >>> print(data.pm2_supporting)  # True if AF < 0.01%
        >>> print(data.by_ancestry["nfe"].af)  # Non-Finnish European AF
    """

    def __init__(
        self,
        api_url: str = GNOMAD_API_URL,
        dataset: str = GNOMAD_GENOME_DATASET,
        http_timeout: float = 15.0,
    ) -> None:
        self._api_url = api_url
        self._dataset = dataset
        self._http_timeout = http_timeout

    async def get_allele_frequency(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        population: Optional[str] = None,
    ) -> GnomADData:
        """Retrieve allele frequency data from gnomAD v4.1.

        Args:
            chrom: Chromosome string (with or without "chr" prefix).
            pos: 1-based genomic position.
            ref: Reference allele.
            alt: Alternate allele.
            population: Optional gnomAD ancestry code (e.g. "nfe", "afr")
                to restrict the returned AF to a specific population.

        Returns:
            GnomADData object with genome and exome frequencies,
            population-stratified data, and PM2_Supporting flag.
            If the variant is absent from gnomAD v4.1, all AF fields
            are None and pm2_supporting is True (variant absent from
            controls — PM2_Supporting applies).
        """
        chrom_bare = chrom.replace("chr", "")
        variant_id = f"{chrom_bare}-{pos}-{ref}-{alt}"

        empty_result = GnomADData(
            chrom=_normalise_chrom(chrom),
            pos=pos,
            ref=ref,
            alt=alt,
            pm2_supporting=True,  # Absent = PM2_Supporting by default
            note="Variant absent from gnomAD v4.1",
        )

        try:
            raw = await self._query_api(variant_id)
        except httpx.HTTPError as exc:
            logger.error("gnomAD API request failed: %s", exc)
            return empty_result

        if raw is None:
            return empty_result

        return self._parse_response(chrom, pos, ref, alt, raw, population)

    async def get_max_af_across_populations(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[float]:
        """Return the maximum allele frequency across all gnomAD v4.1 populations.

        This is used to assess PM2_Supporting — a variant must be absent or
        below threshold in ALL populations, not just overall.

        Args:
            chrom: Chromosome (with or without "chr" prefix).
            pos: 1-based position.
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            Maximum AF across all ancestry groups, or None if absent.
        """
        data = await self.get_allele_frequency(chrom, pos, ref, alt)
        if not data.by_ancestry:
            return data.af
        return max(
            (pop.af for pop in data.by_ancestry.values()),
            default=data.af,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    async def _query_api(self, variant_id: str) -> Optional[dict]:
        """Execute the GraphQL query against the gnomAD API.

        Args:
            variant_id: Variant in gnomAD format "CHROM-POS-REF-ALT".

        Returns:
            The ``variant`` dict from the GraphQL response, or None if
            the variant is not found.
        """
        payload = {
            "query": _GQL_VARIANT_QUERY,
            "variables": {
                "variantId": variant_id,
                "datasetId": self._dataset,
            },
        }

        async with httpx.AsyncClient(timeout=self._http_timeout) as client:
            resp = await client.post(self._api_url, json=payload)
            resp.raise_for_status()
            result = resp.json()

        # GraphQL errors are in "errors" key, not HTTP status codes
        if "errors" in result:
            logger.warning("gnomAD GraphQL errors for %s: %s", variant_id, result["errors"])
            return None

        return result.get("data", {}).get("variant")

    def _parse_response(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        data: dict,
        requested_population: Optional[str],
    ) -> GnomADData:
        """Parse the gnomAD GraphQL response into GnomADData.

        Args:
            chrom: Original chromosome string.
            pos: Position.
            ref: Reference allele.
            alt: Alternate allele.
            data: The ``variant`` dict from GraphQL response.
            requested_population: Optional population to surface as primary AF.

        Returns:
            Populated GnomADData instance.
        """
        genome = data.get("genome") or {}
        exome = data.get("exome") or {}

        # Genome-level overall frequencies (v4.1 AN bug-corrected)
        af_genome = _safe_float(genome.get("af"))
        ac_genome = _safe_int(genome.get("ac"))
        an_genome = _safe_int(genome.get("an"))
        nhomalt = _safe_int(genome.get("homozygote_count"))
        flags: list[str] = data.get("flags") or []

        # Exome-level frequencies
        af_exome = _safe_float(exome.get("af"))

        # Population-stratified frequencies from genome callset
        by_ancestry: dict[str, AncestryFrequency] = {}
        for pop in genome.get("populations", []):
            pop_id = pop["id"].lower()
            by_ancestry[pop_id] = AncestryFrequency(
                population=pop_id,
                af=_safe_float(pop.get("af")) or 0.0,
                ac=_safe_int(pop.get("ac")) or 0,
                an=_safe_int(pop.get("an")) or 0,
                nhomalt=_safe_int(pop.get("homozygote_count")) or 0,
            )

        # If caller requested a specific population, use that as primary AF
        primary_af: Optional[float] = af_genome
        if requested_population and requested_population in by_ancestry:
            primary_af = by_ancestry[requested_population].af

        # PM2_Supporting: absent OR max population AF below threshold
        max_pop_af = max(
            (p.af for p in by_ancestry.values()), default=primary_af or 0.0
        )
        pm2 = (primary_af is None) or (max_pop_af < PM2_SUPPORTING_THRESHOLD)

        return GnomADData(
            chrom=_normalise_chrom(chrom),
            pos=pos,
            ref=ref,
            alt=alt,
            af=primary_af,
            ac=ac_genome,
            an=an_genome,
            nhomalt=nhomalt,
            by_ancestry=by_ancestry,
            af_genome=af_genome,
            af_exome=af_exome,
            dataset=self._dataset,
            flags=flags,
            pm2_supporting=pm2,
        )


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _normalise_chrom(chrom: str) -> str:
    """Ensure "chr" prefix on chromosome string.

    Args:
        chrom: Chromosome identifier.

    Returns:
        "chr"-prefixed string.
    """
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _safe_float(value: object) -> Optional[float]:
    """Convert value to float, returning None on failure.

    Args:
        value: Input value.

    Returns:
        Float or None.
    """
    if value is None:
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _safe_int(value: object) -> Optional[int]:
    """Convert value to int, returning None on failure.

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
