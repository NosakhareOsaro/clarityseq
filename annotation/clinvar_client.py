"""
ClinVar clinical significance lookup client.

ClinVar is NCBI's public archive of reports of relationships among human
variations and phenotypes, with supporting evidence.

References
----------
Landrum MJ, et al. "ClinVar: public archive of interpretations of clinically
relevant variants." Nucleic Acids Research. 2016;44(D1):D862–D868.
PMID:26582918. DOI:10.1093/nar/gkv1222

API
---
ClinVar E-utilities API: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
ClinVar VCV/RCV REST API: https://api.ncbi.nlm.nih.gov/variation/v0/

Star Ratings (review status)
-----------------------------
0 stars — no assertion criteria provided
1 star  — criteria provided, single submitter
2 stars — criteria provided, multiple submitters, no conflicts
3 stars — reviewed by expert panel
4 stars — practice guideline

For clinical reporting, variants with ≥2 stars (expert panel review or
multiple concordant submitters) should be cited as PP5 (Supporting) or
PS4 evidence per the ACMG/AMP framework.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import httpx

logger = logging.getLogger(__name__)

# NCBI E-utilities base URL
_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ClinVar Variation API (v0)
_CLINVAR_VARIATION_API = "https://api.ncbi.nlm.nih.gov/variation/v0"

# Map star count to review status descriptions
_STAR_TO_STATUS: dict[int, str] = {
    0: "no assertion criteria provided",
    1: "criteria provided, single submitter",
    2: "criteria provided, multiple submitters, no conflicts",
    3: "reviewed by expert panel",
    4: "practice guideline",
}

# ClinVar clinical significance strings → normalised category
_SIGNIFICANCE_MAP: dict[str, str] = {
    "Pathogenic": "Pathogenic",
    "Likely pathogenic": "Likely pathogenic",
    "Pathogenic/Likely pathogenic": "Pathogenic/Likely pathogenic",
    "Uncertain significance": "Uncertain significance",
    "Likely benign": "Likely benign",
    "Benign": "Benign",
    "Benign/Likely benign": "Benign/Likely benign",
    "Conflicting classifications of pathogenicity": "Conflicting",
    "not provided": "Not provided",
    "other": "Other",
    "drug response": "Drug response",
    "risk factor": "Risk factor",
    "association": "Association",
    "protective": "Protective",
}


@dataclass
class ClinVarSubmission:
    """A single ClinVar submission from one organisation.

    Attributes:
        submitter: Submitter organisation name.
        classification: Clinical significance reported.
        date_last_evaluated: Date of last evaluation (ISO format).
        review_status: Submission review status string.
        condition: Condition/disease name.
        method: Submission method (e.g. "clinical testing").
    """

    submitter: str
    classification: str
    date_last_evaluated: str | None = None
    review_status: str | None = None
    condition: str | None = None
    method: str | None = None


@dataclass
class ClinVarData:
    """Aggregated ClinVar data for a genomic variant.

    Attributes:
        variation_id: ClinVar Variation ID (VCV accession number).
        rcv_accession: Primary RCV accession (e.g. "RCV000031282").
        star_rating: 0–4 star aggregate review status rating.
        review_status: Human-readable aggregate review status.
        classification: Aggregate clinical significance.
        condition_names: List of associated conditions/diseases.
        submissions: Individual submitter records.
        last_updated: Date of last ClinVar update (ISO format).
        acmg_evidence: Suggested ACMG evidence code (PP5/BP6 based on stars).
    """

    variation_id: str | None = None
    rcv_accession: str | None = None
    star_rating: int = 0
    review_status: str = _STAR_TO_STATUS[0]
    classification: str | None = None
    condition_names: list[str] = field(default_factory=list)
    submissions: list[ClinVarSubmission] = field(default_factory=list)
    last_updated: str | None = None
    acmg_evidence: str | None = None


class ClinVarClient:
    """Client for querying ClinVar via the NCBI variation API.

    Supports lookup by genomic position (GRCh38) using the SPDI format
    (Sequence:Position:Deletion:Insertion) or by ClinVar variation ID.

    Args:
        api_key: NCBI API key for higher rate limits (10 req/s vs 3 req/s).
            Register at: https://www.ncbi.nlm.nih.gov/account/
        http_timeout: HTTP timeout in seconds.

    Example:
        >>> client = ClinVarClient(api_key="your_key_here")
        >>> data = await client.get_clinvar_data("chr17", 43094692, "G", "A")
        >>> if data:
        ...     print(data.classification, data.star_rating)
    """

    def __init__(
        self,
        api_key: str | None = None,
        http_timeout: float = 10.0,
    ) -> None:
        self._api_key = api_key
        self._http_timeout = http_timeout

    async def get_clinvar_data(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> ClinVarData | None:
        """Retrieve ClinVar clinical significance for a variant.

        Queries ClinVar using SPDI notation on the GRCh38 reference sequence.

        Args:
            chrom: Chromosome (with or without "chr" prefix).
            pos: 1-based genomic position.
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            ClinVarData if the variant has ClinVar entries, otherwise None.
        """
        chrom_bare = chrom.replace("chr", "")
        # SPDI format uses 0-based position and RefSeq accessions
        # We use the genomic lookup endpoint instead
        spdi = f"NC_{_chrom_to_refseq(chrom_bare)}:{pos - 1}:{ref}:{alt}"

        try:
            variation_data = await self._lookup_by_spdi(spdi)
            if variation_data is None:
                return None
            return self._parse_variation(variation_data)
        except httpx.HTTPError as exc:
            logger.warning("ClinVar lookup failed for %s:%d: %s", chrom, pos, exc)
            return None

    async def get_by_variation_id(self, variation_id: str) -> ClinVarData | None:
        """Retrieve ClinVar data by variation ID (VCV accession).

        Args:
            variation_id: ClinVar variation ID (numeric or VCV-prefixed).

        Returns:
            ClinVarData or None if not found.
        """
        vid = (
            variation_id.replace("VCV", "").lstrip("0")
            if "VCV" in variation_id
            else variation_id
        )
        url = f"{_CLINVAR_VARIATION_API}/variation/{vid}/clinically_significant_alleles"
        params: dict[str, str] = {}
        if self._api_key:
            params["api_key"] = self._api_key

        try:
            async with httpx.AsyncClient(timeout=self._http_timeout) as client:
                resp = await client.get(url, params=params)
                resp.raise_for_status()
                data = resp.json()
                return self._parse_variation(data)
        except httpx.HTTPError as exc:
            logger.warning("ClinVar VID lookup failed for %s: %s", variation_id, exc)
            return None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    async def _lookup_by_spdi(self, spdi: str) -> dict | None:
        """Look up a variant by SPDI notation.

        Args:
            spdi: SPDI-format variant string.

        Returns:
            Parsed JSON response dict or None.
        """
        url = f"{_CLINVAR_VARIATION_API}/spdi/{spdi}/clinically_significant_alleles"
        params: dict[str, str] = {}
        if self._api_key:
            params["api_key"] = self._api_key

        async with httpx.AsyncClient(timeout=self._http_timeout) as client:
            resp = await client.get(url, params=params)
            if resp.status_code == 404:
                return None  # Variant not in ClinVar
            resp.raise_for_status()
            return resp.json()  # type: ignore[return-value]

    def _parse_variation(self, data: dict) -> ClinVarData | None:
        """Extract ClinVarData from a variation API response.

        Args:
            data: Raw JSON response from the ClinVar API.

        Returns:
            Populated ClinVarData or None if parsing fails.
        """
        if not data:
            return None

        # The API response structure varies by endpoint; handle both forms
        record = data if "classification" in data else data.get("result", {})
        if not record:
            return None

        classification = record.get("germline_classification", {})
        agg_class = classification.get("description") or record.get("classification")
        if isinstance(agg_class, dict):
            agg_class = agg_class.get("description")

        # Normalise classification string
        normalised = _SIGNIFICANCE_MAP.get(agg_class or "", agg_class)

        # Star rating from review status
        review_status_str = (
            classification.get("review_status") or record.get("review_status") or ""
        )
        star_rating = _review_status_to_stars(review_status_str)

        # Suggest ACMG evidence code based on stars and classification
        acmg_evidence = _suggest_acmg_evidence(normalised, star_rating)

        # RCV accessions
        rcv_list = record.get("accession_list", [])
        rcv = rcv_list[0] if rcv_list else record.get("accession")

        # Condition names
        conditions = [
            c.get("name", "") for c in record.get("trait_set", []) if c.get("name")
        ]

        # Submissions
        submissions: list[ClinVarSubmission] = []
        for sub in record.get("submissions", []):
            submissions.append(
                ClinVarSubmission(
                    submitter=sub.get("submitter_name", "Unknown"),
                    classification=sub.get("classification", ""),
                    date_last_evaluated=sub.get("last_evaluated"),
                    review_status=sub.get("review_status"),
                    condition=sub.get("condition_name"),
                    method=sub.get("collection_method"),
                )
            )

        return ClinVarData(
            variation_id=str(record.get("variation_id", "")),
            rcv_accession=rcv,
            star_rating=star_rating,
            review_status=review_status_str or _STAR_TO_STATUS.get(star_rating, ""),
            classification=normalised,
            condition_names=conditions,
            submissions=submissions,
            last_updated=record.get("date_last_updated"),
            acmg_evidence=acmg_evidence,
        )


# ---------------------------------------------------------------------------
# Private utilities
# ---------------------------------------------------------------------------


def _chrom_to_refseq(chrom_bare: str) -> str:
    """Convert a bare chromosome name to a partial RefSeq NC_ accession stub.

    This is a simplified mapping; in production use the full RefSeq accession
    (e.g. NC_000017.11 for chr17 GRCh38).

    Args:
        chrom_bare: Chromosome without "chr" prefix, e.g. "17", "X".

    Returns:
        Partial RefSeq accession stub.
    """
    # GRCh38 RefSeq accessions for chromosomes 1-22, X, Y, MT
    accession_map: dict[str, str] = {
        "1": "000001.11",
        "2": "000002.12",
        "3": "000003.12",
        "4": "000004.12",
        "5": "000005.10",
        "6": "000006.12",
        "7": "000007.14",
        "8": "000008.11",
        "9": "000009.12",
        "10": "000010.11",
        "11": "000011.10",
        "12": "000012.12",
        "13": "000013.11",
        "14": "000014.9",
        "15": "000015.10",
        "16": "000016.10",
        "17": "000017.11",
        "18": "000018.10",
        "19": "000019.10",
        "20": "000020.11",
        "21": "000021.9",
        "22": "000022.11",
        "X": "000023.11",
        "Y": "000024.10",
        "MT": "012920.1",
        "M": "012920.1",
    }
    stub = accession_map.get(chrom_bare.upper(), "000001.11")
    return stub


def _review_status_to_stars(review_status: str) -> int:
    """Convert a ClinVar review status string to a 0–4 star rating.

    Args:
        review_status: ClinVar review status description.

    Returns:
        Integer star rating 0–4.
    """
    rs = review_status.lower()
    if "no assertion" in rs:
        # e.g. "no assertion criteria provided" / "no assertion provided" —
        # must be checked before the "criteria provided" substring match
        # below, which would otherwise incorrectly award 1 star.
        return 0
    if "practice guideline" in rs:
        return 4
    if "expert panel" in rs:
        return 3
    if "multiple submitters" in rs and "no conflicts" in rs:
        return 2
    if "single submitter" in rs or "criteria provided" in rs:
        return 1
    return 0


def _suggest_acmg_evidence(classification: str | None, stars: int) -> str | None:
    """Suggest an ACMG/AMP evidence code based on ClinVar classification.

    Per ACMG/AMP 2015 guidelines and ClinGen SVI 2024:
    - PP5: Reputable source reports P/LP (criteria provided, ≥1 star)
    - BP6: Reputable source reports B/LB (criteria provided, ≥1 star)
    Note: PP5/BP6 are not in the original ACMG 2015 criteria but are widely
    used; ClinGen SVI cautions against over-reliance on these codes.

    Args:
        classification: Normalised ClinVar classification string.
        stars: Star rating (0–4).

    Returns:
        "PP5", "BP6", or None.
    """
    if not classification or stars < 1:
        return None
    if classification in (
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
    ):
        return "PP5"
    if classification in ("Benign", "Likely benign", "Benign/Likely benign"):
        return "BP6"
    return None
