"""
AlphaMissense client for missense pathogenicity scoring.

Reference
---------
Cheng et al. (2023) "Accurate proteome-wide missense variant effect prediction
with AlphaMissense." Science 381, eadg7492. PMID:37703350
DOI:10.1126/science.adg7492

AlphaMissense scores range from 0 (benign) to 1 (pathogenic).

ClinGen SVI 2024 Thresholds
---------------------------
These thresholds were calibrated against ClinVar P/LP and B/LB variants using
the Tavtigian et al. 2020 (PMID:32645316) Bayesian framework:

  score >= 0.564  → PP3 (Supporting Pathogenic)
  score <= 0.340  → BP4 (Supporting Benign)
  0.340 < score < 0.564 → ambiguous (no evidence contribution)

Source: ClinGen SVI Working Group 2024 calibration memo (v2024-05).

Download URL
------------
https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

Tabix Indexing Command
----------------------
bgzip -d -k AlphaMissense_hg38.tsv.gz  # keep original
bgzip AlphaMissense_hg38.tsv
tabix -s 1 -b 2 -e 2 -c '#' AlphaMissense_hg38.tsv.gz

The file uses 1-based coordinates; column layout:
  #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id
  protein_variant  am_pathogenicity  am_class
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

# ClinGen SVI 2024 calibrated thresholds (PMID:37703350 + SVI 2024 memo)
AM_PP3_THRESHOLD: float = 0.564  # >= this → PP3 Supporting Pathogenic
AM_BP4_THRESHOLD: float = 0.340  # <= this → BP4 Supporting Benign

# Public download URL for the hg38 scored TSV
AM_DOWNLOAD_URL: str = (
    "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
)

# Expected column indices in the decompressed TSV (0-based, after '#CHROM' header)
_COL_CHROM = 0
_COL_POS = 1
_COL_REF = 2
_COL_ALT = 3
_COL_SCORE = 8  # am_pathogenicity
_COL_CLASS = 9  # am_class (likely_pathogenic / likely_benign / ambiguous)

# REST API fallback – not an official endpoint; use tabix-indexed file in production
_AM_API_BASE = "https://alphamissense.hegelab.org/api/v1"


@dataclass(frozen=True)
class AlphaMissenseResult:
    """Container for a single AlphaMissense lookup result.

    Attributes:
        chrom: Chromosome (e.g. "chr1").
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
        score: am_pathogenicity score in [0, 1], or None if not found.
        am_class: Raw class string from the TSV ("likely_pathogenic",
            "likely_benign", "ambiguous"), or None.
        evidence_code: ACMG/AMP evidence code ("PP3", "BP4", or "ambiguous").
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    score: Optional[float]
    am_class: Optional[str]
    evidence_code: str


class AlphaMissenseClient:
    """Client for retrieving AlphaMissense scores via tabix or HTTP fallback.

    Prefers local tabix-indexed file (fastest, no network dependency).
    Falls back to the community REST API when the file is absent.

    Args:
        tsv_path: Path to the bgzip-compressed, tabix-indexed TSV file
            (``AlphaMissense_hg38.tsv.gz``).  May be ``None`` to force
            HTTP fallback.
        http_timeout: Timeout in seconds for HTTP requests.

    Example:
        >>> client = AlphaMissenseClient(tsv_path="/data/AlphaMissense_hg38.tsv.gz")
        >>> result = await client.get_am_score("chr17", 43094692, "G", "A")
        >>> print(result.evidence_code)  # "PP3", "BP4", or "ambiguous"
    """

    def __init__(
        self,
        tsv_path: Optional[str | Path] = None,
        http_timeout: float = 10.0,
    ) -> None:
        self._tsv_path: Optional[Path] = Path(tsv_path) if tsv_path else None
        self._http_timeout = http_timeout
        self._pysam_available = self._check_pysam()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def get_am_score(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> AlphaMissenseResult:
        """Look up the AlphaMissense score for a missense variant.

        Args:
            chrom: Chromosome string, with or without "chr" prefix.
            pos: 1-based genomic position.
            ref: Reference nucleotide (single base).
            alt: Alternate nucleotide (single base).

        Returns:
            AlphaMissenseResult with score and ACMG evidence code.
            If the variant is not found, score is None and evidence_code
            is "ambiguous".
        """
        # Normalise chromosome name to "chr"-prefixed form
        chrom = _normalise_chrom(chrom)

        score: Optional[float] = None
        am_class: Optional[str] = None

        if self._tsv_path and self._tsv_path.exists() and self._pysam_available:
            score, am_class = self._tabix_lookup(chrom, pos, ref, alt)
        else:
            score, am_class = await self._http_lookup(chrom, pos, ref, alt)

        evidence_code = classify_am_score(score)
        return AlphaMissenseResult(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            score=score,
            am_class=am_class,
            evidence_code=evidence_code,
        )

    # ------------------------------------------------------------------
    # Tabix lookup
    # ------------------------------------------------------------------

    def _tabix_lookup(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> tuple[Optional[float], Optional[str]]:
        """Query the tabix-indexed TSV for a specific variant.

        Args:
            chrom: Chromosome with "chr" prefix.
            pos: 1-based position.
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            Tuple of (score, am_class), both None if not found.
        """
        try:
            import pysam  # type: ignore[import-untyped]

            tbx = pysam.TabixFile(str(self._tsv_path))
            # tabix fetch uses 0-based half-open coordinates internally
            for row in tbx.fetch(chrom, pos - 1, pos):
                fields = row.split("\t")
                if (
                    len(fields) > _COL_CLASS
                    and fields[_COL_REF] == ref
                    and fields[_COL_ALT] == alt
                ):
                    raw_score = fields[_COL_SCORE]
                    am_class = fields[_COL_CLASS].strip()
                    return float(raw_score), am_class
        except Exception as exc:  # noqa: BLE001
            logger.warning("Tabix lookup failed: %s", exc)
        return None, None

    # ------------------------------------------------------------------
    # HTTP fallback
    # ------------------------------------------------------------------

    async def _http_lookup(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> tuple[Optional[float], Optional[str]]:
        """Retrieve score from community REST API fallback.

        Args:
            chrom: Chromosome with "chr" prefix.
            pos: 1-based position.
            ref: Reference allele.
            alt: Alternate allele.

        Returns:
            Tuple of (score, am_class), both None on failure.
        """
        # Strip "chr" for the API query string
        bare_chrom = chrom.replace("chr", "")
        url = f"{_AM_API_BASE}/variant/{bare_chrom}/{pos}/{ref}/{alt}"
        try:
            async with httpx.AsyncClient(timeout=self._http_timeout) as client:
                resp = await client.get(url)
                resp.raise_for_status()
                data = resp.json()
                score = data.get("am_pathogenicity")
                am_class = data.get("am_class")
                return (float(score) if score is not None else None), am_class
        except httpx.HTTPError as exc:
            logger.warning("AlphaMissense HTTP lookup failed for %s:%d: %s", chrom, pos, exc)
            return None, None

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _check_pysam() -> bool:
        """Return True if pysam is importable."""
        try:
            import pysam  # noqa: F401  # type: ignore[import-untyped]

            return True
        except ImportError:
            logger.info(
                "pysam not available; tabix lookups disabled. "
                "Install with: pip install pysam"
            )
            return False


# ---------------------------------------------------------------------------
# Module-level classification helper
# ---------------------------------------------------------------------------


def classify_am_score(score: Optional[float]) -> str:
    """Map an AlphaMissense score to a ClinGen SVI 2024 ACMG evidence code.

    Thresholds from ClinGen SVI Working Group 2024 calibration:
    - score >= 0.564 → PP3 (Supporting Pathogenic evidence)
    - score <= 0.340 → BP4 (Supporting Benign evidence)
    - 0.340 < score < 0.564 → "ambiguous" (no evidence contribution)

    Args:
        score: AlphaMissense am_pathogenicity score in [0.0, 1.0], or None
            for variants not scored (e.g. synonymous, intronic).

    Returns:
        "PP3", "BP4", or "ambiguous".

    Examples:
        >>> classify_am_score(0.85)
        'PP3'
        >>> classify_am_score(0.20)
        'BP4'
        >>> classify_am_score(0.45)
        'ambiguous'
        >>> classify_am_score(None)
        'ambiguous'
    """
    if score is None:
        return "ambiguous"
    if score >= AM_PP3_THRESHOLD:
        return "PP3"
    if score <= AM_BP4_THRESHOLD:
        return "BP4"
    return "ambiguous"


# ---------------------------------------------------------------------------
# Private utilities
# ---------------------------------------------------------------------------


def _normalise_chrom(chrom: str) -> str:
    """Ensure chromosome string has a "chr" prefix.

    Args:
        chrom: Chromosome identifier, e.g. "1", "chr1", "X".

    Returns:
        Chromosome string with "chr" prefix, e.g. "chr1", "chrX".
    """
    if not chrom.startswith("chr"):
        return f"chr{chrom}"
    return chrom
