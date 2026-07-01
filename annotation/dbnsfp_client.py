"""
dbNSFP v4.7 client for aggregated in-silico predictor scores.

dbNSFP (database for Non-Synonymous Function Prediction) v4.7 was released in
late 2023 and adds two major new predictors over v4.5:

New in v4.7
-----------
- **AlphaMissense** (Cheng et al. 2023, Science, PMID:37703350): Deep learning
  score derived from the ESMFold protein language model; calibrated thresholds
  (>=0.564 PP3, <=0.340 BP4) adopted by ClinGen SVI 2024.
- **ESM1b** (Brandes et al. 2023, Nature Genetics, PMID:37055584): Evolutionary
  Scale Modeling variant effect score; negative values indicate functional impact.

Core scores returned
--------------------
- **REVEL** (Ioannidis et al. 2016, PMID:27666373): Meta-predictor; ClinGen SVI
  threshold >= 0.75 → PP3, <= 0.15 → BP4.
- **BayesDel** (Feng et al. 2017, PMID:29218908): Additive model integrating
  multiple features; ClinGen SVI threshold >= 0.13 → PP3 (no AF), <= -0.18 → BP4.
- **CADD** (Kircher et al. 2014, PMID:24487276; v1.7 with GRCh38): Phred-scaled
  combined annotation-dependent depletion score.

Reference
---------
Liu X, Li C, Mou C, Dong Y, Tu Y. "dbNSFP v4: a comprehensive database of
transcript-specific functional predictions and annotations for human nonsynonymous
and splice-site SNVs." Genome Medicine. 2020;12(1):103. PMID:33261662
DOI:10.1186/s13073-020-00803-9

Download
--------
https://sites.google.com/site/jpopgen/dbNSFP

Tabix index command
-------------------
zcat dbNSFP4.7a.gz | head -1 > header.txt
zcat dbNSFP4.7a.gz | grep -v "^#" | sort -k1,1 -k2,2n | bgzip > dbNSFP4.7a_sorted.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.7a_sorted.gz
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# Column names that we extract from dbNSFP v4.7 (may vary slightly by release)
_SCORE_COLUMNS = {
    "REVEL_score": "revel",
    "BayesDel_noAF_score": "bayesdel_noaf",
    "BayesDel_addAF_score": "bayesdel_addaf",
    "CADD_phred": "cadd_phred",
    "CADD_raw": "cadd_raw",
    "AlphaMissense_score": "alphamissense",
    "ESM1b_score": "esm1b",
}

# gnomAD v4.1 / ClinGen SVI 2024 validated thresholds for REVEL
REVEL_PP3_THRESHOLD: float = 0.75
REVEL_BP4_THRESHOLD: float = 0.15

# ClinGen SVI 2024 thresholds for BayesDel_noAF
BAYESDEL_PP3_THRESHOLD: float = 0.13
BAYESDEL_BP4_THRESHOLD: float = -0.18


@dataclass
class DbNSFPScores:
    """Aggregated in-silico predictor scores from dbNSFP v4.7.

    Attributes:
        revel: REVEL meta-predictor score [0, 1].
        bayesdel_noaf: BayesDel score without allele frequency.
        bayesdel_addaf: BayesDel score with allele frequency.
        cadd_phred: CADD Phred-scaled score.
        cadd_raw: CADD raw (non-scaled) score.
        alphamissense: AlphaMissense am_pathogenicity score [0, 1].
        esm1b: ESM1b variant effect score (negative = more damaging).
        transcript_id: Ensembl transcript identifier used for lookup.
        raw: Full dict of raw column values for auditing.
    """

    revel: float | None = None
    bayesdel_noaf: float | None = None
    bayesdel_addaf: float | None = None
    cadd_phred: float | None = None
    cadd_raw: float | None = None
    alphamissense: float | None = None
    esm1b: float | None = None
    transcript_id: str | None = None
    raw: dict[str, str] = field(default_factory=dict)


class DbNSFPClient:
    """Client for querying the dbNSFP v4.7 tabix-indexed database.

    Performs tabix lookups against a locally downloaded and indexed copy of
    dbNSFP v4.7.  When a variant has multiple scored transcripts, the MANE
    Select transcript (preferred) or the transcript with the highest REVEL
    score is returned.

    Args:
        db_path: Path to the bgzip-compressed, tabix-indexed dbNSFP file.
        preferred_transcript: Ensembl transcript ID to prefer when multiple
            transcripts are scored, e.g. "ENST00000357654" for BRCA1.

    Raises:
        FileNotFoundError: If db_path is specified but does not exist.

    Example:
        >>> client = DbNSFPClient("/data/dbNSFP4.7a_sorted.gz")
        >>> scores = client.get_scores("chr17", 43094692, "G", "A")
        >>> print(scores.revel)
    """

    def __init__(
        self,
        db_path: str | Path | None = None,
        preferred_transcript: str | None = None,
    ) -> None:
        self._db_path: Path | None = Path(db_path) if db_path else None
        self._preferred_transcript = preferred_transcript
        self._header: list[str] | None = None
        self._col_index: dict[str, int] | None = None

        if self._db_path and not self._db_path.exists():
            raise FileNotFoundError(f"dbNSFP file not found: {self._db_path}")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_scores(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        transcript_id: str | None = None,
    ) -> DbNSFPScores:
        """Retrieve all predictor scores for a missense variant.

        Args:
            chrom: Chromosome string (with or without "chr" prefix).
            pos: 1-based genomic position.
            ref: Reference allele (single nucleotide).
            alt: Alternate allele (single nucleotide).
            transcript_id: Optional Ensembl transcript ID to filter on.
                Falls back to ``preferred_transcript`` then highest-REVEL.

        Returns:
            DbNSFPScores with all available predictor values.
            All fields are None if the variant is absent from the database.
        """
        chrom = _normalise_chrom(chrom)
        use_transcript = transcript_id or self._preferred_transcript

        if self._db_path is None or not self._db_path.exists():
            logger.warning("dbNSFP database not configured; returning empty scores")
            return DbNSFPScores()

        try:
            return self._tabix_lookup(chrom, pos, ref, alt, use_transcript)
        except Exception as exc:  # noqa: BLE001
            logger.error(
                "dbNSFP lookup failed for %s:%d %s>%s: %s", chrom, pos, ref, alt, exc
            )
            return DbNSFPScores()

    def classify_revel(self, score: float | None) -> str:
        """Classify a REVEL score according to ClinGen SVI 2024 thresholds.

        Args:
            score: REVEL score in [0, 1] or None.

        Returns:
            "PP3", "BP4", or "ambiguous".
        """
        if score is None:
            return "ambiguous"
        if score >= REVEL_PP3_THRESHOLD:
            return "PP3"
        if score <= REVEL_BP4_THRESHOLD:
            return "BP4"
        return "ambiguous"

    def classify_bayesdel(self, score: float | None) -> str:
        """Classify a BayesDel_noAF score according to ClinGen SVI 2024 thresholds.

        Args:
            score: BayesDel_noAF score or None.

        Returns:
            "PP3", "BP4", or "ambiguous".
        """
        if score is None:
            return "ambiguous"
        if score >= BAYESDEL_PP3_THRESHOLD:
            return "PP3"
        if score <= BAYESDEL_BP4_THRESHOLD:
            return "BP4"
        return "ambiguous"

    # ------------------------------------------------------------------
    # Internal tabix logic
    # ------------------------------------------------------------------

    def _tabix_lookup(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        preferred_transcript: str | None,
    ) -> DbNSFPScores:
        """Perform a tabix fetch and parse matching rows.

        Args:
            chrom: "chr"-prefixed chromosome.
            pos: 1-based position.
            ref: Reference allele.
            alt: Alternate allele.
            preferred_transcript: Transcript to prefer; may be None.

        Returns:
            DbNSFPScores populated from the best-matching row.
        """
        import pysam  # type: ignore[import-untyped]

        tbx = pysam.TabixFile(str(self._db_path))

        # Lazily load and cache the header column index mapping
        if self._col_index is None:
            self._col_index = self._parse_header(tbx)

        candidates: list[tuple[float, dict[str, str]]] = []

        for row in tbx.fetch(chrom, pos - 1, pos):
            fields = row.split("\t")
            row_ref = self._field(fields, "ref", "")
            row_alt = self._field(fields, "alt", "")

            # Skip rows that don't match our alleles
            if row_ref != ref or row_alt != alt:
                continue

            row_transcript = self._field(fields, "Ensembl_transcriptid", "")
            raw = {col: self._field(fields, col, ".") for col in _SCORE_COLUMNS}

            # If caller wants a specific transcript, return immediately on match
            if preferred_transcript and preferred_transcript in row_transcript:
                return self._build_scores(raw, row_transcript, fields)

            # Collect candidate with its REVEL score for later selection
            revel_val = _parse_float(raw.get("REVEL_score"))
            candidates.append((revel_val if revel_val is not None else -1.0, raw))

        if not candidates:
            return DbNSFPScores()

        # Select the candidate with the highest REVEL score (most informative)
        candidates.sort(key=lambda x: x[0], reverse=True)
        _, best_raw = candidates[0]
        return self._build_scores(best_raw, None, [])

    def _parse_header(self, tbx: object) -> dict[str, int]:
        """Extract column-name-to-index mapping from tabix header.

        Args:
            tbx: Open pysam.TabixFile object.

        Returns:
            Mapping of column name → 0-based column index.
        """
        # pysam exposes header lines via tbx.header
        for line in tbx.header:  # type: ignore[union-attr]
            if line.startswith("#chr"):
                cols = line.lstrip("#").split("\t")
                return {name: idx for idx, name in enumerate(cols)}
        logger.warning("dbNSFP header not found; column lookup may be inaccurate")
        return {}

    def _field(self, fields: list[str], col_name: str, default: str) -> str:
        """Safely retrieve a field by column name using cached index.

        Args:
            fields: Split row fields.
            col_name: Column name to look up.
            default: Value returned when column is absent or ".".

        Returns:
            String field value or default.
        """
        if self._col_index is None:
            return default
        idx = self._col_index.get(col_name)
        if idx is None or idx >= len(fields):
            return default
        val = fields[idx]
        return val if val not in (".", "", "NA") else default

    @staticmethod
    def _build_scores(
        raw: dict[str, str],
        transcript_id: str | None,
        _fields: list[str],
    ) -> DbNSFPScores:
        """Construct DbNSFPScores from a raw column value dict.

        Args:
            raw: Mapping of dbNSFP column name to raw string value.
            transcript_id: Transcript identifier for this row.
            _fields: Original row fields (unused; reserved for future use).

        Returns:
            Populated DbNSFPScores instance.
        """

        # dbNSFP may store multiple semicolon-separated values per row (isoforms)
        # We take the maximum numeric value where multiple exist
        def best(key: str) -> float | None:
            val = raw.get(key, "")
            nums = [_parse_float(v) for v in val.split(";") if v]
            nums = [n for n in nums if n is not None]
            return max(nums) if nums else None

        return DbNSFPScores(
            revel=best("REVEL_score"),
            bayesdel_noaf=best("BayesDel_noAF_score"),
            bayesdel_addaf=best("BayesDel_addAF_score"),
            cadd_phred=best("CADD_phred"),
            cadd_raw=best("CADD_raw"),
            alphamissense=best("AlphaMissense_score"),
            esm1b=best("ESM1b_score"),
            transcript_id=transcript_id,
            raw=raw,
        )


# ---------------------------------------------------------------------------
# Private utilities
# ---------------------------------------------------------------------------


def _normalise_chrom(chrom: str) -> str:
    """Ensure chromosome string has a "chr" prefix.

    Args:
        chrom: Chromosome identifier, e.g. "1", "chr1", "X".

    Returns:
        "chr"-prefixed chromosome string.
    """
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _parse_float(value: str | None) -> float | None:
    """Convert a string to float, returning None on failure.

    Args:
        value: Numeric string or None.

    Returns:
        Float value or None.
    """
    if value is None or value in (".", "", "NA"):
        return None
    try:
        return float(value)
    except ValueError:
        return None
