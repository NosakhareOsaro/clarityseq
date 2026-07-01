"""
annotation.tests.test_dbnsfp_client
======================================
pytest tests for the dbNSFP v4.7 in-silico predictor score client.

Tests cover:
    - _normalise_chrom: chr prefix normalisation.
    - _parse_float: string-to-float conversion.
    - DbNSFPClient.classify_revel: ClinGen SVI 2024 REVEL thresholds.
    - DbNSFPClient.classify_bayesdel: ClinGen SVI 2024 BayesDel thresholds.
    - DbNSFPClient.get_scores: no-db fallback, tabix error handling.
    - DbNSFPClient.__init__: FileNotFoundError on missing db.
    - DbNSFPClient._build_scores: multi-value best-score selection.
    - DbNSFPClient._field: column index lookup.
    - DbNSFPScores dataclass.

References:
    Liu X et al. 2020 Genome Medicine PMID:33261662 (dbNSFP v4).
    Ioannidis et al. 2016 PMID:27666373 (REVEL).
    ClinGen SVI 2024 thresholds (REVEL ≥0.75 PP3, ≤0.15 BP4).
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from annotation.dbnsfp_client import (
    BAYESDEL_BP4_THRESHOLD,
    BAYESDEL_PP3_THRESHOLD,
    REVEL_BP4_THRESHOLD,
    REVEL_PP3_THRESHOLD,
    DbNSFPClient,
    DbNSFPScores,
    _normalise_chrom,
    _parse_float,
)

# ---------------------------------------------------------------------------
# Shared fake dbNSFP header/row fixtures for _tabix_lookup / _parse_header
# ---------------------------------------------------------------------------

# Header line as returned via tbx.header (dbNSFP convention: leading '#chr').
_FAKE_HEADER_LINE = (
    "#chr\tpos\tref\talt\tEnsembl_transcriptid\tREVEL_score\t"
    "BayesDel_noAF_score\tBayesDel_addAF_score\tCADD_phred\tCADD_raw\t"
    "AlphaMissense_score\tESM1b_score"
)

# Matching data row: chr17 43094692 G>A on ENST00000357654.
_FAKE_MATCHING_ROW = (
    "17\t43094692\tG\tA\tENST00000357654\t0.85\t0.20\t0.18\t32.5\t5.2\t0.95\t-3.5"
)

# Non-matching alleles (different ref/alt) — should be skipped.
_FAKE_NONMATCHING_ROW = (
    "17\t43094692\tC\tT\tENST00000999999\t0.10\t.\t.\t.\t.\t.\t."
)


def _make_fake_pysam(header_lines, fetch_rows):
    """Build a fake ``pysam`` module whose TabixFile returns canned data."""
    fake_pysam = MagicMock()
    fake_tbx = MagicMock()
    fake_tbx.header = header_lines
    fake_tbx.fetch.return_value = fetch_rows
    fake_pysam.TabixFile.return_value = fake_tbx
    return fake_pysam


# ---------------------------------------------------------------------------
# _normalise_chrom tests
# ---------------------------------------------------------------------------


class TestNormaliseChrom:
    """Tests for _normalise_chrom()."""

    def test_adds_chr_prefix(self) -> None:
        """'1' → 'chr1'."""
        assert _normalise_chrom("1") == "chr1"

    def test_x_gets_prefix(self) -> None:
        """'X' → 'chrX'."""
        assert _normalise_chrom("X") == "chrX"

    def test_already_prefixed_unchanged(self) -> None:
        """'chr17' unchanged."""
        assert _normalise_chrom("chr17") == "chr17"

    def test_chry_unchanged(self) -> None:
        """'chrY' unchanged."""
        assert _normalise_chrom("chrY") == "chrY"


# ---------------------------------------------------------------------------
# _parse_float tests
# ---------------------------------------------------------------------------


class TestParseFloat:
    """Tests for _parse_float()."""

    def test_valid_string_parses(self) -> None:
        """Valid numeric string returns float."""
        assert _parse_float("0.75") == pytest.approx(0.75)

    def test_dot_returns_none(self) -> None:
        """'.' (missing value) returns None."""
        assert _parse_float(".") is None

    def test_empty_string_returns_none(self) -> None:
        """Empty string returns None."""
        assert _parse_float("") is None

    def test_none_input_returns_none(self) -> None:
        """None input returns None."""
        assert _parse_float(None) is None

    def test_na_returns_none(self) -> None:
        """'NA' returns None."""
        assert _parse_float("NA") is None

    def test_negative_value_parses(self) -> None:
        """Negative float string parses correctly."""
        assert _parse_float("-0.25") == pytest.approx(-0.25)

    def test_invalid_string_returns_none(self) -> None:
        """Non-numeric string returns None."""
        assert _parse_float("not_a_number") is None


# ---------------------------------------------------------------------------
# DbNSFPScores dataclass
# ---------------------------------------------------------------------------


class TestDbNSFPScores:
    """Tests for DbNSFPScores dataclass."""

    def test_all_fields_none_by_default(self) -> None:
        """Default DbNSFPScores has all score fields as None."""
        scores = DbNSFPScores()
        assert scores.revel is None
        assert scores.bayesdel_noaf is None
        assert scores.cadd_phred is None
        assert scores.alphamissense is None

    def test_construction_with_values(self) -> None:
        """DbNSFPScores can be constructed with explicit values."""
        scores = DbNSFPScores(
            revel=0.85,
            bayesdel_noaf=0.20,
            cadd_phred=32.5,
            alphamissense=0.95,
        )
        assert scores.revel == pytest.approx(0.85)
        assert scores.cadd_phred == pytest.approx(32.5)


# ---------------------------------------------------------------------------
# DbNSFPClient.__init__ tests
# ---------------------------------------------------------------------------


class TestDbNSFPClientInit:
    """Tests for DbNSFPClient constructor."""

    def test_no_db_path_creates_client(self) -> None:
        """Client can be created without a db_path."""
        client = DbNSFPClient()
        assert client._db_path is None

    def test_missing_db_path_raises_file_not_found(self, tmp_path: Path) -> None:
        """Non-existent db_path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="dbNSFP file not found"):
            DbNSFPClient(db_path=tmp_path / "nonexistent.gz")

    def test_preferred_transcript_stored(self) -> None:
        """preferred_transcript is stored on the client."""
        client = DbNSFPClient(preferred_transcript="ENST00000357654")
        assert client._preferred_transcript == "ENST00000357654"


# ---------------------------------------------------------------------------
# DbNSFPClient.classify_revel tests
# ---------------------------------------------------------------------------


class TestClassifyRevel:
    """Tests for DbNSFPClient.classify_revel() ClinGen SVI 2024 thresholds."""

    def setup_method(self) -> None:
        self.client = DbNSFPClient()

    def test_high_revel_returns_pp3(self) -> None:
        """REVEL ≥ 0.75 → PP3 (pathogenic in silico evidence)."""
        assert self.client.classify_revel(0.90) == "PP3"

    def test_exactly_threshold_pp3(self) -> None:
        """REVEL == 0.75 → PP3."""
        assert self.client.classify_revel(REVEL_PP3_THRESHOLD) == "PP3"

    def test_low_revel_returns_bp4(self) -> None:
        """REVEL ≤ 0.15 → BP4 (benign in silico evidence)."""
        assert self.client.classify_revel(0.10) == "BP4"

    def test_exactly_threshold_bp4(self) -> None:
        """REVEL == 0.15 → BP4."""
        assert self.client.classify_revel(REVEL_BP4_THRESHOLD) == "BP4"

    def test_intermediate_revel_is_ambiguous(self) -> None:
        """REVEL between thresholds → ambiguous."""
        assert self.client.classify_revel(0.50) == "ambiguous"

    def test_none_revel_is_ambiguous(self) -> None:
        """None REVEL → ambiguous."""
        assert self.client.classify_revel(None) == "ambiguous"


# ---------------------------------------------------------------------------
# DbNSFPClient.classify_bayesdel tests
# ---------------------------------------------------------------------------


class TestClassifyBayesDel:
    """Tests for DbNSFPClient.classify_bayesdel() ClinGen SVI 2024 thresholds."""

    def setup_method(self) -> None:
        self.client = DbNSFPClient()

    def test_high_bayesdel_returns_pp3(self) -> None:
        """BayesDel ≥ 0.13 → PP3."""
        assert self.client.classify_bayesdel(0.25) == "PP3"

    def test_exactly_threshold_pp3(self) -> None:
        """BayesDel == 0.13 → PP3."""
        assert self.client.classify_bayesdel(BAYESDEL_PP3_THRESHOLD) == "PP3"

    def test_low_bayesdel_returns_bp4(self) -> None:
        """BayesDel ≤ -0.18 → BP4."""
        assert self.client.classify_bayesdel(-0.30) == "BP4"

    def test_exactly_threshold_bp4(self) -> None:
        """BayesDel == -0.18 → BP4."""
        assert self.client.classify_bayesdel(BAYESDEL_BP4_THRESHOLD) == "BP4"

    def test_intermediate_bayesdel_is_ambiguous(self) -> None:
        """BayesDel between thresholds → ambiguous."""
        assert self.client.classify_bayesdel(0.00) == "ambiguous"

    def test_none_bayesdel_is_ambiguous(self) -> None:
        """None BayesDel → ambiguous."""
        assert self.client.classify_bayesdel(None) == "ambiguous"


# ---------------------------------------------------------------------------
# DbNSFPClient.get_scores tests
# ---------------------------------------------------------------------------


class TestGetScores:
    """Tests for DbNSFPClient.get_scores()."""

    def test_no_db_configured_returns_empty_scores(self) -> None:
        """Client without db_path returns empty DbNSFPScores."""
        client = DbNSFPClient()
        scores = client.get_scores("chr17", 43094692, "G", "A")
        assert isinstance(scores, DbNSFPScores)
        assert scores.revel is None

    def test_tabix_exception_returns_empty_scores(self, tmp_path: Path) -> None:
        """Exception during tabix lookup returns empty DbNSFPScores (no raise)."""
        db_file = tmp_path / "dbnsfp.gz"
        db_file.write_bytes(b"fake")

        client = DbNSFPClient(db_path=db_file)
        with patch.object(client, "_tabix_lookup", side_effect=Exception("tabix error")):
            scores = client.get_scores("chr17", 43094692, "G", "A")

        assert isinstance(scores, DbNSFPScores)
        assert scores.revel is None

    def test_chrom_normalised_before_lookup(self, tmp_path: Path) -> None:
        """Chromosome is normalised to chr-prefixed form before tabix lookup."""
        db_file = tmp_path / "dbnsfp.gz"
        db_file.write_bytes(b"fake")

        client = DbNSFPClient(db_path=db_file)
        with patch.object(
            client, "_tabix_lookup", return_value=DbNSFPScores(revel=0.85)
        ) as mock_lookup:
            client.get_scores("17", 43094692, "G", "A")

        call_chrom = mock_lookup.call_args[0][0]
        assert call_chrom == "chr17"


# ---------------------------------------------------------------------------
# DbNSFPClient._build_scores tests
# ---------------------------------------------------------------------------


class TestBuildScores:
    """Tests for DbNSFPClient._build_scores()."""

    def test_semicolon_separated_selects_max(self) -> None:
        """Semicolon-separated scores → maximum value selected."""
        raw = {
            "REVEL_score": "0.60;0.85;0.70",
            "BayesDel_noAF_score": ".",
            "BayesDel_addAF_score": ".",
            "CADD_phred": ".",
            "CADD_raw": ".",
            "AlphaMissense_score": ".",
            "ESM1b_score": ".",
        }
        scores = DbNSFPClient._build_scores(raw, None, [])
        assert scores.revel == pytest.approx(0.85)

    def test_missing_value_returns_none(self) -> None:
        """'.' value returns None for that score."""
        raw = {
            "REVEL_score": ".",
            "BayesDel_noAF_score": ".",
            "BayesDel_addAF_score": ".",
            "CADD_phred": ".",
            "CADD_raw": ".",
            "AlphaMissense_score": ".",
            "ESM1b_score": ".",
        }
        scores = DbNSFPClient._build_scores(raw, None, [])
        assert scores.revel is None

    def test_transcript_id_set(self) -> None:
        """Transcript ID from argument is stored in DbNSFPScores."""
        raw = {k: "." for k in ["REVEL_score", "BayesDel_noAF_score",
                                   "BayesDel_addAF_score", "CADD_phred",
                                   "CADD_raw", "AlphaMissense_score", "ESM1b_score"]}
        scores = DbNSFPClient._build_scores(raw, "ENST00000357654", [])
        assert scores.transcript_id == "ENST00000357654"

    def test_all_scores_populated(self) -> None:
        """All score fields populated from raw dict."""
        raw = {
            "REVEL_score": "0.85",
            "BayesDel_noAF_score": "0.20",
            "BayesDel_addAF_score": "0.18",
            "CADD_phred": "32.5",
            "CADD_raw": "5.2",
            "AlphaMissense_score": "0.95",
            "ESM1b_score": "-3.5",
        }
        scores = DbNSFPClient._build_scores(raw, None, [])
        assert scores.revel == pytest.approx(0.85)
        assert scores.bayesdel_noaf == pytest.approx(0.20)
        assert scores.cadd_phred == pytest.approx(32.5)
        assert scores.alphamissense == pytest.approx(0.95)
        assert scores.esm1b == pytest.approx(-3.5)


# ---------------------------------------------------------------------------
# DbNSFPClient._field tests
# ---------------------------------------------------------------------------


class TestDbNSFPClientField:
    """Tests for DbNSFPClient._field()."""

    def test_returns_field_by_column_name(self) -> None:
        """Returns field value when col_index is set."""
        client = DbNSFPClient()
        client._col_index = {"REVEL_score": 2}
        fields = ["chr17", "43094692", "0.85"]
        assert client._field(fields, "REVEL_score", ".") == "0.85"

    def test_returns_default_when_col_not_in_index(self) -> None:
        """Returns default when column not in col_index."""
        client = DbNSFPClient()
        client._col_index = {}
        assert client._field(["a", "b"], "REVEL_score", "X") == "X"

    def test_returns_default_when_no_col_index(self) -> None:
        """Returns default when col_index is None."""
        client = DbNSFPClient()
        client._col_index = None
        assert client._field(["a", "b"], "REVEL_score", "default") == "default"

    def test_returns_default_for_dot_value(self) -> None:
        """Returns default when field value is '.'."""
        client = DbNSFPClient()
        client._col_index = {"REVEL_score": 0}
        assert client._field(["."], "REVEL_score", "missing") == "missing"


# ---------------------------------------------------------------------------
# DbNSFPClient._parse_header tests
# ---------------------------------------------------------------------------


class TestParseHeader:
    """Tests for DbNSFPClient._parse_header()."""

    def test_parses_column_index_from_hash_chr_line(self) -> None:
        client = DbNSFPClient()
        fake_tbx = MagicMock()
        fake_tbx.header = [_FAKE_HEADER_LINE]
        col_index = client._parse_header(fake_tbx)
        assert col_index["REVEL_score"] == 5
        assert col_index["chr"] == 0
        assert col_index["Ensembl_transcriptid"] == 4

    def test_no_hash_chr_line_returns_empty_dict(self) -> None:
        """Header without a '#chr'-prefixed line returns an empty mapping."""
        client = DbNSFPClient()
        fake_tbx = MagicMock()
        fake_tbx.header = ["##some other metadata line"]
        col_index = client._parse_header(fake_tbx)
        assert col_index == {}


# ---------------------------------------------------------------------------
# DbNSFPClient._tabix_lookup tests (pysam mocked via sys.modules)
# ---------------------------------------------------------------------------


class TestTabixLookup:
    """Tests for DbNSFPClient._tabix_lookup() with a fake pysam module.

    pysam is not installed in this environment (and the client imports it
    lazily inside the method), so we inject a fake module into
    ``sys.modules`` for the duration of each test.
    """

    def _make_client(self, tmp_path: Path) -> DbNSFPClient:
        db_file = tmp_path / "dbnsfp.gz"
        db_file.write_bytes(b"fake")
        return DbNSFPClient(db_path=db_file)

    def test_preferred_transcript_match_returns_immediately(
        self, tmp_path: Path
    ) -> None:
        """When preferred_transcript matches the row, scores are returned
        without waiting to collect/sort all candidates."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam(
            [_FAKE_HEADER_LINE], [_FAKE_MATCHING_ROW]
        )
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = client._tabix_lookup(
                "chr17", 43094692, "G", "A", "ENST00000357654"
            )

        assert result.revel == pytest.approx(0.85)
        assert result.transcript_id == "ENST00000357654"
        assert result.cadd_phred == pytest.approx(32.5)

    def test_no_preferred_transcript_selects_best_revel_candidate(
        self, tmp_path: Path
    ) -> None:
        """Without a preferred transcript, the row is collected as a
        candidate and the highest-REVEL candidate is returned."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam(
            [_FAKE_HEADER_LINE], [_FAKE_MATCHING_ROW]
        )
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = client._tabix_lookup("chr17", 43094692, "G", "A", None)

        assert result.revel == pytest.approx(0.85)

    def test_mismatched_alleles_are_skipped(self, tmp_path: Path) -> None:
        """Rows with a different ref/alt than requested must be skipped,
        leaving no candidates and empty DbNSFPScores."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam(
            [_FAKE_HEADER_LINE], [_FAKE_NONMATCHING_ROW]
        )
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = client._tabix_lookup("chr17", 43094692, "G", "A", None)

        assert result.revel is None

    def test_no_fetch_results_returns_empty_scores(self, tmp_path: Path) -> None:
        """No rows returned by tabix fetch yields an empty DbNSFPScores."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam([_FAKE_HEADER_LINE], [])
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = client._tabix_lookup("chr17", 43094692, "G", "A", None)

        assert isinstance(result, DbNSFPScores)
        assert result.revel is None

    def test_header_cached_after_first_lookup(self, tmp_path: Path) -> None:
        """_col_index should be populated (cached) after a lookup."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam(
            [_FAKE_HEADER_LINE], [_FAKE_MATCHING_ROW]
        )
        assert client._col_index is None
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            client._tabix_lookup("chr17", 43094692, "G", "A", None)

        assert client._col_index is not None
        assert "REVEL_score" in client._col_index

    def test_get_scores_end_to_end_with_fake_pysam(self, tmp_path: Path) -> None:
        """get_scores() should surface real _tabix_lookup results end-to-end."""
        client = self._make_client(tmp_path)
        fake_pysam = _make_fake_pysam(
            [_FAKE_HEADER_LINE], [_FAKE_MATCHING_ROW]
        )
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = client.get_scores("17", 43094692, "G", "A")

        assert result.revel == pytest.approx(0.85)
        assert result.alphamissense == pytest.approx(0.95)
