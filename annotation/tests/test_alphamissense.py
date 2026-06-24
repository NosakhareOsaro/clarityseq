"""
Tests for the AlphaMissense client and classification thresholds.

ClinGen SVI 2024 thresholds under test:
  score >= 0.564 → PP3
  score <= 0.340 → BP4
  0.340 < score < 0.564 → ambiguous
"""

from __future__ import annotations

import pytest
from unittest.mock import AsyncMock, MagicMock, patch

from annotation.alphamissense_client import (
    AM_BP4_THRESHOLD,
    AM_PP3_THRESHOLD,
    AlphaMissenseClient,
    AlphaMissenseResult,
    classify_am_score,
)


# ---------------------------------------------------------------------------
# classify_am_score — unit tests (no IO)
# ---------------------------------------------------------------------------


class TestClassifyAmScore:
    """Tests for the classify_am_score() module-level function."""

    def test_above_pp3_threshold_returns_pp3(self) -> None:
        """Score at exactly the PP3 threshold should return PP3."""
        assert classify_am_score(AM_PP3_THRESHOLD) == "PP3"

    def test_well_above_pp3_returns_pp3(self) -> None:
        """Score clearly above threshold returns PP3."""
        assert classify_am_score(0.95) == "PP3"

    def test_at_unity_returns_pp3(self) -> None:
        """Maximum possible score of 1.0 should return PP3."""
        assert classify_am_score(1.0) == "PP3"

    def test_just_below_pp3_returns_ambiguous(self) -> None:
        """Score just below PP3 threshold falls in ambiguous zone."""
        assert classify_am_score(AM_PP3_THRESHOLD - 0.001) == "ambiguous"

    def test_midpoint_returns_ambiguous(self) -> None:
        """Midpoint between thresholds is ambiguous."""
        midpoint = (AM_PP3_THRESHOLD + AM_BP4_THRESHOLD) / 2
        assert classify_am_score(midpoint) == "ambiguous"

    def test_at_bp4_threshold_returns_bp4(self) -> None:
        """Score at exactly the BP4 threshold should return BP4."""
        assert classify_am_score(AM_BP4_THRESHOLD) == "BP4"

    def test_below_bp4_threshold_returns_bp4(self) -> None:
        """Score clearly below BP4 threshold returns BP4."""
        assert classify_am_score(0.10) == "BP4"

    def test_at_zero_returns_bp4(self) -> None:
        """Minimum possible score of 0.0 should return BP4."""
        assert classify_am_score(0.0) == "BP4"

    def test_just_above_bp4_returns_ambiguous(self) -> None:
        """Score just above BP4 threshold falls in ambiguous zone."""
        assert classify_am_score(AM_BP4_THRESHOLD + 0.001) == "ambiguous"

    def test_none_returns_ambiguous(self) -> None:
        """None score (variant not scored) should return ambiguous."""
        assert classify_am_score(None) == "ambiguous"

    def test_threshold_boundary_values(self) -> None:
        """Test all four quadrant boundary values explicitly."""
        # Exactly at PP3 boundary
        assert classify_am_score(0.564) == "PP3"
        # Exactly at BP4 boundary
        assert classify_am_score(0.340) == "BP4"
        # Mid-ambiguous
        assert classify_am_score(0.45) == "ambiguous"
        # Clearly pathogenic
        assert classify_am_score(0.999) == "PP3"
        # Clearly benign
        assert classify_am_score(0.001) == "BP4"


# ---------------------------------------------------------------------------
# AlphaMissenseClient — HTTP fallback
# ---------------------------------------------------------------------------


class TestAlphaMissenseClientHTTPFallback:
    """Tests for the HTTP fallback path when no tabix file is configured."""

    @pytest.mark.asyncio
    async def test_http_fallback_returns_result(self) -> None:
        """HTTP fallback should return a valid AlphaMissenseResult."""
        client = AlphaMissenseClient(tsv_path=None)

        mock_response_data = {
            "am_pathogenicity": 0.82,
            "am_class": "likely_pathogenic",
        }

        with patch.object(client, "_http_lookup", new_callable=AsyncMock) as mock_lookup:
            mock_lookup.return_value = (0.82, "likely_pathogenic")
            result = await client.get_am_score("chr17", 43094692, "G", "A")

        assert isinstance(result, AlphaMissenseResult)
        assert result.score == pytest.approx(0.82)
        assert result.evidence_code == "PP3"
        assert result.chrom == "chr17"
        assert result.pos == 43094692

    @pytest.mark.asyncio
    async def test_http_fallback_not_found_returns_ambiguous(self) -> None:
        """Missing variant should yield ambiguous evidence code."""
        client = AlphaMissenseClient(tsv_path=None)

        with patch.object(client, "_http_lookup", new_callable=AsyncMock) as mock_lookup:
            mock_lookup.return_value = (None, None)
            result = await client.get_am_score("chr1", 12345, "A", "C")

        assert result.score is None
        assert result.evidence_code == "ambiguous"

    @pytest.mark.asyncio
    async def test_benign_score_returns_bp4(self) -> None:
        """Score below BP4 threshold should yield BP4 evidence code."""
        client = AlphaMissenseClient(tsv_path=None)

        with patch.object(client, "_http_lookup", new_callable=AsyncMock) as mock_lookup:
            mock_lookup.return_value = (0.15, "likely_benign")
            result = await client.get_am_score("chr1", 99999, "T", "C")

        assert result.evidence_code == "BP4"
        assert result.score == pytest.approx(0.15)

    @pytest.mark.asyncio
    async def test_chrom_prefix_normalised(self) -> None:
        """Client should normalise both '1' and 'chr1' to 'chr1'."""
        client = AlphaMissenseClient(tsv_path=None)

        with patch.object(client, "_http_lookup", new_callable=AsyncMock) as mock_lookup:
            mock_lookup.return_value = (0.70, "likely_pathogenic")
            result = await client.get_am_score("1", 12345, "A", "G")

        # Chromosome should be normalised to chr-prefixed form
        assert result.chrom == "chr1"


# ---------------------------------------------------------------------------
# AlphaMissenseClient — tabix path
# ---------------------------------------------------------------------------


class TestAlphaMissenseClientTabix:
    """Tests for the tabix lookup path."""

    @pytest.mark.asyncio
    async def test_tabix_lookup_called_when_file_exists(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        """Tabix lookup should be attempted when TSV file exists."""
        # Create a dummy file so Path.exists() returns True
        tsv_file = tmp_path / "AlphaMissense_hg38.tsv.gz"
        tsv_file.touch()

        client = AlphaMissenseClient(tsv_path=str(tsv_file))
        client._pysam_available = True  # Pretend pysam is available

        with patch.object(client, "_tabix_lookup") as mock_tabix:
            mock_tabix.return_value = (0.75, "likely_pathogenic")
            result = await client.get_am_score("chr17", 43094692, "G", "A")

        mock_tabix.assert_called_once_with("chr17", 43094692, "G", "A")
        assert result.evidence_code == "PP3"

    @pytest.mark.asyncio
    async def test_http_fallback_used_when_tabix_unavailable(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        """HTTP fallback should be used when pysam is not available."""
        tsv_file = tmp_path / "AlphaMissense_hg38.tsv.gz"
        tsv_file.touch()

        client = AlphaMissenseClient(tsv_path=str(tsv_file))
        client._pysam_available = False  # Simulate pysam not installed

        with patch.object(client, "_http_lookup", new_callable=AsyncMock) as mock_http:
            mock_http.return_value = (0.25, "likely_benign")
            result = await client.get_am_score("chr1", 1000, "A", "T")

        mock_http.assert_called_once()
        assert result.evidence_code == "BP4"


# ---------------------------------------------------------------------------
# Constants sanity check
# ---------------------------------------------------------------------------


class TestThresholdConstants:
    """Verify the ClinGen SVI 2024 threshold constants are correct."""

    def test_pp3_threshold_value(self) -> None:
        """PP3 threshold must be 0.564 (ClinGen SVI 2024)."""
        assert AM_PP3_THRESHOLD == pytest.approx(0.564)

    def test_bp4_threshold_value(self) -> None:
        """BP4 threshold must be 0.340 (ClinGen SVI 2024)."""
        assert AM_BP4_THRESHOLD == pytest.approx(0.340)

    def test_pp3_greater_than_bp4(self) -> None:
        """PP3 threshold must be strictly greater than BP4 threshold."""
        assert AM_PP3_THRESHOLD > AM_BP4_THRESHOLD

    def test_ambiguous_range_is_positive_width(self) -> None:
        """There must be a positive-width ambiguous range between thresholds."""
        ambiguous_width = AM_PP3_THRESHOLD - AM_BP4_THRESHOLD
        assert ambiguous_width > 0, "Ambiguous range must have positive width"
