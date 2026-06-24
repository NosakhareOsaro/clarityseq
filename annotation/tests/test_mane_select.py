"""
Tests for MANE Select transcript utility functions.

Validates:
- is_mane_select() with RefSeq and Ensembl IDs
- get_mane_select_for_gene() lookup
- adjust_pvs1_for_mane() strength downgrade logic (ACGS 2024 §5)
"""

from __future__ import annotations

import pytest

from annotation.mane_select import (
    adjust_pvs1_for_mane,
    get_mane_select_ensembl_for_gene,
    get_mane_select_for_gene,
    is_mane_select,
)


# ---------------------------------------------------------------------------
# Tests for is_mane_select()
# ---------------------------------------------------------------------------


class TestIsManeSelect:
    """Tests for the is_mane_select() function."""

    def test_brca1_refseq_transcript_recognised(self) -> None:
        """BRCA1 MANE Select RefSeq transcript should be recognised."""
        assert is_mane_select("NM_007294.4") is True

    def test_brca2_refseq_transcript_recognised(self) -> None:
        """BRCA2 MANE Select RefSeq transcript should be recognised."""
        assert is_mane_select("NM_000059.4") is True

    def test_brca1_ensembl_transcript_recognised(self) -> None:
        """BRCA1 MANE Select Ensembl transcript should be recognised."""
        assert is_mane_select("ENST00000357654.9") is True

    def test_unknown_refseq_returns_false(self) -> None:
        """Unknown RefSeq ID should return False."""
        assert is_mane_select("NM_999999.1") is False

    def test_unknown_ensembl_returns_false(self) -> None:
        """Unknown Ensembl ID should return False."""
        assert is_mane_select("ENST99999999.1") is False

    def test_version_stripped_for_matching(self) -> None:
        """Version suffix should be stripped for flexible matching."""
        # NM_007294.4 is BRCA1; version .99 should still match the accession base
        assert is_mane_select("NM_007294.99") is True

    def test_empty_string_returns_false(self) -> None:
        """Empty string should return False."""
        assert is_mane_select("") is False

    def test_cftr_transcript_recognised(self) -> None:
        """CFTR MANE Select transcript should be recognised."""
        assert is_mane_select("NM_000492.4") is True

    def test_tp53_transcript_recognised(self) -> None:
        """TP53 MANE Select transcript should be recognised."""
        assert is_mane_select("NM_000546.6") is True

    def test_hbb_ensembl_recognised(self) -> None:
        """HBB Ensembl MANE Select transcript should be recognised."""
        assert is_mane_select("ENST00000335295.4") is True


# ---------------------------------------------------------------------------
# Tests for get_mane_select_for_gene()
# ---------------------------------------------------------------------------


class TestGetManeSelectForGene:
    """Tests for the get_mane_select_for_gene() function."""

    def test_brca1_returns_correct_refseq(self) -> None:
        """BRCA1 should return its known MANE Select RefSeq ID."""
        result = get_mane_select_for_gene("BRCA1")
        assert result == "NM_007294.4"

    def test_brca2_returns_correct_refseq(self) -> None:
        """BRCA2 should return its known MANE Select RefSeq ID."""
        result = get_mane_select_for_gene("BRCA2")
        assert result == "NM_000059.4"

    def test_tp53_returns_correct_refseq(self) -> None:
        """TP53 should return the correct RefSeq transcript."""
        result = get_mane_select_for_gene("TP53")
        assert result == "NM_000546.6"

    def test_cftr_returns_correct_refseq(self) -> None:
        """CFTR should return its MANE Select transcript."""
        result = get_mane_select_for_gene("CFTR")
        assert result == "NM_000492.4"

    def test_unknown_gene_returns_none(self) -> None:
        """Unknown gene symbol should return None."""
        assert get_mane_select_for_gene("FAKEGENE123") is None

    def test_lowercase_gene_symbol_returns_none(self) -> None:
        """Gene symbols are case-sensitive; lowercase should return None."""
        # The registry uses uppercase; lowercase should miss
        result = get_mane_select_for_gene("brca1")
        assert result is None

    def test_dmd_returns_correct_refseq(self) -> None:
        """DMD (Duchenne muscular dystrophy) should return the correct RefSeq."""
        result = get_mane_select_for_gene("DMD")
        assert result == "NM_004006.3"


class TestGetManeSelectEnsemblForGene:
    """Tests for get_mane_select_ensembl_for_gene()."""

    def test_brca1_returns_ensembl_id(self) -> None:
        """BRCA1 should return its Ensembl MANE Select ID."""
        result = get_mane_select_ensembl_for_gene("BRCA1")
        assert result == "ENST00000357654.9"

    def test_unknown_gene_returns_none(self) -> None:
        """Unknown gene should return None."""
        assert get_mane_select_ensembl_for_gene("NOTREAL") is None


# ---------------------------------------------------------------------------
# Tests for adjust_pvs1_for_mane() — ACGS 2024 §5
# ---------------------------------------------------------------------------


class TestAdjustPvs1ForMane:
    """Tests for PVS1 strength adjustment per ACGS 2024 §5."""

    # When is_mane=True, no adjustment should occur
    def test_pvs1_mane_select_unchanged(self) -> None:
        """PVS1 on MANE Select transcript should remain PVS1."""
        assert adjust_pvs1_for_mane("PVS1", is_mane=True) == "PVS1"

    def test_ps1_mane_select_unchanged(self) -> None:
        """PS1 on MANE Select transcript should remain PS1."""
        assert adjust_pvs1_for_mane("PS1", is_mane=True) == "PS1"

    def test_pm1_mane_select_unchanged(self) -> None:
        """PM1 on MANE Select transcript should remain PM1."""
        assert adjust_pvs1_for_mane("PM1", is_mane=True) == "PM1"

    def test_pp1_mane_select_unchanged(self) -> None:
        """PP1 on MANE Select transcript should remain PP1."""
        assert adjust_pvs1_for_mane("PP1", is_mane=True) == "PP1"

    # When is_mane=False, strength should be reduced by one level
    def test_pvs1_non_mane_downgraded_to_ps1(self) -> None:
        """PVS1 on non-MANE Select → PS1 (ACGS 2024 §5)."""
        assert adjust_pvs1_for_mane("PVS1", is_mane=False) == "PS1"

    def test_ps1_non_mane_downgraded_to_pm1(self) -> None:
        """PS1 on non-MANE Select → PM1 (ACGS 2024 §5)."""
        assert adjust_pvs1_for_mane("PS1", is_mane=False) == "PM1"

    def test_pm1_non_mane_downgraded_to_pp1(self) -> None:
        """PM1 on non-MANE Select → PP1 (ACGS 2024 §5)."""
        assert adjust_pvs1_for_mane("PM1", is_mane=False) == "PP1"

    def test_pp1_non_mane_no_contribution(self) -> None:
        """PP1 on non-MANE Select → no_contribution (minimum strength)."""
        assert adjust_pvs1_for_mane("PP1", is_mane=False) == "no_contribution"

    def test_very_strong_downgraded_to_strong(self) -> None:
        """Lowercase 'very_strong' on non-MANE → 'strong'."""
        assert adjust_pvs1_for_mane("very_strong", is_mane=False) == "strong"

    def test_supporting_downgraded_to_no_contribution(self) -> None:
        """Lowercase 'supporting' on non-MANE → 'no_contribution'."""
        assert adjust_pvs1_for_mane("supporting", is_mane=False) == "no_contribution"

    def test_unknown_strength_returned_unchanged(self) -> None:
        """Unknown strength label should be returned unchanged with a warning."""
        result = adjust_pvs1_for_mane("UNKNOWN_STRENGTH", is_mane=False)
        assert result == "UNKNOWN_STRENGTH"

    def test_double_downgrade_requires_two_calls(self) -> None:
        """Two non-MANE adjustments should require two separate calls."""
        first = adjust_pvs1_for_mane("PVS1", is_mane=False)
        second = adjust_pvs1_for_mane(first, is_mane=False)
        assert first == "PS1"
        assert second == "PM1"
