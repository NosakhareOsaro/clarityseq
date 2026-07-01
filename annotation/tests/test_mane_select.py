"""
Tests for MANE Select transcript utility functions.

Validates:
- is_mane_select() with RefSeq and Ensembl IDs
- get_mane_select_for_gene() lookup
- adjust_pvs1_for_mane() strength downgrade logic (ACGS 2024 §5)
"""

from __future__ import annotations

import gzip

from annotation.mane_select import (
    adjust_pvs1_for_mane,
    get_mane_select_ensembl_for_gene,
    get_mane_select_for_gene,
    is_mane_select,
    load_mane_summary,
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


# ---------------------------------------------------------------------------
# Tests for load_mane_summary()
#
# Note: load_mane_summary is decorated with @lru_cache(maxsize=1), so each
# test below must use a distinct summary_path (tmp_path is unique per test)
# to avoid reading a stale cached result from a previous test.
# ---------------------------------------------------------------------------


# Sample rows (tab-separated) matching the documented MANE summary format:
# #NCBI_GeneID Ensembl_Gene GeneSymbol name RefSeq_nuc RefSeq_prot
# Ensembl_nuc Ensembl_prot MANE_status GRCh38_chr chr_start chr_end chr_strand
_MANE_SELECT_ROW = (
    "123\tENSG00000000123\tTESTGENE\tTest gene\tNM_999999.1\tNP_999999.1\t"
    "ENST00000999999.1\tENSP00000999999.1\tMANE Select\tchr1\t100\t200\t+\n"
)
_MANE_PLUS_CLINICAL_ROW = (
    "456\tENSG00000000456\tOTHERGENE\tOther gene\tNM_888888.1\tNP_888888.1\t"
    "ENST00000888888.1\tENSP00000888888.1\tMANE Plus Clinical\tchr2\t300\t400\t-\n"
)
_SHORT_ROW = "789\tENSG1\tSHORTGENE\n"  # Fewer than 9 columns — must be skipped
_COMMENT_LINE = "#NCBI_GeneID\tEnsembl_Gene\tGeneSymbol\n"


class TestLoadManeSummary:
    """Tests for load_mane_summary() file parsing."""

    def test_missing_file_returns_empty_dict(self, tmp_path) -> None:
        """A summary_path that doesn't exist returns an empty registry."""
        missing = tmp_path / "does_not_exist.summary.txt"
        result = load_mane_summary(str(missing))
        assert result == {}

    def test_plain_text_file_parsed(self, tmp_path) -> None:
        """A plain-text (non-gzipped) summary file should be read via open()."""
        summary_file = tmp_path / "MANE.GRCh38.plain.summary.txt"
        summary_file.write_text(
            _COMMENT_LINE + _MANE_SELECT_ROW + _MANE_PLUS_CLINICAL_ROW + _SHORT_ROW
        )

        registry = load_mane_summary(str(summary_file))

        assert registry["TESTGENE"] == ("NM_999999.1", "ENST00000999999.1")

    def test_mane_plus_clinical_not_registered(self, tmp_path) -> None:
        """Only MANE Select (not MANE Plus Clinical) rows are registered."""
        summary_file = tmp_path / "MANE.GRCh38.plusclin.summary.txt"
        summary_file.write_text(
            _COMMENT_LINE + _MANE_SELECT_ROW + _MANE_PLUS_CLINICAL_ROW
        )

        registry = load_mane_summary(str(summary_file))

        assert "OTHERGENE" not in registry
        assert "TESTGENE" in registry

    def test_short_row_skipped(self, tmp_path) -> None:
        """Rows with fewer than 9 tab-separated columns are skipped."""
        summary_file = tmp_path / "MANE.GRCh38.short.summary.txt"
        summary_file.write_text(_COMMENT_LINE + _SHORT_ROW)

        registry = load_mane_summary(str(summary_file))

        assert registry == {}

    def test_gzipped_file_parsed(self, tmp_path) -> None:
        """A .gz summary file should be opened via gzip.open()."""
        summary_file = tmp_path / "MANE.GRCh38.gzipped.summary.txt.gz"
        with gzip.open(summary_file, "wt", encoding="utf-8") as fh:
            fh.write(_COMMENT_LINE)
            fh.write(_MANE_SELECT_ROW)

        registry = load_mane_summary(str(summary_file))

        assert registry["TESTGENE"] == ("NM_999999.1", "ENST00000999999.1")

    def test_oserror_during_read_returns_partial_or_empty_registry(
        self, tmp_path
    ) -> None:
        """An OSError while reading (e.g. path is a directory) is caught and
        a registry (possibly empty) is returned rather than raising."""
        # A directory satisfies path.exists() but cannot be opened as a file,
        # triggering IsADirectoryError (a subclass of OSError) inside the
        # try/except block.
        bad_path = tmp_path / "not_a_file.summary.txt"
        bad_path.mkdir()

        registry = load_mane_summary(str(bad_path))

        assert registry == {}
