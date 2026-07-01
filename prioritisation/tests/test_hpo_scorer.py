"""
prioritisation.tests.test_hpo_scorer
======================================
pytest tests for HPO phenotype scoring of candidate genes.

Tests cover:
    - jaccard_score: set overlap computation.
    - score_genes_by_hpo: ranking of genes by HPO similarity.
    - get_top_genes: threshold filtering and top-N selection.
    - load_hpo_gene_annotations: file loading and parsing.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from prioritisation.hpo_scorer import (
    GeneHPOScore,
    get_top_genes,
    jaccard_score,
    load_hpo_gene_annotations,
    score_genes_by_hpo,
)


# ---------------------------------------------------------------------------
# jaccard_score tests
# ---------------------------------------------------------------------------


class TestJaccardScore:
    """Tests for jaccard_score() set-overlap function."""

    def test_identical_sets_returns_one(self) -> None:
        """Identical sets yield Jaccard index of 1.0."""
        s = {"HP:0001250", "HP:0004322"}
        assert jaccard_score(s, s) == pytest.approx(1.0)

    def test_disjoint_sets_returns_zero(self) -> None:
        """Disjoint sets yield Jaccard index of 0.0."""
        a = {"HP:0001250"}
        b = {"HP:0004322"}
        assert jaccard_score(a, b) == pytest.approx(0.0)

    def test_partial_overlap(self) -> None:
        """Partial overlap: intersection/union = 1/3."""
        a = {"HP:0001250", "HP:0004322"}
        b = {"HP:0001250", "HP:0000823"}
        # intersection={HP:0001250}, union={HP:0001250, HP:0004322, HP:0000823}
        assert jaccard_score(a, b) == pytest.approx(1 / 3)

    def test_empty_sets_returns_zero(self) -> None:
        """Both empty sets return 0.0 (not division by zero)."""
        assert jaccard_score(set(), set()) == pytest.approx(0.0)

    def test_one_empty_set_returns_zero(self) -> None:
        """One empty set yields 0.0 (nothing in common)."""
        assert jaccard_score({"HP:0001250"}, set()) == pytest.approx(0.0)

    def test_superset_relationship(self) -> None:
        """Patient set is a subset of gene terms."""
        patient = {"HP:0001250"}
        gene = {"HP:0001250", "HP:0004322", "HP:0000823"}
        # intersection=1, union=3 → 1/3
        assert jaccard_score(patient, gene) == pytest.approx(1 / 3)


# ---------------------------------------------------------------------------
# score_genes_by_hpo tests
# ---------------------------------------------------------------------------


class TestScoreGenesByHPO:
    """Tests for score_genes_by_hpo() gene ranking function."""

    def _make_gene_map(self) -> dict[str, set[str]]:
        """Return a minimal gene→HPO-terms map for testing."""
        return {
            "BRCA1": {"HP:0001250", "HP:0004322", "HP:0000823"},
            "SCN1A": {"HP:0001250", "HP:0004322"},
            "CFTR": {"HP:0000998", "HP:0002099"},
            "TP53": {"HP:0001250"},
        }

    def test_returns_sorted_by_score_descending(self) -> None:
        """Genes are returned sorted by similarity score, highest first."""
        patient_terms = ["HP:0001250", "HP:0004322"]
        gene_map = self._make_gene_map()
        results = score_genes_by_hpo(patient_terms, gene_map)
        scores = [r.score for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_only_genes_with_score_gt_zero(self) -> None:
        """Genes with zero HPO overlap are excluded."""
        patient_terms = ["HP:0001250"]
        gene_map = self._make_gene_map()
        results = score_genes_by_hpo(patient_terms, gene_map)
        # CFTR has no overlap with patient → should be excluded
        gene_symbols = [r.gene_symbol for r in results]
        assert "CFTR" not in gene_symbols

    def test_empty_patient_terms_returns_empty(self) -> None:
        """Empty patient HPO terms → empty result list."""
        results = score_genes_by_hpo([], {"BRCA1": {"HP:0001250"}})
        assert results == []

    def test_matched_terms_populated(self) -> None:
        """matched_terms contains the intersecting HPO terms."""
        patient_terms = ["HP:0001250", "HP:0004322"]
        gene_map = {"BRCA1": {"HP:0001250", "HP:0004322", "HP:0000823"}}
        results = score_genes_by_hpo(patient_terms, gene_map)
        assert len(results) == 1
        r = results[0]
        assert set(r.matched_terms) == {"HP:0001250", "HP:0004322"}

    def test_result_contains_gene_symbol(self) -> None:
        """GeneHPOScore objects contain the gene_symbol."""
        results = score_genes_by_hpo(
            ["HP:0001250"],
            {"SCN1A": {"HP:0001250"}},
        )
        assert len(results) == 1
        assert results[0].gene_symbol == "SCN1A"

    def test_result_method_field(self) -> None:
        """Default method is 'jaccard'."""
        results = score_genes_by_hpo(
            ["HP:0001250"],
            {"SCN1A": {"HP:0001250"}},
        )
        assert results[0].method == "jaccard"

    def test_bma_method_falls_back_to_jaccard(self) -> None:
        """Unsupported 'bma' method falls back to Jaccard silently."""
        results_jaccard = score_genes_by_hpo(
            ["HP:0001250"], {"SCN1A": {"HP:0001250"}}, method="jaccard"
        )
        results_bma = score_genes_by_hpo(
            ["HP:0001250"], {"SCN1A": {"HP:0001250"}}, method="bma"
        )
        # Scores should be equal since BMA falls back to Jaccard
        assert results_jaccard[0].score == pytest.approx(results_bma[0].score)

    def test_multiple_genes_ranked_correctly(self) -> None:
        """Gene with more HPO overlap ranks higher."""
        patient_terms = ["HP:0001250", "HP:0004322"]
        gene_map = {
            "HIGH": {"HP:0001250", "HP:0004322"},  # exact match → score=1.0
            "LOW": {"HP:0001250"},                  # partial → lower score
        }
        results = score_genes_by_hpo(patient_terms, gene_map)
        assert results[0].gene_symbol == "HIGH"
        assert results[0].score == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# get_top_genes tests
# ---------------------------------------------------------------------------


class TestGetTopGenes:
    """Tests for get_top_genes() filtering function."""

    def _make_gene_map(self) -> dict[str, set[str]]:
        return {
            "BRCA1": {"HP:0001250", "HP:0004322"},
            "SCN1A": {"HP:0001250"},
            "CFTR": {"HP:0000998"},  # no overlap
        }

    def test_returns_at_most_top_n(self) -> None:
        """Result list does not exceed top_n genes."""
        patient = ["HP:0001250", "HP:0004322"]
        results = get_top_genes(patient, self._make_gene_map(), top_n=1)
        assert len(results) <= 1

    def test_min_score_filter_applied(self) -> None:
        """Genes below min_score are excluded."""
        patient = ["HP:0001250", "HP:0004322"]
        gene_map = {
            "HIGH": {"HP:0001250", "HP:0004322"},  # score=1.0
            "LOW": {"HP:0001250", "X", "Y", "Z", "W", "V"},  # low score
        }
        results = get_top_genes(patient, gene_map, min_score=0.9)
        symbols = [r.gene_symbol for r in results]
        assert "HIGH" in symbols

    def test_empty_patient_terms(self) -> None:
        """Empty patient terms → empty result."""
        results = get_top_genes([], self._make_gene_map())
        assert results == []


# ---------------------------------------------------------------------------
# load_hpo_gene_annotations tests
# ---------------------------------------------------------------------------


class TestLoadHPOGeneAnnotations:
    """Tests for load_hpo_gene_annotations() file loader."""

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        """Missing annotations file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="HPO annotations not found"):
            load_hpo_gene_annotations(tmp_path / "nonexistent.txt")

    def test_parses_tsv_format(self, tmp_path: Path) -> None:
        """Parse the HPO phenotype_to_genes.txt format correctly."""
        annotations_file = tmp_path / "phenotype_to_genes.txt"
        annotations_file.write_text(
            "# Header comment\n"
            "HP:0001250\tSeizure\t1131\tBRCA1\n"
            "HP:0004322\tShort stature\t1131\tBRCA1\n"
            "HP:0001250\tSeizure\t6323\tSCN1A\n",
            encoding="utf-8",
        )
        gene_map = load_hpo_gene_annotations(annotations_file)
        assert "BRCA1" in gene_map
        assert "HP:0001250" in gene_map["BRCA1"]
        assert "HP:0004322" in gene_map["BRCA1"]
        assert "SCN1A" in gene_map
        assert "HP:0001250" in gene_map["SCN1A"]

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        """Lines starting with # are skipped."""
        f = tmp_path / "phenotype_to_genes.txt"
        f.write_text("# Comment\n# Another comment\n", encoding="utf-8")
        gene_map = load_hpo_gene_annotations(f)
        assert len(gene_map) == 0

    def test_gene_symbols_uppercased(self, tmp_path: Path) -> None:
        """Gene symbols are uppercased for consistent lookup."""
        f = tmp_path / "phenotype_to_genes.txt"
        f.write_text("HP:0001250\tSeizure\t1131\tbrca1\n", encoding="utf-8")
        gene_map = load_hpo_gene_annotations(f)
        assert "BRCA1" in gene_map

    def test_skips_lines_with_too_few_fields(self, tmp_path: Path) -> None:
        """Lines with fewer than 4 tab-separated fields are skipped."""
        f = tmp_path / "phenotype_to_genes.txt"
        f.write_text(
            "HP:0001250\tSeizure\t1131\n"  # malformed: only 3 fields
            "HP:0004322\tShort stature\t1131\tBRCA1\n",
            encoding="utf-8",
        )
        gene_map = load_hpo_gene_annotations(f)
        assert list(gene_map.keys()) == ["BRCA1"], (
            "Malformed line with <4 fields should be skipped, not raise or "
            "be added as a bogus entry"
        )
        assert gene_map["BRCA1"] == {"HP:0004322"}
