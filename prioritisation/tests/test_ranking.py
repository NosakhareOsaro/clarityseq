"""
prioritisation.tests.test_ranking
====================================
pytest tests for composite variant ranking.

Tests cover:
    - compute_composite_score: ACMG/HPO/inheritance/panel weighting.
    - rank_variants: sorting and rank assignment.
    - _normalise_acmg: classification string normalisation.
"""

from __future__ import annotations

import pytest

from prioritisation.ranking import (
    RankedVariant,
    compute_composite_score,
    rank_variants,
)


# ---------------------------------------------------------------------------
# compute_composite_score tests
# ---------------------------------------------------------------------------


class TestCompositeScore:
    """Tests for compute_composite_score()."""

    def test_pathogenic_in_panel_passes_inheritance(self) -> None:
        """P variant in panel with inheritance pass → high composite score."""
        composite, acmg_s, hpo_s, inh_b, panel_b = compute_composite_score(
            acmg_class="Pathogenic",
            hpo_score=0.8,
            passes_inheritance=True,
            in_panel=True,
        )
        assert composite == pytest.approx(
            0.4 * 1.0 + 0.3 * 0.8 + 0.2 * 1.0 + 0.1 * 1.0
        )
        assert acmg_s == pytest.approx(1.0)
        assert inh_b == pytest.approx(1.0)
        assert panel_b == pytest.approx(1.0)

    def test_benign_not_in_panel_no_inheritance(self) -> None:
        """B variant not in panel, inheritance fail → low composite score."""
        composite, acmg_s, hpo_s, inh_b, panel_b = compute_composite_score(
            acmg_class="Benign",
            hpo_score=0.0,
            passes_inheritance=False,
            in_panel=False,
        )
        assert composite == pytest.approx(0.0)
        assert acmg_s == pytest.approx(0.0)
        assert inh_b == pytest.approx(0.0)
        assert panel_b == pytest.approx(0.0)

    def test_vus_default_score(self) -> None:
        """VUS maps to ACMG score 0.5."""
        _, acmg_s, _, _, _ = compute_composite_score(
            acmg_class="VUS", hpo_score=0.0,
            passes_inheritance=False, in_panel=False,
        )
        assert acmg_s == pytest.approx(0.5)

    def test_likely_pathogenic_score(self) -> None:
        """LP maps to ACMG score 0.8."""
        _, acmg_s, _, _, _ = compute_composite_score(
            acmg_class="Likely_Pathogenic", hpo_score=0.0,
            passes_inheritance=False, in_panel=False,
        )
        assert acmg_s == pytest.approx(0.8)

    def test_normalisation_alias_p(self) -> None:
        """Short alias 'P' is normalised to Pathogenic score 1.0."""
        _, acmg_s, _, _, _ = compute_composite_score(
            acmg_class="P", hpo_score=0.0,
            passes_inheritance=False, in_panel=False,
        )
        assert acmg_s == pytest.approx(1.0)

    def test_normalisation_alias_lp(self) -> None:
        """Short alias 'LP' is normalised to Likely_Pathogenic score 0.8."""
        _, acmg_s, _, _, _ = compute_composite_score(
            acmg_class="LP", hpo_score=0.0,
            passes_inheritance=False, in_panel=False,
        )
        assert acmg_s == pytest.approx(0.8)

    def test_unknown_class_defaults_to_vus(self) -> None:
        """Unknown classification defaults to VUS score 0.5."""
        _, acmg_s, _, _, _ = compute_composite_score(
            acmg_class="UNKNOWN_CLASS", hpo_score=0.0,
            passes_inheritance=False, in_panel=False,
        )
        assert acmg_s == pytest.approx(0.5)

    def test_custom_weights(self) -> None:
        """Custom weights are applied correctly."""
        composite, _, _, _, _ = compute_composite_score(
            acmg_class="Pathogenic",
            hpo_score=1.0,
            passes_inheritance=True,
            in_panel=True,
            acmg_weight=0.25,
            hpo_weight=0.25,
            inheritance_weight=0.25,
            panel_weight=0.25,
        )
        # All scores=1.0, weights sum to 1.0 → composite=1.0
        assert composite == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# rank_variants tests
# ---------------------------------------------------------------------------


class TestRankVariants:
    """Tests for rank_variants() main ranking function."""

    def _make_variant(
        self,
        gene: str = "BRCA1",
        acmg: str = "VUS",
        chrom: str = "chr17",
        pos: int = 43044295,
    ) -> dict:
        return {
            "gene_symbol": gene,
            "acmg_class": acmg,
            "chrom": chrom,
            "pos": pos,
            "ref": "G",
            "alt": "A",
        }

    def test_returns_ranked_variants_list(self) -> None:
        """rank_variants returns a list of RankedVariant objects."""
        variants = [self._make_variant()]
        results = rank_variants(variants)
        assert len(results) == 1
        assert isinstance(results[0], RankedVariant)

    def test_rank_assigned_starting_from_one(self) -> None:
        """First ranked variant has rank=1."""
        variants = [self._make_variant(), self._make_variant("SCN1A", "VUS")]
        results = rank_variants(variants)
        ranks = [r.rank for r in results]
        assert 1 in ranks
        assert sorted(ranks) == list(range(1, len(ranks) + 1))

    def test_pathogenic_ranked_above_vus(self) -> None:
        """Pathogenic variant ranks higher than VUS with same HPO/inheritance."""
        variants = [
            self._make_variant("SCN1A", "VUS"),
            self._make_variant("BRCA1", "Pathogenic"),
        ]
        results = rank_variants(variants)
        assert results[0].gene_symbol == "BRCA1"
        assert results[0].acmg_class == "Pathogenic"

    def test_hpo_scores_applied(self) -> None:
        """Gene with higher HPO score ranks higher when ACMG is equal."""
        variants = [
            self._make_variant("BRCA1", "VUS"),
            self._make_variant("SCN1A", "VUS"),
        ]
        hpo_scores = {"BRCA1": 0.9, "SCN1A": 0.1}
        results = rank_variants(variants, hpo_gene_scores=hpo_scores)
        assert results[0].gene_symbol == "BRCA1"

    def test_panel_genes_ranked_higher(self) -> None:
        """Gene on clinical panel ranks higher than off-panel gene with same ACMG."""
        variants = [
            self._make_variant("BRCA1", "VUS"),
            self._make_variant("UNKNOWN_GENE", "VUS"),
        ]
        panel = {"BRCA1"}
        results = rank_variants(variants, panel_genes=panel)
        assert results[0].gene_symbol == "BRCA1"

    def test_inheritance_filter_bonus(self) -> None:
        """Genes passing inheritance filter rank higher."""
        variants = [
            self._make_variant("GENE_PASS", "VUS"),
            self._make_variant("GENE_FAIL", "VUS"),
        ]
        passing = {"GENE_PASS"}
        results = rank_variants(variants, passing_inheritance_genes=passing)
        assert results[0].gene_symbol == "GENE_PASS"

    def test_empty_variants_returns_empty(self) -> None:
        """Empty input returns empty list."""
        assert rank_variants([]) == []

    def test_evidence_summary_populated(self) -> None:
        """Evidence summary string is non-empty and contains ACMG info."""
        variants = [self._make_variant("BRCA1", "Pathogenic")]
        results = rank_variants(variants)
        assert "ACMG" in results[0].evidence_summary

    def test_composite_score_field_populated(self) -> None:
        """composite_score field is a float between 0 and 1."""
        variants = [self._make_variant()]
        results = rank_variants(variants)
        assert 0.0 <= results[0].composite_score <= 1.0

    def test_extra_keys_stored(self) -> None:
        """Extra variant keys (not standard fields) are stored in extra dict."""
        var = self._make_variant()
        var["custom_field"] = "custom_value"
        results = rank_variants([var])
        assert results[0].extra.get("custom_field") == "custom_value"
