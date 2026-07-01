"""
Tests for the ClinGen SVI Splicing Subgroup evidence framework.

Tests all SpliceAI Δ score branches from Walker et al. 2023 PMID:36898414:
    ≥ 0.5  → PP3 Strong
    ≥ 0.2  → PP3 Moderate
    < 0.1  + synonymous → BP7
    0.1–0.2 → inconclusive (no evidence applied)

Guidelines:
    Walker et al. 2023 Am J Hum Genet 110:1046 PMID:36898414
    ACGS 2024 v1.2 §5 (splicing evidence)
"""

from __future__ import annotations


from bayesacmg.models import EvidenceStrength, VariantInput, VariantType
from bayesacmg.rules.splicing import _resolve_splice_score, rule_splicing_pp3_bp4_bp7


class TestSplicingPP3BP4BP7:
    """Walker et al. 2023 splicing framework tests."""

    def test_spliceai_strong_threshold_pp3_strong(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """SpliceAI Δ ≥ 0.5 → PP3 Strong.

        Variant: BRCA1 c.5277+1G>A (RCV000048342) — canonical splice donor
        SpliceAI Δ = 0.95 → strong splice impact
        Expected: PP3 Strong
        """
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.95,
            pangolin_score=0.92,
        )
        assert result.applies is True
        assert result.rule_id == "PP3"
        assert result.strength == EvidenceStrength.STRONG

    def test_spliceai_exact_strong_boundary(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """SpliceAI Δ exactly 0.5 → PP3 Strong (at boundary)."""
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.5,
            pangolin_score=None,
        )
        assert result.applies is True
        assert result.strength == EvidenceStrength.STRONG

    def test_spliceai_moderate_threshold_pp3_moderate(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """SpliceAI Δ ≥ 0.2 but < 0.5 → PP3 Moderate."""
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.35,
            pangolin_score=0.31,
        )
        assert result.applies is True
        assert result.rule_id == "PP3"
        assert result.strength == EvidenceStrength.MODERATE

    def test_spliceai_synonymous_no_splice_impact_bp7(
        self, synonymous_no_splice_impact: VariantInput
    ) -> None:
        """Synonymous variant + SpliceAI < 0.1 → BP7 (Walker 2023).

        Walker et al. 2023: synonymous variants with SpliceAI < 0.1
        do not affect splicing → BP7 (benign supporting).
        """
        result = rule_splicing_pp3_bp4_bp7(
            variant=synonymous_no_splice_impact,
            spliceai_delta=0.03,  # < 0.1
            pangolin_score=0.02,
        )
        assert result.applies is True
        assert result.rule_id == "BP7"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_spliceai_below_strong_not_synonymous_inconclusive(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """SpliceAI 0.1–0.2 for non-synonymous variant → inconclusive (no evidence)."""
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.15,  # 0.1 ≤ Δ < 0.2 = inconclusive
            pangolin_score=0.12,
        )
        assert (
            result.applies is False
        ), "SpliceAI in 0.1–0.2 range is inconclusive — no PP3/BP4/BP7 evidence"

    def test_pangolin_spliceai_disagreement_uses_conservative(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """When Pangolin and SpliceAI disagree, use the more conservative score.

        Walker 2023: when tools disagree, use the lower Δ score.
        SpliceAI = 0.6 (Strong); Pangolin = 0.25 (Moderate) → use 0.25 → Moderate
        """
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.6,  # Strong
            pangolin_score=0.25,  # Moderate — lower
        )
        # Conservative = lower score = 0.25 = Moderate
        assert result.strength == EvidenceStrength.MODERATE, (
            "When SpliceAI and Pangolin disagree, the lower (more conservative) "
            "score should determine the evidence strength per Walker 2023"
        )
        # Disagreement should be documented
        evidence_str = " ".join(result.evidence_items)
        assert (
            "disagree" in evidence_str.lower() or "pangolin" in evidence_str.lower()
        ), "Disagreement between SpliceAI and Pangolin should be documented in evidence_items"

    def test_no_scores_available_at_all(self) -> None:
        """No SpliceAI or Pangolin score anywhere (args or variant) → applies=False."""
        variant = VariantInput(
            chrom="17",
            pos=41_215_918,
            ref="G",
            alt="A",
            variant_type=VariantType.SPLICE_SITE,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            spliceai_delta=None,
            spliceai_max_delta=None,
            pangolin_score=None,
        )
        result = rule_splicing_pp3_bp4_bp7(variant)
        assert result.applies is False
        assert result.rule_id == "PP3"
        assert "no spliceai or pangolin score" in " ".join(result.evidence_items).lower()

    def test_non_synonymous_no_splice_impact_gives_bp4(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """Non-synonymous variant + score < 0.1 → BP4 (not BP7, since not synonymous)."""
        # canonical_splice_donor is variant_type=SPLICE_SITE, not synonymous
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.02,
            pangolin_score=0.01,
        )
        assert result.applies is True
        assert result.rule_id == "BP4"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN


class TestResolveSpliceScore:
    """Direct tests for the internal _resolve_splice_score() helper."""

    def test_both_scores_none(self) -> None:
        """Both SpliceAI and Pangolin absent → (None, 'none', False)."""
        score, tool, disagreement = _resolve_splice_score(None, None)
        assert score is None
        assert tool == "none"
        assert disagreement is False

    def test_only_pangolin_available(self) -> None:
        """SpliceAI absent, Pangolin present → uses Pangolin, no disagreement flag."""
        score, tool, disagreement = _resolve_splice_score(None, 0.42)
        assert score == 0.42
        assert tool == "Pangolin"
        assert disagreement is False

    def test_only_spliceai_available(self) -> None:
        """Pangolin absent, SpliceAI present → uses SpliceAI, no disagreement flag."""
        score, tool, disagreement = _resolve_splice_score(0.6, None)
        assert score == 0.6
        assert tool == "SpliceAI"
        assert disagreement is False

    def test_both_available_spliceai_lower_or_equal_is_conservative(self) -> None:
        """When SpliceAI <= Pangolin, SpliceAI (the lower score) is chosen as conservative."""
        score, tool, disagreement = _resolve_splice_score(0.3, 0.3)
        assert score == 0.3
        assert tool == "SpliceAI (conservative)"
        # Same categorical bucket (both moderate_impact) → no disagreement
        assert disagreement is False

    def test_both_available_pangolin_lower_is_conservative(self) -> None:
        """When Pangolin < SpliceAI, Pangolin (the lower score) is chosen as conservative."""
        score, tool, disagreement = _resolve_splice_score(0.6, 0.25)
        assert score == 0.25
        assert tool == "Pangolin (conservative)"
        assert disagreement is True  # strong_impact vs moderate_impact
