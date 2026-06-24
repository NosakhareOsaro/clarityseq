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

import pytest

from bayesacmg.models import EvidenceStrength, VariantInput, VariantType
from bayesacmg.rules.splicing import rule_splicing_pp3_bp4_bp7


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
            spliceai_delta=0.03,   # < 0.1
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
            spliceai_delta=0.15,   # 0.1 ≤ Δ < 0.2 = inconclusive
            pangolin_score=0.12,
        )
        assert result.applies is False, (
            "SpliceAI in 0.1–0.2 range is inconclusive — no PP3/BP4/BP7 evidence"
        )

    def test_pangolin_spliceai_disagreement_uses_conservative(
        self, canonical_splice_donor: VariantInput
    ) -> None:
        """When Pangolin and SpliceAI disagree, use the more conservative score.

        Walker 2023: when tools disagree, use the lower Δ score.
        SpliceAI = 0.6 (Strong); Pangolin = 0.25 (Moderate) → use 0.25 → Moderate
        """
        result = rule_splicing_pp3_bp4_bp7(
            variant=canonical_splice_donor,
            spliceai_delta=0.6,    # Strong
            pangolin_score=0.25,   # Moderate — lower
        )
        # Conservative = lower score = 0.25 = Moderate
        assert result.strength == EvidenceStrength.MODERATE, (
            "When SpliceAI and Pangolin disagree, the lower (more conservative) "
            "score should determine the evidence strength per Walker 2023"
        )
        # Disagreement should be documented
        evidence_str = " ".join(result.evidence_items)
        assert "disagree" in evidence_str.lower() or "pangolin" in evidence_str.lower(), (
            "Disagreement between SpliceAI and Pangolin should be documented in evidence_items"
        )
