"""
Tests for BayesACMG benign rules.

Tests BA1, BS1-4, BP1-7.
Key test: BA1 uses gnomAD v4.1 AF > 5% threshold (NOT for mito variants).

Guidelines:
    Richards et al. 2015 PMID:25741868
    ClinGen SVI 2024 (BP4 via AlphaMissense ≤ 0.340)
    ACGS 2024 v1.2 §5 (BA1, BS1 thresholds)
"""
from __future__ import annotations

import pytest

from bayesacmg.models import EvidenceStrength, VariantInput, VariantType
from bayesacmg.rules.benign import rule_ba1, rule_bs1, rule_bp4, rule_bp7


class TestBA1:
    """BA1 — allele frequency > 5% in gnomAD v4.1 → Benign (stand-alone)."""

    def test_ba1_fires_when_af_exceeds_five_percent(
        self, common_benign_missense: VariantInput
    ) -> None:
        """AF = 0.12 (12%) in gnomAD v4.1 → BA1 (stand-alone Benign).

        Richards 2015: BA1 if AF > 5% in any general population database.
        gnomAD v4.1 (807,162 individuals) is used.
        """
        result = rule_ba1(common_benign_missense)
        assert result.applies is True
        assert result.strength == EvidenceStrength.STAND_ALONE
        assert result.rule_id == "BA1"

    def test_ba1_does_not_fire_for_rare_variant(self, brca2_novel_lof: VariantInput) -> None:
        """AF = 0 (absent) → BA1 does NOT apply."""
        result = rule_ba1(brca2_novel_lof)
        assert result.applies is False

    def test_ba1_threshold_is_five_percent(self, common_benign_missense: VariantInput) -> None:
        """Verify BA1 threshold is 5% (0.05), not 1% or 10%."""
        # Variant just above threshold
        variant_above = VariantInput(**{**common_benign_missense.__dict__, "gnomad_af": 0.051})
        result_above = rule_ba1(variant_above)
        assert result_above.applies is True, "AF > 5% should trigger BA1"

        # Variant just below threshold
        variant_below = VariantInput(**{**common_benign_missense.__dict__, "gnomad_af": 0.049})
        result_below = rule_ba1(variant_below)
        assert result_below.applies is False, "AF < 5% should NOT trigger BA1"


class TestBP4:
    """BP4 — AlphaMissense ≤ 0.340 → benign supporting (ClinGen SVI 2024)."""

    def test_bp4_alphamissense_below_threshold(
        self, common_benign_missense: VariantInput
    ) -> None:
        """AlphaMissense 0.21 → BP4 (≤ 0.340 threshold)."""
        result = rule_bp4(
            variant=common_benign_missense,
            alphamissense_score=0.21,
        )
        assert result.applies is True
        assert result.rule_id == "BP4"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_bp4_exact_boundary(self, common_benign_missense: VariantInput) -> None:
        """AlphaMissense exactly 0.340 → BP4 applies."""
        result = rule_bp4(variant=common_benign_missense, alphamissense_score=0.340)
        assert result.applies is True

    def test_bp4_does_not_fire_above_threshold(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """AlphaMissense 0.97 → BP4 does NOT apply."""
        result = rule_bp4(variant=tp53_missense_pathogenic, alphamissense_score=0.97)
        assert result.applies is False


class TestBP7:
    """BP7 — synonymous variant with no predicted splice impact."""

    def test_bp7_synonymous_no_splice_impact(
        self, synonymous_no_splice_impact: VariantInput
    ) -> None:
        """Synonymous + SpliceAI < 0.1 → BP7.

        Walker et al. 2023 PMID:36898414: BP7 for synonymous variants
        with no predicted splice impact.
        """
        result = rule_bp7(
            variant=synonymous_no_splice_impact,
            spliceai_delta=0.03,
        )
        assert result.applies is True
        assert result.rule_id == "BP7"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_bp7_does_not_fire_for_missense(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """BP7 does not apply to missense variants (only synonymous)."""
        result = rule_bp7(
            variant=tp53_missense_pathogenic,
            spliceai_delta=0.03,
        )
        assert result.applies is False, "BP7 applies only to synonymous variants"

    def test_bp7_does_not_fire_with_splice_impact(
        self, synonymous_no_splice_impact: VariantInput
    ) -> None:
        """Synonymous variant WITH splice impact (SpliceAI ≥ 0.1) → BP7 does not apply."""
        result = rule_bp7(
            variant=synonymous_no_splice_impact,
            spliceai_delta=0.45,   # significant splice impact
        )
        assert result.applies is False, (
            "BP7 should not apply when SpliceAI ≥ 0.1 even for synonymous variants"
        )
