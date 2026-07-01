"""
Tests for bayesacmg.combinations — ClinGen SVI 2024 novel evidence combinations
and point-to-classification thresholds.

Covers:
    - classify_by_points() threshold boundaries (P / LP / VUS / LB / B)
    - stand_alone_benign override
    - evaluate_all_combinations() empty vs non-empty
    - classify_variant() end-to-end (stand-alone, novel combo, standard points)
    - compute_total_points()

Guidelines:
    Tavtigian et al. 2020 PMID:32645316 — point thresholds.
    Richards et al. 2015 PMID:25741868 — classification categories.
    ClinGen SVI Working Group 2024 — PVS1+PM2_Supporting=LP novel combination.
"""

from __future__ import annotations

from bayesacmg.combinations import (
    classify_by_points,
    classify_variant,
    compute_total_points,
    evaluate_all_combinations,
    evaluate_pvs1_pm2_supporting,
)
from bayesacmg.models import ACMGRule, EvidenceStrength, VariantInput


def _dummy_variant() -> VariantInput:
    return VariantInput(chrom="1", pos=100, ref="A", alt="G", variant_type="snv")


def _rule(rule_id: str, strength: EvidenceStrength, applies: bool = True) -> ACMGRule:
    return ACMGRule(
        rule_id=rule_id,
        strength=strength,
        evidence_items=[f"{rule_id} evidence"],
        applies=applies,
    )


# ---------------------------------------------------------------------------
# classify_by_points
# ---------------------------------------------------------------------------


class TestClassifyByPoints:
    def test_pathogenic_at_threshold(self) -> None:
        """>=10 pts -> Pathogenic (Tavtigian 2020)."""
        assert classify_by_points(10) == "Pathogenic"

    def test_pathogenic_above_threshold(self) -> None:
        assert classify_by_points(15) == "Pathogenic"

    def test_likely_pathogenic_range(self) -> None:
        """6-9 pts -> Likely Pathogenic."""
        assert classify_by_points(6) == "Likely Pathogenic"
        assert classify_by_points(9) == "Likely Pathogenic"

    def test_vus_positive_range(self) -> None:
        """0-5 pts -> VUS."""
        assert classify_by_points(0) == "VUS"
        assert classify_by_points(5) == "VUS"

    def test_vus_negative_range(self) -> None:
        """-1 to -5 pts -> VUS."""
        assert classify_by_points(-1) == "VUS"
        assert classify_by_points(-5) == "VUS"

    def test_likely_benign_range(self) -> None:
        """-6 to -9 pts -> Likely Benign."""
        assert classify_by_points(-6) == "Likely Benign"
        assert classify_by_points(-9) == "Likely Benign"

    def test_benign_at_threshold(self) -> None:
        """<=-10 pts -> Benign."""
        assert classify_by_points(-10) == "Benign"

    def test_benign_below_threshold(self) -> None:
        assert classify_by_points(-20) == "Benign"

    def test_stand_alone_benign_overrides_points(self) -> None:
        """stand_alone_benign=True forces Benign regardless of total_points (e.g. BA1)."""
        assert classify_by_points(50, stand_alone_benign=True) == "Benign"
        assert classify_by_points(-50, stand_alone_benign=True) == "Benign"


# ---------------------------------------------------------------------------
# evaluate_all_combinations
# ---------------------------------------------------------------------------


class TestEvaluateAllCombinations:
    def test_empty_when_no_combination_applies(self) -> None:
        """Rules without PVS1+PM2_Supporting produce an empty combination list."""
        rules = [_rule("PP3", EvidenceStrength.SUPPORTING)]
        results = evaluate_all_combinations(rules)
        assert results == []

    def test_contains_pvs1_pm2_combo_when_present(self) -> None:
        """PVS1 + PM2 Supporting present -> combination list contains the combo."""
        rules = [
            _rule("PVS1", EvidenceStrength.VERY_STRONG),
            _rule("PM2", EvidenceStrength.SUPPORTING),
        ]
        results = evaluate_all_combinations(rules)
        assert len(results) == 1
        assert results[0].combination_name == "PVS1+PM2_Supporting=LP"
        assert results[0].applies is True

    def test_pm2_mito_variant_also_triggers_combo(self) -> None:
        """PM2_MITO (mito-specific PM2) also participates in the novel combination."""
        rules = [
            _rule("PVS1", EvidenceStrength.VERY_STRONG),
            _rule("PM2_MITO", EvidenceStrength.SUPPORTING),
        ]
        results = evaluate_pvs1_pm2_supporting(rules)
        assert results.applies is True
        assert "PM2_MITO" in results.rules_contributing


# ---------------------------------------------------------------------------
# classify_variant — end-to-end
# ---------------------------------------------------------------------------


class TestClassifyVariant:
    def test_stand_alone_benign_short_circuits_to_benign(self) -> None:
        """A stand-alone rule (e.g. BA1) forces Benign regardless of other rules."""
        variant = _dummy_variant()
        rules = [
            _rule("BA1", EvidenceStrength.STAND_ALONE),
            _rule("PVS1", EvidenceStrength.VERY_STRONG),
        ]
        result = classify_variant(rules, variant)
        assert result.classification == "Benign"
        assert result.stand_alone_benign is True
        assert result.novel_combination is None

    def test_novel_combination_applied_for_pvs1_pm2(self) -> None:
        """PVS1 + PM2_Supporting -> Likely Pathogenic with novel_combination recorded."""
        variant = _dummy_variant()
        rules = [
            _rule("PVS1", EvidenceStrength.VERY_STRONG),
            _rule("PM2", EvidenceStrength.SUPPORTING),
        ]
        result = classify_variant(rules, variant)
        assert result.classification == "Likely Pathogenic"
        assert result.novel_combination == "PVS1+PM2_Supporting=LP"
        assert result.total_points == 9
        assert result.stand_alone_benign is False

    def test_standard_points_no_combination(self) -> None:
        """Simple point-sum classification when no novel combination applies."""
        variant = _dummy_variant()
        rules = [
            _rule("PP3", EvidenceStrength.SUPPORTING),
            _rule("PP4", EvidenceStrength.SUPPORTING),
        ]
        result = classify_variant(rules, variant)
        assert result.total_points == 2
        assert result.classification == "VUS"
        assert result.novel_combination is None

    def test_not_applied_rules_are_separated(self) -> None:
        """rules_applied / rules_not_applied correctly partition on .applies."""
        variant = _dummy_variant()
        rules = [
            _rule("PP3", EvidenceStrength.SUPPORTING, applies=True),
            _rule("PM1", EvidenceStrength.MODERATE, applies=False),
        ]
        result = classify_variant(rules, variant)
        assert len(result.rules_applied) == 1
        assert result.rules_applied[0].rule_id == "PP3"
        assert len(result.rules_not_applied) == 1
        assert result.rules_not_applied[0].rule_id == "PM1"

    def test_variant_preserved_on_result(self) -> None:
        """The returned ClassificationResult carries through the original variant."""
        variant = _dummy_variant()
        result = classify_variant([], variant)
        assert result.variant is variant


# ---------------------------------------------------------------------------
# compute_total_points
# ---------------------------------------------------------------------------


class TestComputeTotalPoints:
    def test_sums_only_applied_rules(self) -> None:
        rules = [
            _rule("PVS1", EvidenceStrength.VERY_STRONG, applies=True),  # +8
            _rule("PM1", EvidenceStrength.MODERATE, applies=False),  # 0 (doesn't apply)
            _rule("PP3", EvidenceStrength.SUPPORTING, applies=True),  # +1
        ]
        assert compute_total_points(rules) == 9

    def test_zero_for_empty_rules(self) -> None:
        assert compute_total_points([]) == 0

    def test_negative_total_for_benign_rules(self) -> None:
        rules = [
            _rule("BS1", EvidenceStrength.STRONG_BENIGN, applies=True),  # -4
            _rule("BP4", EvidenceStrength.SUPPORTING_BENIGN, applies=True),  # -1
        ]
        assert compute_total_points(rules) == -5
