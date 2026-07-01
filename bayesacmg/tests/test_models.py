"""
Tests for bayesacmg.models core data structures.

Covers dataclass behaviour not otherwise exercised through rule modules:
    - ACMGRule.points when applies=False (must return 0)
    - ClassificationResult.posterior_probability (Tavtigian 2020 OddsP formula)
    - ClassificationResult.credible_interval_95 default fallback
    - ClassificationResult.point_summary formatting

Guidelines:
    Tavtigian et al. 2020 PMID:32645316 — point-score / posterior formula.
    Richards et al. 2015 PMID:25741868 — prior probability 0.1.
"""

from __future__ import annotations

from bayesacmg.models import (
    ACMGRule,
    ClassificationResult,
    EvidenceStrength,
    VariantInput,
)


def _dummy_variant() -> VariantInput:
    return VariantInput(chrom="1", pos=100, ref="A", alt="G", variant_type="snv")


class TestACMGRulePoints:
    """ACMGRule.points property."""

    def test_points_zero_when_not_applies(self) -> None:
        """A rule that does not apply must contribute 0 points regardless of strength."""
        rule = ACMGRule(
            rule_id="PM1",
            strength=EvidenceStrength.MODERATE,
            evidence_items=["does not apply"],
            applies=False,
        )
        assert rule.points == 0

    def test_points_zero_when_stand_alone_and_not_applies(self) -> None:
        """STAND_ALONE rule that doesn't apply also returns 0 (not just because STAND_ALONE=0)."""
        rule = ACMGRule(
            rule_id="BA1",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=[],
            applies=False,
        )
        assert rule.points == 0

    def test_points_nonzero_when_applies(self) -> None:
        """A rule that applies returns its strength's point value."""
        rule = ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            evidence_items=["LoF"],
            applies=True,
        )
        assert rule.points == 8


class TestClassificationResultPosteriorProbability:
    """ClassificationResult.posterior_probability — Tavtigian 2020 OddsP formula."""

    def test_zero_points_returns_prior(self) -> None:
        """total_points <= 0 returns the prior probability (0.1) directly."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="VUS",
            total_points=0,
            rules_applied=[],
            rules_not_applied=[],
        )
        assert result.posterior_probability == 0.1

    def test_negative_points_returns_prior(self) -> None:
        """Negative total_points is clamped to 0 and returns the prior."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="Likely_Benign",
            total_points=-7,
            rules_applied=[],
            rules_not_applied=[],
        )
        assert result.posterior_probability == 0.1

    def test_positive_points_uses_odds_formula(self) -> None:
        """total_points=8 (Very Strong) matches the Tavtigian 2020 OddsP formula.

        OddsP = 350^(pts/8); posterior = (prior*odds)/(prior*odds + (1-prior))
        At pts=8: OddsP = 350^1 = 350; prior=0.1
        posterior = (0.1*350)/(0.1*350+0.9) = 35/35.9
        """
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="Likely_Pathogenic",
            total_points=8,
            rules_applied=[],
            rules_not_applied=[],
        )
        expected = (0.1 * 350.0) / (0.1 * 350.0 + 0.9)
        assert abs(result.posterior_probability - expected) < 1e-9

    def test_higher_points_gives_higher_posterior(self) -> None:
        """Monotonic: more points → higher posterior probability of pathogenicity."""
        low = ClassificationResult(
            variant=_dummy_variant(),
            classification="VUS",
            total_points=2,
            rules_applied=[],
            rules_not_applied=[],
        )
        high = ClassificationResult(
            variant=_dummy_variant(),
            classification="Pathogenic",
            total_points=12,
            rules_applied=[],
            rules_not_applied=[],
        )
        assert high.posterior_probability > low.posterior_probability


class TestClassificationResultCredibleInterval95:
    """ClassificationResult.credible_interval_95 — stored bounds or fallback."""

    def test_returns_default_when_bounds_none(self) -> None:
        """When lower/upper are not set, returns the uninformative (0.0, 1.0) fallback."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="VUS",
            total_points=3,
            rules_applied=[],
            rules_not_applied=[],
        )
        assert result.credible_interval_95 == (0.0, 1.0)

    def test_returns_stored_bounds_when_present(self) -> None:
        """When lower/upper are set, returns them as a tuple."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="Pathogenic",
            total_points=12,
            rules_applied=[],
            rules_not_applied=[],
            credible_interval_lower=0.85,
            credible_interval_upper=0.99,
        )
        assert result.credible_interval_95 == (0.85, 0.99)

    def test_returns_default_when_only_one_bound_present(self) -> None:
        """If only one of lower/upper is set, falls back to (0.0, 1.0)."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="VUS",
            total_points=3,
            rules_applied=[],
            rules_not_applied=[],
            credible_interval_lower=0.5,
            credible_interval_upper=None,
        )
        assert result.credible_interval_95 == (0.0, 1.0)


class TestClassificationResultPointSummary:
    """ClassificationResult.point_summary — human-readable rule/point listing."""

    def test_point_summary_lists_applied_rules_and_total(self) -> None:
        """point_summary lists each applied rule with signed points, then the total."""
        pvs1 = ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            evidence_items=["LoF"],
            applies=True,
        )
        pm2 = ACMGRule(
            rule_id="PM2",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Absent from gnomAD"],
            applies=True,
        )
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="Likely_Pathogenic",
            total_points=9,
            rules_applied=[pvs1, pm2],
            rules_not_applied=[],
        )
        summary = result.point_summary
        assert "PVS1" in summary
        assert "PM2" in summary
        assert "+8 pts" in summary
        assert "+1 pts" in summary
        assert "TOTAL: +9 pts" in summary
        assert "Likely_Pathogenic" in summary

    def test_point_summary_with_no_applied_rules(self) -> None:
        """point_summary with no applied rules still shows the TOTAL line."""
        result = ClassificationResult(
            variant=_dummy_variant(),
            classification="VUS",
            total_points=0,
            rules_applied=[],
            rules_not_applied=[],
        )
        summary = result.point_summary
        assert summary.strip() == "TOTAL: +0 pts → VUS"
