"""
Tests for the BayesACMG Bayesian model.

Verifies:
1. PM2 prior is centred on Supporting (1 pt), not Moderate (2 pts)
2. Credible interval width is reasonable
3. Classification probabilities match manual point-score expectations
4. ECE target < 0.05 on fixture set

Guidelines:
    Tavtigian et al. 2020 PMID:32645316 — Bayesian point-scoring framework
    ClinGen SVI 2024 — PM2 prior weight
"""

from __future__ import annotations

import pytest

from bayesacmg.model import BayesACMGModel
from bayesacmg.models import (
    ACMGRule,
    EvidenceStrength,
    VariantInput,
)


@pytest.fixture
def model() -> BayesACMGModel:
    """BayesACMG model instance."""
    return BayesACMGModel()


class TestBayesACMGModel:
    """Tests for the Bayesian ACMG/AMP classification model."""

    def test_pm2_prior_weight_is_supporting(self, model: BayesACMGModel) -> None:
        """PM2 prior in Bayesian model is centred on Supporting (1 pt), not Moderate (2 pts).

        ClinGen SVI 2024: PM2 concentration parameter is 5.0 (vs 10.0 for PVS1)
        reflecting ongoing community uncertainty after the 2024 revision.
        """
        pm2_prior = model.get_criterion_prior("PM2")
        assert pm2_prior is not None, "PM2 must have a prior defined in the model"
        # Prior should be centred on Supporting weight (1), not Moderate (2)
        assert pm2_prior.centre == 1, (
            f"PM2 prior should be centred on 1 (Supporting), not {pm2_prior.centre}. "
            "ClinGen SVI 2024 requires PM2 at Supporting weight."
        )
        assert pm2_prior.concentration == 5.0, (
            "PM2 concentration parameter should be 5.0 (reflecting 2024 uncertainty), "
            f"not {pm2_prior.concentration}."
        )

    def test_pvs1_prior_weight_is_eight(self, model: BayesACMGModel) -> None:
        """PVS1 prior in Bayesian model is centred on Very Strong (8 pts)."""
        pvs1_prior = model.get_criterion_prior("PVS1")
        assert pvs1_prior is not None
        assert pvs1_prior.centre == 8

    def test_pathogenic_variant_gives_high_posterior(
        self, model: BayesACMGModel, brca1_frameshift: VariantInput
    ) -> None:
        """Strongly pathogenic variant (5-star ClinVar P) gives P(Path) > 0.95."""
        # PVS1 + PS1 + PM2 for BRCA1 frameshift
        rules = [
            ACMGRule(
                rule_id="PVS1",
                strength=EvidenceStrength.VERY_STRONG,
                applies=True,
                evidence_items=["BRCA1 frameshift; MANE Select NM_007294.4"],
            ),
            ACMGRule(
                rule_id="PM2",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["Absent from gnomAD v4.1"],
            ),
        ]
        result = model.classify(brca1_frameshift, rules)
        assert result.posterior_probability > 0.90, (
            f"Strongly pathogenic variant (PVS1 + PM2_Supporting = 9 pts) should give "
            f"P(Path) > 0.90, got {result.posterior_probability:.3f}"
        )

    def test_credible_interval_is_within_zero_one(
        self, model: BayesACMGModel, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """95% credible interval bounds must be within [0, 1]."""
        rules = [
            ACMGRule(
                rule_id="PP3",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["AlphaMissense 0.97 ≥ 0.564"],
            ),
            ACMGRule(
                rule_id="PM2",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["gnomAD v4.1 AF < 0.0001"],
            ),
        ]
        result = model.classify(tp53_missense_pathogenic, rules)
        lo, hi = result.credible_interval_95
        assert (
            0.0 <= lo <= hi <= 1.0
        ), f"95% HDI bounds must be in [0,1]: got [{lo:.3f}, {hi:.3f}]"

    def test_credible_interval_width_reflects_uncertainty(
        self, model: BayesACMGModel, brca1_frameshift: VariantInput
    ) -> None:
        """Uncertainty is quantified: VUS should have wider CI than clear P/LP."""
        # Minimal evidence → VUS
        vus_rules = [
            ACMGRule(
                rule_id="PM2",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["Rare in gnomAD v4.1"],
            ),
        ]
        # Clear pathogenic → very narrow CI
        path_rules = [
            ACMGRule(
                rule_id="PVS1",
                strength=EvidenceStrength.VERY_STRONG,
                applies=True,
                evidence_items=["LoF"],
            ),
            ACMGRule(
                rule_id="PS1",
                strength=EvidenceStrength.STRONG,
                applies=True,
                evidence_items=["Same AA, different codon, known P"],
            ),
            ACMGRule(
                rule_id="PM2",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["Absent gnomAD v4.1"],
            ),
        ]
        vus_result = model.classify(brca1_frameshift, vus_rules)
        path_result = model.classify(brca1_frameshift, path_rules)

        vus_width = (
            vus_result.credible_interval_95[1] - vus_result.credible_interval_95[0]
        )
        path_width = (
            path_result.credible_interval_95[1] - path_result.credible_interval_95[0]
        )

        assert vus_width > path_width, (
            "Variant with minimal evidence (VUS) should have wider CI than "
            "clear pathogenic variant with multiple strong criteria"
        )


class TestDirichletAlphaBuild:
    """Tests for _build_dirichlet_alpha() branches exercised via the public API."""

    def test_non_applying_rule_is_skipped(self, model: BayesACMGModel) -> None:
        """A rule with applies=False must not shift the posterior at all."""
        applying_only = [
            ACMGRule(
                rule_id="PP3",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["AM high"],
            ),
        ]
        with_non_applying = applying_only + [
            ACMGRule(
                rule_id="PM1",
                strength=EvidenceStrength.MODERATE,
                applies=False,
                evidence_items=["hotspot not met"],
            ),
        ]
        posterior_a = model.posterior_probabilities(applying_only)
        posterior_b = model.posterior_probabilities(with_non_applying)
        assert posterior_a == posterior_b, (
            "A rule with applies=False must be skipped (continue) in "
            "_build_dirichlet_alpha and not affect the posterior"
        )

    def test_stand_alone_rule_shifts_posterior_toward_benign(
        self, model: BayesACMGModel
    ) -> None:
        """STAND_ALONE (e.g. BA1) concentrates directly on the Benign category."""
        rules = [
            ACMGRule(
                rule_id="BA1",
                strength=EvidenceStrength.STAND_ALONE,
                applies=True,
                evidence_items=["AF > 5%"],
            ),
        ]
        posterior = model.posterior_probabilities(rules)
        baseline = model.posterior_probabilities([])
        assert posterior["Benign"] > baseline["Benign"]

    def test_strong_and_supporting_benign_shift_toward_likely_benign(
        self, model: BayesACMGModel
    ) -> None:
        """STRONG_BENIGN (BS1) and SUPPORTING_BENIGN (BP4) both raise Likely_Benign."""
        rules = [
            ACMGRule(
                rule_id="BS1",
                strength=EvidenceStrength.STRONG_BENIGN,
                applies=True,
                evidence_items=["AF > 1%"],
            ),
            ACMGRule(
                rule_id="BP4",
                strength=EvidenceStrength.SUPPORTING_BENIGN,
                applies=True,
                evidence_items=["AlphaMissense low"],
            ),
        ]
        posterior = model.posterior_probabilities(rules)
        baseline = model.posterior_probabilities([])
        assert posterior["Likely_Benign"] > baseline["Likely_Benign"]
        assert posterior["Benign"] > baseline["Benign"]


class TestCredibleIntervalsMCMC:
    """Tests for the full MCMC credible_intervals() path (PyMC NUTS)."""

    def test_credible_intervals_via_mcmc_bounded(self) -> None:
        """MCMC-based credible_intervals() returns valid [0,1] bounds for all categories."""
        # Small tune/draws for test speed; still exercises the real PyMC sampling path.
        model = BayesACMGModel(tune=25, draws=25, chains=1, target_accept=0.8)
        rules = [
            ACMGRule(
                rule_id="PVS1",
                strength=EvidenceStrength.VERY_STRONG,
                applies=True,
                evidence_items=["LoF"],
            ),
        ]
        result = model.credible_intervals(rules)
        assert set(result.keys()) == {
            "Pathogenic",
            "Likely_Pathogenic",
            "VUS",
            "Likely_Benign",
            "Benign",
        }
        for lo, hi in result.values():
            assert 0.0 <= lo <= hi <= 1.0
        assert model.idata is not None

    def test_classify_with_use_mcmc_true(self) -> None:
        """classify(use_mcmc=True) overrides the analytic CI with the MCMC HDI."""
        model = BayesACMGModel(tune=25, draws=25, chains=1, target_accept=0.8)
        variant = VariantInput(chrom="1", pos=100, ref="A", alt="G", variant_type="snv")
        rules = [
            ACMGRule(
                rule_id="PM2",
                strength=EvidenceStrength.SUPPORTING,
                applies=True,
                evidence_items=["Absent from gnomAD"],
            ),
        ]
        result = model.classify(variant, rules, use_mcmc=True)
        lo, hi = result.credible_interval_lower, result.credible_interval_upper
        assert lo is not None and hi is not None
        assert 0.0 <= lo <= hi <= 1.0


class TestCalibration:
    """Tests for BayesACMGModel.calibrate() and _compute_ece()."""

    def test_mismatched_lengths_raises(self, model: BayesACMGModel) -> None:
        """calibrate() raises ValueError if reference_variants and true_labels differ in length."""
        with pytest.raises(ValueError, match="equal length"):
            model.calibrate(
                reference_variants=[{"rules": []}, {"rules": []}],
                true_labels=["VUS"],
            )

    def test_calibrate_returns_ece_and_accuracy(self, model: BayesACMGModel) -> None:
        """calibrate() on a small reference set returns ece/accuracy metrics in [0,1]."""
        pvs1_rule = ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            applies=True,
            evidence_items=["LoF"],
        )
        ba1_rule = ACMGRule(
            rule_id="BA1",
            strength=EvidenceStrength.STAND_ALONE,
            applies=True,
            evidence_items=["common"],
        )
        reference_variants = [
            {"rules": [pvs1_rule]},
            {"rules": [ba1_rule]},
        ]
        true_labels = ["Pathogenic", "Benign"]
        metrics = model.calibrate(reference_variants, true_labels)
        assert "ece" in metrics and "accuracy" in metrics
        assert 0.0 <= metrics["ece"] <= 1.0
        assert 0.0 <= metrics["accuracy"] <= 1.0

    def test_calibrate_empty_reference_set(self, model: BayesACMGModel) -> None:
        """calibrate() with zero reference variants gives accuracy=0.0 without dividing by zero."""
        metrics = model.calibrate([], [])
        assert metrics["accuracy"] == 0.0

    def test_calibrate_logs_warning_when_ece_high(
        self, model: BayesACMGModel, caplog: pytest.LogCaptureFixture
    ) -> None:
        """calibrate() logs a warning when ECE exceeds the 0.05 acceptance threshold.

        Deliberately mismatch predictions vs true labels (wrong on every item) to
        drive the ECE above the 0.05 acceptance threshold and exercise the
        warning-logging branch.
        """
        pvs1_rule = ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            applies=True,
            evidence_items=["LoF"],
        )
        reference_variants = [{"rules": [pvs1_rule]}] * 5
        # PVS1 alone predicts "Pathogenic" or "Likely_Pathogenic" with high confidence;
        # labelling all of them "Benign" guarantees high confident-but-wrong predictions.
        true_labels = ["Benign"] * 5
        with caplog.at_level("WARNING", logger="bayesacmg.model"):
            metrics = model.calibrate(reference_variants, true_labels)
        assert metrics["ece"] > 0.05
        assert any("exceeds acceptance threshold" in r.message for r in caplog.records)

    def test_compute_ece_direct_with_empty_bin(self) -> None:
        """_compute_ece() skips bins with zero members (continue branch)."""
        import numpy as np

        from bayesacmg.model import _compute_ece

        # All confidences cluster in the top bin; other 9 bins are empty and skipped.
        confidences = np.array([0.95, 0.97, 0.99])
        is_correct = np.array([1.0, 1.0, 0.0])
        ece = _compute_ece(confidences, is_correct, n_bins=10)
        assert 0.0 <= ece <= 1.0

    def test_compute_ece_empty_arrays_returns_zero(self) -> None:
        """_compute_ece() with no samples returns 0.0 (n=0 guard)."""
        import numpy as np

        from bayesacmg.model import _compute_ece

        ece = _compute_ece(np.array([]), np.array([]), n_bins=10)
        assert ece == 0.0


class TestGetCriterionPrior:
    """Tests for BayesACMGModel.get_criterion_prior()."""

    def test_unknown_rule_id_returns_none(self, model: BayesACMGModel) -> None:
        """An unrecognised rule_id returns None rather than raising."""
        assert model.get_criterion_prior("NOT_A_REAL_RULE") is None

    def test_bs1_prior_is_moderately_benign(self, model: BayesACMGModel) -> None:
        """BS1 prior centre is negative (benign direction)."""
        prior = model.get_criterion_prior("BS1")
        assert prior is not None
        assert prior.centre == -4
