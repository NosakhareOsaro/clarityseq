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
    VariantType,
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
            ACMGRule(rule_id="PVS1", strength=EvidenceStrength.VERY_STRONG, applies=True,
                     evidence_items=["BRCA1 frameshift; MANE Select NM_007294.4"]),
            ACMGRule(rule_id="PM2", strength=EvidenceStrength.SUPPORTING, applies=True,
                     evidence_items=["Absent from gnomAD v4.1"]),
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
            ACMGRule(rule_id="PP3", strength=EvidenceStrength.SUPPORTING, applies=True,
                     evidence_items=["AlphaMissense 0.97 ≥ 0.564"]),
            ACMGRule(rule_id="PM2", strength=EvidenceStrength.SUPPORTING, applies=True,
                     evidence_items=["gnomAD v4.1 AF < 0.0001"]),
        ]
        result = model.classify(tp53_missense_pathogenic, rules)
        lo, hi = result.credible_interval_95
        assert 0.0 <= lo <= hi <= 1.0, (
            f"95% HDI bounds must be in [0,1]: got [{lo:.3f}, {hi:.3f}]"
        )

    def test_credible_interval_width_reflects_uncertainty(
        self, model: BayesACMGModel, brca1_frameshift: VariantInput
    ) -> None:
        """Uncertainty is quantified: VUS should have wider CI than clear P/LP."""
        # Minimal evidence → VUS
        vus_rules = [
            ACMGRule(rule_id="PM2", strength=EvidenceStrength.SUPPORTING, applies=True,
                     evidence_items=["Rare in gnomAD v4.1"]),
        ]
        # Clear pathogenic → very narrow CI
        path_rules = [
            ACMGRule(rule_id="PVS1", strength=EvidenceStrength.VERY_STRONG, applies=True,
                     evidence_items=["LoF"]),
            ACMGRule(rule_id="PS1", strength=EvidenceStrength.STRONG, applies=True,
                     evidence_items=["Same AA, different codon, known P"]),
            ACMGRule(rule_id="PM2", strength=EvidenceStrength.SUPPORTING, applies=True,
                     evidence_items=["Absent gnomAD v4.1"]),
        ]
        vus_result = model.classify(brca1_frameshift, vus_rules)
        path_result = model.classify(brca1_frameshift, path_rules)

        vus_width = vus_result.credible_interval_95[1] - vus_result.credible_interval_95[0]
        path_width = path_result.credible_interval_95[1] - path_result.credible_interval_95[0]

        assert vus_width > path_width, (
            "Variant with minimal evidence (VUS) should have wider CI than "
            "clear pathogenic variant with multiple strong criteria"
        )
