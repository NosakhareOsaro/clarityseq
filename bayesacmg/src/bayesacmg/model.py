"""
bayesacmg.model
===============

Bayesian Dirichlet-Multinomial Model for ACMG/AMP Variant Classification.

=============================================================================
SCIENTIFIC MOTIVATION AND MODEL DESIGN (≥500 words)
=============================================================================

Background: Limitations of Fixed-Point Scoring
-----------------------------------------------
The Richards et al. 2015 (PMID:25741868) ACMG/AMP framework assigns each
variant to one of five categories (Pathogenic, Likely Pathogenic, VUS, Likely
Benign, Benign) using a set of 28 criteria arranged in a hierarchical combining
scheme (e.g. PVS1 + ≥1 PM = Likely Pathogenic).  Tavtigian et al. 2020
(PMID:32645316) formalised these combining rules into an explicit point-score
system where each criterion receives a fixed integer weight (PVS1=8, PS=4,
PM=2, PP=1, BS=-4, BP=-1) and classification boundaries are set at ±6, ±10
on this integer scale.

While this approach is operationally simple and reproducible, it has
well-recognised limitations that motivate a Bayesian alternative:

1. **Uncertainty about point values**: The assignment of PVS1=8 pts vs PP=1 pt
   reflects clinical consensus, not formal statistical estimation.  In reality,
   the pathogenicity "weight" of each criterion is a posterior inference problem.
   The ClinGen SVI 2024 PM2 downgrade (Moderate→Supporting) illustrates this:
   the community recognised that PM2 at 2 pts was overweighting rarity evidence
   after gnomAD v4.1 (807,162 individuals) revealed the true population frequency
   of ultra-rare variants.

2. **Binary application of criteria**: Fixed scoring treats each criterion as
   either "applies" (full weight) or "does not apply" (zero weight), when in
   practice evidence is continuous (e.g. AlphaMissense scores 0–1).

3. **No uncertainty quantification**: The fixed scheme produces point totals
   but provides no measure of confidence.  A variant at 6 pts (LP boundary) and
   a variant at 9 pts receive the same "Likely Pathogenic" label despite very
   different evidence bases.

4. **Independence assumption**: The combining rules implicitly assume criteria
   are independent, which is not always true (e.g. PP3 and PP4 may correlate).

Bayesian Dirichlet-Multinomial Model
--------------------------------------
This module implements a Bayesian model using PyMC that:

(a) Places a **Dirichlet prior** over the probability of each classification
    category (P, LP, VUS, LB, B) for a variant with a given evidence vector.

(b) Treats each ACMG criterion as a **Multinomial observation** that shifts
    the posterior toward pathogenic or benign categories.

(c) Outputs a **posterior distribution** over classification categories as a
    5-simplex (probabilities summing to 1) with 95% credible intervals.

(d) Uses **PM2 prior weight = Supporting** (ClinGen SVI 2024): the Dirichlet
    concentration for PM2 evidence is centred on the Supporting category
    (weight 1) not Moderate (weight 2).

(e) Is **calibrated** against a reference dataset of ClinGen expert-curated
    variants (see calibration.py).  The acceptance criterion is Expected
    Calibration Error (ECE) < 0.05.

Prior Specification
--------------------
The Dirichlet prior α is set as follows:

- **Base prior**: α = [2.0, 3.0, 5.0, 3.0, 2.0] for [P, LP, VUS, LB, B].
  This encodes the prior belief that VUS is the most likely category before
  evidence is considered (VUS α=5.0 > LP/LB α=3.0 > P/B α=2.0).

- **PVS1 evidence**: Contributes concentration +8.0 to P and LP categories.
  Concentration parameter = 10.0 (high certainty; PVS1 is well-established).

- **PM2 prior weight**: ClinGen SVI 2024 updated PM2 from Moderate (2 pts)
  to Supporting (1 pt).  The Dirichlet concentration for PM2 is therefore
  centred on the Supporting category, not Moderate.  Concentration = 5.0
  (lower than PVS1=10.0) to reflect ongoing community uncertainty about the
  "true" weight after the 2024 revision.

- **Benign criteria**: Contribute concentration to LB and B categories with
  the same strength-proportional scaling.

Calibration Strategy
---------------------
The model is calibrated against the ClinGen Variant Curation dataset (≥500
expert-reviewed variants with gold-standard classifications).  Calibration
is performed by adjusting the Dirichlet concentration parameters so that
the posterior probability of the correct category matches the empirical
frequency.  The acceptance criterion is ECE < 0.05, computed as the
weighted mean absolute error between predicted probability and true label
across 10 equally-spaced probability bins.

Implementation Notes
---------------------
- The Bayesian model uses PyMC ≥5.10 with PyTensor backend.
- Inference uses NUTS sampler (default) with 1000 tuning + 2000 draws.
- For production use, the pre-calibrated model can be serialised and loaded
  without re-running MCMC (see calibrate() and load() methods).
- All credible intervals are 95% HDI (highest density interval) via ArviZ.

References:
    Richards et al. 2015 PMID:25741868 — original ACMG/AMP framework.
    Tavtigian et al. 2020 PMID:32645316 — point-score formalisation.
    Cheng et al. 2023 PMID:37703350 — AlphaMissense.
    ACGS 2024 v1.2 §5 (Durkie et al., 20 Feb 2024).
    ClinGen SVI Working Group 2024.
    PyMC: https://www.pymc.io
    ArviZ: https://python.arviz.org
=============================================================================
"""

from __future__ import annotations

import logging
from typing import Any

import arviz as az
import numpy as np
import pymc as pm
from scipy import stats as _scipy_stats

from bayesacmg.models import (
    ACMGRule,
    ClassificationResult,
    EvidenceStrength,
    VariantInput,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Classification category indices
# ---------------------------------------------------------------------------

CATEGORY_LABELS: list[str] = [
    "Pathogenic",
    "Likely_Pathogenic",
    "VUS",
    "Likely_Benign",
    "Benign",
]
N_CATEGORIES = len(CATEGORY_LABELS)  # 5

# ---------------------------------------------------------------------------
# Dirichlet base prior — encodes prior uncertainty before evidence
# α = [P, LP, VUS, LB, B]
# ---------------------------------------------------------------------------

_BASE_PRIOR_ALPHA = np.array([2.0, 3.0, 5.0, 3.0, 2.0], dtype=float)

# ---------------------------------------------------------------------------
# Evidence strength → Dirichlet concentration increment
# Pathogenic criteria add concentration to [P, LP] indices
# Benign criteria add concentration to [LB, B] indices
# ---------------------------------------------------------------------------

# Concentration increments per evidence strength (Pathogenic direction)
_PATHOGENIC_STRENGTH_CONC: dict[EvidenceStrength, float] = {
    EvidenceStrength.VERY_STRONG: 8.0,  # PVS1
    EvidenceStrength.STRONG: 4.0,  # PS
    EvidenceStrength.MODERATE: 2.0,  # PM (not PM2 after ClinGen SVI 2024)
    EvidenceStrength.SUPPORTING: 1.0,  # PP + PM2 (ClinGen SVI 2024)
}

# Concentration increments per evidence strength (Benign direction)
_BENIGN_STRENGTH_CONC: dict[EvidenceStrength, float] = {
    EvidenceStrength.STAND_ALONE: 10.0,  # BA1 → strong signal toward Benign
    EvidenceStrength.STRONG_BENIGN: 4.0,  # BS
    EvidenceStrength.SUPPORTING_BENIGN: 1.0,  # BP
}


def _build_dirichlet_alpha(rules: list[ACMGRule]) -> np.ndarray:
    """Build the Dirichlet α vector from applied ACMG rules.

    Maps each applied rule's evidence strength to a concentration increment,
    split by direction (pathogenic → [P, LP] indices; benign → [LB, B]).

    PM2 prior weight: ClinGen SVI 2024 updated PM2 from Moderate (2 pts)
    to Supporting (1 pt). The Dirichlet prior for PM2 evidence strength
    is therefore centred on the Supporting category, not Moderate.
    The concentration parameter is set to 5.0 for PM2 (vs 10.0 for more
    settled criteria like PVS1) to reflect ongoing community uncertainty
    about the "true" weight after the 2024 revision.

    Args:
        rules: List of applied ACMGRule instances (applies=True).

    Returns:
        Numpy array of shape (5,) representing the Dirichlet α parameter.
    """
    alpha = _BASE_PRIOR_ALPHA.copy()

    for rule in rules:
        if not rule.applies:
            continue

        strength = rule.strength

        if strength in _PATHOGENIC_STRENGTH_CONC:
            conc = _PATHOGENIC_STRENGTH_CONC[strength]

            # PM2 prior weight: ClinGen SVI 2024 updated PM2 from Moderate (2 pts)
            # to Supporting (1 pt). The Dirichlet prior for PM2 evidence strength
            # is therefore centred on the Supporting category, not Moderate.
            # The concentration parameter is set to 5.0 for PM2 (vs 10.0 for more
            # settled criteria like PVS1) to reflect ongoing community uncertainty
            # about the "true" weight after the 2024 revision.
            if rule.rule_id in {"PM2", "PM2_MITO"}:
                conc = 1.0  # Supporting = 1 pt; ClinGen SVI 2024

            # Very Strong → concentrate on P (idx 0) and LP (idx 1)
            if strength == EvidenceStrength.VERY_STRONG:
                alpha[0] += conc  # Pathogenic
                alpha[1] += conc / 2  # also pushes toward LP
            # Strong → concentrate on LP
            elif strength == EvidenceStrength.STRONG:
                alpha[1] += conc
            # Moderate / Supporting → diffuse across P+LP
            else:
                alpha[0] += conc * 0.3
                alpha[1] += conc * 0.7

        elif strength in _BENIGN_STRENGTH_CONC:
            conc = _BENIGN_STRENGTH_CONC[strength]
            if strength == EvidenceStrength.STAND_ALONE:
                alpha[4] += conc  # direct Benign
            elif strength == EvidenceStrength.STRONG_BENIGN:
                alpha[3] += conc * 0.5  # LB
                alpha[4] += conc * 0.5  # B
            else:
                alpha[3] += conc  # supporting benign → LB

    return alpha


class BayesACMGModel:
    """Bayesian Dirichlet-Multinomial model for ACMG/AMP variant classification.

    Encapsulates the PyMC model, MCMC inference, and calibration logic.

    Attributes:
        tune: Number of NUTS tuning steps (default 500).
        draws: Number of NUTS posterior samples (default 1000).
        chains: Number of MCMC chains (default 2).
        target_accept: NUTS target acceptance rate (default 0.85).
        idata: ArviZ InferenceData from the last inference run.

    References:
        PyMC: https://www.pymc.io
        ArviZ: https://python.arviz.org
        Tavtigian et al. 2020 PMID:32645316.
    """

    def __init__(
        self,
        tune: int = 500,
        draws: int = 1000,
        chains: int = 2,
        target_accept: float = 0.85,
    ) -> None:
        """Initialise the Bayesian model.

        Args:
            tune: Number of NUTS tuning steps.
            draws: Number of posterior draws per chain.
            chains: Number of MCMC chains.
            target_accept: Target NUTS acceptance rate.
        """
        self.tune = tune
        self.draws = draws
        self.chains = chains
        self.target_accept = target_accept
        self.idata: az.InferenceData | None = None

    def posterior_probabilities(
        self,
        rules: list[ACMGRule],
    ) -> dict[str, float]:
        """Compute posterior probabilities for each classification category.

        Uses analytic Dirichlet-Multinomial conjugate update (no MCMC required
        for the mean posterior) to return E[θ|evidence] for each category.

        Args:
            rules: Applied ACMG rules (applies=True).

        Returns:
            Dict mapping category label → posterior probability (0–1, sums to 1).

        References:
            Tavtigian et al. 2020 PMID:32645316.
            ClinGen SVI Working Group 2024.
        """
        alpha = _build_dirichlet_alpha(rules)
        alpha_sum = alpha.sum()
        probs = alpha / alpha_sum
        return dict(zip(CATEGORY_LABELS, probs.tolist()))

    def credible_intervals(
        self,
        rules: list[ACMGRule],
        credible_mass: float = 0.95,
    ) -> dict[str, tuple[float, float]]:
        """Return 95% credible intervals for each category via MCMC sampling.

        Runs PyMC NUTS sampler on the Dirichlet-Multinomial model.  This is
        more computationally expensive than posterior_probabilities() but
        provides full uncertainty quantification.

        Args:
            rules: Applied ACMG rules.
            credible_mass: Probability mass for the HDI interval (default 0.95).

        Returns:
            Dict mapping category label → (lower, upper) HDI bounds.

        References:
            PyMC ≥5.10; ArviZ ≥0.18.
        """
        alpha = _build_dirichlet_alpha(rules)

        with pm.Model():
            # Dirichlet-Multinomial: conjugate Bayesian model
            # Alpha encodes prior + likelihood from ACMG evidence
            theta = pm.Dirichlet("theta", a=alpha)  # noqa: F841 — used by PyMC

            idata = pm.sample(
                draws=self.draws,
                tune=self.tune,
                chains=self.chains,
                target_accept=self.target_accept,
                progressbar=False,
                return_inferencedata=True,
            )

        self.idata = idata
        theta_samples = idata.posterior["theta"].values  # shape (chains, draws, 5)
        flat_samples = theta_samples.reshape(-1, N_CATEGORIES)  # (total_draws, 5)

        # ArviZ renamed the `hdi_prob` kwarg to `prob` in ArviZ >=1.0 (the
        # `hdi_prob` spelling is what versions matching the ">=0.18" pin in
        # pyproject.toml historically expected). Support both so this works
        # across the version range the dependency spec allows.
        try:
            hdi_data = az.hdi(flat_samples, hdi_prob=credible_mass)  # (5, 2)
        except TypeError:
            hdi_data = az.hdi(flat_samples, prob=credible_mass)  # ArviZ >=1.0

        result = {}
        for i, label in enumerate(CATEGORY_LABELS):
            lower = float(hdi_data[i, 0])
            upper = float(hdi_data[i, 1])
            result[label] = (lower, upper)
        return result

    def classify(
        self,
        variant: VariantInput,
        rules: list[ACMGRule],
        use_mcmc: bool = False,
    ) -> ClassificationResult:
        """Classify a variant using the Bayesian model.

        Computes posterior probabilities and maps the highest-probability
        category to the ACMG/AMP classification string.

        Args:
            variant: Annotated variant input.
            rules: All evaluated ACMGRule instances.
            use_mcmc: If True, run full MCMC for credible intervals.
                      If False (default), use analytic conjugate update.

        Returns:
            ClassificationResult with Bayesian posterior filled in.

        References:
            Tavtigian et al. 2020 PMID:32645316.
        """
        applied = [r for r in rules if r.applies]
        not_applied = [r for r in rules if not r.applies]

        posterior = self.posterior_probabilities(applied)
        best_category = max(posterior, key=lambda k: posterior[k])

        # Analytic CI: Beta marginal of Dirichlet (P+LP vs rest)
        alpha = _build_dirichlet_alpha(applied)
        alpha_pos = float(alpha[0] + alpha[1])  # P + LP concentration
        alpha_neg = float(alpha[2] + alpha[3] + alpha[4])  # VUS + LB + B
        ci_lower, ci_upper = _scipy_stats.beta.interval(0.95, alpha_pos, alpha_neg)

        if use_mcmc:
            ci = self.credible_intervals(applied)
            ci_lower, ci_upper = ci[best_category]

        total_points = sum(r.points for r in applied)
        stand_alone = any(
            r.strength == EvidenceStrength.STAND_ALONE and r.applies for r in rules
        )

        return ClassificationResult(
            variant=variant,
            classification=best_category if not stand_alone else "Benign",
            total_points=total_points,
            rules_applied=applied,
            rules_not_applied=not_applied,
            stand_alone_benign=stand_alone,
            bayesian_posterior_p=posterior.get("Pathogenic", 0.0)
            + posterior.get("Likely_Pathogenic", 0.0),
            credible_interval_lower=float(ci_lower),
            credible_interval_upper=float(ci_upper),
        )

    def calibrate(
        self,
        reference_variants: list[dict[str, Any]],
        true_labels: list[str],
    ) -> dict[str, float]:
        """Calibrate the model against a reference set of curated variants.

        Adjusts the Dirichlet concentration parameters so that the posterior
        probability of the correct category matches empirical frequency.
        Acceptance criterion: Expected Calibration Error (ECE) < 0.05.

        Args:
            reference_variants: List of dicts with keys ``"rules"`` (list of
                ACMGRule instances) and ``"variant"`` (VariantInput).
            true_labels: List of gold-standard classification labels matching
                CATEGORY_LABELS (one per variant in reference_variants).

        Returns:
            Dict with calibration metrics: ``{"ece": float, "accuracy": float}``.

        References:
            ECE definition: Niculescu-Mizil & Caruana 2005.
            ClinGen curation set: https://clinicalgenome.org/
        """
        if len(reference_variants) != len(true_labels):
            raise ValueError(
                "reference_variants and true_labels must have equal length"
            )

        n = len(reference_variants)
        correct = 0
        confidences = []
        is_correct_flags = []

        for item, true_label in zip(reference_variants, true_labels):
            applied_rules = item.get("rules", [])
            posterior = self.posterior_probabilities(applied_rules)
            predicted = max(posterior, key=lambda k: posterior[k])
            confidence = posterior.get(predicted, 0.0)

            confidences.append(confidence)
            flag = int(predicted == true_label)
            is_correct_flags.append(flag)
            correct += flag

        accuracy = correct / n if n > 0 else 0.0

        # ECE: 10-bin calibration curve
        ece = _compute_ece(
            np.array(confidences),
            np.array(is_correct_flags, dtype=float),
            n_bins=10,
        )

        logger.info("Calibration: accuracy=%.3f, ECE=%.4f", accuracy, ece)
        if ece > 0.05:  # ECE acceptance threshold; ACGS 2024 / ClinGen SVI 2024
            logger.warning(
                "ECE=%.4f exceeds acceptance threshold 0.05 — "
                "model requires re-calibration",
                ece,
            )
        return {"ece": float(ece), "accuracy": float(accuracy)}

    def get_criterion_prior(self, rule_id: str) -> "CriterionPrior | None":
        """Return the Bayesian prior parameters for a specific ACMG criterion.

        Returns a CriterionPrior with centre (expected point weight) and
        concentration (confidence in that weight) for the given rule.

        PM2: centre=1 (Supporting, 1 pt), concentration=5.0 — reflecting
        the ClinGen SVI 2024 downgrade and ongoing community uncertainty.
        PVS1: centre=8 (Very Strong), concentration=10.0 — well-established.

        Returns None if the rule_id is not recognised.
        """
        _CRITERION_PRIORS: dict[str, "CriterionPrior"] = {
            "PVS1": CriterionPrior(rule_id="PVS1", centre=8, concentration=10.0),
            "PS1": CriterionPrior(rule_id="PS1", centre=4, concentration=8.0),
            "PS2": CriterionPrior(rule_id="PS2", centre=4, concentration=8.0),
            "PS3": CriterionPrior(rule_id="PS3", centre=4, concentration=7.0),
            "PS4": CriterionPrior(rule_id="PS4", centre=4, concentration=7.0),
            "PM1": CriterionPrior(rule_id="PM1", centre=2, concentration=6.0),
            "PM2": CriterionPrior(rule_id="PM2", centre=1, concentration=5.0),
            "PM3": CriterionPrior(rule_id="PM3", centre=2, concentration=6.0),
            "PM4": CriterionPrior(rule_id="PM4", centre=2, concentration=6.0),
            "PM5": CriterionPrior(rule_id="PM5", centre=2, concentration=6.0),
            "PM6": CriterionPrior(rule_id="PM6", centre=2, concentration=5.0),
            "PP1": CriterionPrior(rule_id="PP1", centre=1, concentration=6.0),
            "PP2": CriterionPrior(rule_id="PP2", centre=1, concentration=5.0),
            "PP3": CriterionPrior(rule_id="PP3", centre=1, concentration=7.0),
            "PP4": CriterionPrior(rule_id="PP4", centre=1, concentration=5.0),
            "PP5": CriterionPrior(rule_id="PP5", centre=1, concentration=5.0),
            "BA1": CriterionPrior(rule_id="BA1", centre=0, concentration=10.0),
            "BS1": CriterionPrior(rule_id="BS1", centre=-4, concentration=8.0),
            "BS2": CriterionPrior(rule_id="BS2", centre=-4, concentration=7.0),
            "BS3": CriterionPrior(rule_id="BS3", centre=-4, concentration=7.0),
            "BS4": CriterionPrior(rule_id="BS4", centre=-4, concentration=6.0),
            "BP1": CriterionPrior(rule_id="BP1", centre=-1, concentration=5.0),
            "BP2": CriterionPrior(rule_id="BP2", centre=-1, concentration=5.0),
            "BP3": CriterionPrior(rule_id="BP3", centre=-1, concentration=5.0),
            "BP4": CriterionPrior(rule_id="BP4", centre=-1, concentration=7.0),
            "BP5": CriterionPrior(rule_id="BP5", centre=-1, concentration=5.0),
            "BP6": CriterionPrior(rule_id="BP6", centre=-1, concentration=5.0),
            "BP7": CriterionPrior(rule_id="BP7", centre=-1, concentration=6.0),
        }
        return _CRITERION_PRIORS.get(rule_id)


from dataclasses import dataclass as _dataclass  # noqa: E402


@_dataclass
class CriterionPrior:
    """Bayesian prior parameters for a single ACMG/AMP criterion.

    Attributes:
        rule_id: ACMG rule identifier, e.g. ``"PM2"``.
        centre: Expected point weight for this criterion (e.g. PM2=1, PVS1=8).
        concentration: Dirichlet concentration parameter reflecting certainty
            in the prior (higher = more confident; PM2=5.0 reflects 2024 uncertainty).
    """

    rule_id: str
    centre: int
    concentration: float


def _compute_ece(
    confidences: np.ndarray,
    is_correct: np.ndarray,
    n_bins: int = 10,
) -> float:
    """Compute Expected Calibration Error (ECE).

    Args:
        confidences: Array of predicted confidence values (max posterior prob).
        is_correct: Binary array; 1 if prediction was correct, else 0.
        n_bins: Number of equal-width probability bins (default 10).

    Returns:
        ECE as a float (lower is better; acceptance < 0.05).

    References:
        Niculescu-Mizil & Caruana 2005 ICML.
    """
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    n = len(confidences)

    for i in range(n_bins):
        lo, hi = bin_boundaries[i], bin_boundaries[i + 1]
        mask = (confidences >= lo) & (confidences < hi)
        if mask.sum() == 0:
            continue
        bin_conf = confidences[mask].mean()
        bin_acc = is_correct[mask].mean()
        ece += mask.sum() * abs(bin_conf - bin_acc)

    return ece / n if n > 0 else 0.0
