# ADR-002: Why PyMC over Stan for the Bayesian ACMG model

**Status:** Accepted  
**Date:** 2026-06-22  
**Deciders:** GenomeForge core team  
**Category:** Bayesian modelling

---

## Context

GenomeForge's BayesACMG module implements a fully Bayesian ACMG/AMP variant
classification model. The model computes posterior probabilities of pathogenicity
given observed evidence (28 ACMG rules, ClinGen SVI 2024 updates, VCEP overrides).

Two leading probabilistic programming frameworks were evaluated:

1. **PyMC v5.x** — Python-native Bayesian modelling library built on PyTensor
2. **Stan** — C++ probabilistic programming language with Python interface (PyStan/CmdStanPy)

---

## Decision

**PyMC v5.x is the Bayesian modelling framework for BayesACMG.**

---

## Rationale

### 1. Python-native ecosystem

GenomeForge is a Python 3.12 project. PyMC integrates seamlessly with:

- **NumPy/SciPy**: for evidence weight matrices
- **pandas**: for ClinVar calibration data loading
- **arviz**: for posterior visualisation and credible interval calculation
- **matplotlib**: for calibration plots in `calibration/`

Stan requires a separate compilation step and a C++ toolchain, complicating
the Docker build and CI environment.

### 2. PyMC 5 has mature NUTS sampling with PyTensor

PyMC 5 uses PyTensor (successor to Theano/Aesara) which:

- Provides automatic differentiation for gradient-based NUTS sampling
- Supports GPU-accelerated sampling via JAX backend (`PYTENSOR_FLAGS=device=cuda`)
- Has stable NUTS (No-U-Turn Sampler) implementation matching Stan's quality

### 3. Easier Docker integration

PyMC installs via pip: `pip install pymc==5.13.1`. No compilation step.

Stan (via CmdStanPy) requires compiling the CmdStan C++ binary in the Docker build,
adding ~15 minutes to build time and requiring `cmake`, `gcc`, etc.

### 4. The BayesACMG model does not require Stan's advanced features

The BayesACMG model uses:
- Beta priors on pathogenicity probability (conjugate to Bernoulli likelihood)
- NUTS sampling for posterior inference
- Bayesian credible intervals for clinical reporting

This is well within PyMC's capabilities. Stan's advantages (better divergence
handling, marginalisation) are relevant for more complex hierarchical models.

### 5. Active development and clinical bioinformatics adoption

PyMC 5 (released 2023) is actively maintained. It is used in:
- VariantInterpretation (Nature 2023)
- Multiple Bayesian ACMG papers (ClinGen Bayesian framework)
- Broad Institute variant interpretation tools

---

## Model architecture

The BayesACMG PyMC model (`bayesacmg/src/bayesacmg/model.py`):

```python
with pm.Model() as acmg_model:
    # Prior: pathogenicity probability for a variant of unknown significance
    # Beta(1, 1) = uniform prior; calibrate against ClinVar LP/LB variants
    p_path = pm.Beta("p_path", alpha=1.0, beta=1.0)

    # Evidence likelihood: each ACMG rule contributes log-odds
    # Weights from Tavtigian et al. 2020 (PMID:32125505)
    # ClinGen SVI 2024 updates applied to PM2, PP3/BP4 (AlphaMissense)
    log_odds = pm.Deterministic("log_odds", compute_log_odds(evidence, p_path))

    # Posterior
    posterior = pm.Deterministic("posterior", pm.math.sigmoid(log_odds))
```

### Calibration

Prior weights are calibrated against ClinVar pathogenic/likely-pathogenic and
benign/likely-benign variants using the `calibration/` module.
Calibration curves are generated and stored in `calibration/results/`.

---

## Consequences

**Positive:**
- Pure Python; no C++ toolchain required
- pip-installable; Docker build is simple
- arviz integration for credible interval reporting (required by ACGS 2024 for
  Bayesian classification reporting)
- GPU sampling available via JAX backend for large cohorts

**Negative:**
- PyMC is generally slower than Stan for complex hierarchical models
- PyTensor compilation step on first run adds ~30s (acceptable for clinical use)

---

## Alternatives Considered

| Framework | Reason Not Chosen |
|-----------|------------------|
| Stan (PyStan/CmdStanPy) | Requires C++ compilation; complex Docker build |
| NumPyro (JAX-based) | Less documentation for clinical bioinformatics; newer and less mature |
| TensorFlow Probability | Heavy dependency; less ergonomic API |
| Variational inference only | Approximation too inaccurate for clinical credible intervals |
| Point estimate (frequentist) | Does not produce the credible intervals required by ACGS 2024 §5.4 |

---

## References

- Tavtigian et al. 2020: Bayesian ACMG framework PMID:32125505
- ClinGen SVI 2024 updated recommendations
- PyMC documentation: <https://www.pymc.io/welcome.html>
- Durkie et al. 2024 (ACGS 2024 v1.2) — Bayesian reporting requirements §5.4
