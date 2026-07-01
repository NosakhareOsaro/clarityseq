"""
bayesacmg — Bayesian ACMG/AMP variant classifier.

Implements all 28 ACMG/AMP criteria following:
- Richards et al. 2015 PMID:25741868 (original framework)
- Tavtigian et al. 2020 PMID:32645316 (point-score system)
- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., Feb 2024)
- ClinGen SVI Working Group 2024 recommendations

Version: 0.1.0
"""

__version__ = "0.1.0"
__author__ = "ClaritySeq"
__license__ = "MIT"

from bayesacmg.models import (
    ACMGRule,
    ClassificationResult,
    EvidenceStrength,
    GeneData,
    TranscriptData,
    VariantInput,
)

__all__ = [
    "__version__",
    "ACMGRule",
    "ClassificationResult",
    "EvidenceStrength",
    "GeneData",
    "TranscriptData",
    "VariantInput",
]
