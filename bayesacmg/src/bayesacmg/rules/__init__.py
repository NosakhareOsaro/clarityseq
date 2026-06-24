"""
bayesacmg.rules
===============

ACMG/AMP rule implementations, split by category.

Submodules:
    pathogenic  — PVS1, PS1-4, PM1-6, PP1-5 (15 pathogenic criteria)
    benign      — BA1, BS1-4, BP1-7 (12 benign criteria)
    splicing    — PP3/BP4/BP7 for splice-impacting variants
                  (Walker et al. 2023 PMID:36898414)
    mito        — Mitochondrial-specific rules (ACGS 2024 §6)

References:
    Richards et al. 2015 PMID:25741868 — original 28-criterion framework
    Tavtigian et al. 2020 PMID:32645316 — point-score system
    ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., 20 Feb 2024)
    ClinGen SVI Working Group 2024 recommendations
"""

from bayesacmg.rules import benign, mito, pathogenic, splicing

__all__ = ["pathogenic", "benign", "splicing", "mito"]
