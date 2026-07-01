"""
annotation — Variant annotation stack for ClaritySeq WGS platform.

Provides unified access to:
- VEP v111 with MANE Select transcript prioritisation (Morales et al. 2022 PMID:35356062)
- AlphaMissense missense pathogenicity scores (Cheng et al. 2023 PMID:37703350)
- dbNSFP v4.7 aggregated in-silico predictors
- gnomAD v4.1 population allele frequencies (April 2024; 807,162 individuals)
- ClinVar clinical significance data
- PanelApp gene panel membership
- ClinGen gene-disease validity curations

All clients are designed to be async-compatible and cache-friendly.

References
----------
- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., Feb 2024)
- ClinGen SVI Working Group 2024 thresholds for computational evidence
"""

from annotation.alphamissense_client import AlphaMissenseClient, classify_am_score
from annotation.clinvar_client import ClinVarClient, ClinVarData
from annotation.gnomad_client import GnomADClient, GnomADData
from annotation.mane_select import (
    adjust_pvs1_for_mane,
    get_mane_select_for_gene,
    is_mane_select,
)
from annotation.vep_runner import VEPRunner

__version__ = "0.1.0"
__all__ = [
    "AlphaMissenseClient",
    "classify_am_score",
    "ClinVarClient",
    "ClinVarData",
    "GnomADClient",
    "GnomADData",
    "VEPRunner",
    "is_mane_select",
    "get_mane_select_for_gene",
    "adjust_pvs1_for_mane",
]
