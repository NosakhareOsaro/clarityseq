"""
multi_ancestry
==============
Multi-ancestry analysis module for GenomeForge.

Handles ancestry inference, population label assignment for admixed samples,
and selection of gnomAD v4.1 ancestry-stratified VQSR training sets.

Submodules:
    somalier_runner   — Run somalier relate + ancestry inference.
    ancestry_assigner — Assign population labels; fallback for admixed.
    vqsr_selector     — Select gnomAD v4.1 ancestry-stratified VQSR training.

References:
    Pedersen et al. 2020 Genome Biology PMID:32620139 (somalier).
    Karczewski et al. 2020 Nature PMID:32461654 (gnomAD v3 ancestry).
    Chen et al. 2024 (gnomAD v4.1 ancestry).
"""
