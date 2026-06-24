"""
prioritisation
==============
Variant prioritisation module for GenomeForge.

Combines HPO phenotype scoring, inheritance mode filtering, Exomiser 14
prioritisation, and composite ACMG/HPO/inheritance ranking.

Submodules:
    hpo_scorer         — HPO phenotype-to-gene scoring.
    inheritance_filter — Filter variants by AD/AR/XL/de_novo mode.
    exomiser_client    — Exomiser 14 REST API client.
    ranking            — Composite ranking (ACMG + HPO + inheritance + panel).

References:
    Köhler et al. 2021 PMID:33264411 (HPO).
    Robinson et al. 2023 Nature Genetics PMID:37604970 (Exomiser).
"""
