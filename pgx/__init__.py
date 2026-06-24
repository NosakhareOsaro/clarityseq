"""
pgx
===
Pharmacogenomics (PGx) module for GenomeForge.

Provides CYP2D6 star allele genotyping (via Cyrius), PharmVar star allele
definitions, and CPIC drug dosing recommendations.

Submodules:
    cyrius_runner    — CYP2D6 star allele calling via Cyrius subprocess.
    pharmvar_client  — PharmVar REST API client for star allele definitions.
    cpic_client      — CPIC REST API client for drug dosing recommendations.

Key references:
    Aliev et al. 2022 NPJ Genomic Medicine PMID:35264608 (Cyrius validation).
    CPIC: https://cpicpgx.org/
    PharmVar: https://www.pharmvar.org/
"""
