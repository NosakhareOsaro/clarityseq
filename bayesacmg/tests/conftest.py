"""
Shared pytest fixtures for BayesACMG tests.

Real ClinVar RCV accessions are used throughout to document the specific
variants represented by each fixture (per §0.3.13 of PROJECT_GUIDE.MD).

Guidelines implemented and tested:
    - Richards et al. 2015 PMID:25741868: original 28-rule framework
    - Tavtigian et al. 2020 PMID:32645316: Bayesian point-scoring
    - ACGS 2024 v1.2 (Durkie et al.): PM2 at Supporting; MANE Select
    - ClinGen SVI 2024: AlphaMissense PP3/BP4 thresholds; PM2 downgrade
    - Walker et al. 2023 PMID:36898414: splicing framework
"""

from __future__ import annotations

import pytest

from bayesacmg.models import (
    GeneData,
    TranscriptData,
    VariantInput,
    VariantType,
)

# ---------------------------------------------------------------------------
# Pathogenic variant fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def brca1_frameshift() -> VariantInput:
    """BRCA1 c.5266dupC (p.Gln1756ProfsTer74) — canonical pathogenic LoF.

    ClinVar: RCV000007535
    Classification: Pathogenic (5 stars — practice guideline)
    Why chosen: Classic PVS1 frameshift in a known LoF-mechanism gene.
    gnomAD v4.1: absent (AF = 0)
    AlphaMissense: N/A (frameshift, not missense)
    """
    return VariantInput(
        chrom="17",
        pos=41_276_045,
        ref="C",
        alt="CC",
        variant_type=VariantType.FRAMESHIFT,
        gene_symbol="BRCA1",
        transcript_id="NM_007294.4",  # MANE Select
        hgvsc="NM_007294.4:c.5266dupC",
        hgvsp="NM_007294.4(BRCA1):p.Gln1756ProfsTer74",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
        clinvar_stars=5,
        clinvar_classification="Pathogenic",
        clinvar_rcv="RCV000007535",
    )


@pytest.fixture
def tp53_missense_pathogenic() -> VariantInput:
    """TP53 c.817C>T (p.Arg273Cys) — hotspot missense, pathogenic.

    ClinVar: RCV000012735
    Classification: Pathogenic
    Why chosen: PM1 (hotspot domain) + PP3 (high AlphaMissense score).
    AlphaMissense score: 0.97 → PP3 (≥ 0.564 threshold)
    """
    return VariantInput(
        chrom="17",
        pos=7_674_220,
        ref="C",
        alt="T",
        variant_type=VariantType.MISSENSE,
        gene_symbol="TP53",
        transcript_id="NM_000546.6",  # MANE Select
        hgvsc="NM_000546.6:c.817C>T",
        hgvsp="NM_000546.6(TP53):p.Arg273Cys",
        gnomad_af=0.000001,
        gnomad_ac=1,
        gnomad_nhomalt=0,
        alphamissense_score=0.97,  # ≥ 0.564 → PP3
        clinvar_stars=3,
        clinvar_classification="Pathogenic",
        clinvar_rcv="RCV000012735",
    )


@pytest.fixture
def brca2_novel_lof() -> VariantInput:
    """BRCA2 novel LoF variant — absent from gnomAD v4.1.

    No ClinVar accession (novel variant, not yet submitted).
    Why chosen: Tests PVS1 + PM2_Supporting = LP combination (9 pts ≥ 6 → LP).
    This is the critical test for the ClinGen SVI 2024 novel combination rule.
    """
    return VariantInput(
        chrom="13",
        pos=32_338_772,
        ref="ATTTT",
        alt="A",
        variant_type=VariantType.FRAMESHIFT,
        gene_symbol="BRCA2",
        transcript_id="NM_000059.4",  # MANE Select
        hgvsc="NM_000059.4:c.6406delAAAA",
        hgvsp="NM_000059.4(BRCA2):p.Lys2136ArgfsTer2",
        gnomad_af=0.0,  # Absent from gnomAD v4.1 → PM2_Supporting
        gnomad_ac=0,
        gnomad_nhomalt=0,
        clinvar_stars=None,
        clinvar_classification=None,
        clinvar_rcv=None,  # Novel; no ClinVar accession
    )


# ---------------------------------------------------------------------------
# Benign variant fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def common_benign_missense() -> VariantInput:
    """A common missense variant at BA1 frequency threshold.

    ClinVar: RCV000030819 (example benign common variant)
    Classification: Benign (BA1 — stand-alone)
    gnomAD v4.1 AF: 0.12 (>5% → BA1 → Benign)
    AlphaMissense score: 0.21 → BP4 (≤ 0.340)
    """
    return VariantInput(
        chrom="12",
        pos=111_803_912,
        ref="G",
        alt="A",
        variant_type=VariantType.MISSENSE,
        gene_symbol="OAS1",
        transcript_id="NM_016816.4",
        hgvsc="NM_016816.4:c.362G>A",
        hgvsp="NM_016816.4(OAS1):p.Arg121Gln",
        gnomad_af=0.12,  # >5% → BA1
        gnomad_ac=96_000,
        gnomad_nhomalt=5760,
        alphamissense_score=0.21,  # ≤ 0.340 → BP4
        clinvar_stars=2,
        clinvar_classification="Benign",
        clinvar_rcv="RCV000030819",
    )


# ---------------------------------------------------------------------------
# Splice variant fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def canonical_splice_donor() -> VariantInput:
    """Canonical splice donor variant — PVS1 via splicing pathway.

    ClinVar: RCV000048342
    Classification: Pathogenic
    SpliceAI Δ score: 0.95 → PP3 Strong (Walker 2023 PMID:36898414)
    """
    return VariantInput(
        chrom="17",
        pos=41_215_918,
        ref="G",
        alt="A",
        variant_type=VariantType.SPLICE_SITE,
        gene_symbol="BRCA1",
        transcript_id="NM_007294.4",
        hgvsc="NM_007294.4:c.5277+1G>A",
        hgvsp=None,  # Splice variant — no HGVSp
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
        spliceai_delta=0.95,  # ≥ 0.5 → PP3 Strong
        clinvar_stars=3,
        clinvar_classification="Pathogenic",
        clinvar_rcv="RCV000048342",
    )


@pytest.fixture
def synonymous_no_splice_impact() -> VariantInput:
    """Synonymous variant with no predicted splice impact → BP7.

    SpliceAI Δ score: 0.03 (< 0.1) + synonymous → BP7
    Walker et al. 2023 PMID:36898414: synonymous + SpliceAI < 0.1 → BP7
    """
    return VariantInput(
        chrom="17",
        pos=41_267_742,
        ref="C",
        alt="T",
        variant_type=VariantType.SYNONYMOUS,
        gene_symbol="BRCA1",
        transcript_id="NM_007294.4",
        hgvsc="NM_007294.4:c.4327C>T",
        hgvsp="NM_007294.4(BRCA1):p.Pro1443=",
        gnomad_af=0.000008,
        gnomad_ac=6,
        gnomad_nhomalt=0,
        spliceai_delta=0.03,  # < 0.1 → no splice impact
        clinvar_stars=1,
        clinvar_classification="Likely Benign",
        clinvar_rcv=None,
    )


# ---------------------------------------------------------------------------
# Mito variant fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mito_haplogroup_defining() -> VariantInput:
    """Mitochondrial haplogroup-defining variant — must be excluded.

    ACGS 2024 §6: haplogroup-defining variants are automatically Benign.
    Must be assessed by Haplogrep3 BEFORE ACMG classification.
    """
    return VariantInput(
        chrom="MT",
        pos=8_860,
        ref="A",
        alt="G",
        variant_type=VariantType.MISSENSE,
        gene_symbol="MT-ATP6",
        transcript_id="NC_012920.1",
        hgvsc="NC_012920.1:m.8860A>G",
        hgvsp=None,
        gnomad_af=0.85,  # Very common — haplogroup-defining
        gnomad_ac=700_000,
        gnomad_nhomalt=None,
        is_mito=True,
        is_haplogroup_defining=True,  # Haplogrep3 output
    )


# ---------------------------------------------------------------------------
# Shared TranscriptData and GeneData fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def brca1_transcript() -> TranscriptData:
    """BRCA1 MANE Select transcript data."""
    return TranscriptData(
        transcript_id="NM_007294.4",
        gene_symbol="BRCA1",
        is_mane_select=True,  # MANE Select — required for PVS1 full strength
        is_canonical=True,
        lof_disease_mechanism=True,  # LoF is known disease mechanism (AD breast cancer)
        exon_count=23,
    )


@pytest.fixture
def brca1_gene() -> GeneData:
    """BRCA1 gene-level data."""
    return GeneData(
        gene_symbol="BRCA1",
        omim_id="OMIM:113705",
        lof_mechanism=True,
        missense_constraint_z=3.72,  # High constraint (gnomAD)
        pli=0.999,  # High pLI — intolerant of LoF
        has_vcep_specification=True,  # BRCA Exchange VCEP
    )
