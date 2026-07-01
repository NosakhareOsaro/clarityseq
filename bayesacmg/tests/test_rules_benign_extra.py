"""
bayesacmg.tests.test_rules_benign_extra
=========================================
Additional pytest tests for the BayesACMG benign rules covering
the applies=True branches not exercised by test_rules_benign.py.

Covers: BS1, BS2, BS3, BS4, BP1, BP2, BP3, BP4 (splice/REVEL), BP5, BP6, BP7.

References:
    Richards et al. 2015 PMID:25741868 (ACMG/AMP criteria).
    ACGS 2024 v1.2 (Durkie et al.) §5 Table 2.
    Cheng et al. 2023 PMID:37703350 (AlphaMissense thresholds).
    Walker et al. 2023 PMID:36898414 (SpliceAI BP7 threshold).
"""

from __future__ import annotations

import pytest

from bayesacmg.models import (
    EvidenceStrength,
    GeneData,
    TranscriptData,
    VariantInput,
    VariantType,
)
from bayesacmg.rules.benign import (
    rule_bs1,
    rule_bs2,
    rule_bs3,
    rule_bs4,
    rule_bp1,
    rule_bp2,
    rule_bp3,
    rule_bp4,
    rule_bp5,
    rule_bp6,
    rule_bp7,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def rare_missense() -> VariantInput:
    """Rare missense variant — AF well below any benign threshold."""
    return VariantInput(
        chrom="7",
        pos=117_548_628,
        ref="G",
        alt="A",
        variant_type=VariantType.MISSENSE,
        gene_symbol="CFTR",
        transcript_id="NM_000492.4",
        gnomad_af=0.0001,
        gnomad_popmax_af=0.0002,
        gnomad_ac=1,
        gnomad_nhomalt=0,
    )


@pytest.fixture
def common_variant_bs1() -> VariantInput:
    """Common variant with AF > 1% — qualifies for BS1."""
    return VariantInput(
        chrom="1",
        pos=11_794_419,
        ref="C",
        alt="T",
        variant_type=VariantType.MISSENSE,
        gene_symbol="MTHFR",
        transcript_id="NM_005957.5",
        gnomad_af=0.035,  # 3.5% > 1% → BS1
        gnomad_popmax_af=0.035,
        gnomad_ac=28_000,
        gnomad_nhomalt=490,
    )


@pytest.fixture
def inframe_variant() -> VariantInput:
    """In-frame insertion variant for BP3 testing."""
    return VariantInput(
        chrom="11",
        pos=108_111_710,
        ref="ATG",
        alt="ATGCTG",
        variant_type="inframe_insertion",
        gene_symbol="ATM",
        transcript_id="NM_000051.4",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
    )


@pytest.fixture
def truncating_only_gene() -> GeneData:
    """Gene where only truncating variants cause disease (BP1 context)."""
    return GeneData(
        gene_symbol="APC",
        lof_is_disease_mechanism=True,
        missense_only_gene=False,  # NOT missense-only → BP1 applies for missense
        gnomad_pli=0.998,
    )


@pytest.fixture
def missense_transcript() -> TranscriptData:
    """Transcript with a missense amino acid change."""
    return TranscriptData(
        transcript_id="NM_000038.6",
        gene_symbol="APC",
        is_mane_select=True,
        is_canonical=True,
        aa_change="p.Gly1438Ser",
        lof_disease_mechanism=True,
    )


@pytest.fixture
def missense_variant_for_bp1() -> VariantInput:
    """Missense variant (snv type) for BP1 testing."""
    return VariantInput(
        chrom="5",
        pos=112_175_770,
        ref="G",
        alt="A",
        variant_type="snv",
        gene_symbol="APC",
        transcript_id="NM_000038.6",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
    )


# ---------------------------------------------------------------------------
# BS1
# ---------------------------------------------------------------------------


class TestBS1:
    """BS1 — allele frequency > 1% → Strong Benign."""

    def test_applies_when_af_exceeds_one_percent(
        self, common_variant_bs1: VariantInput
    ) -> None:
        """AF = 3.5% > 1% → BS1 applies."""
        result = rule_bs1(common_variant_bs1)
        assert result.applies is True
        assert result.rule_id == "BS1"
        assert result.strength == EvidenceStrength.STRONG_BENIGN

    def test_not_applies_when_af_below_threshold(
        self, rare_missense: VariantInput
    ) -> None:
        """AF = 0.02% < 1% → BS1 does not apply."""
        result = rule_bs1(rare_missense)
        assert result.applies is False

    def test_not_applies_when_af_is_none(self) -> None:
        """AF = None → BS1 does not apply."""
        variant = VariantInput(
            chrom="1",
            pos=100,
            ref="A",
            alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="GENE1",
            transcript_id="NM_000001.1",
            gnomad_af=None,
            gnomad_popmax_af=None,
        )
        result = rule_bs1(variant)
        assert result.applies is False

    def test_uses_popmax_af_when_available(self) -> None:
        """BS1 uses gnomad_popmax_af when set (even if gnomad_af is lower)."""
        variant = VariantInput(
            chrom="1",
            pos=100,
            ref="A",
            alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="GENE1",
            transcript_id="NM_000001.1",
            gnomad_af=0.005,  # 0.5% — below threshold
            gnomad_popmax_af=0.02,  # 2% — above threshold → BS1
        )
        result = rule_bs1(variant)
        assert result.applies is True


# ---------------------------------------------------------------------------
# BS2
# ---------------------------------------------------------------------------


class TestBS2:
    """BS2 — observed in healthy adults for dominant/recessive disorder."""

    def test_applies_when_observed_in_healthy_adults(
        self, rare_missense: VariantInput
    ) -> None:
        """Observed in healthy adults → BS2 applies."""
        result = rule_bs2(rare_missense, observed_in_healthy_adults=True)
        assert result.applies is True
        assert result.rule_id == "BS2"
        assert result.strength == EvidenceStrength.STRONG_BENIGN

    def test_applies_when_two_homozygotes(self, rare_missense: VariantInput) -> None:
        """≥2 healthy homozygotes → BS2 applies."""
        result = rule_bs2(
            rare_missense,
            observed_in_healthy_adults=False,
            n_healthy_homozygotes=5,
        )
        assert result.applies is True

    def test_applies_when_gnomad_homalt_two(self) -> None:
        """gnomAD nhomalt ≥ 2 → BS2 applies (recessive context)."""
        variant = VariantInput(
            chrom="1",
            pos=100,
            ref="A",
            alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="CF",
            transcript_id="NM_000492.4",
            gnomad_nhomalt=3,  # ≥2 homozygotes in gnomAD
        )
        result = rule_bs2(variant, observed_in_healthy_adults=False)
        assert result.applies is True

    def test_not_applies_when_no_healthy_observations(
        self, rare_missense: VariantInput
    ) -> None:
        """No healthy observations, no homozygotes → BS2 does not apply."""
        result = rule_bs2(rare_missense, observed_in_healthy_adults=False, n_healthy_homozygotes=0)
        assert result.applies is False

    def test_evidence_includes_both_criteria_when_both_met(
        self, rare_missense: VariantInput
    ) -> None:
        """Evidence items mention both healthy adults and homozygotes."""
        result = rule_bs2(rare_missense, observed_in_healthy_adults=True, n_healthy_homozygotes=4)
        combined = " ".join(result.evidence_items)
        assert "healthy" in combined.lower()
        assert "homozygous" in combined.lower() or "homozygote" in combined.lower()


# ---------------------------------------------------------------------------
# BS3
# ---------------------------------------------------------------------------


class TestBS3:
    """BS3 — well-established functional study shows no damaging effect."""

    def test_applies_when_functional_study_benign(self) -> None:
        """functional_study_result='benign' → BS3 applies."""
        variant = VariantInput(
            chrom="13",
            pos=32_315_480,
            ref="A",
            alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA2",
            transcript_id="NM_000059.4",
            gnomad_af=0.0001,
            functional_study_result="benign",
        )
        result = rule_bs3(variant)
        assert result.applies is True
        assert result.rule_id == "BS3"
        assert result.strength == EvidenceStrength.STRONG_BENIGN

    def test_not_applies_when_functional_study_lof(self) -> None:
        """functional_study_result='loss_of_function' → BS3 does not apply."""
        variant = VariantInput(
            chrom="13",
            pos=32_315_480,
            ref="A",
            alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA2",
            transcript_id="NM_000059.4",
            functional_study_result="loss_of_function",
        )
        result = rule_bs3(variant)
        assert result.applies is False

    def test_not_applies_when_functional_study_none(
        self, rare_missense: VariantInput
    ) -> None:
        """No functional study → BS3 does not apply."""
        result = rule_bs3(rare_missense)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BS4
# ---------------------------------------------------------------------------


class TestBS4:
    """BS4 — lack of segregation in affected family members."""

    def test_applies_when_lack_of_segregation_true(
        self, rare_missense: VariantInput
    ) -> None:
        """lack_of_segregation=True → BS4 applies."""
        result = rule_bs4(rare_missense, lack_of_segregation=True)
        assert result.applies is True
        assert result.rule_id == "BS4"
        assert result.strength == EvidenceStrength.STRONG_BENIGN

    def test_not_applies_when_segregation_supports_pathogenicity(
        self, rare_missense: VariantInput
    ) -> None:
        """lack_of_segregation=False → BS4 does not apply."""
        result = rule_bs4(rare_missense, lack_of_segregation=False)
        assert result.applies is False

    def test_evidence_mentions_segregation(self, rare_missense: VariantInput) -> None:
        """Evidence item mentions lack of segregation."""
        result = rule_bs4(rare_missense, lack_of_segregation=True)
        assert any("segregat" in e.lower() for e in result.evidence_items)


# ---------------------------------------------------------------------------
# BP1
# ---------------------------------------------------------------------------


class TestBP1:
    """BP1 — missense in gene where only truncating variants cause disease."""

    def test_applies_for_missense_in_truncating_only_gene(
        self,
        missense_variant_for_bp1: VariantInput,
        missense_transcript: TranscriptData,
        truncating_only_gene: GeneData,
    ) -> None:
        """Missense in truncating-only gene → BP1 applies."""
        result = rule_bp1(missense_variant_for_bp1, missense_transcript, truncating_only_gene)
        assert result.applies is True
        assert result.rule_id == "BP1"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_not_applies_for_missense_in_missense_only_gene(
        self,
        missense_variant_for_bp1: VariantInput,
        missense_transcript: TranscriptData,
    ) -> None:
        """Missense in missense-only gene → BP1 does not apply."""
        gene = GeneData(
            gene_symbol="TP53",
            lof_is_disease_mechanism=True,
            missense_only_gene=True,  # Only missense causes disease — BP1 does NOT apply
        )
        result = rule_bp1(missense_variant_for_bp1, missense_transcript, gene)
        assert result.applies is False

    def test_not_applies_for_non_missense_variant(
        self,
        truncating_only_gene: GeneData,
        missense_transcript: TranscriptData,
    ) -> None:
        """Frameshift variant → BP1 does not apply (only missense check)."""
        frameshift = VariantInput(
            chrom="5",
            pos=112_175_770,
            ref="G",
            alt="GA",
            variant_type=VariantType.FRAMESHIFT,
            gene_symbol="APC",
            transcript_id="NM_000038.6",
        )
        result = rule_bp1(frameshift, missense_transcript, truncating_only_gene)
        assert result.applies is False

    def test_not_applies_when_no_aa_change(
        self,
        truncating_only_gene: GeneData,
    ) -> None:
        """SNV with no aa_change → BP1 does not apply."""
        variant = VariantInput(
            chrom="5",
            pos=112_175_770,
            ref="G",
            alt="A",
            variant_type="snv",
            gene_symbol="APC",
            transcript_id="NM_000038.6",
        )
        transcript = TranscriptData(
            transcript_id="NM_000038.6",
            gene_symbol="APC",
            aa_change=None,  # No amino acid change
        )
        result = rule_bp1(variant, transcript, truncating_only_gene)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP2
# ---------------------------------------------------------------------------


class TestBP2:
    """BP2 — in trans with pathogenic variant in dominant disorder."""

    def test_applies_when_in_trans_with_pathogenic_dominant(
        self, rare_missense: VariantInput
    ) -> None:
        """In trans with pathogenic in dominant disorder → BP2 applies."""
        result = rule_bp2(rare_missense, in_trans_with_pathogenic_dominant=True)
        assert result.applies is True
        assert result.rule_id == "BP2"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_not_applies_when_not_in_trans(self, rare_missense: VariantInput) -> None:
        """Not in trans with pathogenic dominant → BP2 does not apply."""
        result = rule_bp2(rare_missense, in_trans_with_pathogenic_dominant=False)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP3
# ---------------------------------------------------------------------------


class TestBP3:
    """BP3 — in-frame indel in repeat region without known function."""

    def test_applies_for_inframe_insertion_in_repeat_region(
        self, inframe_variant: VariantInput
    ) -> None:
        """In-frame insertion in repeat region → BP3 applies."""
        transcript = TranscriptData(
            transcript_id="NM_000051.4",
            gene_symbol="ATM",
        )
        result = rule_bp3(inframe_variant, transcript, in_repeat_region=True)
        assert result.applies is True
        assert result.rule_id == "BP3"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_applies_for_inframe_deletion_in_repeat_region(self) -> None:
        """In-frame deletion in repeat region → BP3 applies."""
        variant = VariantInput(
            chrom="11",
            pos=108_111_710,
            ref="ATGCTG",
            alt="ATG",
            variant_type="inframe_deletion",
            gene_symbol="ATM",
            transcript_id="NM_000051.4",
        )
        transcript = TranscriptData(
            transcript_id="NM_000051.4",
            gene_symbol="ATM",
        )
        result = rule_bp3(variant, transcript, in_repeat_region=True)
        assert result.applies is True

    def test_not_applies_when_not_in_repeat_region(
        self, inframe_variant: VariantInput
    ) -> None:
        """In-frame insertion NOT in repeat region → BP3 does not apply."""
        transcript = TranscriptData(
            transcript_id="NM_000051.4",
            gene_symbol="ATM",
        )
        result = rule_bp3(inframe_variant, transcript, in_repeat_region=False)
        assert result.applies is False

    def test_not_applies_for_frameshift(self) -> None:
        """Frameshift (out-of-frame) → BP3 does not apply."""
        variant = VariantInput(
            chrom="11",
            pos=108_111_710,
            ref="G",
            alt="GA",
            variant_type=VariantType.FRAMESHIFT,
            gene_symbol="ATM",
            transcript_id="NM_000051.4",
        )
        transcript = TranscriptData(transcript_id="NM_000051.4", gene_symbol="ATM")
        result = rule_bp3(variant, transcript, in_repeat_region=True)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP4 extra
# ---------------------------------------------------------------------------


class TestBP4Extra:
    """BP4 — additional branches: splice variant and REVEL secondary predictor."""

    def test_not_applies_for_splice_canonical_variant(self) -> None:
        """Splice canonical variant → BP4 defers to splicing.py; applies=False."""
        variant = VariantInput(
            chrom="17",
            pos=41_215_918,
            ref="G",
            alt="A",
            variant_type="splice_canonical",
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
        )
        result = rule_bp4(variant)
        assert result.applies is False
        assert "splicing" in " ".join(result.evidence_items).lower()

    def test_not_applies_for_splice_region_variant(self) -> None:
        """Splice region variant → BP4 defers to splicing.py; applies=False."""
        variant = VariantInput(
            chrom="17",
            pos=41_215_920,
            ref="G",
            alt="A",
            variant_type="splice_region",
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
        )
        result = rule_bp4(variant)
        assert result.applies is False

    def test_applies_via_revel_secondary_when_alphamissense_none(
        self, rare_missense: VariantInput
    ) -> None:
        """REVEL < 0.15, no AlphaMissense → BP4 applies via secondary predictor."""
        result = rule_bp4(rare_missense, alphamissense_score=None, revel_score=0.08)
        assert result.applies is True
        assert result.rule_id == "BP4"
        assert "REVEL" in " ".join(result.evidence_items)

    def test_revel_above_threshold_does_not_apply(
        self, rare_missense: VariantInput
    ) -> None:
        """REVEL = 0.20 (≥ 0.15) → BP4 does not apply via REVEL."""
        result = rule_bp4(rare_missense, alphamissense_score=None, revel_score=0.20)
        assert result.applies is False

    def test_not_applies_when_no_scores_provided(
        self, rare_missense: VariantInput
    ) -> None:
        """No scores provided → BP4 does not apply."""
        result = rule_bp4(rare_missense)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP5
# ---------------------------------------------------------------------------


class TestBP5:
    """BP5 — variant found in case with alternate molecular basis for disease."""

    def test_applies_when_alternate_basis_found(
        self, rare_missense: VariantInput
    ) -> None:
        """Alternate molecular basis found → BP5 applies."""
        result = rule_bp5(rare_missense, alternate_molecular_basis_found=True)
        assert result.applies is True
        assert result.rule_id == "BP5"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_not_applies_when_no_alternate_basis(
        self, rare_missense: VariantInput
    ) -> None:
        """No alternate basis found → BP5 does not apply."""
        result = rule_bp5(rare_missense, alternate_molecular_basis_found=False)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP6
# ---------------------------------------------------------------------------


class TestBP6:
    """BP6 — reputable source reported as benign (ClinVar ≥2 stars + B/LB)."""

    def test_applies_when_clinvar_benign_two_stars(self) -> None:
        """ClinVar Benign with 2 stars → BP6 applies."""
        variant = VariantInput(
            chrom="12",
            pos=111_803_912,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="OAS1",
            transcript_id="NM_016816.4",
            clinvar_stars=2,
            clinvar_classification="Benign",
        )
        result = rule_bp6(variant)
        assert result.applies is True
        assert result.rule_id == "BP6"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_applies_for_likely_benign_three_stars(self) -> None:
        """ClinVar Likely benign with 3 stars → BP6 applies."""
        variant = VariantInput(
            chrom="12",
            pos=111_803_912,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="OAS1",
            transcript_id="NM_016816.4",
            clinvar_stars=3,
            clinvar_classification="Likely benign",
        )
        result = rule_bp6(variant)
        assert result.applies is True

    def test_applies_for_benign_likely_benign_combined(self) -> None:
        """ClinVar Benign/Likely benign with 2 stars → BP6 applies."""
        variant = VariantInput(
            chrom="12",
            pos=111_803_912,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="OAS1",
            transcript_id="NM_016816.4",
            clinvar_stars=2,
            clinvar_classification="Benign/Likely benign",
        )
        result = rule_bp6(variant)
        assert result.applies is True

    def test_not_applies_when_only_one_star(self) -> None:
        """ClinVar Benign with only 1 star → BP6 does not apply (need ≥2)."""
        variant = VariantInput(
            chrom="12",
            pos=111_803_912,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="OAS1",
            transcript_id="NM_016816.4",
            clinvar_stars=1,
            clinvar_classification="Benign",
        )
        result = rule_bp6(variant)
        assert result.applies is False

    def test_not_applies_when_clinvar_pathogenic(self, rare_missense: VariantInput) -> None:
        """ClinVar Pathogenic (not Benign) → BP6 does not apply."""
        variant = VariantInput(
            chrom="17",
            pos=7_674_220,
            ref="C",
            alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TP53",
            transcript_id="NM_000546.6",
            clinvar_stars=3,
            clinvar_classification="Pathogenic",
        )
        result = rule_bp6(variant)
        assert result.applies is False

    def test_not_applies_when_no_clinvar_data(self, rare_missense: VariantInput) -> None:
        """No ClinVar classification → BP6 does not apply."""
        result = rule_bp6(rare_missense)
        assert result.applies is False


# ---------------------------------------------------------------------------
# BP7 extra
# ---------------------------------------------------------------------------


class TestBP7Extra:
    """BP7 — additional branches not in test_rules_benign.py."""

    def test_applies_when_no_splice_impact_flag_set(self) -> None:
        """Synonymous + no_splice_impact=True → BP7 applies (flag-based)."""
        variant = VariantInput(
            chrom="17",
            pos=41_267_742,
            ref="C",
            alt="T",
            variant_type=VariantType.SYNONYMOUS,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
        )
        result = rule_bp7(variant, no_splice_impact=True)
        assert result.applies is True
        assert result.rule_id == "BP7"

    def test_applies_using_spliceai_score_alias(self) -> None:
        """spliceai_score= kwarg (alias for spliceai_delta) → BP7 applies."""
        variant = VariantInput(
            chrom="17",
            pos=41_267_742,
            ref="C",
            alt="T",
            variant_type=VariantType.SYNONYMOUS,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
        )
        result = rule_bp7(variant, spliceai_score=0.04)
        assert result.applies is True

    def test_applies_using_variant_spliceai_delta(self) -> None:
        """spliceai_delta on variant object itself → BP7 resolves no_splice_impact."""
        variant = VariantInput(
            chrom="17",
            pos=41_267_742,
            ref="C",
            alt="T",
            variant_type=VariantType.SYNONYMOUS,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            spliceai_delta=0.05,  # < 0.1 → no splice impact
        )
        result = rule_bp7(variant)
        assert result.applies is True

    def test_not_applies_when_no_splice_impact_false(self) -> None:
        """Synonymous + no_splice_impact=False → BP7 does not apply."""
        variant = VariantInput(
            chrom="17",
            pos=41_267_742,
            ref="C",
            alt="T",
            variant_type=VariantType.SYNONYMOUS,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
        )
        result = rule_bp7(variant, no_splice_impact=False)
        assert result.applies is False
