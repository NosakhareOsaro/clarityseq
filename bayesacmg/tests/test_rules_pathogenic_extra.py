"""
Extended tests for BayesACMG pathogenic rules.

Covers branches not exercised by test_rules_pathogenic.py:
    - PS1 applies=True (same amino acid change)
    - PS2 applies=True (de novo)
    - PS3 applies=True (functional study)
    - PS4 applies=True (case-control significant)
    - PM1 applies=True (hotspot domain)
    - PM2 applies=True (absent from gnomAD, extremely rare)
    - PM3 applies=True (in trans with pathogenic)
    - PM4 applies=True (in-frame indel protein length change)
    - PM5 applies=True (different amino acid same position)
    - PM6 applies=True (assumed de novo)
    - PP1 applies=True (cosegregation)
    - PP2 applies=True (missense in constraint gene)
    - PP3 applies=True (REVEL/AlphaMissense)
    - PP4 applies=True (phenotype matches disease)
    - PP5 applies=True (ClinVar pathogenic assertion)
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
from bayesacmg.rules.pathogenic import (
    rule_bp4_from_alphamissense,
    rule_pm1,
    rule_pm2,
    rule_pm3,
    rule_pm4,
    rule_pm5,
    rule_pm6,
    rule_pp1,
    rule_pp2,
    rule_pp3,
    rule_pp4,
    rule_pp5,
    rule_ps1,
    rule_ps2,
    rule_ps3,
    rule_ps4,
    rule_pvs1,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def de_novo_variant() -> VariantInput:
    """De novo missense variant in a neurodevelopmental gene."""
    return VariantInput(
        chrom="2",
        pos=165_990_000,
        ref="C",
        alt="T",
        variant_type=VariantType.MISSENSE,
        gene_symbol="SCN1A",
        transcript_id="NM_006920.6",
        hgvsc="NM_006920.6:c.4010C>T",
        hgvsp="NM_006920.6(SCN1A):p.Pro1337Leu",
        hgvs_p="p.Pro1337Leu",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
        is_de_novo=True,
        alphamissense_score=0.92,
    )


@pytest.fixture
def in_trans_variant() -> VariantInput:
    """Variant in trans with pathogenic variant (AR disease model)."""
    return VariantInput(
        chrom="7",
        pos=117_548_628,
        ref="G",
        alt="A",
        variant_type=VariantType.MISSENSE,
        gene_symbol="CFTR",
        transcript_id="NM_000492.4",
        hgvsc="NM_000492.4:c.1521_1523delCTT",
        hgvsp="NM_000492.4(CFTR):p.Phe508del",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
        in_trans_pathogenic=True,
    )


@pytest.fixture
def functional_damaging_variant() -> VariantInput:
    """Variant with confirmed functional loss-of-function study."""
    return VariantInput(
        chrom="17",
        pos=7_673_802,
        ref="G",
        alt="A",
        variant_type=VariantType.MISSENSE,
        gene_symbol="TP53",
        transcript_id="NM_000546.6",
        hgvsc="NM_000546.6:c.413G>A",
        hgvsp="NM_000546.6(TP53):p.Arg138His",
        gnomad_af=0.0,
        gnomad_ac=0,
        gnomad_nhomalt=0,
        functional_study_result="loss_of_function",
        alphamissense_score=0.95,
    )


@pytest.fixture
def hotspot_transcript() -> TranscriptData:
    """Transcript with hotspot domain annotation."""
    return TranscriptData(
        transcript_id="NM_000546.6",
        gene_symbol="TP53",
        is_mane_select=True,
        is_canonical=True,
        lof_disease_mechanism=True,
        exon_count=11,
        domain_annotations=["DNA-binding domain (aa 94-292)"],
    )


@pytest.fixture
def tp53_gene() -> GeneData:
    """TP53 gene data with hotspot domain."""
    return GeneData(
        gene_symbol="TP53",
        omim_id="OMIM:191170",
        lof_mechanism=True,
        missense_constraint_z=4.5,
        pli=0.999,
        has_hotspot_domain=True,
        hotspot_domains=["DNA-binding domain"],
        has_vcep_specification=True,
    )


# ---------------------------------------------------------------------------
# PS1 tests
# ---------------------------------------------------------------------------


class TestPS1:
    """Tests for rule_ps1()."""

    def test_applies_when_same_aa_change_pathogenic(self, de_novo_variant: VariantInput) -> None:
        """PS1 applies when same amino acid change is established pathogenic."""
        result = rule_ps1(de_novo_variant, same_aa_change_pathogenic=True)
        assert result.applies is True
        assert result.rule_id == "PS1"
        assert result.strength == EvidenceStrength.STRONG

    def test_does_not_apply_when_no_same_aa_change(self, de_novo_variant: VariantInput) -> None:
        """PS1 does not apply when no same amino acid change is pathogenic."""
        result = rule_ps1(de_novo_variant, same_aa_change_pathogenic=False)
        assert result.applies is False

    def test_does_not_apply_without_hgvsp(self) -> None:
        """PS1 does not apply when variant has no HGVSp annotation."""
        variant = VariantInput(
            chrom="17", pos=1000, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TEST",
            gnomad_af=0.0,
        )
        result = rule_ps1(variant, same_aa_change_pathogenic=True)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PS2 tests
# ---------------------------------------------------------------------------


class TestPS2:
    """Tests for rule_ps2()."""

    def test_applies_when_de_novo(self, de_novo_variant: VariantInput) -> None:
        """PS2 applies when variant is confirmed de novo."""
        result = rule_ps2(de_novo_variant)
        assert result.applies is True
        assert result.rule_id == "PS2"
        assert result.strength == EvidenceStrength.STRONG

    def test_does_not_apply_when_not_de_novo(self, in_trans_variant: VariantInput) -> None:
        """PS2 does not apply when variant is not confirmed de novo."""
        result = rule_ps2(in_trans_variant)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PS3 tests
# ---------------------------------------------------------------------------


class TestPS3:
    """Tests for rule_ps3()."""

    def test_applies_for_loss_of_function_study(
        self, functional_damaging_variant: VariantInput
    ) -> None:
        """PS3 applies when functional study shows loss-of-function."""
        result = rule_ps3(functional_damaging_variant)
        assert result.applies is True
        assert result.rule_id == "PS3"

    def test_applies_for_dominant_negative(self) -> None:
        """PS3 applies for dominant negative functional study result."""
        variant = VariantInput(
            chrom="17", pos=7_673_802, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TP53",
            gnomad_af=0.0,
            functional_study_result="dominant_negative",
        )
        result = rule_ps3(variant)
        assert result.applies is True

    def test_does_not_apply_when_no_study(self, de_novo_variant: VariantInput) -> None:
        """PS3 does not apply when no functional study is available."""
        result = rule_ps3(de_novo_variant)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PS4 tests
# ---------------------------------------------------------------------------


class TestPS4:
    """Tests for rule_ps4()."""

    def test_applies_when_case_control_significant(self, de_novo_variant: VariantInput) -> None:
        """PS4 applies when case-control evidence is statistically significant."""
        result = rule_ps4(de_novo_variant, case_control_or_significant=True)
        assert result.applies is True
        assert result.rule_id == "PS4"

    def test_does_not_apply_when_not_significant(self, de_novo_variant: VariantInput) -> None:
        """PS4 does not apply when case-control evidence is insufficient."""
        result = rule_ps4(de_novo_variant, case_control_or_significant=False)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PM1 tests
# ---------------------------------------------------------------------------


class TestPM1:
    """Tests for rule_pm1()."""

    def test_applies_when_in_hotspot_domain(
        self,
        functional_damaging_variant: VariantInput,
        hotspot_transcript: TranscriptData,
        tp53_gene: GeneData,
    ) -> None:
        """PM1 applies when variant is in a mutational hotspot domain."""
        result = rule_pm1(functional_damaging_variant, hotspot_transcript, tp53_gene)
        assert result.applies is True
        assert result.rule_id == "PM1"
        assert result.strength == EvidenceStrength.MODERATE

    def test_does_not_apply_when_no_hotspot(
        self,
        de_novo_variant: VariantInput,
        hotspot_transcript: TranscriptData,
    ) -> None:
        """PM1 does not apply when gene has no hotspot domain."""
        gene = GeneData(
            gene_symbol="SCN1A",
            omim_id="OMIM:182389",
            has_hotspot_domain=False,
        )
        result = rule_pm1(de_novo_variant, hotspot_transcript, gene)
        assert result.applies is False

    def test_does_not_apply_when_no_domain_annotations(
        self,
        functional_damaging_variant: VariantInput,
        tp53_gene: GeneData,
    ) -> None:
        """PM1 does not apply when transcript has no domain annotations."""
        transcript = TranscriptData(
            transcript_id="NM_000546.6",
            gene_symbol="TP53",
            is_mane_select=True,
            domain_annotations=[],
        )
        result = rule_pm1(functional_damaging_variant, transcript, tp53_gene)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PM2 tests
# ---------------------------------------------------------------------------


class TestPM2Extra:
    """Additional tests for rule_pm2() — extremely rare branch."""

    def test_applies_when_af_is_none(self) -> None:
        """PM2 applies at Supporting when gnomad_af is None (absent)."""
        from bayesacmg.rules.pathogenic import rule_pm2

        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=None,
        )
        result = rule_pm2(variant)
        assert result.applies is True
        assert result.strength == EvidenceStrength.SUPPORTING

    def test_applies_when_very_rare(self) -> None:
        """PM2 applies at Supporting for extremely rare variants (AF < 0.0001)."""
        from bayesacmg.rules.pathogenic import rule_pm2

        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.00005,  # Below PM2 threshold
        )
        result = rule_pm2(variant)
        assert result.applies is True
        assert result.strength == EvidenceStrength.SUPPORTING


# ---------------------------------------------------------------------------
# PM3 tests
# ---------------------------------------------------------------------------


class TestPM3:
    """Tests for rule_pm3()."""

    def test_applies_when_in_trans(self, in_trans_variant: VariantInput) -> None:
        """PM3 applies when variant is in trans with a pathogenic variant."""
        result = rule_pm3(in_trans_variant)
        assert result.applies is True
        assert result.rule_id == "PM3"
        assert result.strength == EvidenceStrength.MODERATE

    def test_does_not_apply_when_not_in_trans(self, de_novo_variant: VariantInput) -> None:
        """PM3 does not apply when variant is not in trans with pathogenic."""
        result = rule_pm3(de_novo_variant)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PM6 tests
# ---------------------------------------------------------------------------


class TestPM6:
    """Tests for rule_pm6()."""

    def test_applies_for_assumed_de_novo(self) -> None:
        """PM6 applies when de novo status is assumed (not confirmed)."""
        variant = VariantInput(
            chrom="2", pos=165990000, ref="C", alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="SCN1A",
            gnomad_af=0.0,
            assumed_de_novo=True,  # Assumed but not confirmed
        )
        result = rule_pm6(variant)
        assert result.applies is True
        assert result.rule_id == "PM6"

    def test_does_not_apply_when_not_assumed_de_novo(self, in_trans_variant: VariantInput) -> None:
        """PM6 does not apply when assumed_de_novo is False."""
        result = rule_pm6(in_trans_variant)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PP1 tests
# ---------------------------------------------------------------------------


class TestPP1:
    """Tests for rule_pp1()."""

    def test_applies_when_cosegregation_confirmed(self, de_novo_variant: VariantInput) -> None:
        """PP1 applies when variant cosegregates with disease in family."""
        result = rule_pp1(de_novo_variant, segregation_supports=True)
        assert result.applies is True
        assert result.rule_id == "PP1"

    def test_does_not_apply_without_cosegregation(self, de_novo_variant: VariantInput) -> None:
        """PP1 does not apply without cosegregation data."""
        result = rule_pp1(de_novo_variant, segregation_supports=False)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PP2 tests
# ---------------------------------------------------------------------------


class TestPP2:
    """Tests for rule_pp2()."""

    def test_applies_for_missense_in_constrained_gene(self) -> None:
        """PP2 applies for missense (snv) in gene with high pLI."""
        variant = VariantInput(
            chrom="2", pos=165990000, ref="C", alt="T",
            variant_type="snv",
            gene_symbol="SCN1A",
            gnomad_af=0.0,
        )
        transcript = TranscriptData(
            transcript_id="NM_006920.6",
            gene_symbol="SCN1A",
            is_mane_select=True,
            aa_change="p.Pro1337Leu",
        )
        gene = GeneData(
            gene_symbol="SCN1A",
            omim_id="OMIM:182389",
            gnomad_pli=0.999,  # High pLI
        )
        result = rule_pp2(variant, transcript, gene)
        assert result.applies is True
        assert result.rule_id == "PP2"

    def test_does_not_apply_for_non_snv(self, in_trans_variant: VariantInput) -> None:
        """PP2 does not apply when variant is not an SNV missense."""
        transcript = TranscriptData(
            transcript_id="NM_000492.4",
            gene_symbol="CFTR",
            is_mane_select=True,
        )
        gene = GeneData(gene_symbol="CFTR", gnomad_pli=0.1)
        result = rule_pp2(in_trans_variant, transcript, gene)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PP3 tests
# ---------------------------------------------------------------------------


class TestPP3Extra:
    """Extended tests for rule_pp3() with different predictors."""

    def test_applies_with_high_alphamissense(
        self, functional_damaging_variant: VariantInput
    ) -> None:
        """PP3 applies when AlphaMissense score is high (≥ 0.564)."""
        result = rule_pp3(functional_damaging_variant, alphamissense_score=0.95, revel_score=None)
        assert result.applies is True
        assert result.rule_id == "PP3"

    def test_applies_with_high_revel(self) -> None:
        """PP3 applies when REVEL score is high (≥ 0.7)."""
        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.0,
        )
        result = rule_pp3(variant, alphamissense_score=None, revel_score=0.85)
        assert result.applies is True


# ---------------------------------------------------------------------------
# PP4 tests
# ---------------------------------------------------------------------------


class TestPP4:
    """Tests for rule_pp4()."""

    def test_applies_when_phenotype_specific_to_disease(
        self, de_novo_variant: VariantInput
    ) -> None:
        """PP4 applies when phenotype is highly specific for the gene's disease."""
        result = rule_pp4(
            de_novo_variant,
            phenotype_highly_specific=True,
        )
        assert result.applies is True
        assert result.rule_id == "PP4"

    def test_does_not_apply_when_nonspecific_phenotype(
        self, de_novo_variant: VariantInput
    ) -> None:
        """PP4 does not apply when phenotype is not specific enough."""
        result = rule_pp4(de_novo_variant, phenotype_highly_specific=False)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PP5 tests
# ---------------------------------------------------------------------------


class TestPP5:
    """Tests for rule_pp5()."""

    def test_applies_when_clinvar_pathogenic(self) -> None:
        """PP5 applies when ClinVar has pathogenic assertion from reputable source."""
        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.0,
            clinvar_classification="Pathogenic",
            clinvar_stars=3,
        )
        result = rule_pp5(variant)
        assert result.applies is True
        assert result.rule_id == "PP5"

    def test_does_not_apply_when_no_clinvar(self, de_novo_variant: VariantInput) -> None:
        """PP5 does not apply when ClinVar has no assertion."""
        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.0,
            clinvar_classification=None,
        )
        result = rule_pp5(variant)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PVS1 additional branches (not covered by TestPVS1 in test_rules_pathogenic.py)
# ---------------------------------------------------------------------------


class TestPVS1Extra:
    """Extra PVS1 branches: LoF-typed variant but no disease mechanism; NMD escape."""

    def test_lof_type_but_not_disease_mechanism(self) -> None:
        """A LoF-typed variant (frameshift) in a gene where LoF is NOT the
        established disease mechanism → PVS1 does not apply.

        This differs from the existing non-LoF-gene test in
        test_rules_pathogenic.py, which uses a missense (non-LoF-type)
        variant and therefore short-circuits earlier in the decision tree.
        """
        variant = VariantInput(
            chrom="11",
            pos=2_181_009,
            ref="C",
            alt="CC",
            variant_type=VariantType.FRAMESHIFT,
            gene_symbol="KCNQ1",
        )
        transcript = TranscriptData(
            transcript_id="NM_000218.3",
            gene_symbol="KCNQ1",
            is_mane_select=True,
            lof_disease_mechanism=False,
        )
        gene = GeneData(
            gene_symbol="KCNQ1",
            lof_is_disease_mechanism=False,
            lof_mechanism=False,  # Gain-of-function disease; LoF not established
        )
        result = rule_pvs1(variant, transcript, gene)
        assert result.applies is False
        assert "not established disease mechanism" in " ".join(result.evidence_items).lower()

    def test_nmd_escape_last_exon_reduces_to_strong(
        self, tp53_gene: GeneData
    ) -> None:
        """MANE Select transcript, but variant is in the last exon (NMD escape)
        → PVS1 reduced from Very Strong to Strong (Abou Tayoun 2018 decision tree)."""
        variant = VariantInput(
            chrom="17",
            pos=7_669_690,
            ref="C",
            alt="CC",
            variant_type=VariantType.FRAMESHIFT,
            gene_symbol="TP53",
        )
        transcript = TranscriptData(
            transcript_id="NM_000546.6",
            gene_symbol="TP53",
            is_mane_select=True,
            lof_disease_mechanism=True,
            is_last_exon=True,
        )
        result = rule_pvs1(variant, transcript, tp53_gene)
        assert result.applies is True
        assert result.strength == EvidenceStrength.STRONG
        assert "nmd" in " ".join(result.evidence_items).lower()

    def test_nmd_escapes_last_exon_rule_flag_reduces_to_strong(
        self, tp53_gene: GeneData
    ) -> None:
        """nmd_escapes_last_exon_rule=True (not last exon itself) also reduces to Strong."""
        variant = VariantInput(
            chrom="17",
            pos=7_669_690,
            ref="C",
            alt="CC",
            variant_type=VariantType.NONSENSE,
            gene_symbol="TP53",
        )
        transcript = TranscriptData(
            transcript_id="NM_000546.6",
            gene_symbol="TP53",
            is_mane_select=True,
            lof_disease_mechanism=True,
            is_last_exon=False,
            nmd_escapes_last_exon_rule=True,  # <55 nt from final junction
        )
        result = rule_pvs1(variant, transcript, tp53_gene)
        assert result.applies is True
        assert result.strength == EvidenceStrength.STRONG


# ---------------------------------------------------------------------------
# PM4 tests
# ---------------------------------------------------------------------------


class TestPM4:
    """Tests for rule_pm4() — in-frame indel / stop-loss protein length change."""

    def test_does_not_apply_for_non_inframe_type(self, de_novo_variant: VariantInput) -> None:
        """Missense variant is not an in-frame indel or stop-loss → PM4 does not apply."""
        transcript = TranscriptData(transcript_id="NM_006920.6", gene_symbol="SCN1A")
        result = rule_pm4(de_novo_variant, transcript)
        assert result.applies is False
        assert result.rule_id == "PM4"

    def test_applies_when_length_change_meets_threshold(self) -> None:
        """In-frame deletion changing protein length by >=10 aa → PM4 applies."""
        variant = VariantInput(
            chrom="7",
            pos=117_548_628,
            ref="ATGATGATGATGATGATGATGATGATG",
            alt="A",
            variant_type="inframe_deletion",
            gene_symbol="CFTR",
        )
        transcript = TranscriptData(
            transcript_id="NM_000492.4",
            gene_symbol="CFTR",
            prot_length_original=1480,
            prot_length_alt=1468,  # delta = 12 aa >= 10
        )
        result = rule_pm4(variant, transcript)
        assert result.applies is True
        assert "12 aa" in " ".join(result.evidence_items)

    def test_applies_provisionally_when_length_data_unavailable(self) -> None:
        """In-frame indel with no protein length data → PM4 applies provisionally."""
        variant = VariantInput(
            chrom="7",
            pos=117_548_628,
            ref="ATG",
            alt="ATGCTG",
            variant_type="inframe_insertion",
            gene_symbol="CFTR",
        )
        transcript = TranscriptData(
            transcript_id="NM_000492.4",
            gene_symbol="CFTR",
            prot_length_original=None,
            prot_length_alt=None,
        )
        result = rule_pm4(variant, transcript)
        assert result.applies is True
        assert "provisionally" in (result.notes or "").lower()

    def test_applies_provisionally_when_length_change_below_threshold(self) -> None:
        """In-frame indel with a small (<10 aa) protein length change still applies (provisional)."""
        variant = VariantInput(
            chrom="7",
            pos=117_548_628,
            ref="ATG",
            alt="ATGCTG",
            variant_type="stop_loss",
            gene_symbol="CFTR",
        )
        transcript = TranscriptData(
            transcript_id="NM_000492.4",
            gene_symbol="CFTR",
            prot_length_original=1480,
            prot_length_alt=1485,  # delta = 5 aa < 10
        )
        result = rule_pm4(variant, transcript)
        assert result.applies is True
        assert "uncertain" in " ".join(result.evidence_items).lower()


# ---------------------------------------------------------------------------
# PM5 tests
# ---------------------------------------------------------------------------


class TestPM5:
    """Tests for rule_pm5() — novel missense at same position as known pathogenic."""

    def test_applies_when_different_aa_pathogenic(self, de_novo_variant: VariantInput) -> None:
        """A pathogenic missense at the same codon (different AA change) → PM5 applies."""
        result = rule_pm5(de_novo_variant, different_aa_at_same_position_pathogenic=True)
        assert result.applies is True
        assert result.rule_id == "PM5"
        assert result.strength == EvidenceStrength.MODERATE

    def test_does_not_apply_when_no_known_pathogenic_at_position(
        self, de_novo_variant: VariantInput
    ) -> None:
        """No known pathogenic missense at the same position → PM5 does not apply."""
        result = rule_pm5(de_novo_variant, different_aa_at_same_position_pathogenic=False)
        assert result.applies is False

    def test_does_not_apply_without_hgvsp(self) -> None:
        """Even if flagged pathogenic-at-position, missing HGVSp → PM5 does not apply."""
        variant = VariantInput(
            chrom="17", pos=1000, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TEST",
            gnomad_af=0.0,
        )
        result = rule_pm5(variant, different_aa_at_same_position_pathogenic=True)
        assert result.applies is False


# ---------------------------------------------------------------------------
# PP1 extra: segregation LOD score included in evidence
# ---------------------------------------------------------------------------


class TestPP1Extra:
    """PP1 evidence message includes LOD score when available."""

    def test_evidence_includes_lod_score_when_present(self) -> None:
        """When segregation_lod is set, the evidence message reports it."""
        variant = VariantInput(
            chrom="2", pos=165_990_000, ref="C", alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="SCN1A",
            gnomad_af=0.0,
            segregation_lod=3.5,
        )
        result = rule_pp1(variant, segregation_supports=True)
        assert result.applies is True
        assert "LOD=3.50" in " ".join(result.evidence_items)


# ---------------------------------------------------------------------------
# PP2 extra: LOEUF-only branch and "not LoF-intolerant" branch
# ---------------------------------------------------------------------------


class TestPP2Extra:
    """PP2 branches not covered by TestPP2: LOEUF-only qualification; no constraint."""

    def test_applies_via_loeuf_only(self) -> None:
        """PP2 applies when LOEUF is low, even if pLI is not (or is unavailable)."""
        variant = VariantInput(
            chrom="2", pos=165_990_000, ref="C", alt="T",
            variant_type="snv",
            gene_symbol="SCN1A",
            gnomad_af=0.0,
        )
        transcript = TranscriptData(
            transcript_id="NM_006920.6",
            gene_symbol="SCN1A",
            aa_change="p.Pro1337Leu",
        )
        gene = GeneData(
            gene_symbol="SCN1A",
            gnomad_pli=None,  # pLI unavailable
            gnomad_loeuf=0.10,  # <= 0.35 threshold → LOEUF-low
        )
        result = rule_pp2(variant, transcript, gene)
        assert result.applies is True
        assert any("LOEUF" in item for item in result.evidence_items)

    def test_does_not_apply_when_gene_not_constrained(self) -> None:
        """Missense in a gene with neither high pLI nor low LOEUF → PP2 does not apply."""
        variant = VariantInput(
            chrom="2", pos=165_990_000, ref="C", alt="T",
            variant_type="snv",
            gene_symbol="OR4F5",  # olfactory receptor — not LoF-intolerant
            gnomad_af=0.0,
        )
        transcript = TranscriptData(
            transcript_id="NM_001005484.2",
            gene_symbol="OR4F5",
            aa_change="p.Leu20Pro",
        )
        gene = GeneData(
            gene_symbol="OR4F5",
            gnomad_pli=0.02,
            gnomad_loeuf=1.5,
        )
        result = rule_pp2(variant, transcript, gene)
        assert result.applies is False
        assert "not classified as lof-intolerant" in " ".join(result.evidence_items).lower()


# ---------------------------------------------------------------------------
# PP3 extra: splice-blocked, BP4-zone, CADD tertiary, no-evidence fallback
# ---------------------------------------------------------------------------


class TestPP3ExtraBranches:
    """PP3 branches not covered elsewhere: splice block, BP4 zone, CADD, fallback."""

    def test_blocked_for_splice_canonical_variant(self) -> None:
        """PP3 defers to splicing.py for splice_canonical variants."""
        variant = VariantInput(
            chrom="17", pos=41_215_918, ref="G", alt="A",
            variant_type="splice_canonical",
            gene_symbol="BRCA1",
        )
        result = rule_pp3(variant, alphamissense_score=0.9, revel_score=0.9)
        assert result.applies is False
        assert "splicing" in " ".join(result.evidence_items).lower()

    def test_blocked_for_splice_region_variant(self) -> None:
        """PP3 defers to splicing.py for splice_region variants."""
        variant = VariantInput(
            chrom="17", pos=41_215_920, ref="G", alt="A",
            variant_type="splice_region",
            gene_symbol="BRCA1",
        )
        result = rule_pp3(variant, alphamissense_score=0.9, revel_score=None)
        assert result.applies is False

    def test_alphamissense_at_or_below_bp4_threshold_does_not_fire_pp3(
        self, common_benign_missense: VariantInput
    ) -> None:
        """AlphaMissense <= 0.340 → PP3 does not apply (BP4 territory instead)."""
        result = rule_pp3(common_benign_missense, alphamissense_score=0.10, revel_score=None)
        assert result.applies is False
        assert "bp4 applies" in " ".join(result.evidence_items).lower()

    def test_applies_via_cadd_tertiary_predictor(self) -> None:
        """CADD PHRED >= 25 used as tertiary predictor when AM/REVEL unavailable."""
        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.0,
            cadd_phred=31.2,
        )
        result = rule_pp3(variant, alphamissense_score=None, revel_score=None)
        assert result.applies is True
        assert result.rule_id == "PP3"
        assert "CADD" in " ".join(result.evidence_items)
        assert "tertiary" in (result.notes or "").lower()

    def test_no_evidence_at_all_does_not_apply(self) -> None:
        """No AlphaMissense, REVEL, or CADD score → PP3 does not apply (fallback)."""
        variant = VariantInput(
            chrom="17", pos=43044295, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            gnomad_af=0.0,
            cadd_phred=None,
        )
        result = rule_pp3(variant, alphamissense_score=None, revel_score=None)
        assert result.applies is False
        assert "no in silico scores" in " ".join(result.evidence_items).lower()


# ---------------------------------------------------------------------------
# rule_bp4_from_alphamissense extra branches
# ---------------------------------------------------------------------------


class TestBP4FromAlphamissenseExtra:
    """Branches of rule_bp4_from_alphamissense() not covered in test_rules_pathogenic.py."""

    def test_falls_back_to_variant_alphamissense_score(self) -> None:
        """When alphamissense_score kwarg is None, falls back to variant.alphamissense_score."""
        variant = VariantInput(
            chrom="12", pos=111_803_912, ref="G", alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="OAS1",
            alphamissense_score=0.15,  # ≤ 0.340 → BP4 applies
        )
        result = rule_bp4_from_alphamissense(variant)  # no explicit score kwarg
        assert result.applies is True
        assert result.rule_id == "BP4"

    def test_does_not_apply_when_variant_score_above_threshold(self) -> None:
        """Falls back to variant.alphamissense_score, which is above threshold → does not apply."""
        variant = VariantInput(
            chrom="17", pos=7_674_220, ref="C", alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TP53",
            alphamissense_score=0.97,  # > 0.340
        )
        result = rule_bp4_from_alphamissense(variant)
        assert result.applies is False

    def test_does_not_apply_when_no_score_anywhere(self) -> None:
        """Neither kwarg nor variant.alphamissense_score set → BP4 does not apply."""
        variant = VariantInput(
            chrom="1", pos=100, ref="A", alt="G",
            variant_type=VariantType.MISSENSE,
            gene_symbol="TEST",
            alphamissense_score=None,
        )
        result = rule_bp4_from_alphamissense(variant)
        assert result.applies is False
        assert result.rule_id == "BP4"
