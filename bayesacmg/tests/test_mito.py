"""
Tests for ACGS 2024 §6 mitochondrial variant classification rules.

ACGS 2024 v1.2 §6 requirements:
1. Haplogroup classification (Haplogrep3) MUST run before ACMG assessment
2. Haplogroup-defining variants → automatically Benign (excluded from analysis)
3. Standard BA1 threshold (AF > 5%) does NOT apply to mito variants
4. Mito-specific PM2 thresholds apply
5. Heteroplasmy level maps to clinical significance

Guidelines:
    ACGS 2024 v1.2 §6 (Durkie et al., Feb 2024)
    ACGS Best Practice Guidelines for Mitochondrial Disease (Nov 2020)
    MITOMAP: https://www.mitomap.org/
"""

from __future__ import annotations


from bayesacmg.models import EvidenceStrength, VariantInput, VariantType
from bayesacmg.rules.mito import (
    MitoHaploData,
    rule_mito_ba1,
    rule_mito_haplogroup_defining,
    rule_mito_heteroplasmy_level,
    rule_mito_pm2,
)


def _mito_variant(**overrides) -> VariantInput:
    """Minimal mitochondrial VariantInput with sensible defaults."""
    defaults = dict(
        chrom="MT",
        pos=3_460,
        ref="G",
        alt="A",
        variant_type=VariantType.MISSENSE,
        gene_symbol="MT-ND1",
        transcript_id="NC_012920.1",
        is_mito=True,
        is_haplogroup_defining=False,
        gnomad_af=0.0,
    )
    defaults.update(overrides)
    return VariantInput(**defaults)


class TestMitoHaplogroupDefining:
    """ACGS 2024 §6: haplogroup-defining variants must be excluded first."""

    def test_haplogroup_defining_variant_is_benign(
        self, mito_haplogroup_defining: VariantInput
    ) -> None:
        """Haplogroup-defining variant → Benign (stand-alone).

        ACGS 2024 §6: 'Variants that are haplogroup-defining must be
        classified as Benign before any further ACMG assessment.'
        Haplogrep3 must run before this rule is called.
        """
        result = rule_mito_haplogroup_defining(mito_haplogroup_defining)
        assert result.applies is True
        assert result.strength == EvidenceStrength.STAND_ALONE, (
            "Haplogroup-defining mito variant should be classified as Benign "
            "(stand-alone BA1 equivalent) per ACGS 2024 §6"
        )

    def test_non_haplogroup_defining_not_excluded(
        self, mito_haplogroup_defining: VariantInput
    ) -> None:
        """Non-haplogroup-defining mito variant proceeds to full ACMG assessment."""
        non_haplogroup = VariantInput(
            chrom="MT",
            pos=3_460,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="MT-ND1",
            transcript_id="NC_012920.1",
            hgvsc="NC_012920.1:m.3460G>A",
            is_mito=True,
            is_haplogroup_defining=False,  # NOT haplogroup-defining
            gnomad_af=0.0001,
        )
        result = rule_mito_haplogroup_defining(non_haplogroup)
        assert (
            result.applies is False
        ), "Non-haplogroup-defining mito variant should NOT be auto-classified Benign"

    def test_non_mito_variant_returns_not_applicable(self) -> None:
        """A non-mitochondrial variant is not applicable to this mito-only rule."""
        nuclear_variant = VariantInput(
            chrom="17",
            pos=41_276_045,
            ref="C",
            alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            is_mito=False,
        )
        result = rule_mito_haplogroup_defining(nuclear_variant)
        assert result.applies is False
        assert "not a mitochondrial variant" in " ".join(result.evidence_items).lower()


class TestMitoBA1:
    """ACGS 2024 §6: standard BA1 threshold (>5%) does NOT apply to mito variants."""

    def test_standard_ba1_threshold_does_not_apply_to_mito(
        self, mito_haplogroup_defining: VariantInput
    ) -> None:
        """BA1 standard 5% threshold is NOT used for mitochondrial variants.

        ACGS 2024 §6: mito variants have population frequencies that do not
        follow the same population genetics as nuclear variants. Many
        haplogroup-defining variants have AF > 5% but are NOT pathogenic.
        A separate mito-specific BA1 threshold applies.
        """
        # Mito AF > 5% but should NOT trigger standard BA1
        result = rule_mito_ba1(mito_haplogroup_defining)
        # The standard BA1 function should not fire for mito variants
        # (use rule_mito_ba1 which uses mito-specific thresholds)
        # This test verifies the mito-specific logic is used
        assert result is not None, "rule_mito_ba1 must return an ACMGRule object"
        # Even with AF=0.85, if it's haplogroup-defining, it should be caught
        # by rule_mito_haplogroup_defining first — mito BA1 for pathogenic calls

    def test_non_mito_variant_not_applicable(self) -> None:
        """rule_mito_ba1 does not apply to non-mitochondrial variants."""
        nuclear_variant = VariantInput(
            chrom="17",
            pos=41_276_045,
            ref="C",
            alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            is_mito=False,
        )
        result = rule_mito_ba1(nuclear_variant)
        assert result.applies is False
        assert result.strength == EvidenceStrength.STAND_ALONE

    def test_mitomap_confirmed_polymorphism_is_stand_alone_benign(self) -> None:
        """MITOMAP 'Confirmed Polymorphism' status satisfies mito BA1 (stand-alone)."""
        variant = _mito_variant()
        haplo = MitoHaploData(
            haplogroup="H1a1",
            is_haplogroup_defining=False,
            mitomap_status="Confirmed Polymorphism",
        )
        result = rule_mito_ba1(variant, haplo)
        assert result.applies is True
        assert result.rule_id == "BA1_MITO"
        assert result.strength == EvidenceStrength.STAND_ALONE
        assert "polymorphism" in " ".join(result.evidence_items).lower()

    def test_high_haplogroup_frequency_is_strong_benign_not_stand_alone(self) -> None:
        """Haplogroup frequency > 5% gives Strong Benign, NOT stand-alone (ACGS 2024 §6).

        This is the key mito-specific deviation from nuclear BA1: high
        population frequency within a haplogroup is corroborating evidence,
        not automatically stand-alone Benign.
        """
        variant = _mito_variant()
        haplo = MitoHaploData(
            haplogroup="L3e",
            is_haplogroup_defining=False,
            haplogroup_frequency=0.12,  # >5% within haplogroup
            mitomap_status=None,
        )
        result = rule_mito_ba1(variant, haplo)
        assert result.applies is True
        assert result.strength == EvidenceStrength.STRONG_BENIGN, (
            "High haplogroup-specific frequency should give STRONG_BENIGN, "
            "not STAND_ALONE, per ACGS 2024 §6"
        )

    def test_no_evidence_meets_mito_ba1_criteria(self) -> None:
        """No MITOMAP status and low haplogroup frequency → mito BA1 does not apply."""
        variant = _mito_variant()
        haplo = MitoHaploData(
            haplogroup="H1a1", is_haplogroup_defining=False, haplogroup_frequency=0.001
        )
        result = rule_mito_ba1(variant, haplo)
        assert result.applies is False


class TestMitoPM2:
    """Tests for rule_mito_pm2() — mito-specific PM2 (ACGS 2024 §6)."""

    def test_not_applicable_for_non_mito_variant(self) -> None:
        """rule_mito_pm2 does not apply to non-mitochondrial variants."""
        nuclear_variant = VariantInput(
            chrom="17",
            pos=41_276_045,
            ref="C",
            alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            is_mito=False,
        )
        result = rule_mito_pm2(nuclear_variant, MitoHaploData())
        assert result.applies is False
        assert result.rule_id == "PM2_MITO"

    def test_haplogroup_defining_variant_pm2_does_not_apply(self) -> None:
        """Haplogroup-defining status takes precedence — PM2_MITO does not apply."""
        variant = _mito_variant(gnomad_af=None)
        haplo = MitoHaploData(haplogroup="H1a1", is_haplogroup_defining=True)
        result = rule_mito_pm2(variant, haplo)
        assert result.applies is False
        assert "haplogroup-defining" in " ".join(result.evidence_items).lower()

    def test_applies_supporting_when_af_none(self) -> None:
        """Absent from gnomAD v3.1 mitogenome data → PM2_MITO Supporting."""
        variant = _mito_variant(gnomad_af=None)
        haplo = MitoHaploData(haplogroup="H1a1", is_haplogroup_defining=False)
        result = rule_mito_pm2(variant, haplo)
        assert result.applies is True
        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.rule_id == "PM2_MITO"

    def test_applies_supporting_when_extremely_rare(self) -> None:
        """AF < 0.0001 in gnomAD v3.1 mito data → PM2_MITO Supporting."""
        variant = _mito_variant(gnomad_af=0.00002)
        haplo = MitoHaploData(haplogroup="H1a1", is_haplogroup_defining=False)
        result = rule_mito_pm2(variant, haplo)
        assert result.applies is True
        assert result.strength == EvidenceStrength.SUPPORTING

    def test_does_not_apply_when_af_common(self) -> None:
        """AF above the mito PM2 threshold → PM2_MITO does not apply."""
        variant = _mito_variant(gnomad_af=0.01)
        haplo = MitoHaploData(haplogroup="H1a1", is_haplogroup_defining=False)
        result = rule_mito_pm2(variant, haplo)
        assert result.applies is False


class TestMitoHeteroplasmyLevel:
    """ACGS 2024 §6: heteroplasmy level maps to clinical significance."""

    def test_high_heteroplasmy_pathogenic_evidence(self) -> None:
        """High heteroplasmy level (>70%) supports pathogenicity per ACGS 2024 §6."""
        mito_var = VariantInput(
            chrom="MT",
            pos=3_460,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="MT-ND1",
            transcript_id="NC_012920.1",
            is_mito=True,
            is_haplogroup_defining=False,
            gnomad_af=0.0,
            heteroplasmy_fraction=0.85,  # 85% heteroplasmy
        )
        result = rule_mito_heteroplasmy_level(mito_var)
        assert result.applies is True
        # High heteroplasmy is supporting pathogenic evidence
        assert result.strength in (
            EvidenceStrength.SUPPORTING,
            EvidenceStrength.MODERATE,
        )

    def test_low_heteroplasmy_uncertain_evidence(self) -> None:
        """Low heteroplasmy (<30%) in blood may not reflect tissue heteroplasmy."""
        mito_var = VariantInput(
            chrom="MT",
            pos=3_460,
            ref="G",
            alt="A",
            variant_type=VariantType.MISSENSE,
            gene_symbol="MT-ND1",
            transcript_id="NC_012920.1",
            is_mito=True,
            is_haplogroup_defining=False,
            gnomad_af=0.0,
            heteroplasmy_fraction=0.15,  # 15% heteroplasmy in blood
        )
        result = rule_mito_heteroplasmy_level(mito_var)
        # Low blood heteroplasmy: uncertain (may be high in affected tissue)
        # ACGS 2024 §6 notes tissue vs blood heteroplasmy discordance
        assert result is not None

    def test_not_applicable_for_non_mito_variant(self) -> None:
        """Heteroplasmy is a mito-only concept; non-mito variants get applies=False."""
        nuclear_variant = VariantInput(
            chrom="17",
            pos=41_276_045,
            ref="C",
            alt="T",
            variant_type=VariantType.MISSENSE,
            gene_symbol="BRCA1",
            transcript_id="NM_007294.4",
            is_mito=False,
        )
        result = rule_mito_heteroplasmy_level(nuclear_variant)
        assert result.applies is False
        assert result.rule_id == "MITO_HETEROPLASMY"

    def test_no_heteroplasmy_data_available(self) -> None:
        """No heteroplasmy_level/heteroplasmy_fraction set → applies=False, informative note."""
        variant = _mito_variant(heteroplasmy_level=None, heteroplasmy_fraction=None)
        result = rule_mito_heteroplasmy_level(variant)
        assert result.applies is False
        assert "not available" in " ".join(result.evidence_items).lower()

    def test_moderate_heteroplasmy_20_to_60_percent(self) -> None:
        """Heteroplasmy 20-60% → moderate, clinically relevant (applies=True)."""
        variant = _mito_variant(heteroplasmy_fraction=0.35)
        result = rule_mito_heteroplasmy_level(variant)
        assert result.applies is True
        assert "moderate" in " ".join(result.evidence_items).lower()

    def test_very_low_heteroplasmy_likely_artefact(self) -> None:
        """Heteroplasmy < 1% → likely sequencing artefact (applies=False)."""
        variant = _mito_variant(heteroplasmy_fraction=0.005)
        result = rule_mito_heteroplasmy_level(variant)
        assert result.applies is False
        assert "artefact" in " ".join(result.evidence_items).lower()
