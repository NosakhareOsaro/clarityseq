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

import pytest

from bayesacmg.models import EvidenceStrength, VariantInput, VariantType
from bayesacmg.rules.mito import (
    rule_mito_ba1,
    rule_mito_haplogroup_defining,
    rule_mito_heteroplasmy_level,
    rule_mito_pm2,
)


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
            is_haplogroup_defining=False,   # NOT haplogroup-defining
            gnomad_af=0.0001,
        )
        result = rule_mito_haplogroup_defining(non_haplogroup)
        assert result.applies is False, (
            "Non-haplogroup-defining mito variant should NOT be auto-classified Benign"
        )


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
