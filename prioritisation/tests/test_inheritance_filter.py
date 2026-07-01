"""
prioritisation.tests.test_inheritance_filter
=============================================
pytest tests for AD/AR/XL/de_novo inheritance mode filters.

Tests cover:
    - filter_ad: het passes, hom-ref fails, hom-alt fails.
    - filter_ar: hom-alt passes, compound het passes, het alone fails.
    - filter_xl: male hemizygous passes, female het passes, autosomal fails.
    - filter_de_novo: absent in parents passes, inherited fails.
    - apply_inheritance_filter: batch filter application.

References:
    Richards et al. 2015 PMID:25741868 (ACMG/AMP criteria PS2, PM3, PM6).
    ACGS 2024 v1.2 §5 Table 2.
"""

from __future__ import annotations

from typing import Any

import pytest

from prioritisation.inheritance_filter import (
    FilterResult,
    VariantRecord,
    apply_inheritance_filter,
    filter_ad,
    filter_ar,
    filter_de_novo,
    filter_xl,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def make_variant(
    chrom: str = "chr17",
    pos: int = 43044295,
    ref: str = "G",
    alt: str = "A",
    gene: str = "BRCA1",
    proband_gt: str = "0/1",
    mother_gt: str | None = None,
    father_gt: str | None = None,
    acmg_class: str = "VUS",
) -> VariantRecord:
    """Create a VariantRecord for testing.

    Args:
        chrom: Chromosome.
        pos: Position.
        ref: Reference allele.
        alt: Alternate allele.
        gene: Gene symbol.
        proband_gt: Proband genotype.
        mother_gt: Mother genotype.
        father_gt: Father genotype.
        acmg_class: ACMG classification.

    Returns:
        VariantRecord for testing.
    """
    return VariantRecord(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene_symbol=gene,
        proband_gt=proband_gt,
        mother_gt=mother_gt,
        father_gt=father_gt,
        acmg_class=acmg_class,
    )


# ---------------------------------------------------------------------------
# AD filter tests
# ---------------------------------------------------------------------------


class TestFilterAD:
    """Tests for Autosomal Dominant (AD) inheritance filter."""

    def test_het_passes_ad_filter(self) -> None:
        """Heterozygous variant passes AD filter.

        AD requires heterozygous in the proband.
        Richards 2015 PMID:25741868.
        """
        v = make_variant(proband_gt="0/1")
        result = filter_ad(v)
        assert result.passes is True, (
            "Het variant should pass AD filter"
        )
        assert result.mode_applied == "AD"

    def test_hom_ref_fails_ad_filter(self) -> None:
        """Homozygous reference fails AD filter."""
        v = make_variant(proband_gt="0/0")
        result = filter_ad(v)
        assert result.passes is False, (
            "Hom-ref (0/0) should fail AD filter"
        )

    def test_hom_alt_fails_ad_filter(self) -> None:
        """Homozygous alternate fails AD filter.

        1/1 is not heterozygous; AD requires het.
        """
        v = make_variant(proband_gt="1/1")
        result = filter_ad(v)
        assert result.passes is False, (
            "Hom-alt (1/1) should fail AD filter (not het)"
        )

    def test_phased_het_passes_ad_filter(self) -> None:
        """Phased heterozygous (0|1) passes AD filter."""
        v = make_variant(proband_gt="0|1")
        result = filter_ad(v)
        assert result.passes is True, (
            "Phased het (0|1) should pass AD filter"
        )

    def test_phased_het_reverse_passes_ad(self) -> None:
        """Phased het (1|0) passes AD filter."""
        v = make_variant(proband_gt="1|0")
        result = filter_ad(v)
        assert result.passes is True


# ---------------------------------------------------------------------------
# AR filter tests
# ---------------------------------------------------------------------------


class TestFilterAR:
    """Tests for Autosomal Recessive (AR) inheritance filter."""

    def test_hom_alt_passes_ar_filter(self) -> None:
        """Homozygous alternate passes AR filter.

        1/1 satisfies AR (homozygous recessive).
        Richards 2015 PMID:25741868 PM3.
        """
        v = make_variant(proband_gt="1/1")
        result = filter_ar(v)
        assert result.passes is True, (
            "Hom-alt (1/1) should pass AR filter"
        )

    def test_het_without_partner_fails_ar_filter(self) -> None:
        """Heterozygous without a partner variant fails AR filter."""
        v = make_variant(proband_gt="0/1")
        result = filter_ar(v, partner_variant=None)
        assert result.passes is False, (
            "Single het variant without partner should fail AR filter"
        )

    def test_compound_het_passes_ar_filter(self) -> None:
        """Two het variants in the same gene pass AR (compound het).

        Compound heterozygous: two different het variants in the same gene.
        PM3 (Moderate) applies when one is in trans with known pathogenic.
        Richards 2015 PMID:25741868.
        """
        v1 = make_variant(
            chrom="chr7", pos=117540000, ref="C", alt="T",
            gene="CFTR", proband_gt="0/1",
        )
        v2 = make_variant(
            chrom="chr7", pos=117548000, ref="G", alt="A",
            gene="CFTR", proband_gt="0/1",
        )
        result = filter_ar(v1, partner_variant=v2)
        assert result.passes is True, (
            "Compound het (two het variants in same gene) should pass AR"
        )

    def test_compound_het_different_gene_fails(self) -> None:
        """Compound het in different genes fails AR filter."""
        v1 = make_variant(gene="CFTR", proband_gt="0/1")
        v2 = make_variant(gene="BRCA1", proband_gt="0/1")
        result = filter_ar(v1, partner_variant=v2)
        assert result.passes is False, (
            "Compound het with partner in different gene should fail AR"
        )

    def test_hom_ref_fails_ar_filter(self) -> None:
        """Homozygous reference fails AR filter."""
        v = make_variant(proband_gt="0/0")
        result = filter_ar(v)
        assert result.passes is False

    def test_phased_hom_alt_passes_ar(self) -> None:
        """Phased homozygous alt (1|1) passes AR filter."""
        v = make_variant(proband_gt="1|1")
        result = filter_ar(v)
        assert result.passes is True


# ---------------------------------------------------------------------------
# XL filter tests
# ---------------------------------------------------------------------------


class TestFilterXL:
    """Tests for X-Linked (XL) inheritance filter."""

    def test_male_chrX_het_passes_xl(self) -> None:
        """Male proband with chrX het variant passes XL filter.

        In males, chrX variants appear hemizygous; het GT is expected.
        Richards 2015 PMID:25741868.
        """
        v = make_variant(chrom="chrX", pos=150000000, gene="MECP2", proband_gt="0/1")
        result = filter_xl(v, proband_sex="MALE")
        assert result.passes is True, (
            "Male chrX het should pass XL filter"
        )

    def test_female_chrX_het_passes_xl(self) -> None:
        """Female proband with chrX het variant passes XL filter."""
        v = make_variant(chrom="chrX", pos=150000000, gene="MECP2", proband_gt="0/1")
        result = filter_xl(v, proband_sex="FEMALE")
        assert result.passes is True, (
            "Female chrX het should pass XL filter (X-linked dominant)"
        )

    def test_female_chrX_hom_alt_passes_xl(self) -> None:
        """Female proband with chrX hom-alt passes XL filter.

        Rare carrier manifesting female: hom-alt on chrX.
        """
        v = make_variant(chrom="chrX", pos=150000000, gene="MECP2", proband_gt="1/1")
        result = filter_xl(v, proband_sex="FEMALE")
        assert result.passes is True

    def test_autosomal_fails_xl_filter(self) -> None:
        """Variant on chr17 fails XL filter (not chrX)."""
        v = make_variant(chrom="chr17", pos=43044295, gene="BRCA1", proband_gt="0/1")
        result = filter_xl(v, proband_sex="FEMALE")
        assert result.passes is False, (
            "Autosomal variant should fail XL filter"
        )

    def test_chrx_hom_ref_male_fails_xl(self) -> None:
        """Male proband with chrX hom-ref fails XL filter."""
        v = make_variant(chrom="chrX", pos=150000000, gene="MECP2", proband_gt="0/0")
        result = filter_xl(v, proband_sex="MALE")
        assert result.passes is False

    def test_xl_unknown_sex_accepts_nonref(self) -> None:
        """Unknown sex proband with chrX non-ref variant passes XL filter."""
        v = make_variant(chrom="chrX", pos=150000000, gene="MECP2", proband_gt="0/1")
        result = filter_xl(v, proband_sex="UNKNOWN_SEX")
        assert result.passes is True


# ---------------------------------------------------------------------------
# de_novo filter tests
# ---------------------------------------------------------------------------


class TestFilterDeNovo:
    """Tests for de novo inheritance filter."""

    def test_absent_both_parents_passes_de_novo(self) -> None:
        """Variant absent in both parents passes de novo filter.

        Richards 2015 PMID:25741868 PS2 (confirmed) / PM6 (assumed).
        """
        v = make_variant(
            proband_gt="0/1",
            mother_gt="0/0",
            father_gt="0/0",
        )
        result = filter_de_novo(v)
        assert result.passes is True, (
            "Variant absent in both parents should pass de_novo filter"
        )

    def test_present_in_mother_fails_de_novo(self) -> None:
        """Variant present in mother fails de novo filter."""
        v = make_variant(
            proband_gt="0/1",
            mother_gt="0/1",  # present in mother
            father_gt="0/0",
        )
        result = filter_de_novo(v)
        assert result.passes is False, (
            "Variant present in mother should fail de_novo filter"
        )

    def test_present_in_father_fails_de_novo(self) -> None:
        """Variant present in father fails de novo filter."""
        v = make_variant(
            proband_gt="0/1",
            mother_gt="0/0",
            father_gt="0/1",  # present in father
        )
        result = filter_de_novo(v)
        assert result.passes is False, (
            "Variant present in father should fail de_novo filter"
        )

    def test_no_parental_data_passes_with_caveat(self) -> None:
        """No parental genotype data passes de novo filter with caveat.

        Without parental data, de novo cannot be excluded.
        Apply PM6 (assumed de novo) in this case.
        """
        v = make_variant(
            proband_gt="0/1",
            mother_gt=None,
            father_gt=None,
        )
        result = filter_de_novo(v)
        assert result.passes is True, (
            "No parental data: de novo cannot be excluded; should pass with PM6 caveat"
        )

    def test_confirmed_de_novo_passes(self) -> None:
        """Confirmed de novo (PS2) passes de novo filter."""
        v = make_variant(proband_gt="0/1", mother_gt="0/0", father_gt="0/0")
        result = filter_de_novo(v, confirmed=True)
        assert result.passes is True
        assert "PS2" in result.reason, (
            "Confirmed de novo should mention PS2 criterion"
        )

    def test_hom_ref_proband_fails_de_novo(self) -> None:
        """Proband with hom-ref genotype fails de novo filter."""
        v = make_variant(proband_gt="0/0")
        result = filter_de_novo(v)
        assert result.passes is False


# ---------------------------------------------------------------------------
# apply_inheritance_filter tests
# ---------------------------------------------------------------------------


class TestApplyInheritanceFilter:
    """Tests for batch inheritance mode filtering."""

    def test_apply_ad_filter_batch(self) -> None:
        """Batch AD filter correctly classifies het and hom variants."""
        variants = [
            make_variant(proband_gt="0/1", gene="GENE_A"),
            make_variant(proband_gt="1/1", gene="GENE_B"),
            make_variant(proband_gt="0/0", gene="GENE_C"),
        ]
        results = apply_inheritance_filter(variants, mode="AD")
        assert len(results) == 3
        assert results[0].passes is True   # het → passes
        assert results[1].passes is False  # hom-alt → fails
        assert results[2].passes is False  # hom-ref → fails

    def test_apply_any_mode_passes_all(self) -> None:
        """Mode 'any' passes all variants regardless of genotype."""
        variants = [
            make_variant(proband_gt="0/0"),
            make_variant(proband_gt="1/1"),
            make_variant(proband_gt="0/1"),
        ]
        results = apply_inheritance_filter(variants, mode="any")
        assert all(r.passes for r in results), (
            "Mode 'any' should pass all variants"
        )

    def test_apply_xl_filter_female(self) -> None:
        """XL batch filter correctly handles female sex on chrX."""
        variants = [
            make_variant(chrom="chrX", proband_gt="0/1"),
            make_variant(chrom="chr1", proband_gt="0/1"),  # not chrX
        ]
        results = apply_inheritance_filter(variants, mode="XL", proband_sex="FEMALE")
        assert results[0].passes is True   # chrX het → passes
        assert results[1].passes is False  # chr1 → fails XL

    def test_apply_de_novo_batch(self) -> None:
        """Batch de novo filter correctly identifies de novo variants."""
        variants = [
            make_variant(
                proband_gt="0/1",
                mother_gt="0/0",
                father_gt="0/0",
                gene="GENE_A",
            ),
            make_variant(
                proband_gt="0/1",
                mother_gt="0/1",  # inherited from mother
                father_gt="0/0",
                gene="GENE_B",
            ),
        ]
        results = apply_inheritance_filter(variants, mode="de_novo")
        assert results[0].passes is True   # absent in parents → de novo
        assert results[1].passes is False  # present in mother → not de novo

    def test_apply_ar_filter_batch_hom_alt_no_compound_het_map(self) -> None:
        """Batch AR filter (no compound_het_gene_map) evaluates hom-alt only."""
        variants = [
            make_variant(proband_gt="1/1", gene="GENE_A"),  # hom-alt → passes
            make_variant(proband_gt="0/1", gene="GENE_B"),  # het alone → fails
            make_variant(proband_gt="0/0", gene="GENE_C"),  # hom-ref → fails
        ]
        results = apply_inheritance_filter(variants, mode="AR")
        assert len(results) == 3
        assert results[0].passes is True
        assert results[0].mode_applied == "AR"
        assert results[1].passes is False
        assert results[2].passes is False

    def test_apply_ar_filter_with_compound_het_gene_map(self) -> None:
        """Batch AR filter finds a partner variant via compound_het_gene_map
        and evaluates compound heterozygosity.
        """
        v1 = make_variant(
            chrom="chr7", pos=117540000, ref="C", alt="T",
            gene="CFTR", proband_gt="0/1",
        )
        v2 = make_variant(
            chrom="chr7", pos=117548000, ref="G", alt="A",
            gene="CFTR", proband_gt="0/1",
        )
        compound_het_gene_map = {"CFTR": [v1, v2]}

        results = apply_inheritance_filter(
            [v1, v2], mode="AR", compound_het_gene_map=compound_het_gene_map
        )
        assert len(results) == 2
        assert results[0].passes is True, "v1 should find v2 as partner and pass AR"
        assert "Compound heterozygous" in results[0].reason
        assert results[1].passes is True, "v2 should find v1 as partner and pass AR"

    def test_apply_ar_filter_compound_het_map_no_other_variant(self) -> None:
        """When the gene maps only to the variant itself, no partner is found
        and the AR filter falls through to the failing branch.
        """
        v1 = make_variant(gene="SOLO_GENE", proband_gt="0/1")
        compound_het_gene_map = {"SOLO_GENE": [v1]}  # no distinct partner

        results = apply_inheritance_filter(
            [v1], mode="AR", compound_het_gene_map=compound_het_gene_map
        )
        assert len(results) == 1
        assert results[0].passes is False, (
            "No distinct partner variant available; AR filter should fail"
        )

    def test_apply_unknown_mode_passes_with_warning(self, caplog) -> None:
        """An unrecognised inheritance mode logs a warning and passes the
        variant through unfiltered.
        """
        import logging

        variants = [make_variant(proband_gt="0/1", gene="GENE_X")]
        with caplog.at_level(logging.WARNING, logger="prioritisation.inheritance_filter"):
            results = apply_inheritance_filter(variants, mode="MITOCHONDRIAL_WEIRD")

        assert len(results) == 1
        assert results[0].passes is True
        assert results[0].mode_applied == "MITOCHONDRIAL_WEIRD"
        assert "Unknown mode" in results[0].reason
        assert any("Unknown inheritance mode" in rec.message for rec in caplog.records)
