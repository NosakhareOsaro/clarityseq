"""
prioritisation.inheritance_filter
===================================
Filter candidate variants by inheritance mode.

Implements filters for:
    - AD  (Autosomal Dominant): heterozygous variants in dominant genes.
    - AR  (Autosomal Recessive): homozygous or compound heterozygous.
    - XL  (X-Linked): hemizygous in males; het/hom in females.
    - de_novo: variants not present in either parent.

Genotype codes (GT field from VCF):
    0/0  — homozygous reference
    0/1  — heterozygous
    1/1  — homozygous alt
    0/2, 1/2 — multiallelic
    ./.  — missing

Compound heterozygous (AR):
    Two different heterozygous variants in the same gene, one on each
    haplotype (in trans), together causing the recessive phenotype.
    Must be in different alleles — two variants on the same haplotype
    do NOT satisfy AR compound het.

References:
    Richards et al. 2015 PMID:25741868 (ACMG/AMP inheritance criteria).
    ACGS 2024 v1.2 §5 Table 2 (PM3, PS2, PM6 criteria).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Genotype constants
# ---------------------------------------------------------------------------

_HET_PATTERNS = frozenset({"0/1", "0|1", "1|0", "0/2", "1/2", "0|2", "1|2"})
_HOM_ALT_PATTERNS = frozenset({"1/1", "1|1", "2/2", "2|2"})
_HOM_REF_PATTERNS = frozenset({"0/0", "0|0"})
_HEMIZYGOUS_PATTERNS = frozenset({"1", "1/."})  # X-linked hemizygous in males


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class VariantRecord:
    """A candidate variant with genotype and gene information.

    Attributes:
        chrom: Chromosome (GRCh38).
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
        gene_symbol: HGNC gene symbol.
        proband_gt: Proband genotype string (e.g. ``"0/1"``).
        mother_gt: Mother genotype (None if not available).
        father_gt: Father genotype (None if not available).
        gene_disease_mode: Gene-level inheritance mode from OMIM/Orphanet
            (``"AD"``, ``"AR"``, ``"XL"``, ``"Mito"``).
        acmg_class: ACMG/AMP classification string.
        extra: Additional variant annotations.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: str
    proband_gt: str
    mother_gt: str | None = None
    father_gt: str | None = None
    gene_disease_mode: str = ""
    acmg_class: str = "VUS"
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class FilterResult:
    """Result of applying an inheritance mode filter.

    Attributes:
        variant: The variant that was evaluated.
        passes: True if the variant passes the filter.
        mode_applied: Inheritance mode filter that was applied.
        reason: Human-readable reason for pass/fail.
    """

    variant: VariantRecord
    passes: bool
    mode_applied: str
    reason: str


# ---------------------------------------------------------------------------
# Genotype helpers
# ---------------------------------------------------------------------------


def _is_het(gt: str) -> bool:
    """Return True if genotype string is heterozygous.

    Args:
        gt: VCF GT field string.

    Returns:
        True if the genotype represents a heterozygous call.
    """
    return gt in _HET_PATTERNS


def _is_hom_alt(gt: str) -> bool:
    """Return True if genotype string is homozygous alternate.

    Args:
        gt: VCF GT field string.

    Returns:
        True if the genotype represents a homozygous alternate call.
    """
    return gt in _HOM_ALT_PATTERNS


def _is_hom_ref(gt: str) -> bool:
    """Return True if genotype string is homozygous reference.

    Args:
        gt: VCF GT field string.

    Returns:
        True if the genotype represents a homozygous reference (wildtype) call.
    """
    return gt in _HOM_REF_PATTERNS


def _is_hemizygous(gt: str) -> bool:
    """Return True if genotype string is hemizygous (X-linked in males).

    Args:
        gt: VCF GT field string.

    Returns:
        True if the genotype represents a hemizygous call.
    """
    return gt in _HEMIZYGOUS_PATTERNS


# ---------------------------------------------------------------------------
# Inheritance mode filters
# ---------------------------------------------------------------------------


def filter_ad(variant: VariantRecord) -> FilterResult:
    """Apply Autosomal Dominant (AD) inheritance filter.

    AD requires the variant to be heterozygous in the proband.
    For confirmed de novo, PS2 (Strong) applies.
    For inherited AD variants, standard dominant disease segregation is expected.

    Args:
        variant: Candidate variant with proband genotype.

    Returns:
        FilterResult with passes=True if the variant is heterozygous.

    References:
        Richards et al. 2015 PMID:25741868 — AD filter criteria.
        ACGS 2024 v1.2 §5 Table 2.
    """
    if not _is_het(variant.proband_gt):
        return FilterResult(
            variant=variant,
            passes=False,
            mode_applied="AD",
            reason=(
                f"Proband GT={variant.proband_gt!r} is not heterozygous; "
                "AD filter requires het genotype."
            ),
        )
    return FilterResult(
        variant=variant,
        passes=True,
        mode_applied="AD",
        reason=f"Proband GT={variant.proband_gt!r} is heterozygous — passes AD filter.",
    )


def filter_ar(
    variant: VariantRecord,
    partner_variant: VariantRecord | None = None,
) -> FilterResult:
    """Apply Autosomal Recessive (AR) inheritance filter.

    AR requires either:
    1. Homozygous alternate in the proband, OR
    2. Compound heterozygous: two different heterozygous variants in the same
       gene, each inherited from one parent (PM3 criterion applies).

    Args:
        variant: Primary candidate variant.
        partner_variant: Second variant in the same gene for compound
            heterozygous evaluation.  None for homozygous evaluation only.

    Returns:
        FilterResult with passes=True if variant satisfies AR criteria.

    References:
        Richards et al. 2015 PMID:25741868 — AR filter criteria, PM3.
        ACGS 2024 v1.2 §5 Table 2.
    """
    # Case 1: Homozygous alternate
    if _is_hom_alt(variant.proband_gt):
        return FilterResult(
            variant=variant,
            passes=True,
            mode_applied="AR",
            reason=(
                f"Proband GT={variant.proband_gt!r} is homozygous alternate — "
                "satisfies AR (homozygous)."
            ),
        )

    # Case 2: Compound heterozygous
    if _is_het(variant.proband_gt) and partner_variant is not None:
        if partner_variant.gene_symbol != variant.gene_symbol:
            return FilterResult(
                variant=variant,
                passes=False,
                mode_applied="AR",
                reason=(
                    f"Partner variant is in different gene "
                    f"({partner_variant.gene_symbol} vs {variant.gene_symbol}); "
                    "compound het requires same gene."
                ),
            )
        if _is_het(partner_variant.proband_gt):
            return FilterResult(
                variant=variant,
                passes=True,
                mode_applied="AR",
                reason=(
                    f"Compound heterozygous: "
                    f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt} + "
                    f"{partner_variant.chrom}:{partner_variant.pos}:"
                    f"{partner_variant.ref}>{partner_variant.alt} "
                    f"(both het in {variant.gene_symbol}) — passes AR compound het filter."
                ),
            )

    return FilterResult(
        variant=variant,
        passes=False,
        mode_applied="AR",
        reason=(
            f"Proband GT={variant.proband_gt!r} does not satisfy AR filter. "
            "Requires hom-alt or compound het with a partner variant."
        ),
    )


def filter_xl(variant: VariantRecord, proband_sex: str = "UNKNOWN_SEX") -> FilterResult:
    """Apply X-Linked (XL) inheritance filter.

    X-linked dominant:
        Female probands: heterozygous on X chromosome.
        Male probands: hemizygous (chrX has only one copy in males).

    X-linked recessive:
        Female probands: homozygous alternate (rare carrier manifest).
        Male probands: hemizygous (disease-causing).

    Args:
        variant: Candidate variant; ``chrom`` must be ``"chrX"`` or ``"X"``.
        proband_sex: Proband biological sex (``"FEMALE"``, ``"MALE"``,
            ``"UNKNOWN_SEX"``).

    Returns:
        FilterResult with passes=True if variant satisfies XL filter.

    References:
        Richards et al. 2015 PMID:25741868 — XL filter criteria.
        ACGS 2024 v1.2 §5 Table 2.
    """
    chrom = variant.chrom.upper().lstrip("CHR")
    if chrom not in ("X", "CHRX"):
        chrom = variant.chrom
        is_x = "X" in chrom.upper()
    else:
        is_x = True

    if not is_x:
        return FilterResult(
            variant=variant,
            passes=False,
            mode_applied="XL",
            reason=(
                f"Variant on chromosome {variant.chrom!r} — "
                "XL filter requires chrX variants."
            ),
        )

    gt = variant.proband_gt

    if proband_sex == "MALE":
        if _is_het(gt) or _is_hom_alt(gt) or _is_hemizygous(gt):
            return FilterResult(
                variant=variant,
                passes=True,
                mode_applied="XL",
                reason=(
                    f"Male proband, chrX variant, GT={gt!r} — "
                    "hemizygous/het X-linked variant passes XL filter."
                ),
            )
    elif proband_sex == "FEMALE":
        if _is_het(gt) or _is_hom_alt(gt):
            return FilterResult(
                variant=variant,
                passes=True,
                mode_applied="XL",
                reason=(
                    f"Female proband, chrX variant, GT={gt!r} — "
                    "het/hom X-linked variant passes XL filter."
                ),
            )
    else:
        # Unknown sex — accept any non-ref genotype
        if not _is_hom_ref(gt):
            return FilterResult(
                variant=variant,
                passes=True,
                mode_applied="XL",
                reason=(
                    f"Unknown sex, chrX variant, GT={gt!r} — "
                    "non-ref X-linked variant passes XL filter."
                ),
            )

    return FilterResult(
        variant=variant,
        passes=False,
        mode_applied="XL",
        reason=(
            f"chrX variant GT={gt!r} (sex={proband_sex!r}) "
            "does not satisfy XL filter."
        ),
    )


def filter_de_novo(
    variant: VariantRecord,
    confirmed: bool = False,
) -> FilterResult:
    """Apply de novo inheritance filter.

    Confirmed de novo (PS2): variant is absent in both parents (parental
    testing performed).
    Assumed de novo (PM6): variant suspected de novo based on clinical
    grounds without parental testing.

    Args:
        variant: Candidate variant with mother_gt and father_gt fields.
        confirmed: If True, the de novo status was confirmed by parental
            testing (PS2 criterion).  If False, checks available parental
            genotypes.

    Returns:
        FilterResult with passes=True if de novo criteria are met.
        ``reason`` indicates whether PS2 (confirmed) or PM6 (assumed) applies.

    References:
        Richards et al. 2015 PMID:25741868 — PS2, PM6.
        ACGS 2024 v1.2 §5 Table 2.
    """
    gt = variant.proband_gt
    if _is_hom_ref(gt):
        return FilterResult(
            variant=variant,
            passes=False,
            mode_applied="de_novo",
            reason=f"Proband GT={gt!r} is hom-ref; cannot be de novo.",
        )

    if confirmed:
        return FilterResult(
            variant=variant,
            passes=True,
            mode_applied="de_novo",
            reason=(
                "De novo status confirmed by parental testing — "
                "PS2 (Strong) criterion applies. "
                "Richards 2015 PMID:25741868."
            ),
        )

    # Check parental genotypes if available
    mother_gt = variant.mother_gt
    father_gt = variant.father_gt

    if mother_gt is None and father_gt is None:
        return FilterResult(
            variant=variant,
            passes=True,  # Cannot exclude de novo without parental data
            mode_applied="de_novo",
            reason=(
                "No parental genotypes available; cannot exclude de novo. "
                "If phenotype strongly suggests de novo, apply PM6 (assumed de novo)."
            ),
        )

    mother_ref = mother_gt is None or _is_hom_ref(mother_gt)
    father_ref = father_gt is None or _is_hom_ref(father_gt)

    if mother_ref and father_ref:
        return FilterResult(
            variant=variant,
            passes=True,
            mode_applied="de_novo",
            reason=(
                f"Variant absent in mother (GT={mother_gt!r}) and "
                f"father (GT={father_gt!r}) — de novo supported. "
                "Apply PM6 if parental testing was not performed, or PS2 if confirmed."
            ),
        )

    inherited_from = []
    if not mother_ref and mother_gt is not None:
        inherited_from.append(f"mother (GT={mother_gt!r})")
    if not father_ref and father_gt is not None:
        inherited_from.append(f"father (GT={father_gt!r})")

    return FilterResult(
        variant=variant,
        passes=False,
        mode_applied="de_novo",
        reason=(
            f"Variant present in {', '.join(inherited_from)} — not de novo."
        ),
    )


def apply_inheritance_filter(
    variants: list[VariantRecord],
    mode: str,
    proband_sex: str = "UNKNOWN_SEX",
    compound_het_gene_map: dict[str, list[VariantRecord]] | None = None,
) -> list[FilterResult]:
    """Apply an inheritance mode filter to a list of candidate variants.

    Args:
        variants: List of VariantRecord objects to filter.
        mode: Inheritance mode to apply: ``"AD"``, ``"AR"``, ``"XL"``,
            ``"de_novo"``, or ``"any"`` (no filter — return all).
        proband_sex: Proband biological sex (used for XL filtering).
        compound_het_gene_map: Dict of gene_symbol → list of het variants
            in that gene.  Required for AR compound-het evaluation.
            If None, AR filter only applies homozygous check.

    Returns:
        List of FilterResult objects for each variant.

    References:
        Richards et al. 2015 PMID:25741868.
        ACGS 2024 v1.2 §5 Table 2.
    """
    results: list[FilterResult] = []
    mode_upper = mode.upper()

    for variant in variants:
        if mode_upper == "ANY":
            results.append(
                FilterResult(
                    variant=variant,
                    passes=True,
                    mode_applied="any",
                    reason="No inheritance filter applied (mode=any).",
                )
            )
        elif mode_upper == "AD":
            results.append(filter_ad(variant))
        elif mode_upper == "AR":
            # Find partner variant for compound het check
            partner: VariantRecord | None = None
            if compound_het_gene_map and variant.gene_symbol in compound_het_gene_map:
                gene_variants = compound_het_gene_map[variant.gene_symbol]
                others = [v for v in gene_variants if v is not variant]
                partner = others[0] if others else None
            results.append(filter_ar(variant, partner_variant=partner))
        elif mode_upper in ("XL", "XLD", "XLR"):
            results.append(filter_xl(variant, proband_sex=proband_sex))
        elif mode_upper == "DE_NOVO":
            results.append(filter_de_novo(variant))
        else:
            logger.warning("Unknown inheritance mode '%s'; skipping filter.", mode)
            results.append(
                FilterResult(
                    variant=variant,
                    passes=True,
                    mode_applied=mode,
                    reason=f"Unknown mode '{mode}'; no filter applied.",
                )
            )

    return results
