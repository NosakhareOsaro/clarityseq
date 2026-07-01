"""
bayesacmg.rules.benign
======================

All 12 benign ACMG/AMP criteria: BA1, BS1-4, BP1-7.

Every rule function:
  - cites Richards et al. 2015 PMID:25741868 for the original criterion
  - cites the ClinGen SVI or ACGS 2024 v1.2 update where applicable
  - has full Google-style docstrings
  - has complete type annotations

References:
    Richards et al. 2015 PMID:25741868
    Cheng et al. 2023 PMID:37703350 (AlphaMissense thresholds)
    ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al.)
    ClinGen SVI Working Group 2024
"""

from __future__ import annotations

from bayesacmg.models import (
    ACMGRule,
    EvidenceStrength,
    GeneData,
    TranscriptData,
    VariantInput,
    VariantType,
)

# ---------------------------------------------------------------------------
# Thresholds — each with inline citation
# ---------------------------------------------------------------------------

_BA1_AF_THRESHOLD = 0.05  # >5 % in any gnomAD population → Benign stand-alone
# Richards 2015 PMID:25741868; NOT for mito — see mito.py
_BS1_AF_THRESHOLD = 0.01  # >1 % → BS1; ACGS 2024 v1.2 §5 Table 2
_AM_BP4_THRESHOLD = 0.340  # AlphaMissense ≤0.340 → BP4; Cheng 2023 PMID:37703350
_AM_PP3_THRESHOLD = 0.564  # AlphaMissense ≥0.564 → PP3; Cheng 2023 PMID:37703350


# ---------------------------------------------------------------------------
# BA1
# ---------------------------------------------------------------------------


def rule_ba1(variant: VariantInput) -> ACMGRule:
    """BA1 — Allele frequency >5% in gnomAD v4.1 → Stand-alone Benign.

    NOTE: This threshold does NOT apply to mitochondrial variants.
    For mito variants, use rules/mito.py rule_mito_ba1() which implements
    the ACGS 2024 §6 mito-specific BA1 thresholds.

    Args:
        variant: Annotated variant; ``gnomad_popmax_af`` or ``gnomad_af``
            is evaluated.  If ``is_mito`` is True, this rule returns False
            and the caller should use rule_mito_ba1() instead.

    Returns:
        ACMGRule with rule_id ``"BA1"`` and strength STAND_ALONE when
        AF >5 % in any gnomAD v4.1 population.  STAND_ALONE triggers
        direct Benign classification without point scoring.

    References:
        Richards et al. 2015 PMID:25741868 — criterion BA1.
        ACGS 2024 v1.2 §5 Table 2.
        ACGS 2024 v1.2 §6 — mito caveat.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
        "gnomAD v4.1 (April 2024, 807,162 individuals)",
    ]
    if variant.is_mito:
        return ACMGRule(
            rule_id="BA1",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=[
                "Mito variant: standard BA1 threshold does not apply — use rule_mito_ba1()"
            ],
            citations=citations,
            applies=False,
            notes="ACGS 2024 §6: standard 5% AF threshold is not valid for mito variants",
        )

    # Use popmax AF if available; fall back to global AF
    af_to_check = variant.gnomad_popmax_af
    if af_to_check is None:
        af_to_check = variant.gnomad_af

    if (
        af_to_check is not None and af_to_check > _BA1_AF_THRESHOLD
    ):  # >5%; Richards 2015
        return ACMGRule(
            rule_id="BA1",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=[
                f"gnomAD v4.1 popmax AF={af_to_check:.4f} > {_BA1_AF_THRESHOLD} "
                "(>5% — stand-alone Benign; Richards 2015 PMID:25741868)"
            ],
            citations=citations,
            applies=True,
        )

    return ACMGRule(
        rule_id="BA1",
        strength=EvidenceStrength.STAND_ALONE,
        evidence_items=[
            f"AF={af_to_check} does not exceed BA1 threshold ({_BA1_AF_THRESHOLD})"
        ],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BS1
# ---------------------------------------------------------------------------


def rule_bs1(variant: VariantInput) -> ACMGRule:
    """BS1 — Allele frequency greater than expected for the disorder.

    Args:
        variant: Annotated variant; ``gnomad_af`` is the primary field.

    Returns:
        ACMGRule with rule_id ``"BS1"`` and strength STRONG_BENIGN (-4 pts).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BS1.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    af = variant.gnomad_popmax_af or variant.gnomad_af
    if af is not None and af > _BS1_AF_THRESHOLD:  # >1%; ACGS 2024 v1.2 §5 Table 2
        return ACMGRule(
            rule_id="BS1",
            strength=EvidenceStrength.STRONG_BENIGN,
            evidence_items=[
                f"gnomAD v4.1 AF={af:.4f} > {_BS1_AF_THRESHOLD} "
                "(greater than expected for disorder)"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BS1",
        strength=EvidenceStrength.STRONG_BENIGN,
        evidence_items=["AF does not exceed BS1 threshold"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BS2
# ---------------------------------------------------------------------------


def rule_bs2(
    variant: VariantInput,
    observed_in_healthy_adults: bool,
    n_healthy_homozygotes: int = 0,
) -> ACMGRule:
    """BS2 — Variant observed in healthy adult individuals for recessive/dominant conditions.

    For recessive disorders: observe homozygous in healthy adults.
    For dominant disorders: observe heterozygous in healthy adults.

    Args:
        variant: Annotated variant; ``gnomad_nhomalt`` provides homozygous count.
        observed_in_healthy_adults: True if the variant has been observed in
            healthy, unaffected adult individuals with full penetrance expected.
        n_healthy_homozygotes: Number of gnomAD homozygotes (for recessive check).

    Returns:
        ACMGRule with rule_id ``"BS2"`` and strength STRONG_BENIGN (-4 pts).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BS2.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    gnomad_homalt = variant.gnomad_nhomalt or 0
    effective_hom = max(n_healthy_homozygotes, gnomad_homalt)
    if (
        observed_in_healthy_adults or effective_hom >= 2
    ):  # ≥2 homozygotes → strong benign signal
        evidence = []
        if observed_in_healthy_adults:
            evidence.append(
                "Observed in healthy unaffected adults with full penetrance expected"
            )
        if effective_hom >= 2:
            evidence.append(
                f"gnomAD v4.1: {effective_hom} homozygous individuals (unaffected)"
            )
        return ACMGRule(
            rule_id="BS2",
            strength=EvidenceStrength.STRONG_BENIGN,
            evidence_items=evidence,
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BS2",
        strength=EvidenceStrength.STRONG_BENIGN,
        evidence_items=["Not observed in healthy adults with full penetrance"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BS3
# ---------------------------------------------------------------------------


def rule_bs3(variant: VariantInput) -> ACMGRule:
    """BS3 — Well-established functional studies show no damaging effect.

    Args:
        variant: Annotated variant; ``functional_study_result`` is checked.

    Returns:
        ACMGRule with rule_id ``"BS3"`` and strength STRONG_BENIGN (-4 pts).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BS3.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if variant.functional_study_result == "benign":
        return ACMGRule(
            rule_id="BS3",
            strength=EvidenceStrength.STRONG_BENIGN,
            evidence_items=[
                "Well-established functional study demonstrates no damaging effect"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BS3",
        strength=EvidenceStrength.STRONG_BENIGN,
        evidence_items=["No benign functional study result available"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BS4
# ---------------------------------------------------------------------------


def rule_bs4(
    variant: VariantInput,
    lack_of_segregation: bool,
) -> ACMGRule:
    """BS4 — Lack of segregation in affected members of a family.

    Args:
        variant: Annotated variant.
        lack_of_segregation: True if the variant does not segregate with
            disease in affected family members where expected.

    Returns:
        ACMGRule with rule_id ``"BS4"`` and strength STRONG_BENIGN (-4 pts).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BS4.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if lack_of_segregation:
        return ACMGRule(
            rule_id="BS4",
            strength=EvidenceStrength.STRONG_BENIGN,
            evidence_items=[
                "Variant does not segregate with disease in affected family members"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BS4",
        strength=EvidenceStrength.STRONG_BENIGN,
        evidence_items=[
            "Lack of segregation data or segregation supports pathogenicity"
        ],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP1
# ---------------------------------------------------------------------------


def rule_bp1(
    variant: VariantInput,
    transcript: TranscriptData,
    gene: GeneData,
) -> ACMGRule:
    """BP1 — Missense variant in gene where only truncating variants cause disease.

    Args:
        variant: Annotated variant.
        transcript: Transcript annotation; ``aa_change`` is checked.
        gene: Gene data; ``missense_only_gene`` is the inverse flag
            (BP1 applies when gene is NOT missense-only but truncating-only).

    Returns:
        ACMGRule with rule_id ``"BP1"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP1.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    is_missense = (
        variant.variant_type == "snv"
        and transcript.aa_change is not None
        and not transcript.aa_change.startswith("p.=")
    )
    # gene.missense_only_gene == True means only missense causes disease — BP1 does NOT apply
    # BP1 applies when the gene exclusively causes disease via truncating variants,
    # and this variant is a missense.
    truncating_only_gene = not gene.missense_only_gene and gene.lof_is_disease_mechanism

    if is_missense and truncating_only_gene:
        return ACMGRule(
            rule_id="BP1",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                f"Missense variant in {gene.gene_symbol} where "
                "only truncating variants cause disease"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP1",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=["BP1 criteria not met"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP2
# ---------------------------------------------------------------------------


def rule_bp2(
    variant: VariantInput,
    in_trans_with_pathogenic_dominant: bool,
) -> ACMGRule:
    """BP2 — Observed in trans with pathogenic variant in dominant disorder.

    Also applies if observed in cis with established pathogenic variant.

    Args:
        variant: Annotated variant.
        in_trans_with_pathogenic_dominant: True if the variant is in trans
            with a pathogenic variant in a dominant disorder.

    Returns:
        ACMGRule with rule_id ``"BP2"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP2.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if in_trans_with_pathogenic_dominant:
        return ACMGRule(
            rule_id="BP2",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                "Observed in trans with pathogenic variant in dominant disorder"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP2",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=[
            "Not observed in trans with pathogenic variant in dominant disorder"
        ],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP3
# ---------------------------------------------------------------------------


def rule_bp3(
    variant: VariantInput,
    transcript: TranscriptData,
    in_repeat_region: bool,
) -> ACMGRule:
    """BP3 — In-frame indel in repeat region without known function.

    Args:
        variant: Annotated variant; ``variant_type`` is checked.
        transcript: Transcript annotation.
        in_repeat_region: True if the variant falls within an annotated
            simple/tandem repeat region without known function.

    Returns:
        ACMGRule with rule_id ``"BP3"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP3.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    in_frame = variant.variant_type in {"inframe_insertion", "inframe_deletion"}
    if in_frame and in_repeat_region:
        return ACMGRule(
            rule_id="BP3",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                "In-frame indel within repeat region without known function"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP3",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=["Not an in-frame indel in a repeat region"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP4
# ---------------------------------------------------------------------------


def rule_bp4(
    variant: VariantInput,
    alphamissense_score: float | None = None,
    revel_score: float | None = None,
) -> ACMGRule:
    """BP4 — Multiple computational evidence suggest no impact on gene/product.

    PRIMARY PREDICTOR: AlphaMissense ≤ 0.340 → BP4.
    Reference: Cheng et al. 2023 Science PMID:37703350.

    This is the benign complement of PP3.  For splice site variants, use
    rules/splicing.py rule_splicing_pp3_bp4_bp7() instead.

    Args:
        variant: Annotated variant.
        alphamissense_score: AlphaMissense score (0–1).
            ≤0.340 → BP4 fires.  ≥0.564 → PP3 fires (see pathogenic.py).
        revel_score: REVEL score (0–1); secondary predictor.

    Returns:
        ACMGRule with rule_id ``"BP4"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP4.
        Cheng et al. 2023 PMID:37703350 — AlphaMissense thresholds.
        ClinGen SVI in silico recommendation (2024).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "Cheng et al. 2023 PMID:37703350 (AlphaMissense)",
        "ClinGen SVI in silico recommendation (2024)",
        "ACGS 2024 v1.2 §5 Table 2",
    ]

    splice_types = {"splice_canonical", "splice_region"}
    if variant.variant_type in splice_types:
        return ACMGRule(
            rule_id="BP4",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                "Splice variant: use splicing.rule_splicing_pp3_bp4_bp7() instead"
            ],
            citations=citations,
            applies=False,
            notes="BP4 for splice variants handled by rules/splicing.py",
        )

    # PRIMARY: AlphaMissense ≤ 0.340 (Cheng 2023 PMID:37703350)
    if alphamissense_score is not None:
        if alphamissense_score <= _AM_BP4_THRESHOLD:  # ≤0.340; Cheng 2023 PMID:37703350
            return ACMGRule(
                rule_id="BP4",
                strength=EvidenceStrength.SUPPORTING_BENIGN,
                evidence_items=[
                    f"AlphaMissense score {alphamissense_score:.3f} ≤ {_AM_BP4_THRESHOLD} "
                    "(likely benign threshold; Cheng 2023 PMID:37703350)"
                ],
                citations=citations,
                applies=True,
            )
        # Score above BP4 threshold — BP4 does not apply
        return ACMGRule(
            rule_id="BP4",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                f"AlphaMissense score {alphamissense_score:.3f} > {_AM_BP4_THRESHOLD} — "
                "BP4 does not apply"
            ],
            citations=citations,
            applies=False,
        )

    # SECONDARY: REVEL < 0.15 (benign range; ClinGen SVI 2024)
    if revel_score is not None and revel_score < 0.15:  # <0.15 benign; ClinGen SVI 2024
        return ACMGRule(
            rule_id="BP4",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                f"REVEL score {revel_score:.3f} < 0.15 (secondary predictor; "
                "AlphaMissense unavailable)"
            ],
            citations=citations,
            applies=True,
            notes="AlphaMissense unavailable; REVEL used as secondary predictor",
        )

    return ACMGRule(
        rule_id="BP4",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=["No in silico scores indicate benign — BP4 does not apply"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP5
# ---------------------------------------------------------------------------


def rule_bp5(
    variant: VariantInput,
    alternate_molecular_basis_found: bool,
) -> ACMGRule:
    """BP5 — Variant found in case with alternate molecular basis for disease.

    Args:
        variant: Annotated variant.
        alternate_molecular_basis_found: True if another variant fully
            explains the patient's phenotype.

    Returns:
        ACMGRule with rule_id ``"BP5"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP5.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if alternate_molecular_basis_found:
        return ACMGRule(
            rule_id="BP5",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                "Alternate molecular basis found that fully explains the phenotype"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP5",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=["No alternate molecular basis identified"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP6
# ---------------------------------------------------------------------------


def rule_bp6(variant: VariantInput) -> ACMGRule:
    """BP6 — Reputable source recently reported as benign.

    Mirrors PP5 but for benign classifications.  Uses ClinVar review status
    stars as proxy; ≥2 stars with a B/LB classification is sufficient.

    Args:
        variant: Annotated variant; ``clinvar_stars`` and
            ``clinvar_classification`` are key fields.

    Returns:
        ACMGRule with rule_id ``"BP6"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP6.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    benign_terms = {"Benign", "Likely benign", "Benign/Likely benign"}
    if (
        variant.clinvar_stars is not None
        and variant.clinvar_stars >= 2
        and variant.clinvar_classification in benign_terms
    ):
        return ACMGRule(
            rule_id="BP6",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                f"ClinVar: {variant.clinvar_classification} "
                f"({variant.clinvar_stars} stars — reputable source)"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP6",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=[
            "ClinVar does not meet BP6 criteria (insufficient stars or not benign)"
        ],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# BP7
# ---------------------------------------------------------------------------


def rule_bp7(
    variant: VariantInput,
    transcript: TranscriptData | None = None,
    no_splice_impact: bool | None = None,
    spliceai_delta: (
        float | None
    ) = None,  # alternative: derive no_splice_impact from score
    spliceai_score: float | None = None,  # alias for spliceai_delta
) -> ACMGRule:
    """BP7 — Synonymous variant with no predicted splice impact.

    For synonymous variants where in silico splicing predictors (SpliceAI,
    Pangolin) show no splice impact.  See also splicing.py where BP7 is
    assigned within the SpliceAI-based framework.

    Args:
        variant: Annotated variant; ``variant_type`` must be ``"synonymous"``.
        transcript: Transcript annotation.
        no_splice_impact: True if in silico splice predictors (SpliceAI Δ < 0.1
            and/or Pangolin low) indicate no splice impact.

    Returns:
        ACMGRule with rule_id ``"BP7"`` and strength SUPPORTING_BENIGN (-1 pt).

    References:
        Richards et al. 2015 PMID:25741868 — criterion BP7.
        Walker et al. 2023 PMID:36898414 — SpliceAI threshold for BP7.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "Walker et al. 2023 PMID:36898414 (SpliceAI BP7 threshold)",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    # Resolve no_splice_impact from SpliceAI score if not explicitly provided
    if no_splice_impact is None:
        delta = (
            spliceai_delta
            or spliceai_score
            or variant.spliceai_delta
            or variant.spliceai_max_delta
        )
        no_splice_impact = delta is not None and delta < 0.1

    if (
        variant.variant_type in ("synonymous", VariantType.SYNONYMOUS)
        and no_splice_impact
    ):
        return ACMGRule(
            rule_id="BP7",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                "Synonymous variant with no predicted splice impact "
                "(SpliceAI Δ < 0.1; Walker 2023 PMID:36898414)"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="BP7",
        strength=EvidenceStrength.SUPPORTING_BENIGN,
        evidence_items=["Not synonymous or potential splice impact predicted"],
        citations=citations,
        applies=False,
    )
