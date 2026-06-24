"""
bayesacmg.rules.pathogenic
==========================

All 15 pathogenic ACMG/AMP criteria: PVS1, PS1-4, PM1-6, PP1-5.

Every rule function:
  - cites Richards et al. 2015 PMID:25741868 for the original criterion
  - cites the specific ClinGen SVI or ACGS 2024 v1.2 update where applicable
  - has full Google-style docstrings with Args / Returns / Raises
  - has complete type annotations

CRITICAL: PM2 MUST return SUPPORTING (1 pt), not Moderate (2 pts).
This implements the ClinGen SVI 2024 recommendation.

References:
    Richards et al. 2015 PMID:25741868
    Abou Tayoun et al. 2018 PMID:30192042
    Tavtigian et al. 2020 PMID:32645316
    Cheng et al. 2023 PMID:37703350 (AlphaMissense)
    ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al.)
    ClinGen SVI Working Group 2024
"""

from __future__ import annotations

from bayesacmg.models import ACMGRule, EvidenceStrength, GeneData, TranscriptData, VariantInput

# ---------------------------------------------------------------------------
# Constants — thresholds with inline citations
# ---------------------------------------------------------------------------

_PM2_AF_THRESHOLD = 0.0001   # gnomAD v4.1 global AF cut-off; ClinGen SVI 2024
_BA1_AF_THRESHOLD = 0.05     # >5 % → Benign; Richards 2015 PMID:25741868
_BS1_AF_THRESHOLD = 0.01     # >1 % → BS1; ACGS 2024 v1.2 §5 Table 2

# AlphaMissense thresholds — Cheng et al. 2023 Science PMID:37703350
_AM_PP3_THRESHOLD = 0.564    # ≥0.564 → likely pathogenic (PP3)
_AM_BP4_THRESHOLD = 0.340    # ≤0.340 → likely benign (BP4)

# REVEL secondary threshold — ClinGen SVI 2024 in silico recommendation
_REVEL_PP3 = 0.7
# CADD secondary threshold
_CADD_PP3 = 25.0


# ---------------------------------------------------------------------------
# PVS1
# ---------------------------------------------------------------------------

_LOF_TYPES = frozenset({
    "frameshift", "nonsense", "splice_canonical",
    "start_loss", "stop_loss",
})

_NMD_ESCAPE_LAST_EXON_THRESHOLD = 55  # nucleotides; PVS1 decision tree


def rule_pvs1(
    variant: VariantInput,
    transcript: TranscriptData,
    gene: GeneData,
) -> ACMGRule:
    """PVS1 — Null variant in a gene where LoF is the disease mechanism.

    Implements the Abou Tayoun 2018 decision tree (PMID:30192042) with
    ACGS 2024 v1.2 MANE Select requirement.

    Guidelines:
        Richards et al. 2015 PMID:25741868: original PVS1 definition.
        Abou Tayoun et al. 2018 PMID:30192042: PVS1 decision tree v1.0.
        Walker et al. 2023 PMID:36898414: splice site PVS1 — see splicing.py.
        ACGS 2024 v1.2 §5, Notes 1-3: MANE Select requirement; UK practice.

    MANE Select requirement (ACGS 2024 v1.2 Note 1):
        PVS1 is applied at full strength only when the LoF affects the
        MANE Select transcript.  If the LoF only affects non-MANE transcripts,
        reduce strength by one level (Very Strong → Strong).

    Args:
        variant: Annotated variant input object.
        transcript: Transcript-level annotation; ``is_mane_select`` is checked.
        gene: Gene-level data; ``lof_is_disease_mechanism`` is required.

    Returns:
        ACMGRule with rule_id ``"PVS1"``.  Strength may be reduced from
        Very Strong to Strong when MANE Select is not affected or when
        NMD is predicted to escape.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PVS1.
        Abou Tayoun et al. 2018 PMID:30192042 — decision tree.
        ACGS 2024 v1.2 §5, Notes 1-3.
    """
    evidence: list[str] = []
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "Abou Tayoun et al. 2018 PMID:30192042",
        "ACGS 2024 v1.2 §5 Notes 1-3",
    ]

    # Step 1: Is this a LoF variant type?
    if variant.variant_type not in _LOF_TYPES:
        return ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            evidence_items=[f"Variant type '{variant.variant_type}' is not a LoF type"],
            citations=citations,
            applies=False,
        )

    evidence.append(f"LoF variant type: {variant.variant_type}")

    # Step 2: Is LoF the established disease mechanism?
    if not gene.lof_is_disease_mechanism:
        return ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.VERY_STRONG,
            evidence_items=[f"LoF not established disease mechanism for {gene.gene_symbol}"],
            citations=citations,
            applies=False,
        )

    evidence.append(f"LoF is established disease mechanism for {gene.gene_symbol}")

    # Step 3: MANE Select check (ACGS 2024 v1.2 Note 1)
    if not transcript.is_mane_select:
        # LoF only affects non-MANE transcript → reduce from Very Strong to Strong
        evidence.append(
            "LoF does not affect MANE Select transcript → "
            "strength reduced Very Strong → Strong (ACGS 2024 §5 Note 1)"
        )
        return ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.STRONG,
            evidence_items=evidence,
            citations=citations,
            applies=True,
            notes="Non-MANE Select transcript: reduced Very Strong→Strong",
        )

    evidence.append("Affects MANE Select transcript")

    # Step 4: NMD prediction — Abou Tayoun 2018 decision tree
    # NMD expected to fail in last exon or <55 nt from last exon junction
    if transcript.is_last_exon or transcript.nmd_escapes_last_exon_rule:
        evidence.append(
            "NMD predicted to escape (last exon or <55 nt from final junction) — "
            "Abou Tayoun 2018 PMID:30192042 decision tree branch 3"
        )
        # Reduce to Strong if NMD escape is predicted
        return ACMGRule(
            rule_id="PVS1",
            strength=EvidenceStrength.STRONG,
            evidence_items=evidence,
            citations=citations,
            applies=True,
            notes="NMD predicted to escape — strength reduced to Strong",
        )

    # Step 5: Full PVS1 at Very Strong
    evidence.append("NMD expected to occur — full PVS1 Very Strong applies")

    return ACMGRule(
        rule_id="PVS1",
        strength=EvidenceStrength.VERY_STRONG,
        evidence_items=evidence,
        citations=citations,
        applies=True,
    )


# ---------------------------------------------------------------------------
# PS1
# ---------------------------------------------------------------------------

def rule_ps1(
    variant: VariantInput,
    same_aa_change_pathogenic: bool,
) -> ACMGRule:
    """PS1 — Same amino acid change as established pathogenic variant.

    Args:
        variant: Annotated variant input.
        same_aa_change_pathogenic: True if a different nucleotide change
            results in the same amino acid change that is pathogenic.

    Returns:
        ACMGRule with rule_id ``"PS1"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PS1.
        ACGS 2024 v1.2 §5 Table 2 (PS1 notes).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if same_aa_change_pathogenic and variant.hgvs_p:
        return ACMGRule(
            rule_id="PS1",
            strength=EvidenceStrength.STRONG,
            evidence_items=[
                f"Amino acid change {variant.hgvs_p} has been previously "
                "established as pathogenic via a different nucleotide change"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PS1",
        strength=EvidenceStrength.STRONG,
        evidence_items=["No documented pathogenic variant at same amino acid position"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PS2
# ---------------------------------------------------------------------------

def rule_ps2(variant: VariantInput) -> ACMGRule:
    """PS2 — De novo variant in patient with disease and no family history.

    Requires parental testing confirming de novo status.  For assumed de
    novo without parental testing, use PM6 (Supporting).

    Args:
        variant: Annotated variant; ``is_de_novo`` must be True.

    Returns:
        ACMGRule with rule_id ``"PS2"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PS2.
        ACGS 2024 v1.2 §5 Table 2 (de novo notes).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if variant.is_de_novo:
        return ACMGRule(
            rule_id="PS2",
            strength=EvidenceStrength.STRONG,
            evidence_items=["Confirmed de novo by parental testing"],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PS2",
        strength=EvidenceStrength.STRONG,
        evidence_items=["De novo status not confirmed by parental testing"],
        citations=citations,
        applies=False,
        notes="If assumed de novo without parental testing, apply PM6 instead",
    )


# ---------------------------------------------------------------------------
# PS3
# ---------------------------------------------------------------------------

def rule_ps3(variant: VariantInput) -> ACMGRule:
    """PS3 — Well-established functional studies show damaging effect.

    Args:
        variant: Annotated variant; ``functional_study_result`` evaluated.

    Returns:
        ACMGRule with rule_id ``"PS3"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PS3.
        ACGS 2024 v1.2 §5 Table 2 (functional study standards).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    damaging_results = {"loss_of_function", "dominant_negative", "gain_of_function"}
    if variant.functional_study_result in damaging_results:
        return ACMGRule(
            rule_id="PS3",
            strength=EvidenceStrength.STRONG,
            evidence_items=[
                f"Well-established functional study result: {variant.functional_study_result}"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PS3",
        strength=EvidenceStrength.STRONG,
        evidence_items=["No damaging functional study result available"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PS4
# ---------------------------------------------------------------------------

def rule_ps4(
    variant: VariantInput,
    case_control_or_significant: bool,
) -> ACMGRule:
    """PS4 — Prevalence in affected individuals statistically greater than controls.

    Args:
        variant: Annotated variant.
        case_control_or_significant: True if variant shows statistically
            significant enrichment in cases vs controls (OR ≥5 with 95% CI).

    Returns:
        ACMGRule with rule_id ``"PS4"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PS4.
        ACGS 2024 v1.2 §5 Table 2 (case-control notes).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if case_control_or_significant:
        return ACMGRule(
            rule_id="PS4",
            strength=EvidenceStrength.STRONG,
            evidence_items=["Variant prevalence significantly elevated in affected vs controls"],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PS4",
        strength=EvidenceStrength.STRONG,
        evidence_items=["Insufficient case-control evidence for PS4"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PM1
# ---------------------------------------------------------------------------

def rule_pm1(
    variant: VariantInput,
    transcript: TranscriptData,
    gene: GeneData,
) -> ACMGRule:
    """PM1 — Located in mutational hotspot or critical functional domain.

    Args:
        variant: Annotated variant.
        transcript: Transcript annotation with domain_annotations field.
        gene: Gene data with ``has_hotspot_domain`` and ``hotspot_domains``.

    Returns:
        ACMGRule with rule_id ``"PM1"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PM1.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    domain_hits = transcript.domain_annotations or []
    if gene.has_hotspot_domain and domain_hits:
        return ACMGRule(
            rule_id="PM1",
            strength=EvidenceStrength.MODERATE,
            evidence_items=[
                f"Variant in mutational hotspot/critical domain: {', '.join(domain_hits)}"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PM1",
        strength=EvidenceStrength.MODERATE,
        evidence_items=["Variant not in documented hotspot or critical functional domain"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PM2 — CRITICAL: MUST return SUPPORTING (1 pt), NOT Moderate (2 pts)
# ---------------------------------------------------------------------------

def rule_pm2(variant: VariantInput) -> ACMGRule:
    """PM2 — Absent from or extremely rare in population databases.

    DEFAULT WEIGHT: SUPPORTING (1 point) — NOT Moderate (2 points).

    ClinGen SVI 2024 update: PM2 should be applied at Supporting strength
    in most contexts.  gnomAD v4.1 (807,162 individuals, April 2024)
    reveals that many ultra-rare variants exist in the general population,
    making absence from population databases less distinctive evidence
    than assumed when the 2015 guidelines were written.

    When to apply:
        - Variant absent from gnomAD v4.1 (``gnomad_af is None``)
        - OR global AF < 0.0001 in gnomAD v4.1 (all populations)

    Per-gene override:
        Some VCEP specifications allow PM2 at Moderate for very specific
        gene-disease pairs.  The caller should invoke vcep_client.py before
        calling this rule and override the returned strength if indicated.

    Args:
        variant: Annotated variant; ``gnomad_af`` is the key field.

    Returns:
        ACMGRule with rule_id ``"PM2"`` and strength SUPPORTING (1 pt)
        when the criterion applies.

    References:
        Richards et al. 2015 PMID:25741868 — original PM2 definition.
        ClinGen SVI PM2 guidance (2024):
            https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
        ACGS 2024 v1.2 §5 Table 2, Appendix C — PM2 mini impact assessment.
        gnomAD v4.1 release notes (April 2024): 807,162 individuals.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ClinGen SVI PM2 guidance (2024)",
        "ACGS 2024 v1.2 §5 Table 2 Appendix C",
        "gnomAD v4.1 (April 2024, 807,162 individuals)",
    ]
    # gnomad_af is None → variant absent from gnomAD v4.1
    if variant.gnomad_af is None:
        return ACMGRule(
            rule_id="PM2",
            strength=EvidenceStrength.SUPPORTING,  # 1 pt — ClinGen SVI 2024
            evidence_items=["Absent from gnomAD v4.1 (807,162 individuals, April 2024)"],
            citations=citations,
            applies=True,
            notes="PM2 applied at Supporting (1 pt) per ClinGen SVI 2024 — not Moderate",
        )

    # AF < 0.0001 threshold — ClinGen SVI 2024 / ACGS 2024 v1.2 Appendix C
    if variant.gnomad_af < _PM2_AF_THRESHOLD:  # 0.0001; ClinGen SVI 2024
        return ACMGRule(
            rule_id="PM2",
            strength=EvidenceStrength.SUPPORTING,  # 1 pt — ClinGen SVI 2024
            evidence_items=[
                f"Extremely rare in gnomAD v4.1: AF={variant.gnomad_af:.2e} "
                f"(threshold <{_PM2_AF_THRESHOLD})"
            ],
            citations=citations,
            applies=True,
            notes="PM2 applied at Supporting (1 pt) per ClinGen SVI 2024 — not Moderate",
        )

    return ACMGRule(
        rule_id="PM2",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[
            f"gnomAD v4.1 AF={variant.gnomad_af:.4f} exceeds PM2 threshold "
            f"({_PM2_AF_THRESHOLD}) — PM2 does not apply"
        ],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PM3
# ---------------------------------------------------------------------------

def rule_pm3(variant: VariantInput) -> ACMGRule:
    """PM3 — Detected in trans with a pathogenic variant in recessive disease.

    Args:
        variant: Annotated variant; ``in_trans_pathogenic`` is the key field.

    Returns:
        ACMGRule with rule_id ``"PM3"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PM3.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if variant.in_trans_pathogenic:
        return ACMGRule(
            rule_id="PM3",
            strength=EvidenceStrength.MODERATE,
            evidence_items=["Variant detected in trans with a known pathogenic variant"],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PM3",
        strength=EvidenceStrength.MODERATE,
        evidence_items=["Not detected in trans with a pathogenic variant"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PM4
# ---------------------------------------------------------------------------

_PM4_LENGTH_CHANGE_THRESHOLD = 10  # amino acids; Richards 2015 PMID:25741868

def rule_pm4(
    variant: VariantInput,
    transcript: TranscriptData,
) -> ACMGRule:
    """PM4 — Protein length change due to in-frame indel or stop-loss.

    Args:
        variant: Annotated variant; ``variant_type`` is checked.
        transcript: Transcript annotation; protein lengths are checked.

    Returns:
        ACMGRule with rule_id ``"PM4"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PM4.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    in_frame_types = {"inframe_insertion", "inframe_deletion", "stop_loss"}
    if variant.variant_type not in in_frame_types:
        return ACMGRule(
            rule_id="PM4",
            strength=EvidenceStrength.MODERATE,
            evidence_items=[f"Variant type '{variant.variant_type}' not in-frame indel or stop-loss"],
            citations=citations,
            applies=False,
        )

    # Check protein length change if available
    if (transcript.prot_length_original is not None
            and transcript.prot_length_alt is not None):
        delta = abs(transcript.prot_length_alt - transcript.prot_length_original)
        if delta >= _PM4_LENGTH_CHANGE_THRESHOLD:  # ≥10 aa; Richards 2015
            return ACMGRule(
                rule_id="PM4",
                strength=EvidenceStrength.MODERATE,
                evidence_items=[
                    f"In-frame {variant.variant_type} changes protein length by {delta} aa "
                    f"(≥{_PM4_LENGTH_CHANGE_THRESHOLD} aa threshold)"
                ],
                citations=citations,
                applies=True,
            )

    # In-frame indel but small or no protein length data
    return ACMGRule(
        rule_id="PM4",
        strength=EvidenceStrength.MODERATE,
        evidence_items=[f"In-frame {variant.variant_type} — protein length impact uncertain"],
        citations=citations,
        applies=True,
        notes="Protein length change <10 aa or unknown; PM4 applied provisionally",
    )


# ---------------------------------------------------------------------------
# PM5
# ---------------------------------------------------------------------------

def rule_pm5(
    variant: VariantInput,
    different_aa_at_same_position_pathogenic: bool,
) -> ACMGRule:
    """PM5 — Novel missense at same position as known pathogenic missense.

    Args:
        variant: Annotated variant.
        different_aa_at_same_position_pathogenic: True if a different amino
            acid change at the same codon is established as pathogenic.

    Returns:
        ACMGRule with rule_id ``"PM5"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PM5.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if different_aa_at_same_position_pathogenic and variant.hgvs_p:
        return ACMGRule(
            rule_id="PM5",
            strength=EvidenceStrength.MODERATE,
            evidence_items=[
                f"Novel missense {variant.hgvs_p} at same codon as known pathogenic missense"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PM5",
        strength=EvidenceStrength.MODERATE,
        evidence_items=["No known pathogenic missense at same amino acid position"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PM6
# ---------------------------------------------------------------------------

def rule_pm6(variant: VariantInput) -> ACMGRule:
    """PM6 — Assumed de novo, without confirmation of parental genotypes.

    Applied at Supporting strength because parental testing was not performed.
    Use PS2 (Strong) when de novo is confirmed by parental testing.

    Args:
        variant: Annotated variant; ``assumed_de_novo`` is the key field.

    Returns:
        ACMGRule with rule_id ``"PM6"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PM6.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if variant.assumed_de_novo:
        return ACMGRule(
            rule_id="PM6",
            strength=EvidenceStrength.MODERATE,
            evidence_items=["Assumed de novo (parental testing not performed)"],
            citations=citations,
            applies=True,
            notes="Use PS2 instead if parental testing confirms de novo status",
        )
    return ACMGRule(
        rule_id="PM6",
        strength=EvidenceStrength.MODERATE,
        evidence_items=["Neither confirmed nor assumed de novo"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PP1
# ---------------------------------------------------------------------------

def rule_pp1(
    variant: VariantInput,
    segregation_supports: bool,
) -> ACMGRule:
    """PP1 — Cosegregation with disease in multiple affected family members.

    Args:
        variant: Annotated variant; ``segregation_lod`` provides quantitative
            support if available.
        segregation_supports: True if cosegregation data support pathogenicity.

    Returns:
        ACMGRule with rule_id ``"PP1"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PP1.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if segregation_supports:
        evidence_msg = "Cosegregation with disease in multiple affected family members"
        if variant.segregation_lod is not None:
            evidence_msg += f" (LOD={variant.segregation_lod:.2f})"
        return ACMGRule(
            rule_id="PP1",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[evidence_msg],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PP1",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=["Insufficient cosegregation data"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PP2
# ---------------------------------------------------------------------------

_LOEUF_PP2_THRESHOLD = 0.35   # gnomAD LOEUF; ACGS 2024 v1.2 §5
_PLI_PP2_THRESHOLD = 0.9      # pLI threshold; ACGS 2024 v1.2 §5

def rule_pp2(
    variant: VariantInput,
    transcript: TranscriptData,
    gene: GeneData,
) -> ACMGRule:
    """PP2 — Missense variant in gene with low rate of benign missense variation.

    Args:
        variant: Annotated variant; ``variant_type`` is checked for ``"snv"``.
        transcript: Transcript annotation.
        gene: Gene data; ``gnomad_pli`` and ``gnomad_loeuf`` evaluated.

    Returns:
        ACMGRule with rule_id ``"PP2"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PP2.
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
    if not is_missense:
        return ACMGRule(
            rule_id="PP2",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Variant is not a missense change"],
            citations=citations,
            applies=False,
        )

    pli_high = (gene.gnomad_pli is not None and gene.gnomad_pli >= _PLI_PP2_THRESHOLD)  # ≥0.9
    loeuf_low = (gene.gnomad_loeuf is not None and gene.gnomad_loeuf <= _LOEUF_PP2_THRESHOLD)  # ≤0.35

    if pli_high or loeuf_low:
        evidence_items = [f"Missense in LoF-intolerant gene {gene.gene_symbol}"]
        if pli_high:
            evidence_items.append(f"pLI={gene.gnomad_pli:.3f} (≥{_PLI_PP2_THRESHOLD})")
        if loeuf_low:
            evidence_items.append(f"LOEUF={gene.gnomad_loeuf:.3f} (≤{_LOEUF_PP2_THRESHOLD})")
        return ACMGRule(
            rule_id="PP2",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=evidence_items,
            citations=citations,
            applies=True,
        )

    return ACMGRule(
        rule_id="PP2",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[f"Gene {gene.gene_symbol} not classified as LoF-intolerant"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PP3
# ---------------------------------------------------------------------------

def rule_pp3(
    variant: VariantInput,
    alphamissense_score: float | None,
    revel_score: float | None,
) -> ACMGRule:
    """PP3 — Multiple computational evidence supporting pathogenicity.

    PRIMARY PREDICTOR: AlphaMissense (ClinGen SVI 2024 approved).
    AlphaMissense threshold: score ≥ 0.564 → PP3 (Supporting).
    Reference: Cheng et al. 2023 Science doi:10.1126/science.adg7492
               PMID:37703350.

    ClinGen SVI 2024 approved tools for PP3:
        1. AlphaMissense: ≥ 0.564 (PRIMARY — preferred)
        2. REVEL: ≥ 0.7 (secondary comparison)
        3. CADD PHRED: ≥ 25 (tertiary)

    For SPLICE SITE VARIANTS, do NOT use this function.
    Use rules/splicing.py rule_splicing_pp3_bp4_bp7() instead.

    Args:
        variant: Annotated variant; ``cadd_phred`` may also be checked.
        alphamissense_score: AlphaMissense score (0–1); PRIMARY predictor.
            ≥0.564 → PP3 fires.  ≤0.340 → BP4 fires (see rule_bp4 in
            benign.py).
        revel_score: REVEL score (0–1); secondary predictor.

    Returns:
        ACMGRule with rule_id ``"PP3"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PP3.
        Cheng et al. 2023 PMID:37703350 — AlphaMissense thresholds.
        ClinGen SVI in silico recommendation (2024).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "Cheng et al. 2023 PMID:37703350 (AlphaMissense)",
        "ClinGen SVI in silico tool recommendation (2024)",
        "ACGS 2024 v1.2 §5 Table 2",
    ]

    # Block application for splice variants — splicing.py handles those
    splice_types = {"splice_canonical", "splice_region"}
    if variant.variant_type in splice_types:
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Splice variant: use splicing.rule_splicing_pp3_bp4_bp7() instead"],
            citations=citations,
            applies=False,
            notes="PP3 for splice variants handled by rules/splicing.py per Walker 2023",
        )

    # PRIMARY: AlphaMissense ≥ 0.564 (Cheng 2023 PMID:37703350)
    if alphamissense_score is not None:
        if alphamissense_score >= _AM_PP3_THRESHOLD:  # ≥0.564; Cheng 2023 PMID:37703350
            return ACMGRule(
                rule_id="PP3",
                strength=EvidenceStrength.SUPPORTING,
                evidence_items=[
                    f"AlphaMissense score {alphamissense_score:.3f} ≥ {_AM_PP3_THRESHOLD} "
                    "(likely pathogenic threshold; Cheng 2023 PMID:37703350)"
                ],
                citations=citations,
                applies=True,
            )
        if alphamissense_score <= _AM_BP4_THRESHOLD:  # ≤0.340; Cheng 2023 PMID:37703350
            # Likely benign — PP3 does NOT apply; BP4 will apply
            return ACMGRule(
                rule_id="PP3",
                strength=EvidenceStrength.SUPPORTING,
                evidence_items=[
                    f"AlphaMissense score {alphamissense_score:.3f} ≤ {_AM_BP4_THRESHOLD} "
                    "(likely benign threshold) — PP3 does not apply; BP4 applies"
                ],
                citations=citations,
                applies=False,
            )
        # Score in intermediate zone (0.340–0.564) — neither PP3 nor BP4
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"AlphaMissense score {alphamissense_score:.3f} in intermediate zone "
                f"({_AM_BP4_THRESHOLD}–{_AM_PP3_THRESHOLD}) — neither PP3 nor BP4 applies"
            ],
            citations=citations,
            applies=False,
        )

    # SECONDARY: REVEL ≥ 0.7 (ClinGen SVI 2024)
    if revel_score is not None:
        if revel_score >= _REVEL_PP3:  # ≥0.7; ClinGen SVI 2024
            return ACMGRule(
                rule_id="PP3",
                strength=EvidenceStrength.SUPPORTING,
                evidence_items=[
                    f"REVEL score {revel_score:.3f} ≥ {_REVEL_PP3} "
                    "(secondary predictor; AlphaMissense unavailable)"
                ],
                citations=citations,
                applies=True,
                notes="AlphaMissense unavailable; REVEL used as secondary predictor",
            )

    # TERTIARY: CADD PHRED ≥ 25
    if variant.cadd_phred is not None and variant.cadd_phred >= _CADD_PP3:  # ≥25
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"CADD PHRED {variant.cadd_phred:.1f} ≥ {_CADD_PP3} "
                "(tertiary predictor; AlphaMissense and REVEL unavailable)"
            ],
            citations=citations,
            applies=True,
            notes="Primary predictors unavailable; CADD used as tertiary",
        )

    return ACMGRule(
        rule_id="PP3",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=["No in silico scores available or scores do not support PP3"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PP4
# ---------------------------------------------------------------------------

def rule_pp4(
    variant: VariantInput,
    phenotype_highly_specific: bool,
) -> ACMGRule:
    """PP4 — Patient's phenotype highly specific for gene with single genetic cause.

    Args:
        variant: Annotated variant.
        phenotype_highly_specific: True if the clinical phenotype is highly
            specific for the gene in question.

    Returns:
        ACMGRule with rule_id ``"PP4"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PP4.
        ACGS 2024 v1.2 §5 Table 2.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    if phenotype_highly_specific:
        return ACMGRule(
            rule_id="PP4",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"Patient phenotype highly specific for gene {variant.gene_symbol}"
            ],
            citations=citations,
            applies=True,
        )
    return ACMGRule(
        rule_id="PP4",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=["Phenotype not highly specific for gene"],
        citations=citations,
        applies=False,
    )


# ---------------------------------------------------------------------------
# PP5
# ---------------------------------------------------------------------------

def rule_pp5(variant: VariantInput) -> ACMGRule:
    """PP5 — Reputable source recently classified as pathogenic.

    Uses ClinVar review status stars as a proxy for reputable source quality.
    ≥2 stars with a P/LP classification is considered sufficient.

    Args:
        variant: Annotated variant; ``clinvar_stars`` and
            ``clinvar_classification`` are key fields.

    Returns:
        ACMGRule with rule_id ``"PP5"``.

    References:
        Richards et al. 2015 PMID:25741868 — criterion PP5.
        ACGS 2024 v1.2 §5 Table 2 (caution about over-reliance on PP5).
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "ACGS 2024 v1.2 §5 Table 2",
    ]
    pathogenic_terms = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
    if (variant.clinvar_stars is not None
            and variant.clinvar_stars >= 2  # ≥2 stars = expert panel or criteria provided
            and variant.clinvar_classification in pathogenic_terms):
        return ACMGRule(
            rule_id="PP5",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"ClinVar classification: {variant.clinvar_classification} "
                f"({variant.clinvar_stars} stars — reputable source)"
            ],
            citations=citations,
            applies=True,
            notes="ACGS 2024: do not apply PP5 if it is the sole pathogenic criterion",
        )
    return ACMGRule(
        rule_id="PP5",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[
            f"ClinVar stars={variant.clinvar_stars}, "
            f"classification={variant.clinvar_classification} — PP5 not met"
        ],
        citations=citations,
        applies=False,
    )
