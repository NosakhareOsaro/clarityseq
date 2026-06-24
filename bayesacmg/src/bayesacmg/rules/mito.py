"""
bayesacmg.rules.mito
====================

Mitochondrial variant-specific ACMG/AMP rules.

Implements ACGS Best Practice Guidelines 2024 v1.2 §6 mito-specific rules.
These rules MUST be evaluated BEFORE the standard nuclear variant rules for
any chrM variant.  haplogroup status (Haplogrep3 output) MUST be checked
before any pathogenicity assessment.

Key differences from nuclear rules (ACGS 2024 §6):
    - BA1: Standard 5% AF threshold does NOT apply to mito variants.
      Population-specific haplogroup frequency data (MITOMAP / gnomAD) is
      used instead.
    - PM2: Different AF thresholds; haplogroup must be checked first.
    - Heteroplasmy level maps to clinical significance category.
    - Haplogroup-defining variants are classified Benign first.

References:
    ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al.) §6
    Richards et al. 2015 PMID:25741868 (original framework)
    ClinGen SVI Working Group 2024
    MITOMAP: A Human Mitochondrial Genome Database. https://www.mitomap.org
    gnomAD v3.1 mitogenome release (mito-specific population data)
"""

from __future__ import annotations

from dataclasses import dataclass

from bayesacmg.models import ACMGRule, EvidenceStrength, VariantInput

# ---------------------------------------------------------------------------
# Mito-specific constants — ACGS 2024 §6 thresholds
# ---------------------------------------------------------------------------

# Mito BA1: per ACGS 2024 §6 — standard 5% does NOT apply.
# Use MITOMAP confirmed polymorphism status or gnomAD mito AF >5% within
# a specific haplogroup as a proxy, but NOT as stand-alone BA1.
_MITO_BA1_HAPLOGROUP_THRESHOLD = 0.05  # >5% within haplogroup → strong benign

# Mito PM2: variants absent from or very rare in gnomAD v3.1 mito genome
_MITO_PM2_AF_THRESHOLD = 0.0001         # <0.01% in gnomAD mito data; ACGS 2024 §6

# Heteroplasmy thresholds — ACGS 2024 §6
_HETEROPLASMY_HIGH = 0.60     # >60% — likely to cause disease
_HETEROPLASMY_MODERATE = 0.20 # 20-60% — moderate level
_HETEROPLASMY_LOW = 0.01      # <1% — likely passenger/artefact


@dataclass
class MitoHaploData:
    """Mitochondrial haplogroup data for a variant.

    Attributes:
        haplogroup: Haplogrep3 haplogroup string (e.g. ``"H1a1"``).
        is_haplogroup_defining: True if the variant is a haplogroup-defining
            polymorphism per Haplogrep3 / PhyloTree Build 17.
        haplogroup_frequency: Frequency of the variant within its haplogroup
            (0–1); from MITOMAP or gnomAD mito data.
        mitomap_status: MITOMAP variant status string, e.g.
            ``"Confirmed Pathogenic"``, ``"Polymorphism"``, or None.
    """

    haplogroup: str = ""
    is_haplogroup_defining: bool = False
    haplogroup_frequency: float | None = None
    mitomap_status: str | None = None


def rule_mito_haplogroup_defining(
    variant: VariantInput,
    haplo_data: MitoHaploData,
) -> ACMGRule:
    """BA1-equivalent for haplogroup-defining mito variants.

    Must be evaluated FIRST for any mitochondrial variant.  A haplogroup-
    defining variant per Haplogrep3 / PhyloTree Build 17 is classified as
    Benign (stand-alone), overriding all pathogenicity evidence.

    Args:
        variant: Annotated mitochondrial variant; ``is_mito`` should be True.
        haplo_data: Haplogroup data from Haplogrep3 analysis.

    Returns:
        ACMGRule with rule_id ``"BA1_MITO_HAPLO"`` and strength STAND_ALONE
        if the variant is haplogroup-defining.  Otherwise applies=False and
        classification proceeds to other rules.

    References:
        ACGS 2024 v1.2 §6 — haplogroup-first classification.
        PhyloTree Build 17 (van Oven 2009).
        Haplogrep3 (Weissensteiner et al. 2021).

    Raises:
        ValueError: If ``variant.is_mito`` is False.
    """
    citations = [
        "ACGS 2024 v1.2 §6",
        "Richards et al. 2015 PMID:25741868",
        "PhyloTree Build 17 (van Oven 2009)",
        "Haplogrep3 (Weissensteiner et al. 2021)",
    ]
    if not variant.is_mito:
        return ACMGRule(
            rule_id="BA1_MITO_HAPLO",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=["Not a mitochondrial variant"],
            citations=citations,
            applies=False,
        )

    if haplo_data.is_haplogroup_defining:
        return ACMGRule(
            rule_id="BA1_MITO_HAPLO",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=[
                f"Haplogroup-defining variant for haplogroup {haplo_data.haplogroup} "
                "per Haplogrep3/PhyloTree Build 17 — classified Benign (ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes=(
                "ACGS 2024 §6: haplogroup-defining variants are classified Benign. "
                "This must be checked BEFORE any pathogenicity assessment."
            ),
        )

    return ACMGRule(
        rule_id="BA1_MITO_HAPLO",
        strength=EvidenceStrength.STAND_ALONE,
        evidence_items=[
            f"Not a haplogroup-defining variant (haplogroup: {haplo_data.haplogroup or 'unknown'})"
        ],
        citations=citations,
        applies=False,
    )


def rule_mito_ba1(
    variant: VariantInput,
    haplo_data: MitoHaploData,
) -> ACMGRule:
    """BA1 for mitochondrial variants — mito-specific thresholds.

    The standard nuclear BA1 threshold (>5% global AF) does NOT apply to
    mitochondrial variants.  Per ACGS 2024 §6:
    - Use haplogroup-specific frequency data from MITOMAP / gnomAD mito.
    - A variant >5% within its haplogroup may support benign classification
      but requires corroborating evidence (not stand-alone).
    - MITOMAP "Confirmed Polymorphism" status provides strong benign evidence.

    Args:
        variant: Annotated mitochondrial variant.
        haplo_data: Haplogroup data including haplogroup-specific frequency
            and MITOMAP status.

    Returns:
        ACMGRule with rule_id ``"BA1_MITO"`` reflecting the mito-specific
        BA1 assessment.

    References:
        ACGS 2024 v1.2 §6 — mito BA1 caveat.
        MITOMAP: https://www.mitomap.org
        gnomAD v3.1 mitogenome data.
    """
    citations = [
        "ACGS 2024 v1.2 §6",
        "Richards et al. 2015 PMID:25741868",
        "MITOMAP: https://www.mitomap.org",
        "gnomAD v3.1 mitogenome data",
    ]
    if not variant.is_mito:
        return ACMGRule(
            rule_id="BA1_MITO",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=["Not a mitochondrial variant"],
            citations=citations,
            applies=False,
        )

    # MITOMAP "Confirmed Polymorphism" → strong benign evidence
    if haplo_data.mitomap_status and "polymorphism" in haplo_data.mitomap_status.lower():
        return ACMGRule(
            rule_id="BA1_MITO",
            strength=EvidenceStrength.STAND_ALONE,
            evidence_items=[
                f"MITOMAP status: '{haplo_data.mitomap_status}' — confirmed benign polymorphism"
            ],
            citations=citations,
            applies=True,
            notes="MITOMAP Confirmed Polymorphism satisfies mito BA1; ACGS 2024 §6",
        )

    # High haplogroup-specific frequency → benign signal (not stand-alone per ACGS 2024 §6)
    if (haplo_data.haplogroup_frequency is not None
            and haplo_data.haplogroup_frequency > _MITO_BA1_HAPLOGROUP_THRESHOLD):
        return ACMGRule(
            rule_id="BA1_MITO",
            strength=EvidenceStrength.STRONG_BENIGN,  # Strong (not STAND_ALONE) per ACGS 2024 §6
            evidence_items=[
                f"Haplogroup-specific frequency {haplo_data.haplogroup_frequency:.3f} > "
                f"{_MITO_BA1_HAPLOGROUP_THRESHOLD} within haplogroup {haplo_data.haplogroup} "
                "(mito-specific BA1; ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes=(
                "ACGS 2024 §6: high haplogroup frequency → Strong Benign (not stand-alone). "
                "Standard nuclear 5% BA1 threshold does NOT apply to mito variants."
            ),
        )

    return ACMGRule(
        rule_id="BA1_MITO",
        strength=EvidenceStrength.STAND_ALONE,
        evidence_items=["Mito BA1 criteria not met"],
        citations=citations,
        applies=False,
    )


def rule_mito_pm2(
    variant: VariantInput,
    haplo_data: MitoHaploData,
) -> ACMGRule:
    """PM2 for mitochondrial variants — mito-specific thresholds.

    Applies PM2 at SUPPORTING weight per ClinGen SVI 2024 recommendation,
    using mito-specific gnomAD v3.1 AF data.  Haplogroup status MUST be
    checked before calling this rule (rule_mito_haplogroup_defining()).

    Mito PM2 differences from nuclear PM2 (ACGS 2024 §6):
    - Use gnomAD v3.1 mitogenome allele frequency, not gnomAD v4.1 nuclear.
    - Check haplogroup-specific frequency if global frequency is low.
    - AF threshold is the same (< 0.0001) but applied to mito-specific data.

    Args:
        variant: Annotated mitochondrial variant.
        haplo_data: Haplogroup data; haplogroup-defining check must have been
            performed prior to calling this function.

    Returns:
        ACMGRule with rule_id ``"PM2_MITO"`` and strength SUPPORTING (1 pt)
        when the criterion applies.

    References:
        ACGS 2024 v1.2 §6.
        ClinGen SVI PM2 guidance (2024).
        gnomAD v3.1 mitogenome release.
    """
    citations = [
        "ACGS 2024 v1.2 §6",
        "Richards et al. 2015 PMID:25741868",
        "ClinGen SVI PM2 guidance (2024)",
        "gnomAD v3.1 mitogenome data",
    ]
    if not variant.is_mito:
        return ACMGRule(
            rule_id="PM2_MITO",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Not a mitochondrial variant"],
            citations=citations,
            applies=False,
        )

    if haplo_data.is_haplogroup_defining:
        return ACMGRule(
            rule_id="PM2_MITO",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Haplogroup-defining variant — PM2 does not apply; classify as Benign"],
            citations=citations,
            applies=False,
            notes="Haplogroup-defining status takes precedence; see rule_mito_haplogroup_defining()",
        )

    af = variant.gnomad_af  # gnomAD v3.1 mito AF expected here
    if af is None:
        return ACMGRule(
            rule_id="PM2_MITO",
            strength=EvidenceStrength.SUPPORTING,  # 1 pt — ClinGen SVI 2024
            evidence_items=[
                "Absent from gnomAD v3.1 mitogenome dataset "
                "(mito-specific PM2; ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes="PM2 at Supporting (1 pt) per ClinGen SVI 2024",
        )

    if af < _MITO_PM2_AF_THRESHOLD:  # <0.0001; ACGS 2024 §6
        return ACMGRule(
            rule_id="PM2_MITO",
            strength=EvidenceStrength.SUPPORTING,  # 1 pt — ClinGen SVI 2024
            evidence_items=[
                f"Extremely rare in gnomAD v3.1 mitogenome: AF={af:.2e} "
                f"(< {_MITO_PM2_AF_THRESHOLD}; mito PM2; ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes="PM2 at Supporting (1 pt) per ClinGen SVI 2024",
        )

    return ACMGRule(
        rule_id="PM2_MITO",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[
            f"gnomAD v3.1 mito AF={af:.4f} does not meet PM2_MITO threshold"
        ],
        citations=citations,
        applies=False,
    )


def rule_mito_heteroplasmy_level(
    variant: VariantInput,
) -> ACMGRule:
    """Map heteroplasmy level to a clinical significance category.

    Heteroplasmy level (proportion of mitochondrial reads carrying the
    alternate allele) provides evidence about the clinical significance of a
    mitochondrial variant.  Pathogenic mito variants often present at high
    heteroplasmy; very low heteroplasmy may represent an artefact.

    ACGS 2024 §6 heteroplasmy categories:
        > 60%  — High heteroplasmy: variant likely to cause cellular dysfunction
        20–60% — Moderate heteroplasmy: clinically relevant
        1–20%  — Low heteroplasmy: uncertain; may not cause disease
        < 1%   — Very low: likely sequencing artefact; de-prioritise

    Args:
        variant: Annotated mitochondrial variant; ``heteroplasmy_level``
            is the key field (proportion 0–1).

    Returns:
        ACMGRule with rule_id ``"MITO_HETEROPLASMY"`` encoding the
        heteroplasmy category as evidence.  This rule provides contextual
        evidence rather than a point contribution; ``applies`` is set to
        True when heteroplasmy is clinically relevant (>1%).

    References:
        ACGS 2024 v1.2 §6 — heteroplasmy level mapping.
        Richards et al. 2015 PMID:25741868 (mito caveats).
    """
    citations = [
        "ACGS 2024 v1.2 §6",
        "Richards et al. 2015 PMID:25741868",
    ]
    if not variant.is_mito:
        return ACMGRule(
            rule_id="MITO_HETEROPLASMY",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Not a mitochondrial variant"],
            citations=citations,
            applies=False,
        )

    level = variant.heteroplasmy_level
    if level is None:
        return ACMGRule(
            rule_id="MITO_HETEROPLASMY",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=["Heteroplasmy level not available"],
            citations=citations,
            applies=False,
            notes="Heteroplasmy data required for full mito classification (ACGS 2024 §6)",
        )

    if level > _HETEROPLASMY_HIGH:        # >60%; ACGS 2024 §6
        return ACMGRule(
            rule_id="MITO_HETEROPLASMY",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"High heteroplasmy: {level:.1%} > {_HETEROPLASMY_HIGH:.0%} — "
                "variant likely to cause cellular dysfunction (ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes="High heteroplasmy supports pathogenic classification",
        )

    if level >= _HETEROPLASMY_MODERATE:   # 20-60%; ACGS 2024 §6
        return ACMGRule(
            rule_id="MITO_HETEROPLASMY",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"Moderate heteroplasmy: {level:.1%} ({_HETEROPLASMY_MODERATE:.0%}–"
                f"{_HETEROPLASMY_HIGH:.0%}) — clinically relevant (ACGS 2024 §6)"
            ],
            citations=citations,
            applies=True,
            notes="Moderate heteroplasmy: clinically relevant but incomplete penetrance possible",
        )

    if level >= _HETEROPLASMY_LOW:        # 1-20%; ACGS 2024 §6
        return ACMGRule(
            rule_id="MITO_HETEROPLASMY",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                f"Low heteroplasmy: {level:.1%} ({_HETEROPLASMY_LOW:.0%}–"
                f"{_HETEROPLASMY_MODERATE:.0%}) — uncertain clinical significance (ACGS 2024 §6)"
            ],
            citations=citations,
            applies=False,
            notes="Low heteroplasmy: may not cause disease; interpret with caution",
        )

    # < 1% — likely artefact
    return ACMGRule(
        rule_id="MITO_HETEROPLASMY",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[
            f"Very low heteroplasmy: {level:.1%} < {_HETEROPLASMY_LOW:.0%} — "
            "likely sequencing artefact; de-prioritise (ACGS 2024 §6)"
        ],
        citations=citations,
        applies=False,
        notes="Very low heteroplasmy (<1%): likely artefact; recommend repeat sequencing",
    )
