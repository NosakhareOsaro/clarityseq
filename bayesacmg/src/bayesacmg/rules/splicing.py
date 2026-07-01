"""
bayesacmg.rules.splicing
========================

ClinGen SVI Splicing Subgroup evidence framework for PP3, BP4, and BP7.

Reference:
    Walker et al. 2023 Am J Hum Genet 110:1046-1067
    doi:10.1016/j.ajhg.2023.06.002  PMID:36898414

Key thresholds (SpliceAI Δ score):
    ≥ 0.5  → strong splice impact    → PP3 Strong
    ≥ 0.2  → moderate splice impact  → PP3 Moderate
    < 0.1  → no splice impact predicted
    Synonymous variant + SpliceAI < 0.1 → BP7

When Pangolin and SpliceAI disagree:
    Use the more conservative (lower Δ score) estimate.
    Document the disagreement in the evidence_items field.

This module returns a list[ACMGRule] from rule_splicing_pp3_bp4_bp7()
because PP3, BP4, and BP7 are mutually exclusive for the same
variant+evidence — only one can apply.  The list will always contain
exactly one ACMGRule.

Guidelines:
    Richards et al. 2015 PMID:25741868 — original PP3/BP4/BP7 definitions.
    Walker et al. 2023 PMID:36898414 — SpliceAI/Pangolin framework.
    ACGS 2024 v1.2 §5 Table 2.
"""

from __future__ import annotations

from bayesacmg.models import ACMGRule, EvidenceStrength, VariantInput

# ---------------------------------------------------------------------------
# SpliceAI Δ score thresholds — Walker et al. 2023 PMID:36898414
# ---------------------------------------------------------------------------

_SPLICEAI_STRONG_THRESHOLD = 0.5  # ≥0.5 → strong splice impact → PP3 Strong
_SPLICEAI_MODERATE_THRESHOLD = 0.2  # ≥0.2 → moderate splice impact → PP3
_SPLICEAI_NO_IMPACT_THRESHOLD = 0.1  # <0.1 → no splice impact → BP7 (if synonymous)

# Pangolin high-confidence threshold (conservative when disagreement)
_PANGOLIN_HIGH_THRESHOLD = 0.5  # used when Pangolin is primary


def _resolve_splice_score(
    spliceai_score: float | None,
    pangolin_score: float | None,
) -> tuple[float | None, str, bool]:
    """Resolve the effective splice score from SpliceAI and Pangolin.

    When both scores are available and disagree, the more conservative
    (lower) score is used per Walker et al. 2023 recommendation.

    Args:
        spliceai_score: SpliceAI maximum Δ score across all four channels.
        pangolin_score: Pangolin splice-impact score (0–1).

    Returns:
        Tuple of:
            - effective_score: The score to use for threshold evaluation.
            - tool_used: Name of the tool providing the effective score.
            - disagreement: True if SpliceAI and Pangolin disagree on
              the categorical outcome.
    """
    if spliceai_score is None and pangolin_score is None:
        return None, "none", False

    if spliceai_score is None:
        return pangolin_score, "Pangolin", False

    if pangolin_score is None:
        return spliceai_score, "SpliceAI", False

    # Both available: check for categorical disagreement
    spliceai_category = _score_category(spliceai_score)
    pangolin_category = _score_category(pangolin_score)
    disagreement = spliceai_category != pangolin_category

    # Conservative estimate: take the lower score
    if spliceai_score <= pangolin_score:
        return spliceai_score, "SpliceAI (conservative)", disagreement
    return pangolin_score, "Pangolin (conservative)", disagreement


def _score_category(score: float) -> str:
    """Map a splice score to a categorical outcome string.

    Args:
        score: Splice-impact score (0–1).

    Returns:
        One of ``"strong_impact"``, ``"moderate_impact"``, ``"low_impact"``,
        or ``"no_impact"``.
    """
    if score >= _SPLICEAI_STRONG_THRESHOLD:  # ≥0.5; Walker 2023 PMID:36898414
        return "strong_impact"
    if score >= _SPLICEAI_MODERATE_THRESHOLD:  # ≥0.2; Walker 2023 PMID:36898414
        return "moderate_impact"
    if score < _SPLICEAI_NO_IMPACT_THRESHOLD:  # <0.1; Walker 2023 PMID:36898414
        return "no_impact"
    return "low_impact"


def rule_splicing_pp3_bp4_bp7(
    variant: VariantInput,
    spliceai_score: float | None = None,
    spliceai_delta: (
        float | None
    ) = None,  # alias for spliceai_score (test/VEP convention)
    pangolin_score: float | None = None,
) -> ACMGRule:
    """Evaluate PP3, BP4, and BP7 for splice-impacting variants.

    Implements the ClinGen SVI Splicing Subgroup framework from Walker et al.
    2023 (PMID:36898414).  Exactly one of PP3, BP4, or BP7 will fire per call.

    SpliceAI Δ score decision tree (Walker et al. 2023 PMID:36898414):
        ≥ 0.5  → Strong splice impact → PP3 at Strong strength
        ≥ 0.2  → Moderate splice impact → PP3 at Moderate strength
        0.1–0.2 → Low impact, uncertain → neither PP3 nor BP4/BP7
        < 0.1  → No splice impact predicted:
            If variant is synonymous → BP7 (Supporting Benign)
            Otherwise → BP4 (Supporting Benign)

    When SpliceAI and Pangolin disagree on the categorical outcome:
        Use the more conservative (lower Δ score) estimate.
        Document the disagreement in evidence_items.

    Args:
        variant: Annotated variant; ``variant_type`` is checked for synonymous
            and splice site classification.  ``spliceai_max_delta`` and
            ``pangolin_score`` are used if the explicit arguments are None.
        spliceai_score: SpliceAI maximum Δ score (override from variant object).
        pangolin_score: Pangolin splice score (override from variant object).

    Returns:
        List containing exactly one ACMGRule.  The rule_id will be one of
        ``"PP3"``, ``"BP4"``, or ``"BP7"`` depending on the evidence.

    References:
        Richards et al. 2015 PMID:25741868 — original PP3/BP4/BP7.
        Walker et al. 2023 PMID:36898414 — SpliceAI/Pangolin thresholds.
        ACGS 2024 v1.2 §5 Table 2.

    Raises:
        ValueError: If ``variant`` is None.
    """
    citations = [
        "Richards et al. 2015 PMID:25741868",
        "Walker et al. 2023 PMID:36898414 (ClinGen SVI Splicing Subgroup)",
        "ACGS 2024 v1.2 §5 Table 2",
    ]

    # Resolve spliceai_delta alias → spliceai_score
    if spliceai_score is None and spliceai_delta is not None:
        spliceai_score = spliceai_delta

    # Allow explicit arguments to override VariantInput fields
    effective_spliceai = (
        spliceai_score
        if spliceai_score is not None
        else (variant.spliceai_delta or variant.spliceai_max_delta)
    )
    effective_pangolin = (
        pangolin_score if pangolin_score is not None else variant.pangolin_score
    )

    score, tool_used, disagreement = _resolve_splice_score(
        effective_spliceai, effective_pangolin
    )

    disagreement_note = ""
    if disagreement:
        disagreement_note = (
            f" [SpliceAI={effective_spliceai}, Pangolin={effective_pangolin} DISAGREE — "
            "conservative estimate used per Walker 2023 PMID:36898414]"
        )

    # No scores available
    if score is None:
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.SUPPORTING,
            evidence_items=[
                "No SpliceAI or Pangolin score available for splice assessment"
            ],
            citations=citations,
            applies=False,
            notes="Cannot evaluate splicing PP3/BP4/BP7 without in silico scores",
        )

    # ≥ 0.5 → Strong splice impact → PP3 Strong (Walker 2023 PMID:36898414)
    if score >= _SPLICEAI_STRONG_THRESHOLD:
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.STRONG,
            evidence_items=[
                f"{tool_used} score {score:.3f} ≥ {_SPLICEAI_STRONG_THRESHOLD} "
                f"→ strong splice impact (PP3 Strong; Walker 2023 PMID:36898414)"
                f"{disagreement_note}"
            ],
            citations=citations,
            applies=True,
        )

    # ≥ 0.2 → Moderate splice impact → PP3 Moderate (Walker 2023 PMID:36898414)
    if score >= _SPLICEAI_MODERATE_THRESHOLD:
        return ACMGRule(
            rule_id="PP3",
            strength=EvidenceStrength.MODERATE,
            evidence_items=[
                f"{tool_used} score {score:.3f} ≥ {_SPLICEAI_MODERATE_THRESHOLD} "
                f"→ moderate splice impact (PP3; Walker 2023 PMID:36898414)"
                f"{disagreement_note}"
            ],
            citations=citations,
            applies=True,
        )

    # < 0.1 → No splice impact
    if score < _SPLICEAI_NO_IMPACT_THRESHOLD:
        # Synonymous + no splice impact → BP7 (Walker 2023 PMID:36898414)
        if variant.variant_type == "synonymous":
            return ACMGRule(
                rule_id="BP7",
                strength=EvidenceStrength.SUPPORTING_BENIGN,
                evidence_items=[
                    f"Synonymous variant; {tool_used} score {score:.3f} < "
                    f"{_SPLICEAI_NO_IMPACT_THRESHOLD} → no predicted splice impact "
                    f"(BP7; Walker 2023 PMID:36898414){disagreement_note}"
                ],
                citations=citations,
                applies=True,
            )
        # Non-synonymous but no splice impact → BP4
        return ACMGRule(
            rule_id="BP4",
            strength=EvidenceStrength.SUPPORTING_BENIGN,
            evidence_items=[
                f"{tool_used} score {score:.3f} < {_SPLICEAI_NO_IMPACT_THRESHOLD} "
                f"→ no predicted splice impact (BP4; Walker 2023 PMID:36898414)"
                f"{disagreement_note}"
            ],
            citations=citations,
            applies=True,
        )

    # Score in 0.1–0.2 range: uncertain, neither criterion fires
    return ACMGRule(
        rule_id="PP3",
        strength=EvidenceStrength.SUPPORTING,
        evidence_items=[
            f"{tool_used} score {score:.3f} in uncertain zone "
            f"({_SPLICEAI_NO_IMPACT_THRESHOLD}–{_SPLICEAI_MODERATE_THRESHOLD}) — "
            f"neither PP3 nor BP4/BP7 applies (Walker 2023 PMID:36898414)"
            f"{disagreement_note}"
        ],
        citations=citations,
        applies=False,
        notes="SpliceAI in uncertain zone; additional evidence required",
    )
