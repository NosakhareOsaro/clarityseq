"""
bayesacmg.combinations
======================

Novel ClinGen SVI 2024 evidence combinations extending Richards 2015.

Critical combination (addresses PM2 downgrade impact):
    PVS1 (Very Strong, 8 pts) + PM2_Supporting (1 pt) = 9 pts
    → Likely Pathogenic (threshold = 6 pts for LP; Bayesian PP ≥ 0.90)

Without this combination, novel LoF variants in LoF-mechanism genes
where the only secondary evidence is rarity would no longer reach LP
after the PM2 downgrade from Moderate (2 pts) to Supporting (1 pt).
Previously (Richards 2015 + PM2 at Moderate):
    PVS1 (8 pts) + PM2_Moderate (2 pts) = 10 pts → LP
After ClinGen SVI 2024 PM2 downgrade without this combination rule:
    PVS1 (8 pts) + PM2_Supporting (1 pt) = 9 pts
The combination PVS1+PM2_Supporting=LP explicitly preserves LP
classification for such variants because 9 pts > 6 pts (LP threshold).

This combination is documented here because:
    1. It is counter-intuitive (a 9-point variant reaches LP, not P).
    2. The PM2 downgrade is the most significant change in ClinGen SVI 2024.
    3. The classification is still correct — 9 pts is firmly LP territory.

Reference: ClinGen SVI Working Group 2024 recommendations
https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
See also: ACGS 2024 v1.2 §5 Table 2 (UK implementation)
          Tavtigian et al. 2020 PMID:32645316 (point-score thresholds)

Guidelines:
    Richards et al. 2015 PMID:25741868 — original classification categories.
    Tavtigian et al. 2020 PMID:32645316 — point thresholds.
    ACGS 2024 v1.2 §5 Table 2 — UK implementation.
    ClinGen SVI Working Group 2024.
"""

from __future__ import annotations

from dataclasses import dataclass

from bayesacmg.models import ACMGRule, ClassificationResult, EvidenceStrength, VariantInput

# ---------------------------------------------------------------------------
# Point thresholds — Tavtigian et al. 2020 PMID:32645316
# ---------------------------------------------------------------------------

_PATHOGENIC_THRESHOLD = 10        # ≥10 pts → Pathogenic
_LIKELY_PATHOGENIC_THRESHOLD = 6  # ≥6 pts → Likely Pathogenic
_LIKELY_BENIGN_THRESHOLD = -6     # ≤-6 pts → Likely Benign
_BENIGN_THRESHOLD = -10           # ≤-10 pts → Benign


@dataclass
class CombinationResult:
    """Result of evaluating a novel ClinGen SVI 2024 combination.

    Attributes:
        combination_name: Short name, e.g. ``"PVS1+PM2_Supporting=LP"``.
        applies: True if this combination applies to the variant.
        resulting_classification: The classification mandated by the combo.
        total_points: Sum of points from contributing rules.
        rules_contributing: List of rule IDs that make up the combination.
        explanation: Human-readable explanation of why the combination applies.
        reference: Literature / guideline reference for the combination.
    """

    combination_name: str
    applies: bool
    resulting_classification: str
    total_points: int
    rules_contributing: list[str]
    explanation: str
    reference: str


def evaluate_pvs1_pm2_supporting(
    rules: list[ACMGRule],
) -> CombinationResult:
    """Evaluate the PVS1 + PM2_Supporting = LP novel combination.

    This is the critical ClinGen SVI 2024 combination that preserves Likely
    Pathogenic classification for LoF variants after the PM2 downgrade from
    Moderate (2 pts) to Supporting (1 pt).

    Point arithmetic:
        PVS1 (Very Strong): +8 pts
        PM2_Supporting:     +1 pt
        Total:              +9 pts
        LP threshold:       ≥6 pts → Likely Pathogenic
        → Classification: Likely Pathogenic (9 pts > 6 pts threshold)

    Why this is explicitly documented:
        In Richards 2015 + PM2 at Moderate: 8 + 2 = 10 pts → LP (borderline P).
        After ClinGen SVI 2024 PM2 downgrade: 8 + 1 = 9 pts → LP (still LP).
        The combination is preserved but at lower total; both reach LP.
        For variants with ONLY PVS1+PM2 as evidence, the 9 pts still comfortably
        exceeds the LP threshold of 6 pts.

    Args:
        rules: List of all ACMGRule instances evaluated for the variant.

    Returns:
        CombinationResult indicating whether this novel combination applies.

    References:
        ClinGen SVI Working Group 2024 recommendations.
        ACGS 2024 v1.2 §5 Table 2.
        Tavtigian et al. 2020 PMID:32645316 — point thresholds.
    """
    pvs1_rules = [r for r in rules if r.rule_id == "PVS1" and r.applies]
    pm2_rules = [r for r in rules
                 if r.rule_id in {"PM2", "PM2_MITO"}
                 and r.applies
                 and r.strength == EvidenceStrength.SUPPORTING]

    if not pvs1_rules or not pm2_rules:
        return CombinationResult(
            combination_name="PVS1+PM2_Supporting=LP",
            applies=False,
            resulting_classification="",
            total_points=0,
            rules_contributing=[],
            explanation="PVS1 and/or PM2_Supporting are not both present",
            reference="ClinGen SVI Working Group 2024",
        )

    pvs1 = pvs1_rules[0]
    pm2 = pm2_rules[0]
    total = pvs1.points + pm2.points  # 8 + 1 = 9 pts

    explanation = (
        f"PVS1 ({pvs1.strength.value}, {pvs1.points} pts) + "
        f"PM2_Supporting ({pm2.strength.value}, {pm2.points} pt) = "
        f"{total} pts → Likely Pathogenic (LP threshold ≥{_LIKELY_PATHOGENIC_THRESHOLD} pts). "
        "ClinGen SVI 2024 explicitly preserves LP classification for LoF variants "
        "where PM2 is the only secondary evidence, after PM2 downgrade from Moderate "
        "(2 pts) to Supporting (1 pt). 9 pts comfortably exceeds the LP threshold of 6 pts."
    )

    return CombinationResult(
        combination_name="PVS1+PM2_Supporting=LP",
        applies=True,
        resulting_classification="Likely_Pathogenic",
        total_points=total,
        rules_contributing=[pvs1.rule_id, pm2.rule_id],
        explanation=explanation,
        reference=(
            "ClinGen SVI Working Group 2024; "
            "ACGS 2024 v1.2 §5 Table 2; "
            "Tavtigian et al. 2020 PMID:32645316"
        ),
    )


def classify_by_points(
    total_points: int,
    stand_alone_benign: bool = False,
) -> str:
    """Convert a total point score to an ACMG/AMP classification category.

    Uses Tavtigian et al. 2020 PMID:32645316 point thresholds with
    ClinGen SVI 2024 guidelines.

    Classification thresholds (Tavtigian 2020 PMID:32645316):
        ≥ 10 pts  → Pathogenic
        6–9 pts   → Likely Pathogenic
        0–5 pts   → VUS (or LB/LB boundary; see note)
        -1 to -5  → VUS (approaching benign)
        -6 to -9  → Likely Benign
        ≤ -10     → Benign

    Args:
        total_points: Sum of all applied rule point values.
        stand_alone_benign: True if BA1 fired → direct Benign regardless of pts.

    Returns:
        Classification string: one of ``"Pathogenic"``, ``"Likely_Pathogenic"``,
        ``"VUS"``, ``"Likely_Benign"``, ``"Benign"``.

    References:
        Tavtigian et al. 2020 PMID:32645316 — point thresholds.
        Richards et al. 2015 PMID:25741868 — classification categories.
    """
    if stand_alone_benign:
        return "Benign"                                   # BA1: direct Benign; Richards 2015

    if total_points >= _PATHOGENIC_THRESHOLD:             # ≥10; Tavtigian 2020 PMID:32645316
        return "Pathogenic"
    if total_points >= _LIKELY_PATHOGENIC_THRESHOLD:      # ≥6; Tavtigian 2020 PMID:32645316
        return "Likely_Pathogenic"
    if total_points <= _BENIGN_THRESHOLD:                 # ≤-10; Tavtigian 2020 PMID:32645316
        return "Benign"
    if total_points <= _LIKELY_BENIGN_THRESHOLD:          # ≤-6; Tavtigian 2020 PMID:32645316
        return "Likely_Benign"
    return "VUS"


def evaluate_all_combinations(
    rules: list[ACMGRule],
) -> list[CombinationResult]:
    """Evaluate all ClinGen SVI 2024 novel combinations for a variant.

    Iterates through all documented novel combinations and returns those
    that apply to the provided set of ACMG rules.

    Args:
        rules: List of all ACMGRule instances evaluated for the variant.

    Returns:
        List of CombinationResult instances for all combinations that apply.
        Empty list if no novel combinations apply.

    References:
        ClinGen SVI Working Group 2024.
        ACGS 2024 v1.2 §5 Table 2.
    """
    results: list[CombinationResult] = []

    pvs1_pm2 = evaluate_pvs1_pm2_supporting(rules)
    if pvs1_pm2.applies:
        results.append(pvs1_pm2)

    return results


def classify_variant(
    rules: list[ACMGRule],
    variant: VariantInput,
) -> ClassificationResult:
    """Perform full classification of a variant using all applied ACMG rules.

    Checks for BA1 stand-alone, evaluates novel combinations, then falls
    back to standard point-sum classification.

    Args:
        rules: List of all ACMGRule instances evaluated for the variant.
        variant: Original variant input.

    Returns:
        ClassificationResult with the final classification, point totals,
        and any novel combination that applied.

    References:
        Richards et al. 2015 PMID:25741868.
        Tavtigian et al. 2020 PMID:32645316.
        ClinGen SVI Working Group 2024.
    """
    applied = [r for r in rules if r.applies]
    not_applied = [r for r in rules if not r.applies]

    # Check for stand-alone Benign (BA1 or mito haplogroup)
    stand_alone = any(
        r.strength == EvidenceStrength.STAND_ALONE and r.applies for r in rules
    )

    total_points = sum(r.points for r in applied)
    novel_combo: str | None = None

    if stand_alone:
        classification = "Benign"
    else:
        combos = evaluate_all_combinations(applied)
        if combos:
            novel_combo = combos[0].combination_name

        classification = classify_by_points(total_points, stand_alone_benign=stand_alone)

    return ClassificationResult(
        variant=variant,
        classification=classification,
        total_points=total_points,
        rules_applied=applied,
        rules_not_applied=not_applied,
        stand_alone_benign=stand_alone,
        novel_combination=novel_combo,
    )
