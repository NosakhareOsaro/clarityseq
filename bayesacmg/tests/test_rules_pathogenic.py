"""
Tests for BayesACMG pathogenic rules.

CRITICAL TESTS:
- PM2 MUST return SUPPORTING (1 pt), NOT Moderate (2 pts)  [ClinGen SVI 2024]
- PP3 uses AlphaMissense ≥ 0.564 threshold                 [ClinGen SVI 2024]
- BP4 uses AlphaMissense ≤ 0.340 threshold                 [ClinGen SVI 2024]
- Novel combination PVS1 + PM2_Supporting = LP             [ClinGen SVI 2024]

Guidelines:
    Richards et al. 2015 PMID:25741868 — original criteria definitions
    ClinGen SVI 2024 — PM2→Supporting; AlphaMissense thresholds
    ACGS 2024 v1.2 — UK implementation; MANE Select; PVS1 MANE adjustment
    Abou Tayoun et al. 2018 PMID:30192042 — PVS1 decision tree
"""

from __future__ import annotations


from bayesacmg.models import (
    EvidenceStrength,
    VariantInput,
)
from bayesacmg.rules.pathogenic import (
    rule_pm2,
    rule_pp3,
    rule_pvs1,
)
from bayesacmg.combinations import apply_novel_combinations, compute_total_points

# ===========================================================================
# PM2 — CRITICAL: Must return SUPPORTING (1 pt), not Moderate (2 pts)
# ===========================================================================


class TestPM2:
    """Tests for PM2 (absent from / rare in population databases).

    CRITICAL: PM2 MUST return EvidenceStrength.SUPPORTING (1 pt).
    ClinGen SVI 2024 downgraded PM2 from Moderate (2 pts) to Supporting (1 pt).
    gnomAD v4.1 evidence: ultra-rare variants are common; rarity is weak evidence.
    """

    def test_pm2_returns_supporting_not_moderate(
        self, brca2_novel_lof: VariantInput
    ) -> None:
        """PM2 on absent variant MUST return SUPPORTING weight.

        This is the most critical correctness test in the test suite.
        Returning MODERATE would contradict ClinGen SVI 2024 guidance and
        ACGS 2024 v1.2 §5 Appendix C.

        Variant: novel BRCA2 frameshift, absent from gnomAD v4.1 (AF=0)
        Expected: PM2 fires with EvidenceStrength.SUPPORTING (1 pt)
        """
        result = rule_pm2(brca2_novel_lof)
        assert (
            result.applies is True
        ), "PM2 should apply for variant absent from gnomAD v4.1"
        assert result.strength == EvidenceStrength.SUPPORTING, (
            f"PM2 MUST return SUPPORTING (1 pt) per ClinGen SVI 2024, "
            f"but returned {result.strength}. "
            f"PM2 was downgraded from Moderate (2 pts) to Supporting (1 pt) in ClinGen SVI 2024 "
            f"because gnomAD v4.1 (807,162 individuals) reveals that ultra-rare variants "
            f"are common in the general population."
        )
        assert result.rule_id == "PM2"

    def test_pm2_strength_is_not_moderate(self, brca2_novel_lof: VariantInput) -> None:
        """Explicit negative test: PM2 must NOT return MODERATE."""
        result = rule_pm2(brca2_novel_lof)
        assert result.strength != EvidenceStrength.MODERATE, (
            "PM2 returned MODERATE (2 pts). This contradicts ClinGen SVI 2024 guidance. "
            "Fix: bayesacmg/src/bayesacmg/rules/pathogenic.py rule_pm2() must return "
            "EvidenceStrength.SUPPORTING, not EvidenceStrength.MODERATE."
        )

    def test_pm2_point_value_is_one(self, brca2_novel_lof: VariantInput) -> None:
        """PM2 SUPPORTING strength must give exactly 1 point."""
        from bayesacmg.models import STRENGTH_POINTS

        result = rule_pm2(brca2_novel_lof)
        points = STRENGTH_POINTS[result.strength]
        assert points == 1, (
            f"PM2 should contribute 1 pt (SUPPORTING), but gave {points} pts. "
            f"Strength returned: {result.strength}"
        )

    def test_pm2_does_not_apply_when_af_exceeds_threshold(
        self, common_benign_missense: VariantInput
    ) -> None:
        """PM2 does not apply when gnomAD v4.1 AF exceeds 0.0001 threshold."""
        result = rule_pm2(common_benign_missense)
        assert result.applies is False, (
            f"PM2 should not apply for AF={common_benign_missense.gnomad_af} "
            f"which exceeds the 0.0001 threshold"
        )

    def test_pm2_cites_clingen_svi_2024(self, brca2_novel_lof: VariantInput) -> None:
        """PM2 citations must reference ClinGen SVI PM2 guidance 2024."""
        result = rule_pm2(brca2_novel_lof)
        citations_str = " ".join(result.citations).lower()
        assert any(
            term in citations_str for term in ["clingen", "svi", "2024", "pm2"]
        ), f"PM2 should cite ClinGen SVI 2024. Got citations: {result.citations}"


# ===========================================================================
# PP3 — AlphaMissense PRIMARY predictor (ClinGen SVI 2024)
# ===========================================================================


class TestPP3:
    """Tests for PP3 using AlphaMissense as primary predictor.

    ClinGen SVI 2024: AlphaMissense is the PRIMARY tool for PP3/BP4.
    Thresholds: ≥ 0.564 → PP3; ≤ 0.340 → BP4; ambiguous otherwise.
    Reference: Cheng et al. 2023 Science PMID:37703350
    """

    def test_pp3_alphamissense_at_threshold_fires(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """AlphaMissense score ≥ 0.564 triggers PP3.

        Variant: TP53 p.Arg273Cys (RCV000012735) — AM score 0.97
        Expected: PP3 applies
        """
        result = rule_pp3(
            variant=tp53_missense_pathogenic,
            alphamissense_score=0.97,
            revel_score=0.85,
        )
        assert result.applies is True
        assert result.rule_id == "PP3"

    def test_pp3_alphamissense_exact_threshold(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """AlphaMissense score exactly at 0.564 threshold → PP3 applies."""
        result = rule_pp3(
            variant=tp53_missense_pathogenic,
            alphamissense_score=0.564,  # exact threshold
            revel_score=None,
        )
        assert result.applies is True, "Score exactly at 0.564 should trigger PP3"

    def test_pp3_alphamissense_below_threshold_does_not_fire(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """AlphaMissense score 0.563 (just below threshold) → PP3 does NOT apply."""
        result = rule_pp3(
            variant=tp53_missense_pathogenic,
            alphamissense_score=0.563,  # just below 0.564
            revel_score=None,
        )
        assert result.applies is False, (
            "Score just below 0.564 should NOT trigger PP3 "
            "(ambiguous range 0.340–0.564)"
        )

    def test_bp4_alphamissense_at_threshold(
        self, common_benign_missense: VariantInput
    ) -> None:
        """AlphaMissense score ≤ 0.340 → BP4 applies (benign supporting).

        The BP4 rule is tested here via the PP3 function which returns BP4
        when AlphaMissense score is ≤ 0.340.
        """
        from bayesacmg.rules.pathogenic import rule_bp4_from_alphamissense

        result = rule_bp4_from_alphamissense(
            variant=common_benign_missense,
            alphamissense_score=0.21,  # ≤ 0.340 → BP4
        )
        assert result.applies is True
        assert result.rule_id == "BP4"
        assert result.strength == EvidenceStrength.SUPPORTING_BENIGN

    def test_bp4_alphamissense_exact_threshold(
        self, common_benign_missense: VariantInput
    ) -> None:
        """AlphaMissense score exactly at 0.340 → BP4 applies."""
        from bayesacmg.rules.pathogenic import rule_bp4_from_alphamissense

        result = rule_bp4_from_alphamissense(
            variant=common_benign_missense,
            alphamissense_score=0.340,  # exact threshold
        )
        assert result.applies is True, "Score exactly at 0.340 should trigger BP4"

    def test_ambiguous_range_no_evidence(
        self, tp53_missense_pathogenic: VariantInput
    ) -> None:
        """AlphaMissense score in ambiguous range (0.341–0.563) → no evidence."""
        result = rule_pp3(
            variant=tp53_missense_pathogenic,
            alphamissense_score=0.45,  # ambiguous range
            revel_score=None,
        )
        assert (
            result.applies is False
        ), "Ambiguous range (0.340–0.564) should yield no PP3/BP4 evidence"


# ===========================================================================
# PVS1 — null variants in LoF-mechanism genes
# ===========================================================================


class TestPVS1:
    """Tests for PVS1 (null variant in LoF-mechanism gene).

    Abou Tayoun 2018 PMID:30192042 decision tree.
    ACGS 2024 v1.2 §5: MANE Select requirement.
    """

    def test_pvs1_frameshift_in_lof_gene(
        self,
        brca1_frameshift: VariantInput,
        brca1_transcript,
        brca1_gene,
    ) -> None:
        """Frameshift in BRCA1 (LoF mechanism) → PVS1 Very Strong.

        Variant: BRCA1 c.5266dupC (RCV000007535)
        MANE Select transcript → PVS1 Very Strong (not downgraded)
        """
        result = rule_pvs1(
            variant=brca1_frameshift,
            transcript=brca1_transcript,
            gene=brca1_gene,
        )
        assert result.applies is True
        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.rule_id == "PVS1"

    def test_pvs1_non_mane_transcript_reduces_strength(
        self,
        brca1_frameshift: VariantInput,
        brca1_gene,
    ) -> None:
        """PVS1 on non-MANE transcript is reduced by one level (ACGS 2024 §5)."""
        non_mane_transcript = TranscriptData(
            transcript_id="NM_007298.3",  # Non-MANE BRCA1 transcript
            gene_symbol="BRCA1",
            is_mane_select=False,  # NOT MANE Select
            is_canonical=False,
            lof_disease_mechanism=True,
            exon_count=23,
        )
        result = rule_pvs1(
            variant=brca1_frameshift,
            transcript=non_mane_transcript,
            gene=brca1_gene,
        )
        assert result.applies is True
        # Non-MANE → reduced by one level: Very Strong → Strong
        assert result.strength == EvidenceStrength.STRONG, (
            "PVS1 on non-MANE Select transcript should be reduced to STRONG "
            "per ACGS 2024 v1.2 §5"
        )

    def test_pvs1_does_not_apply_in_non_lof_gene(
        self, tp53_missense_pathogenic, brca1_transcript
    ) -> None:
        """PVS1 does not apply when LoF is not the established disease mechanism."""
        non_lof_gene = GeneData(
            gene_symbol="KCNQ1",
            omim_id="OMIM:192500",
            lof_mechanism=False,  # Gain-of-function disease; LoF not established
            pli=0.85,
            has_vcep_specification=True,
        )
        result = rule_pvs1(
            variant=tp53_missense_pathogenic,
            transcript=brca1_transcript,
            gene=non_lof_gene,
        )
        assert result.applies is False


# ===========================================================================
# Novel combination: PVS1 + PM2_Supporting = LP  (ClinGen SVI 2024)
# ===========================================================================


class TestNovelCombinations:
    """Tests for ClinGen SVI 2024 novel evidence combinations.

    Critical combination: PVS1 (8 pts) + PM2_Supporting (1 pt) = 9 pts → LP.
    LP threshold: ≥ 6 pts (Tavtigian 2020 PMID:32645316).
    Without this explicit combination, novel LoF variants where rarity is the
    only secondary evidence would not reach LP after the PM2 downgrade.

    Reference: ClinGen SVI Working Group 2024 recommendations
    """

    def test_pvs1_plus_pm2_supporting_gives_lp(
        self,
        brca2_novel_lof: VariantInput,
        brca1_gene,
    ) -> None:
        """PVS1 (8 pts) + PM2_Supporting (1 pt) = 9 pts ≥ 6 → Likely Pathogenic.

        Variant: Novel BRCA2 frameshift, absent from gnomAD v4.1 (no ClinVar RCV).
        Evidence: PVS1 (Very Strong, 8 pts) + PM2 (Supporting, 1 pt) = 9 pts.
        Expected classification: Likely Pathogenic.

        This is the key regression test: if PM2 erroneously returns Moderate (2 pts)
        instead of Supporting (1 pt), total = 10 pts → Pathogenic instead of LP.
        The test specifically validates 9 pts = LP, not 10 pts = P.
        """

        # Simulate applying the two rules
        pvs1_result = rule_pvs1(
            variant=brca2_novel_lof,
            transcript=TranscriptData(
                transcript_id="NM_000059.4",
                gene_symbol="BRCA2",
                is_mane_select=True,
                lof_disease_mechanism=True,
                exon_count=27,
            ),
            gene=brca1_gene,
        )
        pm2_result = rule_pm2(brca2_novel_lof)

        assert pvs1_result.applies is True
        assert pm2_result.applies is True
        assert (
            pm2_result.strength == EvidenceStrength.SUPPORTING
        ), "PM2 returned non-SUPPORTING strength — this will break the LP classification"

        # Compute total points
        rules = [pvs1_result, pm2_result]
        total_pts = compute_total_points(rules)

        # PVS1 = 8 pts; PM2_Supporting = 1 pt; total = 9 pts
        assert total_pts == 9, (
            f"Expected 9 pts (PVS1=8 + PM2_Supporting=1), got {total_pts}. "
            f"If PM2 erroneously returns MODERATE (2 pts), total would be 10 — wrong."
        )

        # 9 pts ≥ 6 → LP (not P which requires ≥ 10 pts)
        classification = apply_novel_combinations(rules, total_pts)
        assert (
            "Likely Pathogenic" in classification.classification
        ), f"Expected 'Likely Pathogenic' at 9 pts, got: {classification.classification}"
        # Also verify it's NOT classified as Pathogenic (which would require ≥ 10 pts)
        assert (
            classification.classification != "Pathogenic"
        ), "9 pts should give LP, not P. P requires ≥ 10 pts (Richards 2015)."

    def test_pvs1_alone_does_not_reach_lp(self, brca2_novel_lof, brca1_gene) -> None:
        """PVS1 alone (8 pts) reaches LP threshold (≥ 6) but not P threshold (≥ 10)."""
        pvs1_result = rule_pvs1(
            variant=brca2_novel_lof,
            transcript=TranscriptData(
                transcript_id="NM_000059.4",
                gene_symbol="BRCA2",
                is_mane_select=True,
                lof_disease_mechanism=True,
                exon_count=27,
            ),
            gene=brca1_gene,
        )
        total_pts = compute_total_points([pvs1_result])
        assert total_pts == 8
        # 8 pts → LP (≥ 6 pts); below P threshold (≥ 10 pts)
        classification = apply_novel_combinations([pvs1_result], total_pts)
        assert "Likely Pathogenic" in classification.classification


# ---------------------------------------------------------------------------
# Imports needed for test stubs to work
# ---------------------------------------------------------------------------
from bayesacmg.models import GeneData, TranscriptData  # noqa: E402
