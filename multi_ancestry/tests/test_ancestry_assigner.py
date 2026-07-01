"""
multi_ancestry.tests.test_ancestry_assigner
=============================================
pytest tests for population ancestry label assignment.

Tests cover:
    - assign_ancestry: single-ancestry and admixed thresholds.
    - assign_ancestry_batch: batch processing.
    - _recommend_vqsr_training: VQSR training set selection.
"""

from __future__ import annotations

import pytest

from multi_ancestry.somalier_runner import SomalierAncestryResult
from multi_ancestry.ancestry_assigner import (
    AncestryAssignment,
    assign_ancestry,
    assign_ancestry_batch,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_somalier_result(
    sample_id: str = "SAMPLE1",
    predicted_ancestry: str = "EUR",
    predicted_ancestry_p: float = 0.95,
    fractions: dict | None = None,
) -> SomalierAncestryResult:
    """Return a minimal SomalierAncestryResult for testing."""
    if fractions is None:
        fractions = {"EUR": predicted_ancestry_p}
    return SomalierAncestryResult(
        sample_id=sample_id,
        predicted_ancestry=predicted_ancestry,
        predicted_ancestry_p=predicted_ancestry_p,
        ancestry_fractions=fractions,
    )


# ---------------------------------------------------------------------------
# assign_ancestry tests
# ---------------------------------------------------------------------------


class TestAssignAncestry:
    """Tests for assign_ancestry() single-sample assignment."""

    def test_eur_high_fraction_assigned_correctly(self) -> None:
        """EUR sample with fraction 0.95 → EUR label, not admixed."""
        result = _make_somalier_result(
            "S1", "EUR", 0.95, {"EUR": 0.95, "AFR": 0.03, "AMR": 0.02}
        )
        assignment = assign_ancestry(result)
        assert assignment.primary_label == "EUR"
        assert assignment.is_admixed is False

    def test_afr_high_fraction(self) -> None:
        """AFR sample with fraction 0.85 → AFR label."""
        result = _make_somalier_result(
            "S2", "AFR", 0.85, {"AFR": 0.85, "EUR": 0.10, "AMR": 0.05}
        )
        assignment = assign_ancestry(result)
        assert assignment.primary_label == "AFR"
        assert assignment.is_admixed is False

    def test_admixed_sample_below_threshold(self) -> None:
        """Sample with max fraction 0.55 < 0.8 → ADMIXED label."""
        result = _make_somalier_result(
            "S3", "EUR", 0.55, {"EUR": 0.55, "AFR": 0.30, "AMR": 0.15}
        )
        assignment = assign_ancestry(result)
        assert assignment.primary_label == "ADMIXED"
        assert assignment.is_admixed is True

    def test_admixed_component_fractions_populated(self) -> None:
        """Admixed sample has component_fractions with populations ≥10%."""
        result = _make_somalier_result(
            "S4", "EUR", 0.55, {"EUR": 0.55, "AFR": 0.30, "AMR": 0.15}
        )
        assignment = assign_ancestry(result)
        # All fractions ≥ 0.10 should appear in components
        assert "EUR" in assignment.component_fractions
        assert "AFR" in assignment.component_fractions
        assert "AMR" in assignment.component_fractions

    def test_single_ancestry_no_component_fractions(self) -> None:
        """Single-ancestry assignment has no component_fractions."""
        result = _make_somalier_result(
            "S5", "SAS", 0.90, {"SAS": 0.90, "EUR": 0.10}
        )
        assignment = assign_ancestry(result)
        assert not assignment.is_admixed
        assert assignment.component_fractions == {}

    def test_sample_id_preserved(self) -> None:
        """Sample ID is preserved in the assignment."""
        result = _make_somalier_result("MY_SAMPLE", "EUR", 0.90, {"EUR": 0.90})
        assignment = assign_ancestry(result)
        assert assignment.sample_id == "MY_SAMPLE"

    def test_vqsr_training_set_assigned(self) -> None:
        """VQSR training set recommendation is non-empty."""
        result = _make_somalier_result("S6", "EAS", 0.85, {"EAS": 0.85})
        assignment = assign_ancestry(result)
        assert assignment.vqsr_training_set != ""

    def test_full_name_populated(self) -> None:
        """Full population name is set for reporting."""
        result = _make_somalier_result("S7", "EUR", 0.90, {"EUR": 0.90})
        assignment = assign_ancestry(result)
        assert "European" in assignment.full_name

    def test_custom_threshold(self) -> None:
        """Custom threshold 0.6 classifies 0.65 fraction as single ancestry."""
        result = _make_somalier_result(
            "S8", "EUR", 0.65, {"EUR": 0.65, "AFR": 0.35}
        )
        assignment = assign_ancestry(result, threshold=0.6)
        assert assignment.primary_label == "EUR"
        assert assignment.is_admixed is False

    def test_empty_fractions_falls_back_to_predicted(self) -> None:
        """Empty fractions dict falls back to predicted ancestry."""
        result = SomalierAncestryResult(
            sample_id="S9",
            predicted_ancestry="AMR",
            predicted_ancestry_p=0.85,
            ancestry_fractions={},
        )
        assignment = assign_ancestry(result)
        assert assignment.primary_label == "AMR"

    def test_unknown_population_maps_to_admixed(self) -> None:
        """Population label not in known set maps to ADMIXED."""
        result = SomalierAncestryResult(
            sample_id="S10",
            predicted_ancestry="UNKNOWN_POP",
            predicted_ancestry_p=0.90,
            ancestry_fractions={"UNKNOWN_POP": 0.90},
        )
        assignment = assign_ancestry(result)
        assert assignment.primary_label == "ADMIXED"


# ---------------------------------------------------------------------------
# assign_ancestry_batch tests
# ---------------------------------------------------------------------------


class TestAssignAncestryBatch:
    """Tests for assign_ancestry_batch() batch function."""

    def test_returns_list_of_assignments(self) -> None:
        """Batch function returns one assignment per sample."""
        results = [
            _make_somalier_result("S1", "EUR", 0.90, {"EUR": 0.90}),
            _make_somalier_result("S2", "AFR", 0.85, {"AFR": 0.85}),
        ]
        assignments = assign_ancestry_batch(results)
        assert len(assignments) == 2
        assert all(isinstance(a, AncestryAssignment) for a in assignments)

    def test_empty_input_returns_empty(self) -> None:
        """Empty input returns empty list."""
        assert assign_ancestry_batch([]) == []

    def test_admixed_count_correct(self) -> None:
        """Admixed samples are correctly identified in batch."""
        results = [
            _make_somalier_result("S1", "EUR", 0.90, {"EUR": 0.90}),
            _make_somalier_result(
                "S2", "EUR", 0.50, {"EUR": 0.50, "AFR": 0.50}
            ),
        ]
        assignments = assign_ancestry_batch(results)
        admixed = [a for a in assignments if a.is_admixed]
        assert len(admixed) == 1
        assert admixed[0].sample_id == "S2"
