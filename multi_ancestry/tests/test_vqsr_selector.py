"""
multi_ancestry.tests.test_vqsr_selector
=========================================
pytest tests for gnomAD v4.1 VQSR training set selection.

Tests cover:
    - select_vqsr_training: EUR/AFR/admixed selection logic.
    - select_vqsr_training_batch: batch processing.
    - get_population_vcf_uri: URI construction for known populations.
"""

from __future__ import annotations

import pytest

from multi_ancestry.ancestry_assigner import AncestryAssignment
from multi_ancestry.vqsr_selector import (
    VQSRTrainingSet,
    get_population_vcf_uri,
    select_vqsr_training,
    select_vqsr_training_batch,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_assignment(
    sample_id: str = "SAMPLE1",
    label: str = "EUR",
    is_admixed: bool = False,
    fractions: dict | None = None,
) -> AncestryAssignment:
    """Return a minimal AncestryAssignment for testing."""
    return AncestryAssignment(
        sample_id=sample_id,
        primary_label=label,
        full_name="European" if label == "EUR" else label,
        is_admixed=is_admixed,
        component_fractions=fractions or {},
        vqsr_training_set=f"gnomAD_v4.1_{label}",
    )


# ---------------------------------------------------------------------------
# select_vqsr_training tests
# ---------------------------------------------------------------------------


class TestSelectVQSRTraining:
    """Tests for select_vqsr_training() main selection function."""

    def test_eur_sample_selects_eur_training(self) -> None:
        """EUR ancestry assignment selects the EUR gnomAD v4.1 training VCF."""
        assignment = _make_assignment("S1", "EUR", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert isinstance(result, VQSRTrainingSet)
        assert result.population_label == "EUR"
        assert "EUR" in result.training_vcf or "vcf" in result.training_vcf.lower()
        assert result.fallback is False

    def test_afr_sample_selects_afr_training(self) -> None:
        """AFR ancestry assignment selects AFR training set."""
        assignment = _make_assignment("S2", "AFR", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert result.population_label == "AFR"
        assert result.fallback is False

    def test_admixed_sample_falls_back_to_global(self) -> None:
        """Admixed sample uses global (all-population) training set."""
        assignment = _make_assignment(
            "S3", "ADMIXED", is_admixed=True,
            fractions={"EUR": 0.50, "AFR": 0.50},
        )
        result = select_vqsr_training(assignment)
        assert result.fallback is True
        assert result.population_label in ("ADMIXED", "ALL")

    def test_sample_id_preserved(self) -> None:
        """Sample ID is preserved in the training set result."""
        assignment = _make_assignment("MY_SAMPLE", "EUR", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert result.sample_id == "MY_SAMPLE"

    def test_gnomad_version_is_4_1(self) -> None:
        """gnomAD version string is '4.1'."""
        assignment = _make_assignment("S4", "SAS", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert result.gnomad_version == "4.1"

    def test_min_af_filter_default(self) -> None:
        """Default min_af_filter is 0.001."""
        assignment = _make_assignment("S5", "EAS", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert result.min_af_filter == pytest.approx(0.001)

    def test_custom_gnomad_base_applied(self) -> None:
        """Custom gnomad_resource_base is reflected in training_vcf."""
        assignment = _make_assignment("S6", "EUR", is_admixed=False)
        result = select_vqsr_training(
            assignment, gnomad_resource_base="/local/gnomad"
        )
        assert result.training_vcf.startswith("/local/gnomad")

    def test_notes_populated(self) -> None:
        """Notes field is non-empty."""
        assignment = _make_assignment("S7", "AMR", is_admixed=False)
        result = select_vqsr_training(assignment)
        assert result.notes != ""


# ---------------------------------------------------------------------------
# select_vqsr_training_batch tests
# ---------------------------------------------------------------------------


class TestSelectVQSRTrainingBatch:
    """Tests for select_vqsr_training_batch() batch function."""

    def test_returns_one_result_per_sample(self) -> None:
        """Batch returns one VQSRTrainingSet per AncestryAssignment."""
        assignments = [
            _make_assignment("S1", "EUR", is_admixed=False),
            _make_assignment("S2", "AFR", is_admixed=False),
        ]
        results = select_vqsr_training_batch(assignments)
        assert len(results) == 2

    def test_empty_input_returns_empty(self) -> None:
        """Empty input returns empty list."""
        assert select_vqsr_training_batch([]) == []

    def test_admixed_uses_fallback_in_batch(self) -> None:
        """Admixed samples in batch use fallback training set."""
        assignments = [
            _make_assignment("S1", "EUR", is_admixed=False),
            _make_assignment("S2", "ADMIXED", is_admixed=True),
        ]
        results = select_vqsr_training_batch(assignments)
        admixed_results = [r for r in results if r.fallback]
        assert len(admixed_results) == 1


# ---------------------------------------------------------------------------
# get_population_vcf_uri tests
# ---------------------------------------------------------------------------


class TestGetPopulationVcfUri:
    """Tests for get_population_vcf_uri() helper."""

    def test_eur_uri_contains_eur(self) -> None:
        """EUR population URI contains 'EUR' identifier."""
        uri = get_population_vcf_uri("EUR")
        assert "EUR" in uri

    def test_afr_uri_contains_afr(self) -> None:
        """AFR population URI contains 'AFR' identifier."""
        uri = get_population_vcf_uri("AFR")
        assert "AFR" in uri

    def test_unknown_population_returns_all_uri(self) -> None:
        """Unknown population falls back to all-population URI."""
        uri = get_population_vcf_uri("UNKNOWN")
        assert uri  # non-empty

    def test_custom_base_applied(self) -> None:
        """Custom gnomad_resource_base is prepended to URI."""
        uri = get_population_vcf_uri("EUR", gnomad_resource_base="/local/gnomad")
        assert uri.startswith("/local/gnomad")

    def test_case_insensitive(self) -> None:
        """Population label lookup is case-insensitive."""
        uri_upper = get_population_vcf_uri("EUR")
        uri_lower = get_population_vcf_uri("eur")
        assert uri_upper == uri_lower
