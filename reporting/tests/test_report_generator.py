"""
reporting.tests.test_report_generator
=======================================
pytest tests for the NHS GMS clinical report generator.

Tests cover:
    - ReportGenerator initialization with VUS review date scheduling.
    - Novel P/LP ClinVar submission flagging.
    - HTML generation (fallback when Jinja2 unavailable).
    - JSON-LD audit trail generation.
    - generate() creates expected output files.

ACGS 2024 v1.2 compliance checks:
    - VUS review date = report date + 2 years (§9).
    - Novel P/LP → pending_clinvar_submission=True.
    - MANE Select transcript in variant table.

References:
    ACGS 2024 v1.2 §5, §9 (Durkie et al., 20 Feb 2024).
    Morales et al. 2022 PMID:35356062 (MANE Select).
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_config():
    """Return a minimal ReportConfig for testing.

    Returns:
        ReportConfig with test patient and sample IDs.
    """
    from reporting.report_generator import ReportConfig
    return ReportConfig(
        patient_id="TEST-PATIENT-001",
        sample_id="LAB-2024-12345",
        referral_indication="Epileptic encephalopathy (HP:0200134)",
        referring_clinician="Dr Test Clinician",
        pipeline_version="1.0.0",
        report_date="2024-12-13",
    )


@pytest.fixture
def pathogenic_variant():
    """Return a Pathogenic variant entry.

    Returns:
        VariantReportEntry for a known BRCA1 pathogenic variant.
    """
    from reporting.report_generator import VariantReportEntry
    return VariantReportEntry(
        gene_symbol="BRCA1",
        mane_select_transcript="NM_007294.4",
        hgvsc="NM_007294.4:c.5266dupC",
        hgvsp="NP_009225.1:p.Gln1756ProfsTer25",
        acmg_class="Pathogenic",
        gnomad_af=None,
        alphamissense_score=None,
        clinvar_id="RCV000048069",
        clinvar_class="Pathogenic",
        inheritance_mode="AD",
        rules_applied=["PVS1", "PM2"],
        is_novel=False,
    )


@pytest.fixture
def vus_variant():
    """Return a VUS variant entry.

    Returns:
        VariantReportEntry classified as VUS.
    """
    from reporting.report_generator import VariantReportEntry
    return VariantReportEntry(
        gene_symbol="SCN1A",
        mane_select_transcript="NM_006920.6",
        hgvsc="NM_006920.6:c.4000G>A",
        hgvsp="NP_008851.3:p.Val1334Met",
        acmg_class="VUS",
        gnomad_af=0.000005,
        alphamissense_score=0.45,  # intermediate zone
        clinvar_id=None,
        clinvar_class=None,
        inheritance_mode="AD",
        rules_applied=["PM2"],
        is_novel=True,
        posterior_p=0.55,
        hdi_lower=0.42,
        hdi_upper=0.68,
    )


@pytest.fixture
def novel_lp_variant():
    """Return a novel Likely Pathogenic variant for ClinVar flag testing.

    Returns:
        VariantReportEntry: Novel LP, not in ClinVar → should be flagged.
    """
    from reporting.report_generator import VariantReportEntry
    return VariantReportEntry(
        gene_symbol="KCNQ2",
        mane_select_transcript="NM_172107.4",
        hgvsc="NM_172107.4:c.838C>T",
        hgvsp="NP_742105.1:p.Arg280Cys",
        acmg_class="Likely_Pathogenic",
        gnomad_af=None,
        alphamissense_score=0.89,  # above PP3 threshold 0.564
        clinvar_id=None,  # novel — not in ClinVar
        clinvar_class=None,
        inheritance_mode="AD",
        rules_applied=["PM2", "PP3"],
        is_novel=True,
    )


# ---------------------------------------------------------------------------
# VUS review date tests (ACGS 2024 §9)
# ---------------------------------------------------------------------------


class TestVUSReviewDates:
    """Tests for VUS review date scheduling (ACGS 2024 §9)."""

    def test_vus_gets_review_date(
        self, sample_config, vus_variant
    ) -> None:
        """VUS variant gets review date = report_date + 2 years.

        ACGS 2024 §9 mandates VUS review within 2 years of classification.
        """
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[vus_variant])
        assert vus_variant.review_date is not None, (
            "VUS variant must have a review date set (ACGS 2024 §9)"
        )
        # Review date should be 2 years from report date
        report_dt = date.fromisoformat(sample_config.report_date)
        review_dt = date.fromisoformat(vus_variant.review_date)
        expected_year = report_dt.year + 2
        assert review_dt.year == expected_year, (
            f"VUS review date year should be {expected_year}, got {review_dt.year}"
        )

    def test_pathogenic_no_review_date(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Pathogenic variants do not get a VUS review date."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        assert pathogenic_variant.review_date is None, (
            "Pathogenic variants should not have a VUS review date"
        )

    def test_vus_review_date_exact_two_years(
        self, vus_variant
    ) -> None:
        """VUS review date is exactly 2 years from the report date."""
        from reporting.report_generator import ReportConfig, ReportGenerator

        config = ReportConfig(
            patient_id="P001",
            sample_id="S001",
            report_date="2024-06-22",
        )
        gen = ReportGenerator(config=config, variants=[vus_variant])
        assert vus_variant.review_date == "2026-06-22", (
            f"VUS review for 2024-06-22 should be 2026-06-22, got {vus_variant.review_date}"
        )


# ---------------------------------------------------------------------------
# Novel P/LP ClinVar flag tests
# ---------------------------------------------------------------------------


class TestClinVarFlags:
    """Tests for novel P/LP ClinVar submission flagging."""

    def test_novel_lp_flagged_for_clinvar(
        self, sample_config, novel_lp_variant
    ) -> None:
        """Novel LP variant without ClinVar ID is flagged for submission.

        ACGS 2024 §5: novel P/LP variants should be submitted to ClinVar.
        """
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[novel_lp_variant])
        assert novel_lp_variant.pending_clinvar_submission is True, (
            "Novel LP without ClinVar ID should be flagged for submission"
        )

    def test_known_pathogenic_not_flagged(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Known Pathogenic variant with ClinVar ID is NOT flagged.

        If clinvar_id is set, the variant is already in ClinVar.
        """
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        assert pathogenic_variant.pending_clinvar_submission is False, (
            "Pathogenic variant already in ClinVar should not be flagged"
        )

    def test_vus_not_flagged_for_clinvar(
        self, sample_config, vus_variant
    ) -> None:
        """VUS variants are not flagged for ClinVar submission.

        Only P/LP variants require ClinVar submission.
        """
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[vus_variant])
        assert vus_variant.pending_clinvar_submission is False, (
            "VUS variants should not be flagged for ClinVar submission"
        )


# ---------------------------------------------------------------------------
# HTML generation tests
# ---------------------------------------------------------------------------


class TestHTMLGeneration:
    """Tests for HTML report generation."""

    def test_html_contains_patient_id(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Generated HTML contains patient ID."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        html = gen._render_html()
        assert "TEST-PATIENT-001" in html, (
            "HTML report must contain the patient ID"
        )

    def test_html_contains_mane_select(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Generated HTML contains MANE Select transcript.

        ACGS 2024 v1.2 requires MANE Select transcript notation.
        Morales et al. 2022 PMID:35356062.
        """
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        html = gen._render_html()
        assert "NM_007294.4" in html, (
            "HTML report must contain MANE Select transcript NM_007294.4"
        )

    def test_html_contains_acgs_citation(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Generated HTML contains ACGS 2024 v1.2 classification citation."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        html = gen._render_html()
        assert "ACGS" in html, (
            "HTML report must cite ACGS guidelines"
        )
        assert "ClinGen SVI" in html, (
            "HTML report must cite ClinGen SVI (PM2=Supporting)"
        )

    def test_html_contains_pm2_supporting_note(
        self, sample_config, pathogenic_variant
    ) -> None:
        """HTML footer notes PM2=Supporting per ClinGen SVI 2024."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        html = gen._render_html()
        assert "Supporting" in html and "PM2" in html, (
            "HTML report footer should note PM2=Supporting (ClinGen SVI 2024)"
        )

    def test_html_contains_gnomad_version(
        self, sample_config, pathogenic_variant
    ) -> None:
        """Generated HTML references gnomAD v4.1."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        html = gen._render_html()
        assert "4.1" in html, "HTML report should reference gnomAD v4.1"

    def test_html_vus_review_date_visible(
        self, sample_config, vus_variant
    ) -> None:
        """HTML report shows VUS review date (ACGS 2024 §9)."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[vus_variant])
        html = gen._render_html()
        assert "2026" in html, (
            "HTML report should show VUS review year (2024+2=2026)"
        )


# ---------------------------------------------------------------------------
# generate() integration test
# ---------------------------------------------------------------------------


class TestGenerateOutput:
    """Tests for the generate() method output files."""

    def test_generate_creates_html_file(
        self,
        tmp_path: Path,
        sample_config,
        pathogenic_variant,
        vus_variant,
    ) -> None:
        """generate() creates an HTML file at the expected path."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(
            config=sample_config,
            variants=[pathogenic_variant, vus_variant],
        )
        outputs = gen.generate(output_dir=tmp_path, generate_pdf=False)

        assert "html" in outputs, "generate() should return html key"
        html_path = outputs["html"]
        assert html_path.exists(), f"HTML file should exist: {html_path}"
        assert html_path.stat().st_size > 100, "HTML file should not be empty"

    def test_generate_creates_audit_jsonld(
        self,
        tmp_path: Path,
        sample_config,
        pathogenic_variant,
    ) -> None:
        """generate() creates a JSON-LD audit trail file."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(
            config=sample_config,
            variants=[pathogenic_variant],
        )
        outputs = gen.generate(output_dir=tmp_path, generate_pdf=False)

        assert "audit" in outputs, "generate() should return audit key"
        audit_path = outputs["audit"]
        assert audit_path.exists(), f"Audit file should exist: {audit_path}"

        with audit_path.open("r") as fh:
            audit_data = json.load(fh)

        assert "@context" in audit_data, "Audit file should be JSON-LD with @context"
        assert "@type" in audit_data

    def test_audit_contains_classification_scheme(
        self,
        tmp_path: Path,
        sample_config,
        pathogenic_variant,
    ) -> None:
        """Audit trail JSON-LD includes classification scheme with PM2=Supporting."""
        from reporting.report_generator import ReportGenerator

        gen = ReportGenerator(config=sample_config, variants=[pathogenic_variant])
        outputs = gen.generate(output_dir=tmp_path, generate_pdf=False)

        with outputs["audit"].open("r") as fh:
            audit = json.load(fh)

        scheme = audit.get("classification_scheme", {})
        pm2_weight = scheme.get("pm2_weight", "")
        assert "Supporting" in pm2_weight, (
            f"Audit should record PM2=Supporting (ClinGen SVI 2024): {pm2_weight}"
        )
