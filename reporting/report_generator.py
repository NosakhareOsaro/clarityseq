"""
reporting.report_generator
===========================
Generates NHS GMS-style clinical reports (HTML + PDF + JSON-LD audit).

ACGS 2024 v1.2 compliance:
1. MANE Select transcript notation (Morales 2022 PMID:35356062)
2. VUS review date scheduling (§9): VUS → review by date+2yr
3. ClinVar submission flag: novel P/LP → pending_clinvar_submissions
4. Classification scheme citation in footer

JSON-LD audit includes: gnomAD v4.1, VEP 111, AlphaMissense date, PM2=Supporting

NHS GMS report sections:
    - Patient and referral information
    - Summary of findings
    - Variant table (ACMG 2024 columns)
    - Mitochondrial section (§6 — haplogroup first)
    - Repeat expansions section
    - Pharmacogenomics section (CYP2D6)
    - Classification scheme and references
    - VUS review schedule

References:
    ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., 20 Feb 2024).
    Morales et al. 2022 Nature Methods PMID:35356062 (MANE Select).
    Richards et al. 2015 PMID:25741868 (ACMG/AMP framework).
    Tavtigian et al. 2020 PMID:32645316 (point-score system).
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Template directory
# ---------------------------------------------------------------------------

_TEMPLATE_DIR = Path(__file__).parent / "templates"

try:
    from jinja2 import Environment, FileSystemLoader, select_autoescape
    _JINJA2_AVAILABLE = True
    _jinja_env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
except ImportError:
    _JINJA2_AVAILABLE = False
    _jinja_env = None  # type: ignore[assignment]
    logger.warning("Jinja2 not installed; HTML report generation will be limited.")


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class VariantReportEntry:
    """A variant entry for the clinical report.

    ACGS 2024 v1.2 report columns:
        Gene | MANE Select Transcript | HGVSc | HGVSp |
        ACMG Class (ACGS 2024) | P(Path) [95% HDI] |
        gnomAD v4.1 AF | AlphaMissense Score | ClinVar |
        Inheritance | Review Date

    Attributes:
        gene_symbol: HGNC gene symbol.
        mane_select_transcript: MANE Select transcript ID with version
            (e.g. ``"NM_007294.4"``).  Morales 2022 PMID:35356062.
        hgvsc: HGVS cDNA notation on MANE Select transcript.
        hgvsp: HGVS protein notation.
        acmg_class: ACMG/AMP classification (ACGS 2024 v1.2).
        posterior_p: Posterior P(pathogenic) from BayesACMG model.
        hdi_lower: Lower bound of 95% HDI for P(pathogenic).
        hdi_upper: Upper bound of 95% HDI for P(pathogenic).
        gnomad_af: gnomAD v4.1 allele frequency.
        alphamissense_score: AlphaMissense missense pathogenicity score.
        clinvar_id: ClinVar accession.
        clinvar_class: ClinVar classification string.
        inheritance_mode: Inferred inheritance mode (AD/AR/XL/Mito).
        review_date: VUS review date (ACGS 2024 §9: VUS → date+2yr).
        rules_applied: List of ACMG rule IDs that applied.
        is_novel: True if this is a novel variant not in ClinVar.
        pending_clinvar_submission: True if novel P/LP pending submission.
    """

    gene_symbol: str
    mane_select_transcript: str
    hgvsc: str
    hgvsp: str
    acmg_class: str
    posterior_p: float | None = None
    hdi_lower: float | None = None
    hdi_upper: float | None = None
    gnomad_af: float | None = None
    alphamissense_score: float | None = None
    clinvar_id: str | None = None
    clinvar_class: str | None = None
    inheritance_mode: str = ""
    review_date: str | None = None      # ISO8601; set for VUS
    rules_applied: list[str] = field(default_factory=list)
    is_novel: bool = False
    pending_clinvar_submission: bool = False


@dataclass
class ReportConfig:
    """Configuration for a clinical report.

    Attributes:
        patient_id: Patient identifier (pseudonymised for report).
        sample_id: Lab sample identifier.
        referral_indication: Clinical referral indication (HPO/ICD-10).
        referring_clinician: Name of referring clinician.
        laboratory: Laboratory name.
        pipeline_version: GenomeForge pipeline version string.
        report_date: Report generation date (ISO8601).
        assembly: Genome assembly (default ``"GRCh38"``).
        acgs_version: ACGS guidelines version (default ``"2024 v1.2"``).
        vep_version: VEP version used for annotation.
        gnomad_version: gnomAD version used.
        include_mito: True to include mitochondrial section (ACGS 2024 §6).
        include_expansions: True to include repeat expansion section.
        include_pgx: True to include pharmacogenomics section.
    """

    patient_id: str
    sample_id: str
    referral_indication: str = ""
    referring_clinician: str = ""
    laboratory: str = "GenomeForge Genomics Laboratory"
    pipeline_version: str = "1.0.0"
    report_date: str = ""
    assembly: str = "GRCh38"
    acgs_version: str = "2024 v1.2"
    vep_version: str = "111"
    gnomad_version: str = "4.1"
    include_mito: bool = True
    include_expansions: bool = True
    include_pgx: bool = True

    def __post_init__(self) -> None:
        """Set report_date to now if not provided."""
        if not self.report_date:
            self.report_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")


# ---------------------------------------------------------------------------
# ReportGenerator
# ---------------------------------------------------------------------------


class ReportGenerator:
    """NHS GMS-style clinical variant report generator.

    Generates HTML, PDF, and JSON-LD audit trail for a set of variant
    report entries per ACGS Best Practice Guidelines 2024 v1.2.

    ACGS 2024 v1.2 compliance features:
    - MANE Select transcript used for all variant descriptions.
    - VUS variants have a review date scheduled 2 years from report date.
    - Novel P/LP variants are flagged for ClinVar submission.
    - JSON-LD audit trail records all data sources and tool versions.
    - Mitochondrial variants: haplogroup listed first per §6.

    Attributes:
        config: ReportConfig for this report.
        variants: List of VariantReportEntry objects to report.
        mito_variants: Mitochondrial variant entries (ACGS 2024 §6).
        expansion_data: Repeat expansion data dict.
        pgx_data: Pharmacogenomics data dict (CYP2D6 result).
    """

    def __init__(
        self,
        config: ReportConfig,
        variants: list[VariantReportEntry] | None = None,
        mito_variants: list[VariantReportEntry] | None = None,
        expansion_data: dict[str, Any] | None = None,
        pgx_data: dict[str, Any] | None = None,
    ) -> None:
        """Initialise the report generator.

        Args:
            config: ReportConfig with patient and pipeline metadata.
            variants: Nuclear variant report entries.
            mito_variants: Mitochondrial variant entries (ACGS 2024 §6).
            expansion_data: Repeat expansion results dict.
            pgx_data: Pharmacogenomics data (CYP2D6 diplotype, phenotype).
        """
        self.config = config
        self.variants: list[VariantReportEntry] = variants or []
        self.mito_variants: list[VariantReportEntry] = mito_variants or []
        self.expansion_data: dict[str, Any] = expansion_data or {}
        self.pgx_data: dict[str, Any] = pgx_data or {}

        # Apply ACGS 2024 §9: VUS review date (date + 2 years)
        self._schedule_vus_reviews()

        # Flag novel P/LP for ClinVar submission
        self._flag_novel_pathogenic()

    def _schedule_vus_reviews(self) -> None:
        """Schedule VUS review dates per ACGS 2024 §9.

        VUS variants must be reviewed within 2 years of classification.
        Sets review_date = report_date + 2 years for all VUS entries.

        Returns:
            None.

        References:
            ACGS 2024 v1.2 §9 — VUS review schedule.
        """
        from datetime import date

        try:
            report_dt = datetime.strptime(self.config.report_date, "%Y-%m-%d").date()
        except ValueError:
            report_dt = date.today()

        for v in self.variants + self.mito_variants:
            if v.acmg_class in ("VUS", "Uncertain_Significance") and not v.review_date:
                review_dt = report_dt.replace(year=report_dt.year + 2)
                v.review_date = review_dt.isoformat()
                logger.debug(
                    "ACGS 2024 §9: VUS %s:%s review scheduled %s",
                    v.gene_symbol, v.hgvsc, v.review_date,
                )

    def _flag_novel_pathogenic(self) -> None:
        """Flag novel P/LP variants for pending ClinVar submission.

        Novel Pathogenic or Likely Pathogenic variants not in ClinVar
        should be submitted.  Flags pending_clinvar_submission=True.

        Returns:
            None.

        References:
            ACGS 2024 v1.2 §5 — ClinVar submission guidance.
        """
        for v in self.variants + self.mito_variants:
            if (
                v.acmg_class in ("Pathogenic", "Likely_Pathogenic")
                and v.is_novel
                and not v.clinvar_id
            ):
                v.pending_clinvar_submission = True
                logger.debug(
                    "Novel P/LP flagged for ClinVar submission: %s %s",
                    v.gene_symbol, v.hgvsc,
                )

    def _render_html(self) -> str:
        """Render the clinical report HTML using Jinja2 templates.

        Returns:
            HTML string with all report sections.

        References:
            ACGS 2024 v1.2 — NHS GMS report structure.
        """
        if not _JINJA2_AVAILABLE or _jinja_env is None:
            return self._render_html_fallback()

        try:
            template = _jinja_env.get_template("base.html.j2")
            return template.render(
                config=self.config,
                variants=self.variants,
                mito_variants=self.mito_variants,
                expansion_data=self.expansion_data,
                pgx_data=self.pgx_data,
                report_date=self.config.report_date,
                acgs_version=self.config.acgs_version,
                gnomad_version=self.config.gnomad_version,
                vep_version=self.config.vep_version,
            )
        except Exception as exc:
            logger.warning("Jinja2 template rendering failed: %s. Using fallback.", exc)
            return self._render_html_fallback()

    def _render_html_fallback(self) -> str:
        """Minimal HTML fallback when Jinja2 or templates are unavailable.

        Returns:
            Minimal valid HTML report string.
        """
        p_count = len(self.variants)
        vus_count = sum(1 for v in self.variants if v.acmg_class == "VUS")
        pending_clinvar = sum(1 for v in self.variants if v.pending_clinvar_submission)

        var_rows = ""
        for v in self.variants:
            review_cell = v.review_date or ""
            af_display = f"{v.gnomad_af:.2e}" if v.gnomad_af is not None else "absent"
            am_display = f"{v.alphamissense_score:.3f}" if v.alphamissense_score is not None else "N/A"
            p_path = ""
            if v.posterior_p is not None:
                p_path = f"{v.posterior_p:.3f}"
                if v.hdi_lower is not None and v.hdi_upper is not None:
                    p_path += f" [{v.hdi_lower:.3f}–{v.hdi_upper:.3f}]"
            var_rows += (
                f"<tr>"
                f"<td>{v.gene_symbol}</td>"
                f"<td>{v.mane_select_transcript}</td>"
                f"<td>{v.hgvsc}</td>"
                f"<td>{v.hgvsp}</td>"
                f"<td><strong>{v.acmg_class}</strong></td>"
                f"<td>{p_path}</td>"
                f"<td>{af_display}</td>"
                f"<td>{am_display}</td>"
                f"<td>{v.clinvar_id or ''}</td>"
                f"<td>{v.inheritance_mode}</td>"
                f"<td>{review_cell}</td>"
                f"</tr>\n"
            )

        return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GenomeForge Clinical Report — {self.config.patient_id}</title>
<style>
body {{ font-family: Arial, sans-serif; font-size: 11pt; margin: 2cm; }}
h1 {{ color: #003087; }}
h2 {{ color: #003087; border-bottom: 1px solid #003087; }}
table {{ border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 9pt; }}
th {{ background: #003087; color: white; padding: 6px; text-align: left; }}
td {{ border: 1px solid #ccc; padding: 5px; }}
tr:nth-child(even) {{ background: #f5f5f5; }}
.pathogenic {{ color: #c00; font-weight: bold; }}
.vus {{ color: #f90; }}
.benign {{ color: #060; }}
.pending-clinvar {{ background: #ffe0b2; }}
.footer {{ font-size: 8pt; color: #666; border-top: 1px solid #ccc; margin-top: 2em; padding-top: 1em; }}
</style>
</head>
<body>
<h1>NHS GMS WGS Clinical Report</h1>
<table>
<tr><th>Patient ID</th><td>{self.config.patient_id}</td>
    <th>Sample ID</th><td>{self.config.sample_id}</td></tr>
<tr><th>Report Date</th><td>{self.config.report_date}</td>
    <th>Laboratory</th><td>{self.config.laboratory}</td></tr>
<tr><th>Referral</th><td colspan="3">{self.config.referral_indication}</td></tr>
<tr><th>Assembly</th><td>{self.config.assembly}</td>
    <th>Pipeline</th><td>GenomeForge v{self.config.pipeline_version}</td></tr>
</table>

<h2>Summary</h2>
<p>
Total variants reported: <strong>{p_count}</strong> |
VUS: <strong>{vus_count}</strong> |
Pending ClinVar submissions: <strong>{pending_clinvar}</strong>
</p>

<h2>Variants (ACGS {self.config.acgs_version})</h2>
<table>
<tr>
<th>Gene</th>
<th>MANE Select Transcript</th>
<th>HGVSc</th>
<th>HGVSp</th>
<th>ACMG Class (ACGS {self.config.acgs_version})</th>
<th>P(Path) [95% HDI]</th>
<th>gnomAD v{self.config.gnomad_version} AF</th>
<th>AlphaMissense Score</th>
<th>ClinVar</th>
<th>Inheritance</th>
<th>Review Date</th>
</tr>
{var_rows}
</table>

<div class="footer">
<p>
<strong>Classification scheme:</strong>
Richards et al. 2015 (PMID:25741868); Tavtigian et al. 2020 (PMID:32645316);
ACGS Best Practice Guidelines {self.config.acgs_version} (Durkie et al.);
ClinGen SVI Working Group 2024. PM2=Supporting per ClinGen SVI 2024.
AlphaMissense thresholds: ≥0.564→PP3, ≤0.340→BP4 (Cheng et al. 2023 PMID:37703350).
MANE Select transcript notation per Morales et al. 2022 (PMID:35356062).
</p>
<p>
<strong>Data sources:</strong>
VEP {self.config.vep_version} | gnomAD v{self.config.gnomad_version} (April 2024, 807,162 individuals) |
AlphaMissense (Cheng et al. 2023) | ClinVar (accessed {self.config.report_date}).
</p>
<p>
VUS review dates scheduled per ACGS 2024 §9 (2-year review mandate).
Novel P/LP variants flagged for ClinVar submission per ACGS 2024 §5.
</p>
</div>
</body>
</html>"""

    def _write_json_ld_audit(self, output_path: Path) -> None:
        """Write JSON-LD audit trail for the report.

        Records all data sources, tool versions, classification rules,
        and PM2 evidence weight used in this report.

        Args:
            output_path: Path to write the JSON-LD audit file.

        Returns:
            None.

        References:
            ACGS 2024 v1.2 — audit trail requirements.
            JSON-LD: https://json-ld.org/
        """
        from reporting.audit_logger import write_audit_log

        audit_data: dict[str, Any] = {
            "report_date": self.config.report_date,
            "patient_id": self.config.patient_id,
            "sample_id": self.config.sample_id,
            "pipeline_version": self.config.pipeline_version,
            "assembly": self.config.assembly,
            "tools": {
                "vep": self.config.vep_version,
                "gnomad": self.config.gnomad_version,
                "alphamissense": "Cheng et al. 2023 PMID:37703350",
            },
            "classification_scheme": {
                "framework": f"ACGS {self.config.acgs_version}",
                "point_system": "Tavtigian et al. 2020 PMID:32645316",
                "pm2_weight": "Supporting (1 pt) — ClinGen SVI 2024",
                "pp3_bp4_primary": "AlphaMissense (Cheng 2023): ≥0.564→PP3, ≤0.340→BP4",
                "mane_select": "Morales et al. 2022 PMID:35356062",
            },
            "variants": [
                {
                    "gene": v.gene_symbol,
                    "transcript": v.mane_select_transcript,
                    "hgvsc": v.hgvsc,
                    "hgvsp": v.hgvsp,
                    "acmg_class": v.acmg_class,
                    "rules_applied": v.rules_applied,
                    "gnomad_af": v.gnomad_af,
                    "alphamissense_score": v.alphamissense_score,
                    "posterior_p": v.posterior_p,
                    "vus_review_date": v.review_date,
                    "pending_clinvar_submission": v.pending_clinvar_submission,
                }
                for v in self.variants + self.mito_variants
            ],
        }

        write_audit_log(audit_data, output_path)

    def generate(
        self,
        output_dir: Path,
        generate_pdf: bool = True,
    ) -> dict[str, Path]:
        """Generate the clinical report (HTML + PDF + JSON-LD audit).

        Creates three output files:
        1. ``{sample_id}_report.html`` — HTML report (NHS GMS style).
        2. ``{sample_id}_report.pdf`` — PDF generated from HTML via WeasyPrint.
        3. ``{sample_id}_audit.jsonld`` — JSON-LD audit trail.

        Args:
            output_dir: Directory to write all output files.
            generate_pdf: If True (default), generate PDF via WeasyPrint.
                          Set False in environments where WeasyPrint is unavailable.

        Returns:
            Dict mapping ``"html"``, ``"pdf"`` (if generated), ``"audit"``
            to the corresponding output file Path objects.

        References:
            ACGS 2024 v1.2 — NHS GMS report format.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        sample_id = self.config.sample_id.replace("/", "_")
        outputs: dict[str, Path] = {}

        # Generate HTML
        html_path = output_dir / f"{sample_id}_report.html"
        html_content = self._render_html()
        html_path.write_text(html_content, encoding="utf-8")
        outputs["html"] = html_path
        logger.info("HTML report written: %s", html_path)

        # Generate PDF
        if generate_pdf:
            pdf_path = output_dir / f"{sample_id}_report.pdf"
            try:
                from reporting.pdf_renderer import render_pdf
                render_pdf(html_path, pdf_path)
                outputs["pdf"] = pdf_path
                logger.info("PDF report written: %s", pdf_path)
            except Exception as exc:
                logger.warning("PDF generation failed: %s. HTML report still available.", exc)

        # Write JSON-LD audit trail
        audit_path = output_dir / f"{sample_id}_audit.jsonld"
        self._write_json_ld_audit(audit_path)
        outputs["audit"] = audit_path
        logger.info("JSON-LD audit trail written: %s", audit_path)

        return outputs
