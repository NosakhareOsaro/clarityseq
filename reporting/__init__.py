"""
reporting
=========
NHS GMS-style clinical report generation module for ClaritySeq.

Generates HTML + PDF + JSON-LD audit trail clinical reports compliant
with ACGS Best Practice Guidelines 2024 v1.2.

Submodules:
    report_generator — NHS GMS-style report generator (HTML + PDF + JSON-LD).
    pdf_renderer     — WeasyPrint PDF generation from HTML.
    audit_logger     — JSON-LD audit trail writer.

ACGS 2024 v1.2 compliance:
    §5 Table 2: MANE Select transcript notation.
    §9: VUS review date scheduling (date + 2 years).
    Novel P/LP: pending_clinvar_submissions flag.
    §6: Mitochondrial section with haplogroup first.
"""
