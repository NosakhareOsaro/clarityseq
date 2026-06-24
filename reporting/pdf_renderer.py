"""
reporting.pdf_renderer
========================
WeasyPrint PDF generation from HTML clinical reports.

WeasyPrint converts HTML/CSS to PDF.  It is the preferred PDF generator
for NHS GMS-style reports because:
    - Full CSS3 support including custom fonts (NHS fonts).
    - Page breaks, headers, and footers via CSS @page rules.
    - No JavaScript dependency (unlike headless Chrome approaches).
    - Pure Python — no external binary required.

Installation:
    pip install weasyprint
    On Ubuntu: apt-get install libpango-1.0-0 libharfbuzz0b libpangoft2-1.0-0

WeasyPrint version compatibility:
    Requires WeasyPrint ≥ 60.0 for full CSS3 support.
    IMPORTANT: WeasyPrint ≥ 54 changed the HTML() API.

References:
    WeasyPrint docs: https://doc.courtbouillon.org/weasyprint/
    NHS Identity guidelines: https://www.england.nhs.uk/nhsidentity/
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def render_pdf(html_path: Path, pdf_path: Path) -> Path:
    """Render an HTML report file to PDF using WeasyPrint.

    Args:
        html_path: Path to the source HTML report file.
        pdf_path: Destination path for the generated PDF.

    Returns:
        Path to the generated PDF file.

    Raises:
        ImportError: If WeasyPrint is not installed.
        FileNotFoundError: If html_path does not exist.
        RuntimeError: If WeasyPrint fails to generate the PDF.

    References:
        WeasyPrint ≥ 60.0: https://doc.courtbouillon.org/weasyprint/
    """
    if not html_path.exists():
        raise FileNotFoundError(f"HTML source not found: {html_path}")

    try:
        import weasyprint  # type: ignore[import]
    except ImportError as exc:
        raise ImportError(
            "WeasyPrint is required for PDF generation. "
            "Install with: pip install weasyprint\n"
            "On Ubuntu: apt-get install libpango-1.0-0 libharfbuzz0b"
        ) from exc

    try:
        # WeasyPrint ≥ 54: HTML() class accepts filename or url
        doc = weasyprint.HTML(filename=str(html_path))
        doc.write_pdf(str(pdf_path))
        logger.info("PDF generated: %s (%d bytes)", pdf_path, pdf_path.stat().st_size)
        return pdf_path
    except Exception as exc:
        raise RuntimeError(f"WeasyPrint PDF generation failed: {exc}") from exc


def render_pdf_from_string(html_string: str, pdf_path: Path) -> Path:
    """Render an HTML string to PDF using WeasyPrint.

    Args:
        html_string: HTML content as a string.
        pdf_path: Destination path for the generated PDF.

    Returns:
        Path to the generated PDF file.

    Raises:
        ImportError: If WeasyPrint is not installed.
        RuntimeError: If WeasyPrint fails to generate the PDF.

    References:
        WeasyPrint ≥ 60.0: https://doc.courtbouillon.org/weasyprint/
    """
    try:
        import weasyprint  # type: ignore[import]
    except ImportError as exc:
        raise ImportError(
            "WeasyPrint is required for PDF generation. "
            "Install with: pip install weasyprint"
        ) from exc

    try:
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        doc = weasyprint.HTML(string=html_string)
        doc.write_pdf(str(pdf_path))
        logger.info(
            "PDF generated from string: %s (%d bytes)",
            pdf_path,
            pdf_path.stat().st_size,
        )
        return pdf_path
    except Exception as exc:
        raise RuntimeError(f"WeasyPrint PDF generation failed: {exc}") from exc
