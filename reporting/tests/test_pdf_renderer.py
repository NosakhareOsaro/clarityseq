"""
reporting.tests.test_pdf_renderer
====================================
pytest tests for WeasyPrint PDF generation.

Tests cover:
    - render_pdf: file-based HTML to PDF conversion.
    - render_pdf_from_string: string-based HTML to PDF conversion.
    - Error handling: missing HTML file, ImportError, WeasyPrint failures.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from reporting.pdf_renderer import render_pdf, render_pdf_from_string


# ---------------------------------------------------------------------------
# render_pdf tests
# ---------------------------------------------------------------------------


class TestRenderPDF:
    """Tests for render_pdf() file-based conversion."""

    def test_missing_html_raises_file_not_found(self, tmp_path: Path) -> None:
        """Missing HTML source raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="HTML source not found"):
            render_pdf(
                html_path=tmp_path / "nonexistent.html",
                pdf_path=tmp_path / "output.pdf",
            )

    def test_weasyprint_not_installed_raises_import_error(
        self, tmp_path: Path
    ) -> None:
        """Missing WeasyPrint package raises ImportError with install hint."""
        html_path = tmp_path / "report.html"
        html_path.write_text("<html><body>Report</body></html>")

        with (
            patch("builtins.__import__", side_effect=ImportError("No module named 'weasyprint'")),
            pytest.raises(ImportError, match="WeasyPrint"),
        ):
            render_pdf(html_path, tmp_path / "out.pdf")

    def test_successful_render_returns_pdf_path(self, tmp_path: Path) -> None:
        """Successful render returns the pdf_path."""
        html_path = tmp_path / "report.html"
        html_path.write_text("<html><body>Report</body></html>")
        pdf_path = tmp_path / "report.pdf"

        # Create the output file so stat() works
        pdf_path.write_bytes(b"%PDF-1.4 fake pdf content")

        mock_doc = MagicMock()
        mock_html_class = MagicMock(return_value=mock_doc)
        mock_weasyprint = MagicMock()
        mock_weasyprint.HTML = mock_html_class

        with patch.dict("sys.modules", {"weasyprint": mock_weasyprint}):
            result = render_pdf(html_path, pdf_path)

        assert result == pdf_path
        mock_html_class.assert_called_once_with(filename=str(html_path))
        mock_doc.write_pdf.assert_called_once_with(str(pdf_path))

    def test_weasyprint_failure_raises_runtime_error(self, tmp_path: Path) -> None:
        """WeasyPrint runtime error is wrapped in RuntimeError."""
        html_path = tmp_path / "report.html"
        html_path.write_text("<html><body>Report</body></html>")
        pdf_path = tmp_path / "out.pdf"

        mock_doc = MagicMock()
        mock_doc.write_pdf.side_effect = Exception("PDF generation failed")
        mock_html_class = MagicMock(return_value=mock_doc)
        mock_weasyprint = MagicMock()
        mock_weasyprint.HTML = mock_html_class

        with (
            patch.dict("sys.modules", {"weasyprint": mock_weasyprint}),
            pytest.raises(RuntimeError, match="WeasyPrint PDF generation failed"),
        ):
            render_pdf(html_path, pdf_path)


# ---------------------------------------------------------------------------
# render_pdf_from_string tests
# ---------------------------------------------------------------------------


class TestRenderPDFFromString:
    """Tests for render_pdf_from_string() string-based conversion."""

    def test_weasyprint_not_installed_raises_import_error(
        self, tmp_path: Path
    ) -> None:
        """Missing WeasyPrint raises ImportError."""
        with (
            patch("builtins.__import__", side_effect=ImportError("No module named 'weasyprint'")),
            pytest.raises(ImportError, match="WeasyPrint"),
        ):
            render_pdf_from_string(
                "<html><body>Report</body></html>",
                tmp_path / "out.pdf",
            )

    def test_successful_string_render(self, tmp_path: Path) -> None:
        """Successful string render returns the pdf_path."""
        pdf_path = tmp_path / "out.pdf"
        pdf_path.write_bytes(b"%PDF-1.4 fake pdf content")

        html_string = "<html><body>Report</body></html>"

        mock_doc = MagicMock()
        mock_html_class = MagicMock(return_value=mock_doc)
        mock_weasyprint = MagicMock()
        mock_weasyprint.HTML = mock_html_class

        with patch.dict("sys.modules", {"weasyprint": mock_weasyprint}):
            result = render_pdf_from_string(html_string, pdf_path)

        assert result == pdf_path
        mock_html_class.assert_called_once_with(string=html_string)
        mock_doc.write_pdf.assert_called_once_with(str(pdf_path))

    def test_weasyprint_failure_raises_runtime_error(self, tmp_path: Path) -> None:
        """WeasyPrint failure in string render is wrapped in RuntimeError."""
        mock_doc = MagicMock()
        mock_doc.write_pdf.side_effect = Exception("Render error")
        mock_html_class = MagicMock(return_value=mock_doc)
        mock_weasyprint = MagicMock()
        mock_weasyprint.HTML = mock_html_class

        with (
            patch.dict("sys.modules", {"weasyprint": mock_weasyprint}),
            pytest.raises(RuntimeError, match="WeasyPrint PDF generation failed"),
        ):
            render_pdf_from_string("<html></html>", tmp_path / "out.pdf")

    def test_output_dir_created_if_missing(self, tmp_path: Path) -> None:
        """Parent directory of pdf_path is created if it doesn't exist."""
        pdf_path = tmp_path / "nested" / "dir" / "out.pdf"
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.write_bytes(b"%PDF-1.4 fake")

        mock_doc = MagicMock()
        mock_weasyprint = MagicMock()
        mock_weasyprint.HTML.return_value = mock_doc

        with patch.dict("sys.modules", {"weasyprint": mock_weasyprint}):
            render_pdf_from_string("<html></html>", pdf_path)

        assert pdf_path.parent.exists()
