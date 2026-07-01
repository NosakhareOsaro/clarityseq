"""
prioritisation.tests.test_exomiser_client
==========================================
pytest tests for the Exomiser 14 REST API client.

Tests cover:
    - _build_analysis_request: analysis settings construction.
    - _parse_exomiser_results: result parsing and sorting.
    - ExomiserResult dataclass construction.
    - run_exomiser: integration with mocked HTTP.

References:
    Robinson et al. 2023 Nature Genetics PMID:37604970 (Exomiser).
    Jacobsen et al. 2022 PMID:35705716 (Phenopackets v2).
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from prioritisation.exomiser_client import (
    ExomiserResult,
    _build_analysis_request,
    _parse_exomiser_results,
    _poll_results,
    _submit_analysis,
    run_exomiser,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_PHENOPACKET = {
    "id": "PP-001",
    "subject": {"id": "PATIENT-001", "sex": "FEMALE"},
    "phenotypicFeatures": [
        {"type": {"id": "HP:0001250", "label": "Seizures"}},
    ],
    "metaData": {"created": "2024-01-01T00:00:00Z", "createdBy": "test"},
}

MOCK_EXOMISER_RESULTS = {
    "status": "COMPLETED",
    "results": [
        {
            "rank": 1,
            "geneSymbol": "CDKL5",
            "combinedScore": 0.98,
            "phenotypeScore": 0.95,
            "variantScore": 0.90,
            "acmgClassification": "PATHOGENIC",
            "modeOfInheritance": "X_RECESSIVE",
        },
        {
            "rank": 2,
            "geneSymbol": "BRCA1",
            "combinedScore": 0.72,
            "phenotypeScore": 0.65,
            "variantScore": 0.80,
            "acmgClassification": "VUS",
            "modeOfInheritance": "AUTOSOMAL_DOMINANT",
        },
    ],
}


# ---------------------------------------------------------------------------
# ExomiserResult dataclass
# ---------------------------------------------------------------------------


class TestExomiserResult:
    """Tests for ExomiserResult dataclass."""

    def test_construction(self) -> None:
        """ExomiserResult can be constructed with required fields."""
        result = ExomiserResult(
            rank=1,
            gene_symbol="CDKL5",
            combined_score=0.98,
            phenotype_score=0.95,
            variant_score=0.90,
        )
        assert result.rank == 1
        assert result.gene_symbol == "CDKL5"
        assert result.combined_score == pytest.approx(0.98)

    def test_default_acmg_class(self) -> None:
        """Default acmg_class is 'VUS'."""
        result = ExomiserResult(
            rank=1,
            gene_symbol="GENE",
            combined_score=0.5,
            phenotype_score=0.4,
            variant_score=0.6,
        )
        assert result.acmg_class == "VUS"

    def test_default_contributing_variants_empty(self) -> None:
        """Default contributing_variants is empty list."""
        result = ExomiserResult(
            rank=1,
            gene_symbol="GENE",
            combined_score=0.5,
            phenotype_score=0.4,
            variant_score=0.6,
        )
        assert result.contributing_variants == []


# ---------------------------------------------------------------------------
# _build_analysis_request tests
# ---------------------------------------------------------------------------


class TestBuildAnalysisRequest:
    """Tests for _build_analysis_request()."""

    def test_phenopacket_included(self) -> None:
        """Phenopacket is included in the analysis request."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        assert req["phenopacket"] == MINIMAL_PHENOPACKET

    def test_default_inheritance_modes(self) -> None:
        """All 5 default inheritance modes are included."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        modes = req["inheritanceModes"]
        assert "AUTOSOMAL_DOMINANT" in modes
        assert "AUTOSOMAL_RECESSIVE" in modes
        assert "X_DOMINANT" in modes
        assert "X_RECESSIVE" in modes
        assert "MITOCHONDRIAL" in modes

    def test_custom_inheritance_modes(self) -> None:
        """Custom inheritance modes replace defaults."""
        req = _build_analysis_request(
            MINIMAL_PHENOPACKET,
            inheritance_modes=["AUTOSOMAL_DOMINANT"],
        )
        assert list(req["inheritanceModes"].keys()) == ["AUTOSOMAL_DOMINANT"]

    def test_vcf_path_included_when_provided(self) -> None:
        """VCF path is included in request when specified."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET, vcf_path="/data/sample.vcf.gz")
        assert req["vcf"] == "/data/sample.vcf.gz"

    def test_vcf_path_absent_by_default(self) -> None:
        """VCF path key is absent when not specified."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        assert "vcf" not in req

    def test_gene_panel_added_to_steps(self) -> None:
        """Gene panel step is inserted into the analysis steps."""
        req = _build_analysis_request(
            MINIMAL_PHENOPACKET,
            gene_panel=["CDKL5", "KCNQ2"],
        )
        step_types = [list(s.keys())[0] for s in req["steps"]]
        assert "genePanel" in step_types

    def test_gene_panel_absent_by_default(self) -> None:
        """No genePanel step when gene_panel is None."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        step_types = [list(s.keys())[0] for s in req["steps"]]
        assert "genePanel" not in step_types

    def test_custom_min_frequency(self) -> None:
        """Custom min_frequency is applied to frequencyFilter step."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET, min_frequency=0.5)
        freq_steps = [s for s in req["steps"] if "frequencyFilter" in s]
        assert freq_steps[0]["frequencyFilter"]["maxFrequency"] == pytest.approx(0.5)

    def test_analysis_mode_is_pass_only(self) -> None:
        """Default analysisMode is PASS_ONLY."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        assert req["analysisMode"] == "PASS_ONLY"

    def test_gnomad_frequency_sources_included(self) -> None:
        """gnomAD frequency sources are in the request."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        assert "GNOMAD_E_NFE" in req["frequencySources"]

    def test_pathogenicity_sources_include_clinvar(self) -> None:
        """ClinVar is included in pathogenicity sources."""
        req = _build_analysis_request(MINIMAL_PHENOPACKET)
        assert "CLINVAR" in req["pathogenicitySources"]


# ---------------------------------------------------------------------------
# _parse_exomiser_results tests
# ---------------------------------------------------------------------------


class TestParseExomiserResults:
    """Tests for _parse_exomiser_results()."""

    def test_returns_list_of_exomiser_results(self) -> None:
        """Returns a list of ExomiserResult objects."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        assert len(results) == 2
        assert all(isinstance(r, ExomiserResult) for r in results)

    def test_sorted_by_rank(self) -> None:
        """Results are sorted by rank (ascending)."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        ranks = [r.rank for r in results]
        assert ranks == sorted(ranks)

    def test_gene_symbol_extracted(self) -> None:
        """Gene symbol is extracted from geneSymbol field."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        assert results[0].gene_symbol == "CDKL5"

    def test_combined_score_extracted(self) -> None:
        """Combined score is extracted as float."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        assert results[0].combined_score == pytest.approx(0.98)

    def test_acmg_class_extracted(self) -> None:
        """ACMG classification is extracted."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        assert results[0].acmg_class == "PATHOGENIC"

    def test_inheritance_mode_extracted(self) -> None:
        """Inheritance mode is extracted."""
        results = _parse_exomiser_results(MOCK_EXOMISER_RESULTS)
        assert results[0].inheritance_mode == "X_RECESSIVE"

    def test_empty_results_returns_empty_list(self) -> None:
        """Empty results dict returns empty list."""
        results = _parse_exomiser_results({"status": "COMPLETED", "results": []})
        assert results == []

    def test_genes_key_used_as_fallback(self) -> None:
        """'genes' key is used as fallback when 'results' is absent."""
        data = {
            "status": "COMPLETED",
            "genes": [
                {
                    "rank": 1,
                    "geneSymbol": "BRCA2",
                    "combinedScore": 0.85,
                    "phenotypeScore": 0.80,
                    "variantScore": 0.90,
                }
            ],
        }
        results = _parse_exomiser_results(data)
        assert len(results) == 1
        assert results[0].gene_symbol == "BRCA2"

    def test_rank_defaults_to_position_when_absent(self) -> None:
        """When rank key absent, falls back to 1-based position."""
        data = {
            "results": [
                {
                    "geneSymbol": "GENE",
                    "combinedScore": 0.5,
                    "phenotypeScore": 0.4,
                    "variantScore": 0.6,
                }
            ]
        }
        results = _parse_exomiser_results(data)
        assert results[0].rank == 1


# ---------------------------------------------------------------------------
# run_exomiser integration test
# ---------------------------------------------------------------------------


class TestRunExomiser:
    """Integration tests for run_exomiser() with mocked HTTP."""

    def test_returns_ranked_results(self) -> None:
        """run_exomiser returns ranked ExomiserResult list on success."""
        with (
            patch("prioritisation.exomiser_client._submit_analysis", return_value="job-123"),
            patch(
                "prioritisation.exomiser_client._poll_results",
                return_value=MOCK_EXOMISER_RESULTS,
            ),
        ):
            results = run_exomiser(MINIMAL_PHENOPACKET)

        assert len(results) == 2
        assert results[0].gene_symbol == "CDKL5"

    def test_gene_panel_passed_through(self) -> None:
        """Gene panel is included in the analysis request."""
        with (
            patch(
                "prioritisation.exomiser_client._submit_analysis", return_value="job-456"
            ) as mock_submit,
            patch(
                "prioritisation.exomiser_client._poll_results",
                return_value={"status": "COMPLETED", "results": []},
            ),
        ):
            run_exomiser(MINIMAL_PHENOPACKET, gene_panel=["CDKL5"])

        call_args = mock_submit.call_args[0][0]
        step_types = [list(s.keys())[0] for s in call_args["steps"]]
        assert "genePanel" in step_types


# ---------------------------------------------------------------------------
# _submit_analysis tests (mocked httpx.post)
# ---------------------------------------------------------------------------


class TestSubmitAnalysis:
    """Tests for _submit_analysis() REST submission with mocked HTTP."""

    def test_returns_job_id_from_id_key(self) -> None:
        """job_id is extracted from the 'id' key in the response."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"id": "job-abc-123"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.post", return_value=mock_response) as mock_post:
            job_id = _submit_analysis({"phenopacket": {}})

        assert job_id == "job-abc-123"
        mock_post.assert_called_once()

    def test_falls_back_to_jobid_key(self) -> None:
        """job_id falls back to 'jobId' key when 'id' is absent."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"jobId": "job-xyz-789"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.post", return_value=mock_response):
            job_id = _submit_analysis({"phenopacket": {}})

        assert job_id == "job-xyz-789"

    def test_raises_runtime_error_when_no_job_id(self, monkeypatch) -> None:
        """RuntimeError is raised when the API response has no job ID.

        The function is decorated with @retry(stop_after_attempt(3)); patch
        time.sleep to avoid real backoff delays during the retries.
        """
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        mock_response = MagicMock()
        mock_response.json.return_value = {}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.post", return_value=mock_response):
            with pytest.raises(RuntimeError, match="no job ID"):
                _submit_analysis({"phenopacket": {}})

    def test_includes_authorization_header_when_api_key_set(self, monkeypatch) -> None:
        """Authorization header is included when EXOMISER_API_KEY is configured."""
        monkeypatch.setattr("prioritisation.exomiser_client._EXOMISER_API_KEY", "secret-key")
        mock_response = MagicMock()
        mock_response.json.return_value = {"id": "job-1"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.post", return_value=mock_response) as mock_post:
            _submit_analysis({"phenopacket": {}})

        _, kwargs = mock_post.call_args
        assert kwargs["headers"]["Authorization"] == "Bearer secret-key"

    def test_no_authorization_header_when_api_key_unset(self, monkeypatch) -> None:
        """Authorization header is absent when EXOMISER_API_KEY is not set."""
        monkeypatch.setattr("prioritisation.exomiser_client._EXOMISER_API_KEY", None)
        mock_response = MagicMock()
        mock_response.json.return_value = {"id": "job-1"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.post", return_value=mock_response) as mock_post:
            _submit_analysis({"phenopacket": {}})

        _, kwargs = mock_post.call_args
        assert "Authorization" not in kwargs["headers"]


# ---------------------------------------------------------------------------
# _poll_results tests (mocked httpx.get)
# ---------------------------------------------------------------------------


class TestPollResults:
    """Tests for _poll_results() polling loop with mocked HTTP."""

    def test_returns_data_when_completed_on_first_attempt(self, monkeypatch) -> None:
        """Returns the result dict as soon as status is COMPLETED."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        mock_response = MagicMock()
        mock_response.json.return_value = {"status": "COMPLETED", "results": []}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.get", return_value=mock_response):
            data = _poll_results("job-1")

        assert data["status"] == "COMPLETED"

    def test_raises_runtime_error_on_failed_status(self, monkeypatch) -> None:
        """RuntimeError is raised when status is FAILED."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        mock_response = MagicMock()
        mock_response.json.return_value = {"status": "FAILED", "message": "bad phenopacket"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.get", return_value=mock_response):
            with pytest.raises(RuntimeError, match="failed"):
                _poll_results("job-2")

    def test_raises_runtime_error_on_error_status(self, monkeypatch) -> None:
        """RuntimeError is raised when status is ERROR."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        mock_response = MagicMock()
        mock_response.json.return_value = {"status": "ERROR", "message": "internal error"}
        mock_response.raise_for_status = MagicMock()

        with patch("prioritisation.exomiser_client.httpx.get", return_value=mock_response):
            with pytest.raises(RuntimeError, match="failed"):
                _poll_results("job-3")

    def test_polls_multiple_times_until_completed(self, monkeypatch) -> None:
        """Polls repeatedly while status is pending, then returns on COMPLETED."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)

        running_response = MagicMock()
        running_response.json.return_value = {"status": "RUNNING"}
        running_response.raise_for_status = MagicMock()

        completed_response = MagicMock()
        completed_response.json.return_value = {"status": "COMPLETED", "results": [{"rank": 1}]}
        completed_response.raise_for_status = MagicMock()

        with patch(
            "prioritisation.exomiser_client.httpx.get",
            side_effect=[running_response, running_response, completed_response],
        ) as mock_get:
            data = _poll_results("job-4")

        assert data["status"] == "COMPLETED"
        assert mock_get.call_count == 3

    def test_raises_timeout_after_max_attempts(self, monkeypatch) -> None:
        """RuntimeError with 'timed out' is raised after _MAX_POLL_ATTEMPTS."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        monkeypatch.setattr("prioritisation.exomiser_client._MAX_POLL_ATTEMPTS", 2)

        pending_response = MagicMock()
        pending_response.json.return_value = {"status": "RUNNING"}
        pending_response.raise_for_status = MagicMock()

        with patch(
            "prioritisation.exomiser_client.httpx.get", return_value=pending_response
        ):
            with pytest.raises(RuntimeError, match="timed out"):
                _poll_results("job-5")

    def test_includes_authorization_header_when_api_key_set(self, monkeypatch) -> None:
        """Authorization header is included in poll requests when API key is set."""
        monkeypatch.setattr("time.sleep", lambda *_a, **_kw: None)
        monkeypatch.setattr("prioritisation.exomiser_client._EXOMISER_API_KEY", "poll-key")

        mock_response = MagicMock()
        mock_response.json.return_value = {"status": "COMPLETED", "results": []}
        mock_response.raise_for_status = MagicMock()

        with patch(
            "prioritisation.exomiser_client.httpx.get", return_value=mock_response
        ) as mock_get:
            _poll_results("job-6")

        _, kwargs = mock_get.call_args
        assert kwargs["headers"]["Authorization"] == "Bearer poll-key"
