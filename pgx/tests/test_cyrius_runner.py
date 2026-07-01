"""
pgx.tests.test_cyrius_runner
==============================
pytest unit tests for the CYP2D6 Cyrius runner module.

Tests cover:
    - phenotype_from_activity: all phenotype boundary conditions.
    - activity_score_from_diplotype: diplotype parsing.
    - parse_cyrius_output: JSON parsing into CYP2D6Result.
    - run_cyrius: mocked subprocess, error handling.

References:
    CPIC CYP2D6 activity score bins (2023 update):
        0.0        → PM (Poor Metaboliser)
        0.25–1.0   → IM (Intermediate Metaboliser)
        1.25–2.5   → NM (Normal Metaboliser)
        > 2.5      → UM (Ultrarapid Metaboliser)
    Twesigomwe et al. 2022 npj Genomic Medicine PMID:35513406 (Cyrius validation)
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from pgx.cyrius_runner import (
    CYP2D6Result,
    MetaboliserPhenotype,
    activity_score_from_diplotype,
    genotype_sample,
    parse_cyrius_output,
    phenotype_from_activity,
    run_cyrius,
)


# ---------------------------------------------------------------------------
# phenotype_from_activity tests
# ---------------------------------------------------------------------------


class TestPhenotypeFromActivity:
    """Tests for pgx.cyrius_runner.phenotype_from_activity().

    CPIC activity score bins (2023 update):
        0.0        → PM (activity_score == 0)
        0.25–1.0   → IM
        1.25–2.5   → NM
        > 2.5      → UM
    """

    def test_pm_zero_activity(self) -> None:
        """Activity score 0.0 → PM (Poor Metaboliser)."""
        assert phenotype_from_activity(0.0) == MetaboliserPhenotype.POOR

    def test_im_low_boundary(self) -> None:
        """Activity score 0.25 → IM (Intermediate Metaboliser)."""
        assert phenotype_from_activity(0.25) == MetaboliserPhenotype.INTERMEDIATE

    def test_im_typical(self) -> None:
        """Activity score 0.75 → IM."""
        assert phenotype_from_activity(0.75) == MetaboliserPhenotype.INTERMEDIATE

    def test_im_high_boundary(self) -> None:
        """Activity score 1.0 → IM (upper bound of IM range)."""
        assert phenotype_from_activity(1.0) == MetaboliserPhenotype.INTERMEDIATE

    def test_nm_lower_boundary(self) -> None:
        """Activity score 1.25 → NM (Normal Metaboliser)."""
        assert phenotype_from_activity(1.25) == MetaboliserPhenotype.NORMAL

    def test_nm_typical(self) -> None:
        """Activity score 2.0 → NM (typical *1/*1 diplotype)."""
        assert phenotype_from_activity(2.0) == MetaboliserPhenotype.NORMAL

    def test_nm_upper_boundary(self) -> None:
        """Activity score 2.5 → NM (upper bound of NM range)."""
        assert phenotype_from_activity(2.5) == MetaboliserPhenotype.NORMAL

    def test_um_just_above_nm(self) -> None:
        """Activity score 3.0 → UM (gene duplication *1xN)."""
        assert phenotype_from_activity(3.0) == MetaboliserPhenotype.ULTRARAPID

    def test_negative_returns_indeterminate(self) -> None:
        """Negative activity score (indeterminate diplotype) → Indeterminate."""
        assert phenotype_from_activity(-1.0) == MetaboliserPhenotype.INDETERMINATE

    def test_gap_between_im_and_nm_returns_indeterminate(self) -> None:
        """Activity score 1.1 falls in the gap between IM (≤1.0) and NM (≥1.25)."""
        result = phenotype_from_activity(1.1)
        assert result == MetaboliserPhenotype.INDETERMINATE


# ---------------------------------------------------------------------------
# activity_score_from_diplotype tests
# ---------------------------------------------------------------------------


class TestActivityScoreFromDiplotype:
    """Tests for pgx.cyrius_runner.activity_score_from_diplotype()."""

    def test_star1_star1_score(self) -> None:
        """*1/*1 diplotype has activity score 2.0 (1.0 + 1.0)."""
        assert activity_score_from_diplotype("*1/*1") == pytest.approx(2.0)

    def test_star4_star4_score(self) -> None:
        """*4/*4 (two no-function alleles) → activity score 0.0."""
        assert activity_score_from_diplotype("*4/*4") == pytest.approx(0.0)

    def test_star1_star4_score(self) -> None:
        """*1/*4 → activity score 1.0 (1.0 + 0.0)."""
        assert activity_score_from_diplotype("*1/*4") == pytest.approx(1.0)

    def test_star10_star17_score(self) -> None:
        """*10/*17 → 0.25 + 0.5 = 0.75 (reduced function alleles)."""
        assert activity_score_from_diplotype("*10/*17") == pytest.approx(0.75)

    def test_duplication_star1xn_star4(self) -> None:
        """*1xN/*4 → 2.0 + 0.0 = 2.0 (duplication counts as 2 copies)."""
        assert activity_score_from_diplotype("*1xN/*4") == pytest.approx(2.0)

    def test_indeterminate_diplotype_returns_negative(self) -> None:
        """Indeterminate diplotype string returns -1.0."""
        assert activity_score_from_diplotype("Indeterminate") == -1.0

    def test_empty_diplotype_returns_negative(self) -> None:
        """Empty string returns -1.0."""
        assert activity_score_from_diplotype("") == -1.0

    def test_unknown_allele_defaults_to_zero(self) -> None:
        """Unknown allele (*99) defaults to 0.0 (conservative)."""
        score = activity_score_from_diplotype("*1/*99")
        assert score == pytest.approx(1.0)  # *1=1.0, *99=0.0 (default)

    def test_single_allele_no_slash_returns_negative(self) -> None:
        """A diplotype string without a '/' separator returns -1.0 (invalid format)."""
        assert activity_score_from_diplotype("*1") == -1.0

    def test_three_alleles_returns_negative(self) -> None:
        """A diplotype string with more than two alleles returns -1.0 (invalid format)."""
        assert activity_score_from_diplotype("*1/*2/*3") == -1.0


# ---------------------------------------------------------------------------
# parse_cyrius_output tests
# ---------------------------------------------------------------------------


class TestParseCyriusOutput:
    """Tests for pgx.cyrius_runner.parse_cyrius_output()."""

    def _write_cyrius_json(self, tmp_path: Path, **fields: object) -> Path:
        """Write a Cyrius output JSON file and return its path."""
        data = {
            "Sample": fields.get("Sample", "SAMPLE1"),
            "Genotype": fields.get("Genotype", "*1/*4"),
            "Filter": fields.get("Filter", "PASS"),
            "Copy_Number": fields.get("Copy_Number", 2),
        }
        p = tmp_path / "cyrius.json"
        p.write_text(json.dumps(data))
        return p

    def test_parse_nm_diplotype(self, tmp_path: Path) -> None:
        """*1/*1 parses to NM phenotype."""
        json_path = self._write_cyrius_json(tmp_path, Genotype="*1/*1", Copy_Number=2)
        result = parse_cyrius_output(json_path, "SAMPLE1")
        assert isinstance(result, CYP2D6Result)
        assert result.diplotype == "*1/*1"
        assert result.phenotype == MetaboliserPhenotype.NORMAL
        assert result.activity_score == pytest.approx(2.0)
        assert result.sample_id == "SAMPLE1"

    def test_parse_pm_diplotype(self, tmp_path: Path) -> None:
        """*4/*4 parses to PM phenotype."""
        json_path = self._write_cyrius_json(tmp_path, Genotype="*4/*4", Copy_Number=2)
        result = parse_cyrius_output(json_path, "PATIENT")
        assert result.phenotype == MetaboliserPhenotype.POOR
        assert result.activity_score == pytest.approx(0.0)
        assert result.cyrius_filter == "PASS"

    def test_parse_indeterminate_diplotype(self, tmp_path: Path) -> None:
        """Indeterminate diplotype parses correctly."""
        json_path = self._write_cyrius_json(
            tmp_path, Genotype="Indeterminate", Filter="FAIL", Copy_Number=2
        )
        result = parse_cyrius_output(json_path, "S1")
        assert result.phenotype == MetaboliserPhenotype.INDETERMINATE
        assert result.cyrius_filter == "FAIL"

    def test_parse_copy_number(self, tmp_path: Path) -> None:
        """Copy number is parsed from JSON."""
        json_path = self._write_cyrius_json(tmp_path, Copy_Number=3)
        result = parse_cyrius_output(json_path, "S1")
        assert result.gene_copy_number == 3

    def test_raw_output_stored(self, tmp_path: Path) -> None:
        """Raw Cyrius JSON is stored in result."""
        json_path = self._write_cyrius_json(tmp_path, Genotype="*1/*2")
        result = parse_cyrius_output(json_path, "S1")
        assert result.raw_cyrius_output.get("Genotype") == "*1/*2"

    def test_activity_score_clamped_to_zero(self, tmp_path: Path) -> None:
        """Negative activity scores (indeterminate) are clamped to 0.0 in result."""
        json_path = self._write_cyrius_json(tmp_path, Genotype="Indeterminate")
        result = parse_cyrius_output(json_path, "S1")
        assert result.activity_score >= 0.0


# ---------------------------------------------------------------------------
# run_cyrius tests (mocked subprocess)
# ---------------------------------------------------------------------------


class TestRunCyrius:
    """Tests for pgx.cyrius_runner.run_cyrius() with mocked subprocess."""

    def test_run_cyrius_returns_output_json_path(self, tmp_path: Path) -> None:
        """run_cyrius returns the path to the output JSON file."""
        bam_path = tmp_path / "sample.bam"
        bam_path.touch()
        ref_fasta = tmp_path / "ref.fa"
        ref_fasta.touch()
        output_dir = tmp_path / "output"

        expected_json = output_dir / "SAMPLE1.json"

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("pathlib.Path.exists", return_value=True),
            patch("pathlib.Path.mkdir"),
        ):
            result_path = run_cyrius(bam_path, ref_fasta, output_dir, "SAMPLE1")

        assert result_path == expected_json

    def test_run_cyrius_raises_on_nonzero_returncode(self, tmp_path: Path) -> None:
        """run_cyrius raises subprocess.CalledProcessError if Cyrius fails."""
        bam_path = tmp_path / "sample.bam"
        bam_path.touch()

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "Cyrius: error: reference genome not found."

        with (
            patch("subprocess.run", side_effect=Exception("Cyrius failed")),
            patch("pathlib.Path.mkdir"),
            patch("pathlib.Path.exists", return_value=False),
            pytest.raises(Exception),
        ):
            run_cyrius(bam_path, tmp_path / "ref.fa", tmp_path / "out", "SAMPLE1")

    def test_run_cyrius_skips_if_output_exists(self, tmp_path: Path) -> None:
        """run_cyrius skips execution if output JSON already exists."""
        bam_path = tmp_path / "sample.bam"
        bam_path.touch()
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        existing_json = output_dir / "SAMPLE1.json"
        existing_json.touch()

        mock_run = MagicMock()
        with patch("subprocess.run", mock_run):
            result = run_cyrius(bam_path, tmp_path / "ref.fa", output_dir, "SAMPLE1")

        mock_run.assert_not_called()
        assert result == existing_json

    def test_run_cyrius_raises_if_output_missing_after_run(
        self, tmp_path: Path
    ) -> None:
        """run_cyrius raises RuntimeError if output JSON not found after successful run."""
        bam_path = tmp_path / "sample.bam"
        bam_path.touch()
        output_dir = tmp_path / "output"

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("pathlib.Path.mkdir"),
            pytest.raises(RuntimeError, match="output JSON not found"),
        ):
            run_cyrius(bam_path, tmp_path / "ref.fa", output_dir, "SAMPLE1")

    def test_run_cyrius_returns_path_after_successful_run(
        self, tmp_path: Path
    ) -> None:
        """run_cyrius runs Cyrius via subprocess and returns the produced output path."""
        bam_path = tmp_path / "sample.bam"
        bam_path.touch()
        ref_fasta = tmp_path / "ref.fa"
        ref_fasta.touch()
        output_dir = tmp_path / "output"
        expected_json = output_dir / "SAMPLE1.json"

        def fake_run(cmd, check, capture_output):
            # Simulate Cyrius writing its output JSON as a side effect.
            output_dir.mkdir(parents=True, exist_ok=True)
            expected_json.write_text("{}")
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=fake_run) as mock_run:
            result_path = run_cyrius(bam_path, ref_fasta, output_dir, "SAMPLE1")

        assert result_path == expected_json
        assert expected_json.exists()
        mock_run.assert_called_once()


# ---------------------------------------------------------------------------
# genotype_sample tests (mocked run_cyrius + parse_cyrius_output)
# ---------------------------------------------------------------------------


class TestGenotypeSample:
    """Tests for pgx.cyrius_runner.genotype_sample() end-to-end orchestration."""

    def test_genotype_sample_calls_run_then_parse(self, tmp_path: Path) -> None:
        """genotype_sample runs Cyrius then parses its output, returning the result."""
        bam = tmp_path / "sample.bam"
        ref = tmp_path / "ref.fa"
        output_dir = tmp_path / "out"
        cyrius_json = tmp_path / "cyrius.json"

        expected_result = CYP2D6Result(
            sample_id="SAMPLE1",
            diplotype="*1/*4",
            activity_score=1.0,
            phenotype=MetaboliserPhenotype.INTERMEDIATE,
            gene_copy_number=2,
            cyrius_filter="PASS",
        )

        with (
            patch("pgx.cyrius_runner.run_cyrius", return_value=cyrius_json) as mock_run,
            patch(
                "pgx.cyrius_runner.parse_cyrius_output", return_value=expected_result
            ) as mock_parse,
        ):
            result = genotype_sample(bam, ref, output_dir, "SAMPLE1")

        assert result is expected_result
        mock_run.assert_called_once_with(bam, ref, output_dir, "SAMPLE1")
        mock_parse.assert_called_once_with(cyrius_json, "SAMPLE1")
