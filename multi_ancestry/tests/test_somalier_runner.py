"""
multi_ancestry.tests.test_somalier_runner
==========================================
pytest tests for the somalier runner module.

Tests cover:
    - run_somalier_extract: subprocess mocking and error handling.
    - run_somalier_relate: empty input guard and subprocess mocking.
    - _parse_somalier_relate_output: TSV parsing.
    - _parse_somalier_ancestry_output: ancestry TSV parsing.
    - SomalierAncestryResult and SomalierRelatednessResult dataclasses.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from multi_ancestry.somalier_runner import (
    SomalierAncestryResult,
    SomalierRelatednessResult,
    _parse_somalier_ancestry_output,
    _parse_somalier_relate_output,
    run_somalier_ancestry,
    run_somalier_extract,
    run_somalier_relate,
)


# ---------------------------------------------------------------------------
# run_somalier_extract tests
# ---------------------------------------------------------------------------


class TestRunSomalierExtract:
    """Tests for run_somalier_extract()."""

    def test_missing_bam_raises_file_not_found(self, tmp_path: Path) -> None:
        """Missing BAM raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            run_somalier_extract(
                bam_path=tmp_path / "nonexistent.bam",
                ref_fasta=tmp_path / "ref.fa",
                sites_vcf=tmp_path / "sites.vcf.gz",
                output_dir=tmp_path / "out",
            )

    def test_missing_ref_raises_file_not_found(self, tmp_path: Path) -> None:
        """Missing reference FASTA raises FileNotFoundError."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        with pytest.raises(FileNotFoundError):
            run_somalier_extract(
                bam_path=bam,
                ref_fasta=tmp_path / "nonexistent.fa",
                sites_vcf=tmp_path / "sites.vcf.gz",
                output_dir=tmp_path / "out",
            )

    def test_successful_extract_returns_somalier_file(self, tmp_path: Path) -> None:
        """Successful somalier extract returns path to .somalier file."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        sites = tmp_path / "sites.vcf.gz"
        sites.touch()
        output_dir = tmp_path / "out"

        # Create the expected output file
        output_dir.mkdir()
        expected_output = output_dir / "sample.somalier"
        expected_output.touch()

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with patch("subprocess.run", return_value=mock_proc):
            result = run_somalier_extract(bam, ref, sites, output_dir)

        assert result == expected_output

    def test_nonzero_returncode_raises_runtime_error(self, tmp_path: Path) -> None:
        """Non-zero returncode from somalier raises RuntimeError."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        sites = tmp_path / "sites.vcf.gz"
        sites.touch()

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "somalier: reference not found"

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("pathlib.Path.mkdir"),
            pytest.raises(RuntimeError, match="somalier extract failed"),
        ):
            run_somalier_extract(bam, ref, sites, tmp_path / "out")

    def test_missing_output_after_success_raises_runtime_error(
        self, tmp_path: Path
    ) -> None:
        """Successful subprocess.run but no .somalier file produced raises RuntimeError."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()
        sites = tmp_path / "sites.vcf.gz"
        sites.touch()
        output_dir = tmp_path / "out"

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with (
            patch("subprocess.run", return_value=mock_proc),
            pytest.raises(RuntimeError, match="did not produce expected output"),
        ):
            run_somalier_extract(bam, ref, sites, output_dir)


# ---------------------------------------------------------------------------
# run_somalier_relate tests
# ---------------------------------------------------------------------------


class TestRunSomalierRelate:
    """Tests for run_somalier_relate()."""

    def test_empty_files_raises_value_error(self, tmp_path: Path) -> None:
        """Empty somalier_files list raises ValueError."""
        with pytest.raises(ValueError, match="No somalier files"):
            run_somalier_relate([], tmp_path / "out")

    def test_successful_relate(self, tmp_path: Path) -> None:
        """Successful somalier relate parses TSV output."""
        somalier_files = [tmp_path / "s1.somalier"]
        output_dir = tmp_path / "out"
        output_dir.mkdir()

        # Create mock output TSV files
        (output_dir / "somalier.samples.tsv").write_text(
            "#sample_id\tphenotype\n"
            "S1\t1\n",
        )
        (output_dir / "somalier.pairs.tsv").write_text(
            "sample_a\tsample_b\trelatedness\tibs0\n"
            "S1\tS2\t0.5\t10\n",
        )

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with patch("subprocess.run", return_value=mock_proc):
            result = run_somalier_relate(somalier_files, output_dir)

        assert "samples" in result
        assert "pairs" in result

    def test_ped_file_adds_ped_flag(self, tmp_path: Path) -> None:
        """A ped_file that exists adds a --ped flag to the somalier relate command."""
        somalier_files = [tmp_path / "s1.somalier"]
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        ped_file = tmp_path / "cohort.ped"
        ped_file.write_text("FAM1\tS1\t0\t0\t1\t2\n")

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with patch("subprocess.run", return_value=mock_proc) as mock_run:
            run_somalier_relate(somalier_files, output_dir, ped_file=ped_file)

        called_cmd = mock_run.call_args[0][0]
        assert "--ped" in called_cmd
        assert str(ped_file) in called_cmd

    def test_nonexistent_ped_file_omits_ped_flag(self, tmp_path: Path) -> None:
        """A ped_file argument that does not exist on disk is not added to the command."""
        somalier_files = [tmp_path / "s1.somalier"]
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        ped_file = tmp_path / "nonexistent.ped"

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with patch("subprocess.run", return_value=mock_proc) as mock_run:
            run_somalier_relate(somalier_files, output_dir, ped_file=ped_file)

        called_cmd = mock_run.call_args[0][0]
        assert "--ped" not in called_cmd

    def test_nonzero_returncode_raises_runtime_error(self, tmp_path: Path) -> None:
        """Non-zero returncode from somalier relate raises RuntimeError."""
        somalier_files = [tmp_path / "s1.somalier"]
        output_dir = tmp_path / "out"

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "somalier: relate failed unexpectedly"

        with (
            patch("subprocess.run", return_value=mock_proc),
            pytest.raises(RuntimeError, match="somalier relate failed"),
        ):
            run_somalier_relate(somalier_files, output_dir)


# ---------------------------------------------------------------------------
# run_somalier_ancestry tests
# ---------------------------------------------------------------------------


class TestRunSomalierAncestry:
    """Tests for run_somalier_ancestry()."""

    def test_empty_files_raises_value_error(self, tmp_path: Path) -> None:
        """Empty somalier_files list raises ValueError."""
        with pytest.raises(ValueError, match="No somalier files"):
            run_somalier_ancestry([], tmp_path / "ref_panel", tmp_path / "out")

    def test_nonzero_returncode_raises_runtime_error(self, tmp_path: Path) -> None:
        """Non-zero returncode from somalier ancestry raises RuntimeError."""
        somalier_files = [tmp_path / "s1.somalier"]
        ref_panel_dir = tmp_path / "ref_panel"
        ref_panel_dir.mkdir()

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "somalier: ancestry inference failed"

        with (
            patch("subprocess.run", return_value=mock_proc),
            pytest.raises(RuntimeError, match="somalier ancestry failed"),
        ):
            run_somalier_ancestry(somalier_files, ref_panel_dir, tmp_path / "out")

    def test_successful_run_parses_ancestry_output(self, tmp_path: Path) -> None:
        """Successful somalier ancestry run parses the produced TSV into results."""
        somalier_files = [tmp_path / "s1.somalier"]
        ref_panel_dir = tmp_path / "ref_panel"
        ref_panel_dir.mkdir()
        (ref_panel_dir / "ref1.somalier").touch()
        output_dir = tmp_path / "out"

        def fake_run(*args, **kwargs):
            # Simulate somalier ancestry writing its output TSV as a side effect.
            output_dir.mkdir(parents=True, exist_ok=True)
            tsv = output_dir / "ancestry.somalier-ancestry.tsv"
            tsv.write_text(
                "sample_id\tpredicted_ancestry\tgiven_ancestry\tEUR\tAFR\tAMR\tEAS\tMID\tSAS\tPC1\tPC2\tPC3\n"
                "S1\tEUR\t0.9\t0.9\t0.05\t0.02\t0.01\t0.01\t0.01\t1.1\t0.2\t0.1\n"
            )
            return MagicMock(returncode=0)

        with patch("subprocess.run", side_effect=fake_run):
            results = run_somalier_ancestry(somalier_files, ref_panel_dir, output_dir)

        assert len(results) == 1
        assert results[0].sample_id == "S1"
        assert results[0].predicted_ancestry == "EUR"

    def test_reference_panel_labels_used_when_present(self, tmp_path: Path) -> None:
        """The --labels arg points at the reference panel labels file when it exists."""
        somalier_files = [tmp_path / "s1.somalier"]
        ref_panel_dir = tmp_path / "ref_panel"
        ref_panel_dir.mkdir()
        labels_file = ref_panel_dir / "1kg+hgdp.somalier-ancestry.tsv"
        labels_file.write_text("sample_id\tancestry\n")
        output_dir = tmp_path / "out"

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with patch("subprocess.run", return_value=mock_proc) as mock_run:
            run_somalier_ancestry(somalier_files, ref_panel_dir, output_dir)

        called_cmd = mock_run.call_args[0][0]
        assert str(labels_file) in called_cmd


# ---------------------------------------------------------------------------
# _parse_somalier_relate_output tests
# ---------------------------------------------------------------------------


class TestParseSomalierRelateOutput:
    """Tests for _parse_somalier_relate_output()."""

    def test_returns_empty_when_no_files(self, tmp_path: Path) -> None:
        """Returns empty dicts when no TSV files exist."""
        result = _parse_somalier_relate_output(tmp_path)
        assert result["samples"] == []
        assert result["pairs"] == []

    def test_parses_samples_tsv(self, tmp_path: Path) -> None:
        """Parses the samples TSV correctly."""
        (tmp_path / "somalier.samples.tsv").write_text(
            "#sample_id\tphenotype\n"
            "SAMPLE1\t1\n"
            "SAMPLE2\t2\n",
        )
        result = _parse_somalier_relate_output(tmp_path)
        assert len(result["samples"]) == 2

    def test_parses_pairs_tsv(self, tmp_path: Path) -> None:
        """Parses the pairs TSV with relatedness values."""
        (tmp_path / "somalier.pairs.tsv").write_text(
            "sample_a\tsample_b\trelatedness\tibs0\n"
            "S1\tS2\t0.5\t5\n",
        )
        result = _parse_somalier_relate_output(tmp_path)
        assert len(result["pairs"]) == 1
        assert result["pairs"][0]["relatedness"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# _parse_somalier_ancestry_output tests
# ---------------------------------------------------------------------------


class TestParseSomalierAncestryOutput:
    """Tests for _parse_somalier_ancestry_output()."""

    def test_missing_file_returns_empty(self, tmp_path: Path) -> None:
        """Missing ancestry TSV returns empty list (no error)."""
        result = _parse_somalier_ancestry_output(tmp_path / "nonexistent.tsv")
        assert result == []

    def test_parses_ancestry_fractions(self, tmp_path: Path) -> None:
        """Parses ancestry fraction columns from TSV."""
        tsv = tmp_path / "ancestry.somalier-ancestry.tsv"
        tsv.write_text(
            "sample_id\tpredicted_ancestry\tgiven_ancestry\tEUR\tAFR\tAMR\tEAS\tMID\tSAS\tPC1\tPC2\tPC3\n"
            "SAMPLE1\tEUR\t0.90\t0.90\t0.05\t0.02\t0.01\t0.01\t0.01\t1.2\t0.5\t0.3\n",
        )
        results = _parse_somalier_ancestry_output(tsv)
        assert len(results) == 1
        r = results[0]
        assert r.sample_id == "SAMPLE1"
        assert r.predicted_ancestry == "EUR"
        assert r.ancestry_fractions.get("EUR") == pytest.approx(0.90)
        assert r.pc1 == pytest.approx(1.2)

    def test_blank_lines_skipped(self, tmp_path: Path) -> None:
        """Blank lines in the TSV are skipped without disrupting parsing."""
        tsv = tmp_path / "ancestry.somalier-ancestry.tsv"
        tsv.write_text(
            "sample_id\tpredicted_ancestry\tgiven_ancestry\tEUR\n"
            "\n"
            "SAMPLE1\tEUR\t0.9\t0.9\n"
            "\n",
        )
        results = _parse_somalier_ancestry_output(tsv)
        assert len(results) == 1
        assert results[0].sample_id == "SAMPLE1"


# ---------------------------------------------------------------------------
# Data model tests
# ---------------------------------------------------------------------------


class TestSomalierDataModels:
    """Tests for SomalierAncestryResult and SomalierRelatednessResult."""

    def test_ancestry_result_construction(self) -> None:
        """SomalierAncestryResult can be constructed with all fields."""
        result = SomalierAncestryResult(
            sample_id="S1",
            predicted_ancestry="EUR",
            predicted_ancestry_p=0.92,
            ancestry_fractions={"EUR": 0.92, "AFR": 0.08},
            is_admixed=False,
            pc1=1.0,
            pc2=0.5,
            pc3=0.1,
        )
        assert result.sample_id == "S1"
        assert result.predicted_ancestry == "EUR"
        assert result.is_admixed is False

    def test_relatedness_result_construction(self) -> None:
        """SomalierRelatednessResult can be constructed correctly."""
        result = SomalierRelatednessResult(
            sample_a="S1",
            sample_b="S2",
            relatedness=0.5,
            ibs0=5,
            ibs2=100,
            n_sites=17_000,
            expected_relationship="parent-child",
        )
        assert result.relatedness == pytest.approx(0.5)
        assert result.expected_relationship == "parent-child"
