"""
Tests for VEPRunner with mocked subprocess calls.

Validates:
- VEP command construction (pick_order, MANE flags, plugin args)
- JSON output parsing into AnnotatedVariant objects
- MANE Select transcript detection
- SpliceAI score extraction
- gnomAD AF parsing
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Generator
from unittest.mock import MagicMock, patch

import pytest

from annotation.vep_runner import (
    VEP_PICK_ORDER,
    AnnotatedVariant,
    VEPRunner,
    _extract_spliceai_max,
    _safe_float,
    _safe_int,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def runner() -> VEPRunner:
    """Return a VEPRunner with default settings and no real paths."""
    return VEPRunner(vep_binary="vep")


@pytest.fixture
def runner_with_paths(tmp_path: Path) -> VEPRunner:
    """Return a VEPRunner with mock paths that exist."""
    am_tsv = tmp_path / "AlphaMissense_hg38.tsv.gz"
    dbnsfp = tmp_path / "dbNSFP4.7a_sorted.gz"
    spliceai_snv = tmp_path / "spliceai_scores.masked.snv.hg38.vcf.gz"
    spliceai_indel = tmp_path / "spliceai_scores.masked.indel.hg38.vcf.gz"
    for f in [am_tsv, dbnsfp, spliceai_snv, spliceai_indel]:
        f.touch()

    return VEPRunner(
        vep_binary="vep",
        cache_dir=str(tmp_path / "cache"),
        plugin_dir=str(tmp_path / "plugins"),
        alphamissense_tsv=str(am_tsv),
        dbnsfp_db=str(dbnsfp),
        spliceai_snv=str(spliceai_snv),
        spliceai_indel=str(spliceai_indel),
    )


# ---------------------------------------------------------------------------
# Mock VEP JSON output
# ---------------------------------------------------------------------------

MOCK_VEP_RECORD_MANE_SELECT = {
    "seq_region_name": "17",
    "start": 43094692,
    "end": 43094692,
    "allele_string": "G/A",
    "strand": 1,
    "transcript_consequences": [
        {
            "transcript_id": "ENST00000357654",
            "gene_id": "ENSG00000012048",
            "gene_symbol": "BRCA1",
            "consequence_terms": ["missense_variant"],
            "impact": "MODERATE",
            "hgvsc": "ENST00000357654.9:c.5266G>A",
            "hgvsp": "ENSP00000350283.3:p.Glu1756Lys",
            "mane_select": "NM_007294.4",
            "mane_plus_clinical": None,
            "canonical": 1,
            "pick": 1,
            "alphamissense_score": "0.82",
            "revel_score": "0.789",
            "cadd_phred": "28.5",
            "gnomADg_AF": "0.000005",
            "gnomADg_AC": "1",
            "gnomADg_AN": "200000",
            "SpliceAI_pred_DS_AG": "0.01",
            "SpliceAI_pred_DS_AL": "0.02",
            "SpliceAI_pred_DS_DG": "0.01",
            "SpliceAI_pred_DS_DL": "0.03",
        }
    ],
    "colocated_variants": [
        {
            "var_synonyms": {
                "ClinVar": ["VCV000012345"]
            }
        }
    ],
    "extras": {},
}

MOCK_VEP_RECORD_NO_MANE = {
    "seq_region_name": "1",
    "start": 12345678,
    "end": 12345678,
    "allele_string": "A/C",
    "strand": 1,
    "transcript_consequences": [
        {
            "transcript_id": "ENST00000999999",
            "gene_id": "ENSG00000999999",
            "gene_symbol": "TESTGENE",
            "consequence_terms": ["synonymous_variant"],
            "impact": "LOW",
            "hgvsc": "ENST00000999999.1:c.100A>C",
            "hgvsp": None,
            "mane_select": None,
            "mane_plus_clinical": None,
            "canonical": 1,
            "pick": 1,
            "alphamissense_score": None,
            "revel_score": None,
            "cadd_phred": "5.2",
            "gnomADg_AF": "0.25",
            "gnomADg_AC": "187500",
            "gnomADg_AN": "750000",
        }
    ],
    "colocated_variants": [],
    "extras": {},
}


# ---------------------------------------------------------------------------
# Tests for command construction
# ---------------------------------------------------------------------------


class TestVEPCommandConstruction:
    """Tests for _build_command() to verify correct VEP flags."""

    def test_pick_order_in_command(self, runner: VEPRunner, tmp_path: Path) -> None:
        """VEP pick_order flag must include mane_select first."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner._build_command(vcf, out)
        cmd_str = " ".join(cmd)
        assert VEP_PICK_ORDER in cmd_str
        assert "mane_select" in cmd_str

    def test_mane_flag_present(self, runner: VEPRunner, tmp_path: Path) -> None:
        """The --mane flag must be in the command for MANE annotation."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner._build_command(vcf, out)
        assert "--mane" in cmd

    def test_json_flag_present(self, runner: VEPRunner, tmp_path: Path) -> None:
        """The --json flag must be present for machine-readable output."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner._build_command(vcf, out)
        assert "--json" in cmd

    def test_grch38_assembly(self, runner: VEPRunner, tmp_path: Path) -> None:
        """GRCh38 assembly must be specified."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner._build_command(vcf, out)
        assert "GRCh38" in cmd

    def test_alphamissense_plugin_included_when_file_exists(
        self,
        runner_with_paths: VEPRunner,
        tmp_path: Path,
    ) -> None:
        """AlphaMissense plugin should be included when TSV file exists."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner_with_paths._build_command(vcf, out)
        cmd_str = " ".join(cmd)
        assert "AlphaMissense" in cmd_str

    def test_dbnsfp_plugin_included_when_file_exists(
        self,
        runner_with_paths: VEPRunner,
        tmp_path: Path,
    ) -> None:
        """dbNSFP plugin should be included when database file exists."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner_with_paths._build_command(vcf, out)
        cmd_str = " ".join(cmd)
        assert "dbNSFP" in cmd_str

    def test_spliceai_plugin_included_when_file_exists(
        self,
        runner_with_paths: VEPRunner,
        tmp_path: Path,
    ) -> None:
        """SpliceAI plugin should be included when score files exist."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner_with_paths._build_command(vcf, out)
        cmd_str = " ".join(cmd)
        assert "SpliceAI" in cmd_str

    def test_pangolin_plugin_always_included(
        self, runner: VEPRunner, tmp_path: Path
    ) -> None:
        """Pangolin splice scoring plugin should always be requested."""
        vcf = tmp_path / "in.vcf"
        vcf.touch()
        out = tmp_path / "out.json"
        cmd = runner._build_command(vcf, out)
        cmd_str = " ".join(cmd)
        assert "Pangolin" in cmd_str


# ---------------------------------------------------------------------------
# Tests for JSON parsing
# ---------------------------------------------------------------------------


class TestVEPJsonParsing:
    """Tests for _parse_record() and _parse_vep_json()."""

    def test_mane_select_variant_parsed_correctly(self, runner: VEPRunner) -> None:
        """MANE Select variant should be detected and annotated correctly."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)

        assert variant is not None
        assert variant.is_mane_select is True
        assert variant.gene_symbol == "BRCA1"
        assert variant.chrom == "chr17"
        assert variant.pos == 43094692
        assert variant.ref == "G"
        assert variant.alt == "A"

    def test_hgvsc_parsed(self, runner: VEPRunner) -> None:
        """HGVSc notation should be extracted from the picked consequence."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert "c.5266G>A" in (variant.hgvsc or "")

    def test_hgvsp_parsed(self, runner: VEPRunner) -> None:
        """HGVSp notation should be extracted."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert "Glu1756Lys" in (variant.hgvsp or "")

    def test_alphamissense_score_parsed(self, runner: VEPRunner) -> None:
        """AlphaMissense score should be parsed to float."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert variant.alphamissense_score == pytest.approx(0.82)

    def test_revel_score_parsed(self, runner: VEPRunner) -> None:
        """REVEL score should be parsed to float."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert variant.revel_score == pytest.approx(0.789)

    def test_gnomad_af_parsed(self, runner: VEPRunner) -> None:
        """gnomAD v4.1 AF should be parsed to float."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert variant.gnomad_af == pytest.approx(0.000005)

    def test_non_mane_variant_parsed(self, runner: VEPRunner) -> None:
        """Non-MANE variant should have is_mane_select=False."""
        variant = runner._parse_record(MOCK_VEP_RECORD_NO_MANE)
        assert variant is not None
        assert variant.is_mane_select is False
        assert variant.gene_symbol == "TESTGENE"

    def test_chrom_prefix_added(self, runner: VEPRunner) -> None:
        """Chromosomes without 'chr' prefix should be normalised."""
        # Both records use numeric chromosome strings from VEP JSON
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert variant.chrom.startswith("chr")

    def test_consequence_terms_list_populated(self, runner: VEPRunner) -> None:
        """Consequence terms list should be populated."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert "missense_variant" in variant.consequence_terms

    def test_impact_parsed(self, runner: VEPRunner) -> None:
        """VEP impact should be set."""
        variant = runner._parse_record(MOCK_VEP_RECORD_MANE_SELECT)
        assert variant is not None
        assert variant.impact == "MODERATE"

    def test_empty_record_returns_none(self, runner: VEPRunner) -> None:
        """Record with no transcript consequences should return None."""
        empty_record = {
            "seq_region_name": "1",
            "start": 1,
            "allele_string": "A/T",
            "transcript_consequences": [],
        }
        result = runner._parse_record(empty_record)
        assert result is None

    def test_parse_vep_json_reads_file(self, runner: VEPRunner, tmp_path: Path) -> None:
        """_parse_vep_json should parse a real JSON file."""
        json_path = tmp_path / "output.json"
        with open(json_path, "w") as fh:
            fh.write(json.dumps(MOCK_VEP_RECORD_MANE_SELECT) + "\n")
            fh.write(json.dumps(MOCK_VEP_RECORD_NO_MANE) + "\n")

        variants = runner._parse_vep_json(json_path)
        assert len(variants) == 2
        assert variants[0].gene_symbol == "BRCA1"
        assert variants[1].gene_symbol == "TESTGENE"


# ---------------------------------------------------------------------------
# Tests for run_vep subprocess mocking
# ---------------------------------------------------------------------------


class TestRunVEP:
    """Tests for run_vep() using mocked subprocess."""

    def test_run_vep_raises_on_missing_vcf(
        self, runner: VEPRunner, tmp_path: Path
    ) -> None:
        """run_vep should raise FileNotFoundError for non-existent input."""
        with pytest.raises(FileNotFoundError, match="Input VCF not found"):
            runner.run_vep(tmp_path / "nonexistent.vcf")

    def test_run_vep_calls_subprocess(
        self, runner: VEPRunner, tmp_path: Path
    ) -> None:
        """run_vep should call subprocess.run with the VEP command."""
        vcf_path = tmp_path / "input.vcf"
        vcf_path.touch()
        out_path = tmp_path / "output.json"

        # Write a mock VEP output file
        with open(out_path, "w") as fh:
            fh.write(json.dumps(MOCK_VEP_RECORD_MANE_SELECT) + "\n")

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result) as mock_run, \
             patch.object(runner, "_parse_vep_json", return_value=[]) as mock_parse:
            runner.run_vep(vcf_path, out_path)

        mock_run.assert_called_once()
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert "vep" in cmd[0]

    def test_run_vep_raises_on_nonzero_exit(
        self, runner: VEPRunner, tmp_path: Path
    ) -> None:
        """run_vep should raise CalledProcessError on VEP failure."""
        vcf_path = tmp_path / "input.vcf"
        vcf_path.touch()

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "Error: cache not found"
        mock_result.stdout = ""

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(subprocess.CalledProcessError):
                runner.run_vep(vcf_path, tmp_path / "out.json")


# ---------------------------------------------------------------------------
# Tests for utility functions
# ---------------------------------------------------------------------------


class TestUtilities:
    """Tests for module-level utility functions."""

    def test_safe_float_converts_string(self) -> None:
        """_safe_float should convert numeric string to float."""
        assert _safe_float("0.82") == pytest.approx(0.82)

    def test_safe_float_returns_none_for_dot(self) -> None:
        """_safe_float should return None for VEP's missing value '.'."""
        assert _safe_float(".") is None

    def test_safe_float_returns_none_for_none(self) -> None:
        """_safe_float should return None for Python None."""
        assert _safe_float(None) is None

    def test_safe_int_converts_string(self) -> None:
        """_safe_int should convert numeric string to int."""
        assert _safe_int("200000") == 200000

    def test_safe_int_returns_none_for_invalid(self) -> None:
        """_safe_int should return None for non-numeric input."""
        assert _safe_int("not_a_number") is None

    def test_extract_spliceai_max_from_tc(self) -> None:
        """_extract_spliceai_max should return max delta score."""
        tc = {
            "SpliceAI_pred_DS_AG": "0.01",
            "SpliceAI_pred_DS_AL": "0.02",
            "SpliceAI_pred_DS_DG": "0.01",
            "SpliceAI_pred_DS_DL": "0.55",
        }
        result = _extract_spliceai_max(tc, {})
        assert result == pytest.approx(0.55)

    def test_extract_spliceai_max_returns_none_when_missing(self) -> None:
        """_extract_spliceai_max should return None when no scores present."""
        result = _extract_spliceai_max({}, {})
        assert result is None
