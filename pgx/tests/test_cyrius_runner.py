"""
pgx.tests.test_cyrius_runner
==============================
pytest unit tests for the CYP2D6 Cyrius runner module.

Tests cover:
    - classify_phenotype: all four phenotype boundary conditions.
    - run_cyrius: mocked subprocess, JSON parsing, error handling.
    - CYP2D6Result dataclass construction.

References:
    CPIC CYP2D6 guideline activity score thresholds (2022 update):
        ≥2.25 → UM, 1.25-2.25 → NM, 0.25-1.25 → IM, 0 → PM.
    Aliev et al. 2022 NPJ Genomic Medicine PMID:35264608 (Cyrius).
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest


# ---------------------------------------------------------------------------
# classify_phenotype tests
# ---------------------------------------------------------------------------


class TestClassifyPhenotype:
    """Tests for pgx.cyrius_runner.classify_phenotype().

    CPIC activity score thresholds (2022 update):
        ≥ 2.25 → UM (Ultrarapid Metaboliser)
        1.25 – 2.25 → NM (Normal Metaboliser)
        0.25 – 1.25 (exclusive) → IM (Intermediate Metaboliser)
        0.0 → PM (Poor Metaboliser)
    """

    def test_pm_zero_activity(self) -> None:
        """Activity score 0.0 → PM (Poor Metaboliser).

        PM patients cannot convert codeine to morphine (active metabolite).
        CPIC: avoid codeine and tramadol.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(0.0) == "PM", (
            "Activity score 0.0 must map to PM (CPIC 2022 guideline)"
        )

    def test_im_low_boundary(self) -> None:
        """Activity score 0.25 → IM (Intermediate Metaboliser).

        IM patients have reduced enzyme activity.
        CPIC: reduce dose for prodrugs.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(0.25) == "IM", (
            "Activity score 0.25 must map to IM (CPIC 2022 guideline)"
        )

    def test_im_high_boundary(self) -> None:
        """Activity score 1.24 → IM (just below NM threshold 1.25).

        The NM lower boundary is 1.25 (inclusive).
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(1.24) == "IM", (
            "Activity score 1.24 must map to IM (below NM threshold 1.25)"
        )

    def test_nm_lower_boundary(self) -> None:
        """Activity score 1.25 → NM (Normal Metaboliser).

        1.25 is the inclusive lower boundary for NM per CPIC 2022.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(1.25) == "NM", (
            "Activity score 1.25 must map to NM (CPIC 2022 lower NM boundary)"
        )

    def test_nm_typical(self) -> None:
        """Activity score 2.0 → NM (typical *1/*1 diplotype).

        *1/*1 diplotype has activity score 1.0 + 1.0 = 2.0 → NM.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(2.0) == "NM", (
            "Activity score 2.0 (*1/*1) must map to NM"
        )

    def test_nm_upper_boundary(self) -> None:
        """Activity score 2.25 → UM (upper boundary of NM is <2.25).

        2.25 is the UM lower boundary (inclusive) per CPIC 2022.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(2.25) == "UM", (
            "Activity score 2.25 must map to UM (CPIC 2022 UM lower boundary)"
        )

    def test_um_high_activity(self) -> None:
        """Activity score 3.0 → UM (gene duplication *1xN).

        *1/*1xN diplotype has activity 1.0 + 2.0 = 3.0 → UM.
        UM patients may reach toxic morphine levels from codeine.
        """
        from pgx.cyrius_runner import classify_phenotype

        assert classify_phenotype(3.0) == "UM", (
            "Activity score 3.0 (*1/*1xN duplication) must map to UM"
        )

    def test_pm_very_low_nonzero(self) -> None:
        """Activity score > 0 but < 0.25 → IM not PM.

        Values above 0 but below 0.25 are IM per CPIC 2022 thresholds.
        """
        from pgx.cyrius_runner import classify_phenotype

        # 0.01 > 0, so it returns IM (IM threshold is > 0)
        result = classify_phenotype(0.01)
        assert result == "IM", (
            "Activity score 0.01 (>0) must map to IM per CPIC 2022"
        )


# ---------------------------------------------------------------------------
# run_cyrius tests (mocked subprocess)
# ---------------------------------------------------------------------------


class TestRunCyrius:
    """Tests for pgx.cyrius_runner.run_cyrius() with mocked subprocess."""

    def _make_cyrius_json(
        self,
        diplotype: str = "*1/*4",
        activity_score: float = 1.0,
    ) -> str:
        """Create a minimal Cyrius output JSON string for testing.

        Args:
            diplotype: Diplotype string (e.g. ``"*1/*4"``).
            activity_score: CPIC activity score.

        Returns:
            JSON string mimicking Cyrius output format.
        """
        return json.dumps({
            "diplotype": diplotype,
            "activity_score": activity_score,
            "confidence": 0.99,
            "copy_number": 2,
        })

    def test_run_cyrius_nm_diplotype(self, tmp_path: Path) -> None:
        """run_cyrius with *1/*1 diplotype returns NM phenotype.

        *1/*1 has activity score 2.0 (1.0 + 1.0) → NM.
        """
        from pgx.cyrius_runner import CYP2D6Result, run_cyrius

        bam_path = tmp_path / "sample.bam"
        bam_path.touch()
        ref_fasta = tmp_path / "ref.fa"
        ref_fasta.touch()
        output_dir = tmp_path / "output"

        cyrius_output = self._make_cyrius_json("*1/*1", activity_score=2.0)
        output_file = output_dir / "sample_cyrius.json"

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("builtins.open", mock_open(read_data=cyrius_output)),
            patch("pathlib.Path.mkdir"),
        ):
            result = run_cyrius(bam_path, ref_fasta, output_dir)

        assert isinstance(result, CYP2D6Result)
        assert result.diplotype == "*1/*1"
        assert result.metaboliser_phenotype == "NM"
        assert result.activity_score == 2.0
        assert result.sample_id == "sample"

    def test_run_cyrius_pm_diplotype(self, tmp_path: Path) -> None:
        """run_cyrius with *4/*4 diplotype returns PM phenotype.

        *4/*4 (two no-function alleles) has activity score 0.0 → PM.
        PM patients must avoid codeine and tramadol (CPIC Level A).
        """
        from pgx.cyrius_runner import CYP2D6Result, run_cyrius

        bam_path = tmp_path / "patient.bam"
        bam_path.touch()
        ref_fasta = tmp_path / "ref.fa"
        ref_fasta.touch()
        output_dir = tmp_path / "output"

        cyrius_output = self._make_cyrius_json("*4/*4", activity_score=0.0)

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("builtins.open", mock_open(read_data=cyrius_output)),
            patch("pathlib.Path.mkdir"),
        ):
            result = run_cyrius(bam_path, ref_fasta, output_dir)

        assert result.metaboliser_phenotype == "PM"
        assert result.activity_score == 0.0
        assert "codeine" in result.cpic_recommendation.lower(), (
            "PM recommendation should mention codeine"
        )
        assert "tramadol" in result.cpic_recommendation.lower(), (
            "PM recommendation should mention tramadol"
        )

    def test_run_cyrius_raises_on_nonzero_returncode(self, tmp_path: Path) -> None:
        """run_cyrius raises RuntimeError if Cyrius subprocess fails.

        Checks that subprocess failures are propagated as RuntimeError
        with the stderr message included.
        """
        from pgx.cyrius_runner import run_cyrius

        bam_path = tmp_path / "sample.bam"
        bam_path.touch()

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "Cyrius: error: reference genome not found."

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("pathlib.Path.mkdir"),
            pytest.raises(RuntimeError, match="Cyrius failed"),
        ):
            run_cyrius(bam_path, tmp_path / "ref.fa", tmp_path / "out")

    def test_run_cyrius_star_alleles_split(self, tmp_path: Path) -> None:
        """run_cyrius correctly splits diplotype into star_alleles list."""
        from pgx.cyrius_runner import run_cyrius

        bam_path = tmp_path / "s.bam"
        bam_path.touch()
        cyrius_output = self._make_cyrius_json("*1/*41", activity_score=1.5)

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("builtins.open", mock_open(read_data=cyrius_output)),
            patch("pathlib.Path.mkdir"),
        ):
            result = run_cyrius(bam_path, tmp_path / "ref.fa", tmp_path / "out")

        assert result.star_alleles == ["*1", "*41"], (
            f"star_alleles should be ['*1', '*41'], got {result.star_alleles}"
        )

    def test_run_cyrius_im_diplotype(self, tmp_path: Path) -> None:
        """run_cyrius with *1/*41 diplotype returns IM phenotype.

        *1/*41 has activity score 1.5 → NM (1.25-2.25).
        *41 activity = 0.5 (decreased function); *1 = 1.0.
        Total = 1.5 → NM.
        """
        from pgx.cyrius_runner import run_cyrius

        bam_path = tmp_path / "s.bam"
        bam_path.touch()
        cyrius_output = self._make_cyrius_json("*1/*41", activity_score=1.5)

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stderr = ""

        with (
            patch("subprocess.run", return_value=mock_proc),
            patch("builtins.open", mock_open(read_data=cyrius_output)),
            patch("pathlib.Path.mkdir"),
        ):
            result = run_cyrius(bam_path, tmp_path / "ref.fa", tmp_path / "out")

        assert result.metaboliser_phenotype == "NM"


# ---------------------------------------------------------------------------
# CPIC recommendation tests
# ---------------------------------------------------------------------------


class TestCPICRecommendation:
    """Tests for _get_cpic_recommendation() helper."""

    def test_nm_recommendation_contains_standard_dosing(self) -> None:
        """NM recommendation mentions standard dosing."""
        from pgx.cyrius_runner import _get_cpic_recommendation

        rec = _get_cpic_recommendation("NM")
        assert "standard" in rec.lower() or "label" in rec.lower(), (
            f"NM recommendation should mention standard dosing: {rec}"
        )

    def test_pm_recommendation_mentions_codeine(self) -> None:
        """PM recommendation mentions codeine avoidance.

        CPIC Level A: PM patients must avoid codeine (PMID:30447227).
        """
        from pgx.cyrius_runner import _get_cpic_recommendation

        rec = _get_cpic_recommendation("PM")
        assert "codeine" in rec.lower(), (
            f"PM recommendation must mention codeine: {rec}"
        )

    def test_um_recommendation_warns_toxicity(self) -> None:
        """UM recommendation warns about toxicity risk.

        UM patients risk morphine toxicity from standard codeine doses.
        """
        from pgx.cyrius_runner import _get_cpic_recommendation

        rec = _get_cpic_recommendation("UM")
        assert any(
            word in rec.lower()
            for word in ("toxic", "alternative", "avoid")
        ), f"UM recommendation should warn about toxicity: {rec}"

    def test_unknown_phenotype_returns_fallback(self) -> None:
        """Unknown phenotype returns CPIC guidelines fallback."""
        from pgx.cyrius_runner import _get_cpic_recommendation

        rec = _get_cpic_recommendation("UNKNOWN")
        assert "cpic" in rec.lower(), (
            "Unknown phenotype should return CPIC fallback message"
        )
