"""
multi_ancestry.somalier_runner
================================
Run somalier relate and ancestry inference on a cohort of samples.

somalier uses fast site-based genotyping at ~17,000 informative positions
across the genome to:
    1. Check relatedness between samples (identify unexpected duplicates).
    2. Infer genetic ancestry via PCA projection against reference populations.
    3. Identify sex chromosomal anomalies (X/Y copy number).

somalier ancestry:
    Compares sample PCs against 1000 Genomes + HGDP reference panel.
    Reports predicted population: AFR, AMR, EAS, EUR, MID, SAS, and
    admixed fractions for multi-ancestry individuals.

Reference sites:
    GRCh38: https://github.com/brentp/somalier/releases — sites.hg38.vcf.gz

References:
    Pedersen et al. 2020 Genome Biology PMID:32620139 (somalier).
    somalier docs: https://github.com/brentp/somalier
"""

from __future__ import annotations

import json
import logging
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class SomalierAncestryResult:
    """Ancestry inference result from somalier.

    Attributes:
        sample_id: Sample identifier.
        predicted_ancestry: Top predicted ancestry label
            (``"AFR"``, ``"AMR"``, ``"EAS"``, ``"EUR"``, ``"MID"``, ``"SAS"``).
        predicted_ancestry_p: Probability of predicted ancestry (0–1).
        ancestry_fractions: Dict of population label → fraction (sum ≈ 1).
        is_admixed: True if no single ancestry fraction ≥ 0.8.
        pc1: First principal component value.
        pc2: Second principal component value.
        pc3: Third principal component value.
    """

    sample_id: str
    predicted_ancestry: str
    predicted_ancestry_p: float
    ancestry_fractions: dict[str, float] = field(default_factory=dict)
    is_admixed: bool = False
    pc1: float = 0.0
    pc2: float = 0.0
    pc3: float = 0.0


@dataclass
class SomalierRelatednessResult:
    """Pairwise relatedness result from somalier relate.

    Attributes:
        sample_a: First sample ID.
        sample_b: Second sample ID.
        relatedness: Kinship coefficient (0=unrelated, 0.5=parent-child).
        ibs0: IBS0 count (incompatible with parent-child).
        ibs2: IBS2 count.
        n_sites: Number of informative sites used.
        expected_relationship: Expected relationship string
            (``"unrelated"``, ``"parent-child"``, ``"sibling"``, etc.).
    """

    sample_a: str
    sample_b: str
    relatedness: float
    ibs0: int
    ibs2: int
    n_sites: int
    expected_relationship: str = "unrelated"


# ---------------------------------------------------------------------------
# Runner functions
# ---------------------------------------------------------------------------


def run_somalier_extract(
    bam_path: Path,
    ref_fasta: Path,
    sites_vcf: Path,
    output_dir: Path,
) -> Path:
    """Extract somalier sites from a BAM file.

    Args:
        bam_path: Path to coordinate-sorted, indexed BAM file.
        ref_fasta: Path to GRCh38 reference FASTA.
        sites_vcf: Path to somalier sites VCF (sites.hg38.vcf.gz).
        output_dir: Directory for somalier extract output.

    Returns:
        Path to the generated ``.somalier`` file.

    Raises:
        RuntimeError: If somalier extract fails.
        FileNotFoundError: If input files do not exist.
    """
    for f in [bam_path, ref_fasta, sites_vcf]:
        if not f.exists():
            raise FileNotFoundError(f"Required file not found: {f}")

    output_dir.mkdir(parents=True, exist_ok=True)
    sample_id = bam_path.stem

    cmd = [
        "somalier", "extract",
        "--sites", str(sites_vcf),
        "--fasta", str(ref_fasta),
        "--out-dir", str(output_dir),
        str(bam_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"somalier extract failed for {sample_id}: {result.stderr}")

    somalier_file = output_dir / f"{sample_id}.somalier"
    if not somalier_file.exists():
        raise RuntimeError(
            f"somalier extract did not produce expected output: {somalier_file}"
        )

    logger.info("somalier extract complete: %s", somalier_file)
    return somalier_file


def run_somalier_relate(
    somalier_files: list[Path],
    output_dir: Path,
    ped_file: Path | None = None,
) -> dict[str, Any]:
    """Run somalier relate on a cohort of extracted somalier files.

    Args:
        somalier_files: List of ``.somalier`` files from somalier extract.
        output_dir: Directory for somalier relate output.
        ped_file: Optional PED file for expected relatedness validation.

    Returns:
        Dict with ``"samples"`` and ``"pairs"`` keys from somalier TSV output.

    Raises:
        RuntimeError: If somalier relate fails.
        ValueError: If somalier_files is empty.
    """
    if not somalier_files:
        raise ValueError("No somalier files provided to run_somalier_relate().")

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "somalier", "relate",
        "--output-prefix", str(output_dir / "somalier"),
    ]
    if ped_file and ped_file.exists():
        cmd += ["--ped", str(ped_file)]
    cmd += [str(f) for f in somalier_files]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"somalier relate failed: {result.stderr}")

    logger.info("somalier relate complete in %s", output_dir)
    return _parse_somalier_relate_output(output_dir)


def run_somalier_ancestry(
    somalier_files: list[Path],
    reference_panel_dir: Path,
    output_dir: Path,
) -> list[SomalierAncestryResult]:
    """Run somalier ancestry inference on a cohort.

    Args:
        somalier_files: List of ``.somalier`` files.
        reference_panel_dir: Directory containing 1000G+HGDP reference panel
            somalier files and labels.
        output_dir: Directory for ancestry output.

    Returns:
        List of SomalierAncestryResult objects, one per sample.

    Raises:
        RuntimeError: If somalier ancestry fails.
    """
    if not somalier_files:
        raise ValueError("No somalier files provided to run_somalier_ancestry().")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Find reference labels file
    labels_file = reference_panel_dir / "1kg+hgdp.somalier-ancestry.tsv"
    ref_somalier_glob = list(reference_panel_dir.glob("*.somalier"))

    cmd = [
        "somalier", "ancestry",
        "--labels", str(labels_file) if labels_file.exists() else "",
        "--output-prefix", str(output_dir / "ancestry"),
    ]
    cmd += [str(f) for f in ref_somalier_glob]
    cmd += ["++"] + [str(f) for f in somalier_files]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"somalier ancestry failed: {result.stderr}")

    return _parse_somalier_ancestry_output(output_dir / "ancestry.somalier-ancestry.tsv")


def _parse_somalier_relate_output(output_dir: Path) -> dict[str, Any]:
    """Parse somalier relate TSV output files.

    Args:
        output_dir: Directory containing somalier relate output TSV files.

    Returns:
        Dict with ``"samples"`` and ``"pairs"`` lists.
    """
    result: dict[str, Any] = {"samples": [], "pairs": []}

    samples_file = output_dir / "somalier.samples.tsv"
    pairs_file = output_dir / "somalier.pairs.tsv"

    if samples_file.exists():
        with samples_file.open() as fh:
            header = None
            for line in fh:
                if header is None:
                    header = line.strip().split("\t")
                    continue
                values = line.strip().split("\t")
                if len(values) >= len(header or []):
                    result["samples"].append(dict(zip(header or [], values)))

    if pairs_file.exists():
        with pairs_file.open() as fh:
            header = None
            for line in fh:
                if header is None:
                    header = line.strip().split("\t")
                    continue
                values = line.strip().split("\t")
                if len(values) >= 4:
                    result["pairs"].append({
                        "sample_a": values[0],
                        "sample_b": values[1],
                        "relatedness": float(values[2]) if values[2] else 0.0,
                        "ibs0": int(values[3]) if values[3].isdigit() else 0,
                    })

    return result


def _parse_somalier_ancestry_output(tsv_path: Path) -> list[SomalierAncestryResult]:
    """Parse somalier ancestry TSV output.

    Args:
        tsv_path: Path to somalier ancestry output TSV.

    Returns:
        List of SomalierAncestryResult objects.
    """
    results: list[SomalierAncestryResult] = []

    if not tsv_path.exists():
        logger.warning("somalier ancestry output not found: %s", tsv_path)
        return results

    with tsv_path.open() as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if header is None:
                header = line.split("\t")
                continue
            if not header:
                continue
            values = line.split("\t")
            row = dict(zip(header, values))

            sample_id = row.get("sample_id", row.get("#sample_id", ""))
            predicted = row.get("predicted_ancestry", "")
            predicted_p = float(row.get("given_ancestry", row.get("predicted_ancestry_p", 0.0)) or 0.0)

            # Extract ancestry fractions
            pop_labels = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
            fractions = {
                pop: float(row.get(pop, 0.0) or 0.0)
                for pop in pop_labels
                if pop in row
            }

            max_fraction = max(fractions.values()) if fractions else 0.0
            is_admixed = max_fraction < 0.8

            results.append(SomalierAncestryResult(
                sample_id=sample_id,
                predicted_ancestry=predicted,
                predicted_ancestry_p=predicted_p,
                ancestry_fractions=fractions,
                is_admixed=is_admixed,
                pc1=float(row.get("PC1", 0.0) or 0.0),
                pc2=float(row.get("PC2", 0.0) or 0.0),
                pc3=float(row.get("PC3", 0.0) or 0.0),
            ))

    logger.info("Parsed somalier ancestry for %d samples", len(results))
    return results
