"""
MANE Select transcript utilities for clinical variant reporting.

Background
----------
MANE (Matched Annotation from NCBI and EMBL-EBI) is a collaboration between
NCBI and EMBL-EBI that defines a single default transcript per human protein-
coding gene that is:
  - Identical between RefSeq and Ensembl/GENCODE annotation
  - Present on the GRCh38/hg38 reference assembly
  - Supported by strong experimental evidence

Two tiers exist:
  - **MANE Select**: one transcript per gene; the default for clinical reporting
  - **MANE Plus Clinical**: additional transcripts relevant to disease that
    differ from MANE Select

Reference
---------
Morales J, et al. "A joint NCBI and EMBL-EBI transcript set for clinical
genomics and research." Nature. 2022;604:310–315. PMID:35356062.
DOI:10.1038/s41586-022-04558-8

ACGS 2024 §5 — PVS1 Strength Adjustment
-----------------------------------------
Per ACGS Best Practice Guidelines 2024 v1.2 §5, PVS1 strength should be
reduced by one level when the variant is NOT on the MANE Select transcript,
because clinical significance of alternative transcripts is less certain:

  PVS1 (Very Strong) → PS1 (Strong)
  PS1  (Strong)      → PM1 (Moderate)
  PM1  (Moderate)    → PP1 (Supporting)
  PP1  (Supporting)  → no contribution

This function is used by the bayesacmg rule engine to apply this adjustment
automatically.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Mapping of ACMG evidence strength → one-level-reduced strength
# Used by adjust_pvs1_for_mane() per ACGS 2024 §5
_STRENGTH_DOWNGRADE: dict[str, str] = {
    "PVS1": "PS1",
    "PS1": "PM1",
    "PM1": "PP1",
    "PP1": "no_contribution",
    # Handle lowercase/alternative forms
    "very_strong": "strong",
    "strong": "moderate",
    "moderate": "supporting",
    "supporting": "no_contribution",
}

# Hardcoded MANE Select transcript registry for key disease genes.
# In production this should be loaded from the MANE summary file:
#   ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/
#   MANE.GRCh38.v1.3.summary.txt.gz
# Format: gene_symbol → (RefSeq transcript ID, Ensembl transcript ID)
_MANE_SELECT_REGISTRY: dict[str, tuple[str, str]] = {
    "BRCA1": ("NM_007294.4", "ENST00000357654.9"),
    "BRCA2": ("NM_000059.4", "ENST00000380152.8"),
    "TP53": ("NM_000546.6", "ENST00000269305.9"),
    "PTEN": ("NM_000314.8", "ENST00000371953.8"),
    "MLH1": ("NM_000249.4", "ENST00000231790.8"),
    "MSH2": ("NM_000251.3", "ENST00000233146.7"),
    "MSH6": ("NM_000179.3", "ENST00000234420.10"),
    "PMS2": ("NM_000535.7", "ENST00000265849.12"),
    "APC": ("NM_000038.6", "ENST00000257430.9"),
    "RB1": ("NM_000321.3", "ENST00000267163.8"),
    "VHL": ("NM_000551.4", "ENST00000256474.3"),
    "NF1": ("NM_000267.3", "ENST00000358273.9"),
    "NF2": ("NM_000268.4", "ENST00000338641.8"),
    "CFTR": ("NM_000492.4", "ENST00000003084.11"),
    "HTT": ("NM_002111.8", "ENST00000355072.9"),
    "FMR1": ("NM_002024.6", "ENST00000370475.8"),
    "DMD": ("NM_004006.3", "ENST00000357033.9"),
    "LDLR": ("NM_000527.5", "ENST00000558518.6"),
    "PCSK9": ("NM_174936.4", "ENST00000302118.5"),
    "APOB": ("NM_000384.3", "ENST00000233242.4"),
    "KCNQ1": ("NM_000218.3", "ENST00000155840.11"),
    "KCNH2": ("NM_000238.4", "ENST00000262186.5"),
    "SCN5A": ("NM_198056.3", "ENST00000413689.6"),
    "MYH7": ("NM_000257.4", "ENST00000355349.4"),
    "MYBPC3": ("NM_000256.3", "ENST00000545968.6"),
    "PKD1": ("NM_001009944.3", "ENST00000262304.9"),
    "PKD2": ("NM_000297.4", "ENST00000237596.7"),
    "HBB": ("NM_000518.5", "ENST00000335295.4"),
    "HBA1": ("NM_000558.5", "ENST00000320868.10"),
    "HBA2": ("NM_000517.6", "ENST00000251595.9"),
    "G6PD": ("NM_001042351.3", "ENST00000393562.8"),
    "HEXA": ("NM_000520.6", "ENST00000268124.9"),
    "GBA": ("NM_000157.4", "ENST00000368373.8"),
    "LRRK2": ("NM_198578.4", "ENST00000298910.12"),
    "SNCA": ("NM_000345.4", "ENST00000394991.8"),
    "APP": ("NM_000484.4", "ENST00000346798.8"),
    "PSEN1": ("NM_000021.4", "ENST00000324501.9"),
    "PSEN2": ("NM_000447.3", "ENST00000366783.5"),
    "MECP2": ("NM_004992.4", "ENST00000303391.9"),
    "FBN1": ("NM_000138.5", "ENST00000316623.10"),
    "TGFBR1": ("NM_004612.4", "ENST00000374994.6"),
    "TGFBR2": ("NM_001024847.3", "ENST00000359013.4"),
}


def is_mane_select(transcript_id: str) -> bool:
    """Determine whether a transcript ID corresponds to a MANE Select transcript.

    Checks both RefSeq (NM_…) and Ensembl (ENST…) identifiers against the
    bundled MANE Select registry.  Version suffixes (e.g. ".4") are ignored
    to allow minor version flexibility.

    Args:
        transcript_id: RefSeq (e.g. "NM_007294.4") or Ensembl (e.g.
            "ENST00000357654.9") transcript identifier.

    Returns:
        True if the transcript is a known MANE Select transcript.

    Examples:
        >>> is_mane_select("NM_007294.4")
        True
        >>> is_mane_select("ENST00000357654.9")
        True
        >>> is_mane_select("NM_999999.1")
        False
    """
    # Strip version number for fuzzy matching
    bare_id = transcript_id.split(".")[0]

    for refseq_id, ensembl_id in _MANE_SELECT_REGISTRY.values():
        if bare_id == refseq_id.split(".")[0] or bare_id == ensembl_id.split(".")[0]:
            return True

    # Exact match (with version) as fallback
    for refseq_id, ensembl_id in _MANE_SELECT_REGISTRY.values():
        if transcript_id == refseq_id or transcript_id == ensembl_id:
            return True

    return False


def get_mane_select_for_gene(gene_symbol: str) -> Optional[str]:
    """Return the MANE Select RefSeq transcript ID for a gene symbol.

    Args:
        gene_symbol: HGNC gene symbol (case-sensitive), e.g. "BRCA1".

    Returns:
        RefSeq transcript ID (e.g. "NM_007294.4") if the gene is in the
        registry, otherwise None.

    Examples:
        >>> get_mane_select_for_gene("BRCA1")
        'NM_007294.4'
        >>> get_mane_select_for_gene("UNKNOWN")
        None
    """
    entry = _MANE_SELECT_REGISTRY.get(gene_symbol)
    return entry[0] if entry else None


def get_mane_select_ensembl_for_gene(gene_symbol: str) -> Optional[str]:
    """Return the MANE Select Ensembl transcript ID for a gene symbol.

    Args:
        gene_symbol: HGNC gene symbol, e.g. "BRCA2".

    Returns:
        Ensembl transcript ID (e.g. "ENST00000380152.8") or None.

    Examples:
        >>> get_mane_select_ensembl_for_gene("BRCA2")
        'ENST00000380152.8'
    """
    entry = _MANE_SELECT_REGISTRY.get(gene_symbol)
    return entry[1] if entry else None


def adjust_pvs1_for_mane(pvs1_strength: str, is_mane: bool) -> str:
    """Reduce PVS1 strength by one level when variant is not on MANE Select.

    Per ACGS 2024 §5, PVS1 strength applied to a non-MANE Select transcript
    is uncertain because the clinical relevance of that transcript has not been
    independently validated.  The strength is therefore reduced by one level:

    ============  ==============  ==================
    Input         is_mane=True    is_mane=False
    ============  ==============  ==================
    PVS1          PVS1            PS1
    PS1           PS1             PM1
    PM1           PM1             PP1
    PP1           PP1             no_contribution
    ============  ==============  ==================

    Args:
        pvs1_strength: ACMG strength label, one of "PVS1", "PS1", "PM1",
            "PP1" (case-sensitive). Also accepts lowercase "very_strong",
            "strong", "moderate", "supporting".
        is_mane: True if the variant consequence falls on the MANE Select
            transcript for this gene.

    Returns:
        Possibly-adjusted strength label.  The input is returned unchanged
        when ``is_mane`` is True.  "no_contribution" is returned when the
        input is "PP1" (supporting) and is_mane is False.

    Examples:
        >>> adjust_pvs1_for_mane("PVS1", True)
        'PVS1'
        >>> adjust_pvs1_for_mane("PVS1", False)
        'PS1'
        >>> adjust_pvs1_for_mane("PP1", False)
        'no_contribution'
    """
    if is_mane:
        # No adjustment needed — full strength applies
        return pvs1_strength

    downgraded = _STRENGTH_DOWNGRADE.get(pvs1_strength)
    if downgraded is None:
        logger.warning(
            "Unknown PVS1 strength label '%s'; returning unchanged", pvs1_strength
        )
        return pvs1_strength

    logger.debug(
        "PVS1 strength adjusted %s → %s (non-MANE Select, ACGS 2024 §5)",
        pvs1_strength,
        downgraded,
    )
    return downgraded


@lru_cache(maxsize=1)
def load_mane_summary(summary_path: str) -> dict[str, tuple[str, str]]:
    """Load MANE Select mapping from the official NCBI summary file.

    The summary file is available from:
    ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/
    MANE.GRCh38.v1.3.summary.txt.gz

    Expected columns (tab-separated, gzipped):
    ``#NCBI_GeneID  Ensembl_Gene  GeneSymbol  name  RefSeq_nuc  RefSeq_prot
    Ensembl_nuc  Ensembl_prot  MANE_status  GRCh38_chr  chr_start  chr_end
    chr_strand``

    Args:
        summary_path: Path to the (optionally gzipped) MANE summary file.

    Returns:
        Mapping of gene_symbol → (refseq_transcript_id, ensembl_transcript_id).
        Returns an empty dict if the file cannot be read.
    """
    import gzip

    registry: dict[str, tuple[str, str]] = {}
    path = Path(summary_path)

    if not path.exists():
        logger.error("MANE summary file not found: %s", summary_path)
        return registry

    opener = gzip.open if path.suffix == ".gz" else open

    try:
        with opener(summary_path, "rt", encoding="utf-8") as fh:  # type: ignore[call-overload]
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                gene_symbol = parts[2]
                refseq_nuc = parts[4]
                ensembl_nuc = parts[6]
                mane_status = parts[8]
                # Only register MANE Select (not Plus Clinical)
                if mane_status == "MANE Select":
                    registry[gene_symbol] = (refseq_nuc, ensembl_nuc)
    except OSError as exc:
        logger.error("Failed to load MANE summary: %s", exc)

    logger.info("Loaded %d MANE Select transcripts from %s", len(registry), summary_path)
    return registry
