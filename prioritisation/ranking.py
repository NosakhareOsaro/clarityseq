"""
prioritisation.ranking
========================
Composite variant ranking combining ACMG classification, HPO phenotype
similarity, inheritance mode, and gene panel membership.

Ranking formula:
    composite_score = (
        acmg_weight * acmg_score
        + hpo_weight * hpo_score
        + inheritance_weight * inheritance_bonus
        + panel_weight * panel_bonus
    )

Score components:
    acmg_score: Normalised ACMG/AMP classification score (P=1.0, LP=0.8,
        VUS=0.5, LB=0.2, B=0.0).
    hpo_score: HPO phenotypic similarity score (0–1); Jaccard index or
        Exomiser combined score.
    inheritance_bonus: 1.0 if variant passes inheritance mode filter; 0.0 if not.
    panel_bonus: 1.0 if gene is on the clinical gene panel; 0.0 if not.

Default weights (tuned against clinical datasets):
    acmg_weight=0.4, hpo_weight=0.3, inheritance_weight=0.2, panel_weight=0.1

References:
    Richards et al. 2015 PMID:25741868 (ACMG/AMP classification).
    Robinson et al. 2023 PMID:37604970 (Exomiser composite scoring).
    ACGS 2024 v1.2 §5 Table 2.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# ACMG classification → score mapping
# ---------------------------------------------------------------------------

_ACMG_SCORES: dict[str, float] = {
    "Pathogenic": 1.0,
    "Likely_Pathogenic": 0.8,
    "VUS": 0.5,
    "Likely_Benign": 0.2,
    "Benign": 0.0,
}

# Aliases for different ACMG classification string formats
_ACMG_ALIASES: dict[str, str] = {
    "P": "Pathogenic",
    "PATHOGENIC": "Pathogenic",
    "LP": "Likely_Pathogenic",
    "LIKELY_PATHOGENIC": "Likely_Pathogenic",
    "LIKELY PATHOGENIC": "Likely_Pathogenic",
    "VUS": "VUS",
    "UNCERTAIN SIGNIFICANCE": "VUS",
    "VARIANT OF UNCERTAIN SIGNIFICANCE": "VUS",
    "LB": "Likely_Benign",
    "LIKELY_BENIGN": "Likely_Benign",
    "LIKELY BENIGN": "Likely_Benign",
    "B": "Benign",
    "BENIGN": "Benign",
}


def _normalise_acmg(acmg_class: str) -> str:
    """Normalise an ACMG classification string to a canonical form.

    Args:
        acmg_class: ACMG classification string in any format.

    Returns:
        Canonical form (one of Pathogenic, Likely_Pathogenic, VUS,
        Likely_Benign, Benign).
    """
    upper = acmg_class.upper().strip()
    return _ACMG_ALIASES.get(upper, "VUS")


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class RankedVariant:
    """A candidate variant with composite prioritisation score.

    Attributes:
        rank: Overall rank among all candidates (1 = highest priority).
        gene_symbol: HGNC gene symbol.
        chrom: Chromosome.
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
        acmg_class: Normalised ACMG classification.
        acmg_score: ACMG numeric score (0–1).
        hpo_score: HPO phenotypic similarity score (0–1).
        inheritance_bonus: 1.0 if passes inheritance filter; 0.0 otherwise.
        panel_bonus: 1.0 if gene is in clinical gene panel; 0.0 otherwise.
        composite_score: Weighted composite score (0–1).
        evidence_summary: Human-readable evidence summary.
        extra: Additional annotations dict.
    """

    rank: int
    gene_symbol: str
    chrom: str
    pos: int
    ref: str
    alt: str
    acmg_class: str
    acmg_score: float
    hpo_score: float
    inheritance_bonus: float
    panel_bonus: float
    composite_score: float
    evidence_summary: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Ranking function
# ---------------------------------------------------------------------------


def compute_composite_score(
    acmg_class: str,
    hpo_score: float,
    passes_inheritance: bool,
    in_panel: bool,
    acmg_weight: float = 0.4,
    hpo_weight: float = 0.3,
    inheritance_weight: float = 0.2,
    panel_weight: float = 0.1,
) -> tuple[float, float, float, float, float]:
    """Compute the composite prioritisation score for a variant.

    Args:
        acmg_class: ACMG/AMP classification string.
        hpo_score: HPO phenotypic similarity score (0–1).
        passes_inheritance: True if variant passes inheritance mode filter.
        in_panel: True if gene is on the clinical gene panel.
        acmg_weight: Weight for ACMG score (default 0.4).
        hpo_weight: Weight for HPO score (default 0.3).
        inheritance_weight: Weight for inheritance filter (default 0.2).
        panel_weight: Weight for panel membership (default 0.1).

    Returns:
        Tuple of (composite_score, acmg_score, hpo_score,
            inheritance_bonus, panel_bonus).
    """
    canonical = _normalise_acmg(acmg_class)
    acmg_score = _ACMG_SCORES.get(canonical, 0.5)
    inheritance_bonus = 1.0 if passes_inheritance else 0.0
    panel_bonus = 1.0 if in_panel else 0.0

    composite = (
        acmg_weight * acmg_score
        + hpo_weight * hpo_score
        + inheritance_weight * inheritance_bonus
        + panel_weight * panel_bonus
    )
    return composite, acmg_score, hpo_score, inheritance_bonus, panel_bonus


def rank_variants(
    variants: list[dict[str, Any]],
    hpo_gene_scores: dict[str, float] | None = None,
    passing_inheritance_genes: set[str] | None = None,
    panel_genes: set[str] | None = None,
    acmg_weight: float = 0.4,
    hpo_weight: float = 0.3,
    inheritance_weight: float = 0.2,
    panel_weight: float = 0.1,
) -> list[RankedVariant]:
    """Rank candidate variants by composite score.

    Args:
        variants: List of variant dicts with keys:
            ``chrom``, ``pos``, ``ref``, ``alt``, ``gene_symbol``,
            ``acmg_class`` (plus any extra keys in ``extra``).
        hpo_gene_scores: Dict of gene_symbol → HPO similarity score (0–1).
            None uses 0.0 for all genes.
        passing_inheritance_genes: Set of gene symbols that pass the
            inheritance mode filter.  None skips inheritance filtering.
        panel_genes: Set of gene symbols on the clinical panel.
            None treats all genes as not on panel (panel_bonus=0.0).
        acmg_weight: Weight for ACMG score component (default 0.4).
        hpo_weight: Weight for HPO score component (default 0.3).
        inheritance_weight: Weight for inheritance filter (default 0.2).
        panel_weight: Weight for panel membership (default 0.1).

    Returns:
        List of RankedVariant objects sorted by composite_score descending.
        Rank 1 = highest priority candidate.

    References:
        Richards et al. 2015 PMID:25741868.
        Robinson et al. 2023 PMID:37604970.
        ACGS 2024 v1.2 §5.
    """
    hpo_scores = hpo_gene_scores or {}
    inheritance_set = passing_inheritance_genes or set()
    panels = panel_genes or set()

    ranked: list[RankedVariant] = []

    for var in variants:
        gene = var.get("gene_symbol", "")
        acmg_class = var.get("acmg_class", "VUS")

        hpo = hpo_scores.get(gene, 0.0)
        passes_inh = gene in inheritance_set if passing_inheritance_genes is not None else True
        in_panel = gene in panels if panel_genes is not None else False

        (composite, acmg_s, hpo_s, inh_b, panel_b) = compute_composite_score(
            acmg_class=acmg_class,
            hpo_score=hpo,
            passes_inheritance=passes_inh,
            in_panel=in_panel,
            acmg_weight=acmg_weight,
            hpo_weight=hpo_weight,
            inheritance_weight=inheritance_weight,
            panel_weight=panel_weight,
        )

        evidence_parts = [
            f"ACMG={acmg_class}({acmg_s:.2f})",
            f"HPO={hpo_s:.3f}",
            f"Inheritance={'pass' if inh_b > 0 else 'fail'}",
            f"Panel={'yes' if panel_b > 0 else 'no'}",
        ]

        ranked.append(
            RankedVariant(
                rank=0,  # placeholder; assigned after sorting
                gene_symbol=gene,
                chrom=var.get("chrom", ""),
                pos=int(var.get("pos", 0)),
                ref=var.get("ref", ""),
                alt=var.get("alt", ""),
                acmg_class=acmg_class,
                acmg_score=acmg_s,
                hpo_score=hpo_s,
                inheritance_bonus=inh_b,
                panel_bonus=panel_b,
                composite_score=composite,
                evidence_summary=" | ".join(evidence_parts),
                extra={k: v for k, v in var.items()
                       if k not in {"chrom", "pos", "ref", "alt", "gene_symbol", "acmg_class"}},
            )
        )

    # Sort by composite score (descending); ties broken by ACMG score then HPO
    ranked.sort(
        key=lambda r: (r.composite_score, r.acmg_score, r.hpo_score),
        reverse=True,
    )

    # Assign ranks
    for i, rv in enumerate(ranked):
        rv.rank = i + 1

    logger.info("Ranked %d candidate variants.", len(ranked))
    return ranked
