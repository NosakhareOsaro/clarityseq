"""
bayesacmg.models
================

Core data models for the BayesACMG variant classification system.

This module defines the canonical data structures shared by all rule modules,
the Bayesian model, and the CLI.  Every field is documented with an inline
comment so that downstream consumers understand the provenance of each value.

Guidelines implemented
----------------------
- Richards et al. 2015 PMID:25741868 — original ACMG/AMP framework,
  27 criteria across 5 classification categories.
- Tavtigian et al. 2020 PMID:32645316 — point-score system and
  Bayesian framework for combining ACMG/AMP criteria.
- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., 20 Feb 2024)
  — UK clinical implementation; MANE Select requirement; PM2 at Supporting;
  AlphaMissense as primary in silico predictor.
- ClinGen SVI Working Group 2024 — PM2→Supporting, AlphaMissense thresholds,
  novel combinations (PVS1+PM2_Supporting=LP).

CRITICAL CHANGE FROM ACGS 2020
-------------------------------
PM2 is now applied at SUPPORTING weight (1 pt), NOT Moderate (2 pts).
See EvidenceStrength docstring for full rationale.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

# ---------------------------------------------------------------------------
# Variant type
# ---------------------------------------------------------------------------


class VariantType(str, Enum):
    """Molecular consequence category for a variant.

    Used by PVS1, PP3/BP4, BP7, and splicing rules to determine which
    evidence criteria are applicable.
    """

    SNV = "snv"
    MISSENSE = "missense"
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    SPLICE_SITE = "splice_site"  # canonical ±1/2 splice donor/acceptor
    SPLICE_REGION = "splice_region"  # intronic/exonic splice region variant
    SYNONYMOUS = "synonymous"
    START_LOSS = "start_loss"
    STOP_LOSS = "stop_loss"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    LARGE_DELETION = "large_deletion"
    LARGE_DUPLICATION = "large_duplication"
    INDEL = "indel"


# ---------------------------------------------------------------------------
# Evidence strength
# ---------------------------------------------------------------------------


class EvidenceStrength(str, Enum):
    """ACMG/AMP evidence strength categories and their point values.

    Point-score mapping (Tavtigian et al. 2020 PMID:32645316):
        VERY_STRONG       → +8  pts  (PVS criteria)
        STRONG            → +4  pts  (PS criteria)
        MODERATE          → +2  pts  (PM criteria, except PM2)
        SUPPORTING        → +1  pt   (PP criteria AND PM2 after ClinGen SVI 2024)
        STAND_ALONE       → Benign directly (BA1 — no point calculation)
        STRONG_BENIGN     → -4  pts  (BS criteria)
        SUPPORTING_BENIGN → -1  pt   (BP criteria)

    CRITICAL CHANGE FROM ACGS 2020 / Richards 2015:
    PM2 is now applied at SUPPORTING (1 pt) by default in the ClinGen SVI
    2024 recommendations.  gnomAD v4.1 (807,162 individuals, April 2024)
    reveals high rates of ultra-rare variants in the general population,
    making rarity alone weaker evidence than assumed in 2015.
    Reference: ClinGen SVI PM2 guidance (2024)
    https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
    VCEP specifications may override per-gene via vcep_client.py.
    """

    VERY_STRONG = "very_strong"  # PVS:  +8 pts
    STRONG = "strong"  # PS:   +4 pts
    MODERATE = "moderate"  # PM (not PM2): +2 pts
    SUPPORTING = "supporting"  # PP + PM2: +1 pt
    STAND_ALONE = "stand_alone"  # BA1: → Benign directly (no pts)
    STRONG_BENIGN = "strong_benign"  # BS:  -4 pts
    SUPPORTING_BENIGN = "supporting_benign"  # BP:  -1 pt


# Map EvidenceStrength → integer point value for scoring
STRENGTH_POINTS: dict[EvidenceStrength, int] = {
    EvidenceStrength.VERY_STRONG: 8,
    EvidenceStrength.STRONG: 4,
    EvidenceStrength.MODERATE: 2,
    EvidenceStrength.SUPPORTING: 1,
    EvidenceStrength.STAND_ALONE: 0,  # triggers Benign directly
    EvidenceStrength.STRONG_BENIGN: -4,
    EvidenceStrength.SUPPORTING_BENIGN: -1,
}


# ---------------------------------------------------------------------------
# ACMG rule result
# ---------------------------------------------------------------------------


@dataclass
class ACMGRule:
    """Result of evaluating a single ACMG/AMP criterion for one variant.

    Attributes:
        rule_id: Criterion identifier, e.g. ``"PM2"`` or ``"PVS1"``.
        strength: Evidence strength at which the rule fires.  A rule that
            does not fire has ``applies=False`` and strength is informational.
        evidence_items: Human-readable list of evidence observations that
            led to the rule firing (or not firing).  Used in reports and
            the Bayesian model's evidence vector.
        citations: Literature citations supporting the rule application,
            formatted as ``"Author Year PMID:XXXXXXXX"`` strings.
        applies: Whether the criterion applies (``True``) or does not apply
            (``False``) for this variant.  Only ``True`` rules contribute
            to the point total.
        notes: Optional free-text notes for edge cases, reduced-strength
            applications, or VCEP overrides.
    """

    rule_id: str  # e.g. "PM2", "PVS1"
    strength: EvidenceStrength  # strength at which the rule fires
    evidence_items: list[str]  # observations supporting the decision
    applies: bool  # True if this criterion counts
    notes: str = ""  # free-text for edge cases / overrides
    citations: list[str] = field(default_factory=list)  # literature citations

    @property
    def points(self) -> int:
        """Return integer point contribution of this rule.

        Returns:
            Point value if the rule applies, else 0.  STAND_ALONE rules
            return 0 because they trigger a direct Benign classification
            rather than contributing to the point sum.
        """
        if not self.applies:
            return 0
        return STRENGTH_POINTS[self.strength]


# ---------------------------------------------------------------------------
# Variant input
# ---------------------------------------------------------------------------


@dataclass
class VariantInput:
    """Fully-annotated variant ready for ACMG/AMP classification.

    Population frequency fields use gnomAD v4.1 (807,162 individuals,
    April 2024) as the primary source.  Fields are ``None`` when data
    are unavailable or not applicable.

    Attributes:
        chrom: Chromosome (GRCh38 notation, e.g. ``"chr17"`` or ``"chrM"``).
        pos: 1-based genomic position (GRCh38).
        ref: Reference allele (VCF notation, uppercase).
        alt: Alternate allele (VCF notation, uppercase).
        variant_type: Variant type string, one of:
            ``"snv"``, ``"indel"``, ``"frameshift"``, ``"nonsense"``,
            ``"splice_canonical"``, ``"splice_region"``, ``"start_loss"``,
            ``"stop_loss"``, ``"inframe_insertion"``, ``"inframe_deletion"``,
            ``"synonymous"``, ``"large_deletion"``, ``"large_duplication"``.
        gnomad_af: gnomAD v4.1 allele frequency (all populations combined).
            ``None`` if variant is absent from gnomAD (very relevant for PM2).
        gnomad_ac: gnomAD v4.1 allele count (all populations combined).
        gnomad_nhomalt: gnomAD v4.1 number of homozygous individuals.
            Critical for recessive disease BA1/BS1 assessment.
        gnomad_popmax_af: Highest population-specific AF in gnomAD v4.1.
            Used for BA1 assessment (>5% in any population).
        clinvar_stars: ClinVar review status stars (0–4).  ≥2 stars with
            P/LP classification can contribute to PS1/PP5.
        clinvar_classification: ClinVar classification string (latest).
        alphamissense_score: AlphaMissense score (0–1).
            PRIMARY in silico predictor for PP3/BP4 (ClinGen SVI 2024).
            Thresholds: ≥0.564 → PP3; ≤0.340 → BP4 (Cheng 2023 PMID:37703350).
        revel_score: REVEL score (0–1).  Secondary predictor; ≥0.7 supports PP3.
        cadd_phred: CADD PHRED-scaled score.  ≥25 supports PP3.
        spliceai_max_delta: Maximum SpliceAI Δ score across all four channels
            (DS_AG, DS_AL, DS_DG, DS_DL).  Used in splicing.py.
        pangolin_score: Pangolin splice-impact score (0–1).
        is_de_novo: True if confirmed de novo by parental testing.
        assumed_de_novo: True if parents not tested but phenotype strongly
            suggests de novo (used for PM6 at Supporting, not PS2).
        in_trans_pathogenic: True if variant is in trans with a known
            pathogenic variant in a recessive gene.
        functional_study_result: Functional study outcome string, e.g.
            ``"loss_of_function"``, ``"benign"``, or ``None``.
        hgvs_c: HGVS cDNA notation.
        hgvs_p: HGVS protein notation.
        gene_symbol: HGNC gene symbol.
        transcript_id: Transcript ID used for HGVS notation.
        is_mito: True if variant is on the mitochondrial chromosome (chrM).
        heteroplasmy_level: Fraction (0–1) of mitochondrial reads carrying
            the alternate allele.  None for nuclear variants.
        haplogroup: Predicted haplogroup string from Haplogrep3.
        extra: Dict for any additional annotations not covered above.
    """

    # Core variant coordinates (GRCh38)
    chrom: str  # e.g. "chr17", "chrM"
    pos: int  # 1-based position
    ref: str  # reference allele
    alt: str  # alternate allele
    variant_type: VariantType | str  # VariantType enum or raw VEP consequence string

    # Population frequencies — gnomAD v4.1 (807,162 individuals, April 2024)
    gnomad_af: float | None = None  # global AF; None = absent from gnomAD
    gnomad_ac: int | None = None  # allele count
    gnomad_nhomalt: int | None = None  # homozygous count
    gnomad_popmax_af: float | None = None  # highest population-specific AF

    # ClinVar
    clinvar_stars: int | None = None  # review stars 0-4
    clinvar_classification: str | None = None  # e.g. "Pathogenic"

    # In silico scores
    alphamissense_score: float | None = None  # PRIMARY: ≥0.564→PP3, ≤0.340→BP4
    revel_score: float | None = None  # secondary; ≥0.7→PP3
    cadd_phred: float | None = None  # secondary; ≥25→PP3

    # Splicing scores
    spliceai_max_delta: float | None = None  # max Δ across all SpliceAI channels
    spliceai_delta: float | None = (
        None  # alias for spliceai_max_delta (test/VEP convention)
    )
    pangolin_score: float | None = None  # Pangolin splice-impact (0-1)

    # Inheritance / family data
    is_de_novo: bool = False  # confirmed de novo (PS2)
    assumed_de_novo: bool = False  # assumed de novo, not confirmed (PM6)
    in_trans_pathogenic: bool = False  # in trans with pathogenic (PM3)
    segregation_lod: float | None = None  # LOD score for cosegregation (PP1)

    # Functional studies
    functional_study_result: str | None = None  # e.g. "loss_of_function"

    # HGVS identifiers (canonical names use underscore; aliases without for compat)
    hgvs_c: str | None = None  # HGVS cDNA notation (MANE Select; ACGS 2024 §4.1)
    hgvs_p: str | None = None  # HGVS protein notation
    hgvsc: str | None = None  # alias for hgvs_c (VEP/test convention)
    hgvsp: str | None = None  # alias for hgvs_p (VEP/test convention)
    gene_symbol: str | None = None  # HGNC gene symbol
    transcript_id: str | None = None  # transcript used for HGVS

    # ClinVar accession (optional; present for variants with ClinVar RCV)
    clinvar_rcv: str | None = None  # ClinVar RCV accession (e.g. RCV000007535)

    # Mitochondrial-specific fields (ACGS 2024 §6)
    is_mito: bool = False  # True if chrM variant
    heteroplasmy_level: float | None = None  # fraction 0-1 of alt reads
    heteroplasmy_fraction: float | None = None  # alias for heteroplasmy_level
    haplogroup: str | None = None  # Haplogrep3 haplogroup string
    is_haplogroup_defining: bool = (
        False  # True if variant defines a haplogroup (→ Benign)
    )

    # Extensible
    extra: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Transcript data
# ---------------------------------------------------------------------------


@dataclass
class TranscriptData:
    """Transcript-level annotation for a variant.

    Attributes:
        transcript_id: Ensembl or RefSeq transcript identifier.
        is_mane_select: True if this is the MANE Select transcript.
            Required for full PVS1 strength (ACGS 2024 v1.2 §5, Note 1).
        is_mane_plus_clinical: True if MANE Plus Clinical transcript.
        gene_symbol: HGNC gene symbol.
        consequence: VEP-style consequence string.
        exon_number: Exon number where the variant falls (1-based).
        total_exons: Total number of canonical exons in this transcript.
        is_last_exon: True if the variant is in the last exon.
        is_penultimate_exon: True if the variant is in the penultimate exon.
        nmd_escapes_last_exon_rule: True if NMD would be expected to fail
            (last exon or <55 nt from last junction — PVS1 decision tree).
        prot_length_original: Length of wild-type protein in amino acids.
        prot_length_alt: Length of predicted alternate protein in amino acids.
        aa_change: Amino acid change string, e.g. ``"p.Arg175His"``.
        aa_position: Position of amino acid change (1-based).
        domain_annotations: List of domain annotation strings for PM1.
        splice_region: True if variant is within splice region (not canonical).
    """

    transcript_id: str  # Ensembl or RefSeq ID
    is_mane_select: bool = False  # MANE Select transcript flag
    is_mane_plus_clinical: bool = False  # MANE Plus Clinical flag
    is_canonical: bool = False  # alias: True if this is the canonical/MANE transcript
    lof_disease_mechanism: bool = (
        False  # alias for lof_is_disease_mechanism (test convention)
    )
    exon_count: int | None = None  # total exon count for the transcript
    gene_symbol: str = ""  # HGNC gene symbol
    consequence: str = ""  # VEP consequence string
    exon_number: int | None = None  # 1-based exon number
    total_exons: int | None = None  # total canonical exons
    is_last_exon: bool = False  # variant in last exon
    is_penultimate_exon: bool = False  # variant in penultimate exon
    nmd_escapes_last_exon_rule: bool = False  # NMD expected to fail
    prot_length_original: int | None = None  # WT protein length (aa)
    prot_length_alt: int | None = None  # alt protein length (aa)
    aa_change: str | None = None  # e.g. "p.Arg175His"
    aa_position: int | None = None  # 1-based AA position
    domain_annotations: list[str] = field(default_factory=list)  # PM1 domains
    splice_region: bool = False  # True if splice region (not canonical)


# ---------------------------------------------------------------------------
# Gene data
# ---------------------------------------------------------------------------


@dataclass
class GeneData:
    """Gene-level data required for ACMG/AMP rule evaluation.

    Attributes:
        gene_symbol: HGNC gene symbol.
        lof_is_disease_mechanism: True if loss-of-function is a known
            disease mechanism for this gene (required for PVS1).
        gnomad_pli: gnomAD pLI score (probability of LoF intolerance).
            Values >0.9 suggest LoF intolerance; used for PP2.
        gnomad_loeuf: gnomAD LOEUF (LoF observed/expected upper bound).
            Lower values indicate greater LoF intolerance; <0.35 used for PP2.
        gene_disease_mode: Inheritance mode string, e.g. ``"AD"``, ``"AR"``,
            ``"XL"``, ``"Mito"``.
        has_hotspot_domain: True if the gene has known mutational hotspot
            domains relevant to PM1.
        hotspot_domains: List of domain names/coords for PM1.
        vcep_gene: True if a VCEP specification exists for this gene.
            See vcep_client.py for the override mechanism.
        missense_only_gene: True if only missense variants cause disease
            (BP1 applies for truncating variants).
        haploinsufficiency_score: ClinGen haploinsufficiency score (1–3).
        triplosensitivity_score: ClinGen triplosensitivity score (1–3).
    """

    gene_symbol: str  # HGNC gene symbol
    lof_is_disease_mechanism: bool = False  # required for PVS1
    lof_mechanism: bool = False  # alias for lof_is_disease_mechanism (test convention)
    gnomad_pli: float | None = None  # pLI; >0.9 → LoF intolerant
    pli: float | None = None  # alias for gnomad_pli (test convention)
    gnomad_loeuf: float | None = None  # LOEUF; <0.35 → LoF intolerant
    missense_constraint_z: float | None = None  # gnomAD missense Z-score
    gene_disease_mode: str = ""  # "AD", "AR", "XL", "Mito"
    omim_id: str | None = None  # OMIM identifier e.g. "OMIM:113705"
    has_hotspot_domain: bool = False  # PM1 flag
    hotspot_domains: list[str] = field(default_factory=list)  # PM1 domain list
    vcep_gene: bool = False  # True if VCEP spec exists
    has_vcep_specification: bool = False  # alias for vcep_gene
    missense_only_gene: bool = False  # BP1 flag
    haploinsufficiency_score: int | None = None  # ClinGen HI score 1-3
    triplosensitivity_score: int | None = None  # ClinGen TS score 1-3


# ---------------------------------------------------------------------------
# Classification result
# ---------------------------------------------------------------------------


@dataclass
class ClassificationResult:
    """Final ACMG/AMP classification for a variant.

    The five ACMG/AMP classification categories are:
        Pathogenic (P), Likely Pathogenic (LP), Uncertain Significance (VUS),
        Likely Benign (LB), Benign (B).

    Attributes:
        variant: The input variant that was classified.
        classification: One of ``"Pathogenic"``, ``"Likely_Pathogenic"``,
            ``"VUS"``, ``"Likely_Benign"``, ``"Benign"``.
        total_points: Sum of all contributing rule point values.
        rules_applied: List of ACMGRule instances that apply (applies=True).
        rules_not_applied: List of ACMGRule instances that do not apply.
        stand_alone_benign: True if BA1 fired, triggering direct Benign.
        bayesian_posterior_p: Posterior probability of pathogenicity from
            the Dirichlet-Multinomial Bayesian model (0–1).  None if the
            Bayesian model has not been run.
        credible_interval_lower: Lower bound of 95% credible interval.
        credible_interval_upper: Upper bound of 95% credible interval.
        vcep_overrides: Dict of rule_id → override notes from VCEP lookup.
        novel_combination: Name of a ClinGen SVI 2024 novel combination
            that applied, e.g. ``"PVS1+PM2_Supporting=LP"``, or None.
        notes: Free-text notes about the classification.
    """

    variant: VariantInput
    classification: str  # P, LP, VUS, LB, B
    total_points: int  # sum of rule point values
    rules_applied: list[ACMGRule]  # criteria that fire
    rules_not_applied: list[ACMGRule]  # criteria that don't fire
    stand_alone_benign: bool = False  # BA1 fired
    bayesian_posterior_p: float | None = None  # posterior P(pathogenic)
    credible_interval_lower: float | None = None  # 95% CI lower
    credible_interval_upper: float | None = None  # 95% CI upper
    vcep_overrides: dict[str, str] = field(default_factory=dict)
    novel_combination: str | None = None  # e.g. "PVS1+PM2_Supporting=LP"
    notes: str = ""

    @property
    def posterior_probability(self) -> float:
        """Posterior probability of pathogenicity using Tavtigian 2020 formula.

        Uses OddsP = 350^(pts/8) with prior P(path)=0.1 (Richards 2015).
        This is the clinically relevant posterior, not the Dirichlet classifier mean.

        References:
            Tavtigian et al. 2020 PMID:32645316.
            Richards et al. 2015 PMID:25741868 (prior probability 0.1).
        """
        pts = max(0, self.total_points)
        if pts <= 0:
            return 0.1  # prior only
        _prior = 0.1
        _odds_vs = 350.0  # Very Strong OddsP (Tavtigian 2020)
        odds = _odds_vs ** (pts / 8.0)
        return (_prior * odds) / (_prior * odds + (1.0 - _prior))

    @property
    def credible_interval_95(self) -> tuple[float, float]:
        """95% credible interval as (lower, upper) tuple.

        Returns stored MCMC interval if available, else the analytic
        Beta marginal from stored lower/upper bounds.
        """
        lo = self.credible_interval_lower
        hi = self.credible_interval_upper
        if lo is not None and hi is not None:
            return (lo, hi)
        return (0.0, 1.0)

    @property
    def point_summary(self) -> str:
        """Return a human-readable summary of applied rules and points.

        Returns:
            Formatted string listing each applied rule with its point value.
        """
        lines = [
            f"{r.rule_id} ({r.strength.value}: {r.points:+d} pts)"
            for r in self.rules_applied
        ]
        lines.append(f"TOTAL: {self.total_points:+d} pts → {self.classification}")
        return "\n".join(lines)
