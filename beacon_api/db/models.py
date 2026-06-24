"""
beacon_api.db.models
====================
SQLAlchemy 2.x ORM models for the Beacon v2.1.1 PostgreSQL database.

Tables:
    beacon_variants   — Genomic variants with VRS v2.0 identifiers.
    beacon_individuals — Individual metadata (phenotype, disease, sex).
    beacon_datasets   — Dataset registry.

All tables use PostgreSQL-specific types where beneficial:
    - JSONB for complex nested data (phenotypic features, diseases).
    - UUID primary keys for global uniqueness.
    - Indexes on commonly queried columns (chrom, pos, gene_symbol, vrs_id).

References:
    SQLAlchemy 2.x docs: https://docs.sqlalchemy.org/en/20/
    asyncpg: https://magicstack.github.io/asyncpg/
    GA4GH Beacon v2.1.1 data model.
"""

from __future__ import annotations

import uuid
from datetime import datetime

from sqlalchemy import (
    BigInteger,
    Boolean,
    DateTime,
    Float,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import JSONB, UUID
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from sqlalchemy.sql import func


class Base(DeclarativeBase):
    """SQLAlchemy 2.x declarative base for all Beacon ORM models."""

    pass


class BeaconVariant(Base):
    """ORM model for genomic variants in the Beacon database.

    Stores variant-level data including VRS v2.0 identifiers, coordinates,
    population frequencies (gnomAD v4.1), and ACMG classification.

    Attributes:
        id: Internal UUID primary key.
        vrs_id: GA4GH VRS v2.0 computed identifier (24-char digest with prefix).
        chrom: Chromosome in GRCh38 notation (e.g. ``"chr17"``).
        pos: 1-based genomic position.
        ref: Reference allele (VCF notation).
        alt: Alternate allele (VCF notation).
        variant_type: Variant type (SNP, INDEL, MNV, etc.).
        gene_symbol: HGNC gene symbol.
        hgvsc: HGVS cDNA notation on MANE Select transcript.
        hgvsp: HGVS protein notation.
        consequence: VEP consequence string.
        gnomad_af: gnomAD v4.1 global allele frequency.
        gnomad_ac: gnomAD v4.1 allele count.
        gnomad_nhomalt: gnomAD v4.1 homozygous individual count.
        gnomad_popmax_af: gnomAD v4.1 maximum population-specific AF.
        clinvar_id: ClinVar RCV or VCV accession.
        clinvar_classification: ClinVar classification string.
        acmg_class: ACMG/AMP classification (P/LP/VUS/LB/B).
        bayesian_posterior_p: Posterior P(pathogenic) from BayesACMG.
        dataset_id: Foreign key to beacon_datasets.id.
        extra: JSONB for additional annotations.
        created_at: Row creation timestamp.
        updated_at: Last update timestamp.
    """

    __tablename__ = "beacon_variants"

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4,
        comment="Internal UUID primary key.",
    )
    vrs_id: Mapped[str] = mapped_column(
        String(40),
        nullable=False,
        index=True,
        comment="GA4GH VRS v2.0 computed identifier (ga4gh:VA.<24-char-digest>).",
    )
    chrom: Mapped[str] = mapped_column(
        String(10),
        nullable=False,
        comment="Chromosome in GRCh38 notation (e.g. chr17).",
    )
    pos: Mapped[int] = mapped_column(
        BigInteger,
        nullable=False,
        comment="1-based genomic position (GRCh38).",
    )
    ref: Mapped[str] = mapped_column(
        Text,
        nullable=False,
        comment="Reference allele (VCF notation, uppercase).",
    )
    alt: Mapped[str] = mapped_column(
        Text,
        nullable=False,
        comment="Alternate allele (VCF notation, uppercase).",
    )
    variant_type: Mapped[str | None] = mapped_column(
        String(30),
        nullable=True,
        comment="Variant type: SNP, INDEL, MNV, CNV, etc.",
    )
    gene_symbol: Mapped[str | None] = mapped_column(
        String(50),
        nullable=True,
        index=True,
        comment="HGNC gene symbol.",
    )
    hgvsc: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
        comment="HGVS cDNA notation on MANE Select transcript.",
    )
    hgvsp: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
        comment="HGVS protein notation.",
    )
    consequence: Mapped[str | None] = mapped_column(
        String(60),
        nullable=True,
        comment="VEP 111 consequence string (e.g. missense_variant).",
    )
    gnomad_af: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
        comment="gnomAD v4.1 global allele frequency (April 2024, 807,162 individuals).",
    )
    gnomad_ac: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
        comment="gnomAD v4.1 allele count.",
    )
    gnomad_nhomalt: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
        comment="gnomAD v4.1 number of homozygous individuals.",
    )
    gnomad_popmax_af: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
        comment="gnomAD v4.1 maximum population-specific allele frequency.",
    )
    clinvar_id: Mapped[str | None] = mapped_column(
        String(30),
        nullable=True,
        comment="ClinVar RCV or VCV accession number.",
    )
    clinvar_classification: Mapped[str | None] = mapped_column(
        String(40),
        nullable=True,
        comment="ClinVar clinical significance string.",
    )
    acmg_class: Mapped[str | None] = mapped_column(
        String(20),
        nullable=True,
        index=True,
        comment="ACMG/AMP classification: Pathogenic, Likely_Pathogenic, VUS, Likely_Benign, Benign.",
    )
    bayesian_posterior_p: Mapped[float | None] = mapped_column(
        Float,
        nullable=True,
        comment="BayesACMG posterior P(pathogenic) [0-1].",
    )
    dataset_id: Mapped[str | None] = mapped_column(
        String(80),
        nullable=True,
        comment="Dataset identifier (foreign key to beacon_datasets).",
    )
    extra: Mapped[dict | None] = mapped_column(
        JSONB,
        nullable=True,
        comment="JSONB for additional annotations (AlphaMissense, SpliceAI, etc.).",
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
        comment="Row creation timestamp.",
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        onupdate=func.now(),
        nullable=False,
        comment="Last update timestamp.",
    )

    __table_args__ = (
        # Composite index for coordinate-based queries (Beacon /g_variants)
        Index("ix_beacon_variants_coords", "chrom", "pos", "ref", "alt"),
        UniqueConstraint("chrom", "pos", "ref", "alt", name="uq_beacon_variants_coords"),
    )


class BeaconIndividual(Base):
    """ORM model for individual-level data in the Beacon database.

    Individual records are access-controlled via GA4GH Passport.
    Phenotypic features and diseases are stored as JSONB arrays.

    Attributes:
        id: Internal UUID primary key.
        individual_id: External individual identifier (pseudonymised).
        sex: Biological sex (FEMALE/MALE/UNKNOWN_SEX).
        ethnicity: Self-reported ethnicity string.
        ethnicity_id: Ontology ID for ethnicity (e.g. HANCESTRO term).
        phenotypic_features: JSONB array of HPO term objects.
        diseases: JSONB array of disease objects (OMIM/Orphanet).
        dataset_id: Dataset identifier.
        created_at: Row creation timestamp.
    """

    __tablename__ = "beacon_individuals"

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4,
        comment="Internal UUID primary key.",
    )
    individual_id: Mapped[str] = mapped_column(
        String(80),
        nullable=False,
        unique=True,
        index=True,
        comment="External pseudonymised individual identifier.",
    )
    sex: Mapped[str | None] = mapped_column(
        String(20),
        nullable=True,
        comment="Biological sex: FEMALE, MALE, UNKNOWN_SEX.",
    )
    ethnicity: Mapped[str | None] = mapped_column(
        String(100),
        nullable=True,
        comment="Self-reported ethnicity free text.",
    )
    ethnicity_id: Mapped[str | None] = mapped_column(
        String(30),
        nullable=True,
        comment="HANCESTRO ontology identifier for ethnicity.",
    )
    phenotypic_features: Mapped[list | None] = mapped_column(
        JSONB,
        nullable=True,
        comment="JSONB array of HPO phenotype objects {hpo_id, label, excluded}.",
    )
    diseases: Mapped[list | None] = mapped_column(
        JSONB,
        nullable=True,
        comment="JSONB array of disease objects {omim_id, label, stage}.",
    )
    dataset_id: Mapped[str | None] = mapped_column(
        String(80),
        nullable=True,
        comment="Dataset identifier.",
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
        comment="Row creation timestamp.",
    )

    __table_args__ = (
        Index("ix_beacon_individuals_phenotypes", "phenotypic_features", postgresql_using="gin"),
        Index("ix_beacon_individuals_diseases", "diseases", postgresql_using="gin"),
    )


class BeaconDataset(Base):
    """ORM model for Beacon dataset registry.

    Attributes:
        id: Dataset identifier string (primary key).
        name: Human-readable dataset name.
        description: Dataset description.
        assembly_id: Genome assembly (e.g. ``"GRCh38"``).
        variant_count: Cached count of variants in this dataset.
        sample_count: Cached count of individuals.
        is_public: True if dataset is publicly accessible.
        created_at: Dataset creation timestamp.
        updated_at: Last update timestamp.
    """

    __tablename__ = "beacon_datasets"

    id: Mapped[str] = mapped_column(
        String(80),
        primary_key=True,
        comment="Dataset identifier (primary key).",
    )
    name: Mapped[str] = mapped_column(
        String(200),
        nullable=False,
        comment="Human-readable dataset name.",
    )
    description: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
        comment="Dataset description.",
    )
    assembly_id: Mapped[str] = mapped_column(
        String(20),
        nullable=False,
        default="GRCh38",
        comment="Genome assembly identifier (e.g. GRCh38).",
    )
    variant_count: Mapped[int | None] = mapped_column(
        BigInteger,
        nullable=True,
        comment="Cached total variant count.",
    )
    sample_count: Mapped[int | None] = mapped_column(
        BigInteger,
        nullable=True,
        comment="Cached total individual count.",
    )
    is_public: Mapped[bool] = mapped_column(
        Boolean,
        nullable=False,
        default=False,
        comment="True if dataset is publicly accessible (no Passport required).",
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
        comment="Dataset creation timestamp.",
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        onupdate=func.now(),
        nullable=False,
        comment="Last update timestamp.",
    )
