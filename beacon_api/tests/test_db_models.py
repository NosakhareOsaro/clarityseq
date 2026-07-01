"""
beacon_api.tests.test_db_models
================================
pytest tests for the SQLAlchemy 2.x ORM models in beacon_api.db.models.

These are declarative models with no custom business logic, so tests focus
on:
    - Import / class-definition coverage (table names, columns, indexes).
    - Object construction with declarative __init__ kwargs.
    - __table_args__ (composite indexes / unique constraints) are present.
"""

from __future__ import annotations

import uuid

from sqlalchemy import Index, UniqueConstraint

from beacon_api.db.models import Base, BeaconDataset, BeaconIndividual, BeaconVariant


class TestBase:
    """Tests for the shared declarative Base."""

    def test_base_is_declarative_base(self) -> None:
        """Base subclasses DeclarativeBase and has a metadata registry."""
        assert hasattr(Base, "metadata")
        assert hasattr(Base, "registry")


class TestBeaconVariant:
    """Tests for the BeaconVariant ORM model."""

    def test_tablename(self) -> None:
        """__tablename__ is 'beacon_variants'."""
        assert BeaconVariant.__tablename__ == "beacon_variants"

    def test_construction_with_required_fields(self) -> None:
        """BeaconVariant can be constructed with required fields only."""
        variant = BeaconVariant(
            vrs_id="ga4gh:VA.abcdefghijklmnopqrstuvwx",
            chrom="chr17",
            pos=43044295,
            ref="G",
            alt="A",
        )
        assert variant.vrs_id == "ga4gh:VA.abcdefghijklmnopqrstuvwx"
        assert variant.chrom == "chr17"
        assert variant.pos == 43044295
        assert variant.ref == "G"
        assert variant.alt == "A"

    def test_construction_with_all_optional_fields(self) -> None:
        """BeaconVariant accepts all optional annotation fields."""
        variant_id = uuid.uuid4()
        variant = BeaconVariant(
            id=variant_id,
            vrs_id="ga4gh:VA.abcdefghijklmnopqrstuvwx",
            chrom="chr17",
            pos=43044295,
            ref="G",
            alt="A",
            variant_type="SNP",
            gene_symbol="BRCA1",
            hgvsc="c.68_69delAG",
            hgvsp="p.Glu23fs",
            consequence="missense_variant",
            gnomad_af=0.000012,
            gnomad_ac=10,
            gnomad_nhomalt=0,
            gnomad_popmax_af=0.00002,
            clinvar_id="VCV000012345",
            clinvar_classification="Pathogenic",
            acmg_class="Pathogenic",
            bayesian_posterior_p=0.999,
            dataset_id="clarityseq.wgs.grch38",
            extra={"alphamissense": 0.9},
        )
        assert variant.id == variant_id
        assert variant.gene_symbol == "BRCA1"
        assert variant.acmg_class == "Pathogenic"
        assert variant.extra == {"alphamissense": 0.9}

    def test_table_args_include_index_and_unique_constraint(self) -> None:
        """__table_args__ has the composite coordinate index and uniqueness."""
        index_names = {
            arg.name for arg in BeaconVariant.__table_args__ if isinstance(arg, Index)
        }
        unique_names = {
            arg.name
            for arg in BeaconVariant.__table_args__
            if isinstance(arg, UniqueConstraint)
        }
        assert "ix_beacon_variants_coords" in index_names
        assert "uq_beacon_variants_coords" in unique_names

    def test_column_names_present(self) -> None:
        """The mapped table exposes the expected column names."""
        columns = {c.name for c in BeaconVariant.__table__.columns}
        expected = {
            "id",
            "vrs_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_type",
            "gene_symbol",
            "hgvsc",
            "hgvsp",
            "consequence",
            "gnomad_af",
            "gnomad_ac",
            "gnomad_nhomalt",
            "gnomad_popmax_af",
            "clinvar_id",
            "clinvar_classification",
            "acmg_class",
            "bayesian_posterior_p",
            "dataset_id",
            "extra",
            "created_at",
            "updated_at",
        }
        assert expected.issubset(columns)


class TestBeaconIndividual:
    """Tests for the BeaconIndividual ORM model."""

    def test_tablename(self) -> None:
        """__tablename__ is 'beacon_individuals'."""
        assert BeaconIndividual.__tablename__ == "beacon_individuals"

    def test_construction_with_required_fields(self) -> None:
        """BeaconIndividual can be constructed with the required individual_id."""
        individual = BeaconIndividual(individual_id="IND-00001")
        assert individual.individual_id == "IND-00001"

    def test_construction_with_optional_fields(self) -> None:
        """BeaconIndividual accepts phenotype/disease JSONB fields."""
        individual = BeaconIndividual(
            individual_id="IND-00002",
            sex="FEMALE",
            ethnicity="Not specified",
            ethnicity_id="HANCESTRO:0004",
            phenotypic_features=[{"hpo_id": "HP:0001250", "label": "Seizure"}],
            diseases=[{"omim_id": "OMIM:114480", "label": "Breast cancer"}],
            dataset_id="clarityseq.wgs.grch38",
        )
        assert individual.sex == "FEMALE"
        assert individual.phenotypic_features[0]["hpo_id"] == "HP:0001250"
        assert individual.diseases[0]["omim_id"] == "OMIM:114480"

    def test_table_args_include_gin_indexes(self) -> None:
        """__table_args__ has GIN indexes on the JSONB columns."""
        index_names = {
            arg.name for arg in BeaconIndividual.__table_args__ if isinstance(arg, Index)
        }
        assert "ix_beacon_individuals_phenotypes" in index_names
        assert "ix_beacon_individuals_diseases" in index_names


class TestBeaconDataset:
    """Tests for the BeaconDataset ORM model."""

    def test_tablename(self) -> None:
        """__tablename__ is 'beacon_datasets'."""
        assert BeaconDataset.__tablename__ == "beacon_datasets"

    def test_construction_with_required_fields(self) -> None:
        """BeaconDataset can be constructed with id and name."""
        dataset = BeaconDataset(id="clarityseq.wgs.grch38", name="ClaritySeq WGS GRCh38")
        assert dataset.id == "clarityseq.wgs.grch38"
        assert dataset.name == "ClaritySeq WGS GRCh38"

    def test_construction_with_all_fields(self) -> None:
        """BeaconDataset accepts full metadata including counts and visibility."""
        dataset = BeaconDataset(
            id="clarityseq.wgs.grch38",
            name="ClaritySeq WGS GRCh38",
            description="Whole genome sequencing cohort",
            assembly_id="GRCh38",
            variant_count=1000,
            sample_count=50,
            is_public=True,
        )
        assert dataset.description == "Whole genome sequencing cohort"
        assert dataset.assembly_id == "GRCh38"
        assert dataset.variant_count == 1000
        assert dataset.sample_count == 50
        assert dataset.is_public is True

    def test_column_default_metadata(self) -> None:
        """assembly_id and is_public columns declare Python-side defaults."""
        columns = {c.name: c for c in BeaconDataset.__table__.columns}
        assert columns["assembly_id"].default.arg == "GRCh38"
        assert columns["is_public"].default.arg is False
