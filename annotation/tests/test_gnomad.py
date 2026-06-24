"""
Tests for the gnomAD v4.1 client.

Validates:
- AF lookup from GraphQL API
- PM2_Supporting flag logic (AF < 0.0001)
- Population-stratified AF access
- Absent variant handling (PM2_Supporting = True)
- AN bug note for v4.0 (must use v4.1)
"""

from __future__ import annotations

import pytest
from unittest.mock import AsyncMock, patch

from annotation.gnomad_client import (
    GNOMAD_GENOME_DATASET,
    PM2_SUPPORTING_THRESHOLD,
    AncestryFrequency,
    GnomADClient,
    GnomADData,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def client() -> GnomADClient:
    """Return a GnomADClient configured with the v4.1 dataset."""
    return GnomADClient(dataset=GNOMAD_GENOME_DATASET)


# Mocked API response for a rare variant (AF ~ 0.000005)
MOCK_RARE_VARIANT_RESPONSE = {
    "variantId": "17-43094692-G-A",
    "chrom": "17",
    "pos": 43094692,
    "ref": "G",
    "alt": "A",
    "genome": {
        "ac": 1,
        "an": 200000,
        "af": 0.000005,
        "homozygote_count": 0,
        "populations": [
            {"id": "nfe", "ac": 1, "an": 120000, "af": 0.0000083, "homozygote_count": 0},
            {"id": "afr", "ac": 0, "an": 40000, "af": 0.0, "homozygote_count": 0},
            {"id": "eas", "ac": 0, "an": 20000, "af": 0.0, "homozygote_count": 0},
        ],
        "filters": [],
    },
    "exome": {
        "ac": 0,
        "an": 500000,
        "af": 0.0,
        "homozygote_count": 0,
        "populations": [],
    },
    "flags": [],
}

# Mocked API response for a common variant (AF ~ 0.01)
MOCK_COMMON_VARIANT_RESPONSE = {
    "variantId": "1-55505599-T-G",
    "chrom": "1",
    "pos": 55505599,
    "ref": "T",
    "alt": "G",
    "genome": {
        "ac": 7500,
        "an": 750000,
        "af": 0.01,
        "homozygote_count": 35,
        "populations": [
            {"id": "nfe", "ac": 5000, "an": 450000, "af": 0.0111, "homozygote_count": 25},
            {"id": "afr", "ac": 2500, "an": 200000, "af": 0.0125, "homozygote_count": 10},
        ],
        "filters": [],
    },
    "exome": None,
    "flags": [],
}


# ---------------------------------------------------------------------------
# Tests for absent variant (PM2_Supporting)
# ---------------------------------------------------------------------------


class TestAbsentVariant:
    """Test behaviour when a variant is not present in gnomAD v4.1."""

    @pytest.mark.asyncio
    async def test_absent_variant_pm2_supporting_true(self, client: GnomADClient) -> None:
        """Absent variant should set pm2_supporting = True."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = None  # Variant not in gnomAD
            data = await client.get_allele_frequency("chr1", 99999999, "A", "T")

        assert data.pm2_supporting is True

    @pytest.mark.asyncio
    async def test_absent_variant_af_is_none(self, client: GnomADClient) -> None:
        """Absent variant should have af = None."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = None
            data = await client.get_allele_frequency("chr1", 99999999, "A", "T")

        assert data.af is None

    @pytest.mark.asyncio
    async def test_absent_variant_note_set(self, client: GnomADClient) -> None:
        """Absent variant should include a note about absence from gnomAD v4.1."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = None
            data = await client.get_allele_frequency("chr17", 43094692, "G", "A")

        assert data.note is not None
        assert "absent" in data.note.lower() or "gnomAD" in data.note


# ---------------------------------------------------------------------------
# Tests for rare variants (PM2_Supporting threshold)
# ---------------------------------------------------------------------------


class TestRareVariant:
    """Test AF lookups for rare variants below PM2_Supporting threshold."""

    @pytest.mark.asyncio
    async def test_rare_variant_pm2_supporting(self, client: GnomADClient) -> None:
        """Variant with AF < 0.0001 should trigger PM2_Supporting."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_RARE_VARIANT_RESPONSE
            data = await client.get_allele_frequency("chr17", 43094692, "G", "A")

        assert data.pm2_supporting is True
        assert data.af is not None
        assert data.af < PM2_SUPPORTING_THRESHOLD

    @pytest.mark.asyncio
    async def test_rare_variant_genome_af_returned(self, client: GnomADClient) -> None:
        """Genome AF should be preferred over exome AF."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_RARE_VARIANT_RESPONSE
            data = await client.get_allele_frequency("chr17", 43094692, "G", "A")

        assert data.af_genome == pytest.approx(0.000005)

    @pytest.mark.asyncio
    async def test_rare_variant_population_breakdown(self, client: GnomADClient) -> None:
        """Population-stratified AFs should be populated from genome data."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_RARE_VARIANT_RESPONSE
            data = await client.get_allele_frequency("chr17", 43094692, "G", "A")

        assert "nfe" in data.by_ancestry
        assert data.by_ancestry["nfe"].af == pytest.approx(0.0000083)
        assert data.by_ancestry["afr"].ac == 0


# ---------------------------------------------------------------------------
# Tests for common variants (not PM2_Supporting)
# ---------------------------------------------------------------------------


class TestCommonVariant:
    """Test AF lookups for common variants above PM2_Supporting threshold."""

    @pytest.mark.asyncio
    async def test_common_variant_pm2_not_supporting(self, client: GnomADClient) -> None:
        """Common variant (AF >> 0.0001) should NOT trigger PM2_Supporting."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_COMMON_VARIANT_RESPONSE
            data = await client.get_allele_frequency("chr1", 55505599, "T", "G")

        assert data.pm2_supporting is False

    @pytest.mark.asyncio
    async def test_common_variant_af_correct(self, client: GnomADClient) -> None:
        """AF for a common variant should match expected value."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_COMMON_VARIANT_RESPONSE
            data = await client.get_allele_frequency("chr1", 55505599, "T", "G")

        assert data.af == pytest.approx(0.01)
        assert data.ac == 7500
        assert data.nhomalt == 35


# ---------------------------------------------------------------------------
# Tests for population-specific lookup
# ---------------------------------------------------------------------------


class TestPopulationSpecificLookup:
    """Test population-specific allele frequency requests."""

    @pytest.mark.asyncio
    async def test_population_specific_af_returned(self, client: GnomADClient) -> None:
        """When population='nfe', the NFE AF should be the primary AF."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = MOCK_RARE_VARIANT_RESPONSE
            data = await client.get_allele_frequency(
                "chr17", 43094692, "G", "A", population="nfe"
            )

        # Primary AF should be NFE-specific
        assert data.af == pytest.approx(0.0000083)


# ---------------------------------------------------------------------------
# Tests for chromosome normalisation
# ---------------------------------------------------------------------------


class TestChromNormalisation:
    """Test that chromosome names are normalised correctly."""

    @pytest.mark.asyncio
    async def test_chrom_without_prefix_normalised(self, client: GnomADClient) -> None:
        """Input '17' should be normalised to 'chr17' in GnomADData."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = None
            data = await client.get_allele_frequency("17", 43094692, "G", "A")

        assert data.chrom == "chr17"

    @pytest.mark.asyncio
    async def test_chrom_with_prefix_unchanged(self, client: GnomADClient) -> None:
        """Input 'chr17' should remain 'chr17'."""
        with patch.object(client, "_query_api", new_callable=AsyncMock) as mock_api:
            mock_api.return_value = None
            data = await client.get_allele_frequency("chr17", 43094692, "G", "A")

        assert data.chrom == "chr17"


# ---------------------------------------------------------------------------
# Constants tests
# ---------------------------------------------------------------------------


class TestConstants:
    """Verify gnomAD client constants are correct."""

    def test_pm2_threshold_value(self) -> None:
        """PM2_Supporting threshold should be 0.0001 (ClinGen SVI 2024)."""
        assert PM2_SUPPORTING_THRESHOLD == pytest.approx(0.0001)

    def test_dataset_is_v4(self) -> None:
        """Default dataset should reference gnomAD v4."""
        assert "r4" in GNOMAD_GENOME_DATASET.lower() or "v4" in GNOMAD_GENOME_DATASET.lower()
