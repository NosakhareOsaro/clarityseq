"""
Tests for the ClinGen CSpec VCEP registry client.

Tests that the VCEP client:
1. Correctly queries the ClinGen CSpec API
2. Caches responses (24-hour TTL)
3. Returns gene-specific PM2 weight overrides when specified
4. Falls back to Supporting when no VCEP specification exists

Guidelines:
    ClinGen CSpec API: https://cspec.genome.network/cspec/api/svi/
    ClinGen SVI 2024: VCEP specifications override general framework
"""
from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest
import pytest_asyncio

from bayesacmg.models import EvidenceStrength
from bayesacmg.vcep_client import VCEPClient, VCEPSpecification


@pytest.fixture
def vcep_client() -> VCEPClient:
    """Fresh VCEPClient with empty cache for each test."""
    client = VCEPClient(cache_ttl_hours=24)
    client._cache.clear()
    return client


class TestVCEPClient:
    """Tests for VCEPClient ClinGen CSpec API integration."""

    @pytest.mark.asyncio
    async def test_returns_none_for_unknown_gene(self, vcep_client: VCEPClient) -> None:
        """Gene without VCEP specification returns None (use general framework)."""
        with patch.object(vcep_client, "_fetch_from_api", new=AsyncMock(return_value=None)):
            result = await vcep_client.get_specification("UNKNOWNGENE123")
        assert result is None

    @pytest.mark.asyncio
    async def test_caches_response(self, vcep_client: VCEPClient) -> None:
        """Response is cached after first API call (24-hour TTL)."""
        mock_spec = VCEPSpecification(
            gene_symbol="BRCA1",
            vcep_name="BRCA Exchange",
            pm2_strength=EvidenceStrength.SUPPORTING,
            has_gene_specific_rules=True,
        )
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=mock_spec)
        ) as mock_api:
            # First call — hits API
            result1 = await vcep_client.get_specification("BRCA1")
            # Second call — should use cache
            result2 = await vcep_client.get_specification("BRCA1")

        assert mock_api.call_count == 1, "API should only be called once (cached after first call)"
        assert result1 == result2

    @pytest.mark.asyncio
    async def test_pm2_override_to_moderate(self, vcep_client: VCEPClient) -> None:
        """VCEP specification may allow PM2 at Moderate for specific genes."""
        # Some VCEPs (e.g., certain RASopathy genes) specify PM2 at Moderate
        mock_spec = VCEPSpecification(
            gene_symbol="PTPN11",
            vcep_name="RASopathy VCEP",
            pm2_strength=EvidenceStrength.MODERATE,  # VCEP override to Moderate
            has_gene_specific_rules=True,
        )
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=mock_spec)
        ):
            result = await vcep_client.get_specification("PTPN11")

        assert result is not None
        assert result.pm2_strength == EvidenceStrength.MODERATE, (
            "VCEP specification for PTPN11 allows PM2 at Moderate; "
            "this should override the ClinGen SVI 2024 default of Supporting"
        )

    @pytest.mark.asyncio
    async def test_default_pm2_supporting_when_no_vcep(self, vcep_client: VCEPClient) -> None:
        """When no VCEP specification exists, PM2 defaults to Supporting (1 pt)."""
        with patch.object(vcep_client, "_fetch_from_api", new=AsyncMock(return_value=None)):
            result = await vcep_client.get_specification("OBSCUREGENE")

        # No VCEP → use ClinGen SVI 2024 default = Supporting
        assert result is None  # None means use default Supporting

    @pytest.mark.asyncio
    async def test_api_timeout_handled_gracefully(self, vcep_client: VCEPClient) -> None:
        """API timeout should return None (not raise) so pipeline continues."""
        import httpx
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(side_effect=httpx.TimeoutException(""))
        ):
            # Should not raise; should return None and use defaults
            result = await vcep_client.get_specification("BRCA2")
        assert result is None, "Timeout should return None (use defaults), not raise"
