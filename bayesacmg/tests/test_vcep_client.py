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

import time
from unittest.mock import AsyncMock, MagicMock, patch

import httpx
import pytest

import bayesacmg.vcep_client as vcep_client_module
from bayesacmg.models import EvidenceStrength
from bayesacmg.vcep_client import (
    VCEPClient,
    VCEPSpec,
    VCEPSpecification,
    _fetch_spec_from_api,
    _parse_spec,
    clear_cache,
    get_vcep_spec,
    get_vcep_spec_sync,
)


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
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=None)
        ):
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

        assert (
            mock_api.call_count == 1
        ), "API should only be called once (cached after first call)"
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
    async def test_default_pm2_supporting_when_no_vcep(
        self, vcep_client: VCEPClient
    ) -> None:
        """When no VCEP specification exists, PM2 defaults to Supporting (1 pt)."""
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=None)
        ):
            result = await vcep_client.get_specification("OBSCUREGENE")

        # No VCEP → use ClinGen SVI 2024 default = Supporting
        assert result is None  # None means use default Supporting

    @pytest.mark.asyncio
    async def test_api_timeout_handled_gracefully(
        self, vcep_client: VCEPClient
    ) -> None:
        """API timeout should return None (not raise) so pipeline continues."""
        import httpx

        with patch.object(
            vcep_client,
            "_fetch_from_api",
            new=AsyncMock(side_effect=httpx.TimeoutException("")),
        ):
            # Should not raise; should return None and use defaults
            result = await vcep_client.get_specification("BRCA2")
        assert result is None, "Timeout should return None (use defaults), not raise"


# ---------------------------------------------------------------------------
# VCEPSpec.is_expired
# ---------------------------------------------------------------------------


class TestVCEPSpecIsExpired:
    """Tests for VCEPSpec.is_expired property (24-hour TTL)."""

    def test_not_expired_when_recently_retrieved(self) -> None:
        spec = VCEPSpec(gene_symbol="BRCA1", retrieved_at=time.time())
        assert spec.is_expired is False

    def test_expired_when_older_than_ttl(self) -> None:
        """retrieved_at more than 24h (86,400s) in the past is expired."""
        spec = VCEPSpec(gene_symbol="BRCA1", retrieved_at=time.time() - 90_000)
        assert spec.is_expired is True


# ---------------------------------------------------------------------------
# _parse_spec — module-level parsing helper
# ---------------------------------------------------------------------------


class TestParseSpec:
    """Tests for _parse_spec() — CSpec API response → VCEPSpec."""

    def test_empty_response_returns_default_spec(self) -> None:
        """No VCEP specification for gene → default spec, pm2_weight='supporting'."""
        spec = _parse_spec("OBSCUREGENE", [])
        assert spec.gene_symbol == "OBSCUREGENE"
        assert spec.vcep_name == ""
        assert spec.pm2_weight == "supporting"

    def test_response_with_pm2_moderate_override(self) -> None:
        """A criteria entry with PM2 at Moderate strength sets pm2_weight='moderate'."""
        raw = [
            {
                "vcep_name": "RASopathy VCEP",
                "vcep_id": "42",
                "criteria": [
                    {"id": "PM2", "strength": "Moderate"},
                    {"id": "PP3", "strength": "Supporting"},
                ],
            }
        ]
        spec = _parse_spec("PTPN11", raw)
        assert spec.vcep_name == "RASopathy VCEP"
        assert spec.vcep_id == "42"
        assert spec.pm2_weight == "moderate"
        assert spec.raw_spec == raw[0]

    def test_response_without_pm2_criterion_defaults_to_supporting(self) -> None:
        """VCEP spec exists but has no PM2 entry → pm2_weight stays 'supporting'."""
        raw = [
            {
                "vcep_name": "BRCA Exchange",
                "vcep_id": "7",
                "criteria": [{"id": "PS3", "strength": "Strong"}],
            }
        ]
        spec = _parse_spec("BRCA1", raw)
        assert spec.vcep_name == "BRCA Exchange"
        assert spec.pm2_weight == "supporting"

    def test_response_with_pm2_strong_is_not_moderate(self) -> None:
        """PM2 present but not at 'moderate' strength → stays at default 'supporting'."""
        raw = [
            {
                "vcep_name": "Some VCEP",
                "criteria": [{"id": "PM2", "strength": "Supporting"}],
            }
        ]
        spec = _parse_spec("GENE1", raw)
        assert spec.pm2_weight == "supporting"


# ---------------------------------------------------------------------------
# _fetch_spec_from_api — real network call path (httpx client mocked)
# ---------------------------------------------------------------------------


class TestFetchSpecFromApi:
    """Tests for _fetch_spec_from_api() with a mocked httpx.AsyncClient."""

    @pytest.mark.asyncio
    async def test_returns_list_response_directly(self) -> None:
        """When the API returns a JSON list, it is returned as-is."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=[{"vcep_name": "BRCA Exchange"}])
        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        result = await _fetch_spec_from_api("BRCA1", mock_client)
        assert result == [{"vcep_name": "BRCA Exchange"}]
        mock_client.get.assert_awaited_once()

    @pytest.mark.asyncio
    async def test_returns_results_key_from_dict_response(self) -> None:
        """When the API returns a dict with a 'results' key, that list is extracted."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value={"results": [{"vcep_name": "X"}]})
        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        result = await _fetch_spec_from_api("GENEX", mock_client)
        assert result == [{"vcep_name": "X"}]

    @pytest.mark.asyncio
    async def test_dict_response_without_results_key_returns_empty_list(self) -> None:
        """Dict response without a 'results' key falls back to an empty list."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value={"other_key": []})
        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        result = await _fetch_spec_from_api("GENEX", mock_client)
        assert result == []

    @pytest.mark.asyncio
    async def test_raises_on_http_status_error(self) -> None:
        """HTTP error status propagates (retried internally by tenacity, then reraised)."""
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock(
            side_effect=httpx.HTTPStatusError(
                "500", request=MagicMock(), response=MagicMock()
            )
        )
        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        with pytest.raises(httpx.HTTPStatusError):
            await _fetch_spec_from_api("BRCA1", mock_client)


# ---------------------------------------------------------------------------
# get_vcep_spec — module-level cached async accessor
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clear_module_cache():
    """Ensure the module-level _spec_cache never leaks between tests."""
    vcep_client_module._spec_cache.clear()
    yield
    vcep_client_module._spec_cache.clear()


class TestGetVcepSpec:
    """Tests for the module-level get_vcep_spec() function."""

    @pytest.mark.asyncio
    async def test_cache_hit_avoids_api_call(self) -> None:
        """A fresh cached spec is returned without calling the API."""
        cached_spec = VCEPSpec(gene_symbol="BRCA1", vcep_name="Cached", retrieved_at=time.time())
        vcep_client_module._spec_cache["BRCA1"] = cached_spec

        with patch(
            "bayesacmg.vcep_client._fetch_spec_from_api",
            new=AsyncMock(side_effect=AssertionError("API should not be called")),
        ):
            result = await get_vcep_spec("BRCA1")
        assert result is cached_spec

    @pytest.mark.asyncio
    async def test_cache_miss_calls_api_and_caches(self) -> None:
        """On cache miss, the API is queried, parsed, and the result cached."""
        raw = [{"vcep_name": "ENIGMA", "vcep_id": "1", "criteria": []}]
        with patch(
            "bayesacmg.vcep_client._fetch_spec_from_api",
            new=AsyncMock(return_value=raw),
        ) as mock_fetch:
            result = await get_vcep_spec("BRCA2", client=AsyncMock(spec=httpx.AsyncClient))
        mock_fetch.assert_awaited_once()
        assert result.vcep_name == "ENIGMA"
        assert vcep_client_module._spec_cache["BRCA2"] is result

    @pytest.mark.asyncio
    async def test_api_error_returns_default_spec(self) -> None:
        """Any API error is caught and a default (fail-safe) VCEPSpec is returned."""
        with patch(
            "bayesacmg.vcep_client._fetch_spec_from_api",
            new=AsyncMock(side_effect=httpx.RequestError("network down")),
        ):
            result = await get_vcep_spec("BRCA1", client=AsyncMock(spec=httpx.AsyncClient))
        assert result.vcep_name == ""
        assert result.pm2_weight == "supporting"

    @pytest.mark.asyncio
    async def test_creates_own_client_when_none_provided(self) -> None:
        """When no client is passed, get_vcep_spec creates and closes its own."""
        raw = [{"vcep_name": "Own Client VCEP"}]
        with patch(
            "bayesacmg.vcep_client._fetch_spec_from_api",
            new=AsyncMock(return_value=raw),
        ):
            result = await get_vcep_spec("SOMEGENE")
        assert result.vcep_name == "Own Client VCEP"

    @pytest.mark.asyncio
    async def test_expired_cache_entry_triggers_refetch(self) -> None:
        """An expired cached spec is not reused; the API is queried again."""
        stale_spec = VCEPSpec(
            gene_symbol="BRCA1", vcep_name="Stale", retrieved_at=time.time() - 100_000
        )
        vcep_client_module._spec_cache["BRCA1"] = stale_spec
        raw = [{"vcep_name": "Fresh"}]
        with patch(
            "bayesacmg.vcep_client._fetch_spec_from_api",
            new=AsyncMock(return_value=raw),
        ):
            result = await get_vcep_spec("BRCA1", client=AsyncMock(spec=httpx.AsyncClient))
        assert result.vcep_name == "Fresh"


# ---------------------------------------------------------------------------
# get_vcep_spec_sync / clear_cache — module-level sync wrapper & cache clear
# ---------------------------------------------------------------------------


class TestGetVcepSpecSync:
    def test_returns_result_from_async_function(self) -> None:
        """get_vcep_spec_sync() runs the coroutine via asyncio.run and returns its result."""
        expected = VCEPSpec(gene_symbol="BRCA1", vcep_name="Sync Result")
        with patch(
            "bayesacmg.vcep_client.get_vcep_spec",
            new=AsyncMock(return_value=expected),
        ):
            result = get_vcep_spec_sync("BRCA1")
        assert result is expected


class TestClearCacheFunction:
    def test_clear_cache_empties_module_cache(self) -> None:
        """clear_cache() empties the module-level _spec_cache dict."""
        vcep_client_module._spec_cache["BRCA1"] = VCEPSpec(gene_symbol="BRCA1")
        assert len(vcep_client_module._spec_cache) == 1
        clear_cache()
        assert vcep_client_module._spec_cache == {}


# ---------------------------------------------------------------------------
# VCEPClient — real (unmocked) internal methods
# ---------------------------------------------------------------------------


class TestVCEPClientRealFetchFromApi:
    """Tests for VCEPClient._fetch_from_api() real implementation (not mocked)."""

    @pytest.mark.asyncio
    async def test_returns_spec_when_vcep_name_present(self, vcep_client: VCEPClient) -> None:
        """A spec with a non-empty vcep_name is returned as-is."""
        spec = VCEPSpec(gene_symbol="BRCA1", vcep_name="ENIGMA")
        with patch(
            "bayesacmg.vcep_client.get_vcep_spec", new=AsyncMock(return_value=spec)
        ):
            result = await vcep_client._fetch_from_api("BRCA1")
        assert result is spec

    @pytest.mark.asyncio
    async def test_returns_none_when_vcep_name_empty(self, vcep_client: VCEPClient) -> None:
        """A default spec (empty vcep_name) is normalised to None (no VCEP spec)."""
        spec = VCEPSpec(gene_symbol="OBSCUREGENE", vcep_name="")
        with patch(
            "bayesacmg.vcep_client.get_vcep_spec", new=AsyncMock(return_value=spec)
        ):
            result = await vcep_client._fetch_from_api("OBSCUREGENE")
        assert result is None

    @pytest.mark.asyncio
    async def test_returns_none_on_exception(self, vcep_client: VCEPClient) -> None:
        """Any exception from get_vcep_spec is swallowed; returns None."""
        with patch(
            "bayesacmg.vcep_client.get_vcep_spec",
            new=AsyncMock(side_effect=RuntimeError("boom")),
        ):
            result = await vcep_client._fetch_from_api("BRCA1")
        assert result is None


class TestVCEPClientGetSpecSync:
    """Tests for VCEPClient.get_spec_sync()."""

    def test_returns_spec_synchronously(self, vcep_client: VCEPClient) -> None:
        spec = VCEPSpecification(gene_symbol="BRCA1", vcep_name="ENIGMA")
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=spec)
        ):
            result = vcep_client.get_spec_sync("BRCA1")
        assert result is spec


class TestVCEPClientClear:
    """Tests for VCEPClient.clear() — instance-level cache reset."""

    @pytest.mark.asyncio
    async def test_clear_empties_instance_cache(self, vcep_client: VCEPClient) -> None:
        mock_spec = VCEPSpecification(gene_symbol="BRCA1", vcep_name="ENIGMA")
        with patch.object(
            vcep_client, "_fetch_from_api", new=AsyncMock(return_value=mock_spec)
        ):
            await vcep_client.get_specification("BRCA1")
        assert len(vcep_client._cache) == 1
        vcep_client.clear()
        assert vcep_client._cache == {}
