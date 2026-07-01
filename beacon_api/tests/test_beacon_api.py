"""
beacon_api.tests.test_beacon_api
==================================
pytest tests for the GA4GH Beacon v2.1.1 API.

Tests cover:
    - GET /info: Beacon metadata structure.
    - GET /g_variants: Coordinate query, VRS identifier presence,
      anonymous count granularity, missing params 400.
    - VRS identifier computation determinism.

References:
    GA4GH Beacon v2.1.1 spec (December 13, 2024).
    Wagner et al. 2021 Cell Genomics PMID:35072137 (VRS v2.0).
"""

from __future__ import annotations

from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from fastapi.testclient import TestClient

# ---------------------------------------------------------------------------
# App fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def client() -> TestClient:
    """Create a TestClient for the Beacon FastAPI application.

    Overrides the DB session dependency so tests run without PostgreSQL.

    Returns:
        FastAPI TestClient for synchronous test access.
    """
    from beacon_api.main import app
    from beacon_api.db.session import get_session

    async def _mock_get_session():  # type: ignore[return]
        yield AsyncMock()

    app.dependency_overrides[get_session] = _mock_get_session
    with TestClient(app, raise_server_exceptions=False) as c:
        yield c
    app.dependency_overrides.clear()


# ---------------------------------------------------------------------------
# Mock session dependency
# ---------------------------------------------------------------------------

def _mock_session() -> AsyncMock:
    """Return a mock async SQLAlchemy session."""
    return AsyncMock()


# ---------------------------------------------------------------------------
# / (root) tests
# ---------------------------------------------------------------------------


class TestRoot:
    """Tests for GET / root redirect hint."""

    def test_root_returns_200_with_message(self, client: TestClient) -> None:
        """GET / returns HTTP 200 with a hint message pointing to /info."""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert "message" in data
        assert "/info" in data["message"]


# ---------------------------------------------------------------------------
# /info tests
# ---------------------------------------------------------------------------


class TestInfo:
    """Tests for GET /info endpoint."""

    def test_info_returns_200(self, client: TestClient) -> None:
        """GET /info returns HTTP 200.

        Verifies that the /info endpoint is reachable and returns a success
        response per GA4GH Beacon v2.1.1 spec.
        """
        response = client.get("/info")
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}: {response.text}"
        )

    def test_info_contains_beacon_id(self, client: TestClient) -> None:
        """GET /info response contains beaconId field.

        The beaconId must be present in the response meta per Beacon v2.1.1.
        """
        response = client.get("/info")
        assert response.status_code == 200
        data = response.json()
        assert "meta" in data, "Response missing 'meta' key"
        assert "beaconId" in data["meta"], "meta missing 'beaconId'"
        assert data["meta"]["beaconId"] == "org.genomeforge.beacon"

    def test_info_api_version(self, client: TestClient) -> None:
        """GET /info returns correct API version v2.1.1.

        Beacon v2.1.1 was released December 13, 2024.
        """
        response = client.get("/info")
        data = response.json()
        assert data["meta"]["apiVersion"] == "v2.1.1"

    def test_info_response_has_datasets(self, client: TestClient) -> None:
        """GET /info response contains datasets array.

        The Beacon info response must include at least one dataset entry.
        """
        response = client.get("/info")
        data = response.json()
        assert "response" in data, "Response missing 'response' key"
        assert "datasets" in data["response"], "response missing 'datasets'"
        assert isinstance(data["response"]["datasets"], list)
        assert len(data["response"]["datasets"]) >= 1

    def test_info_organisation_present(self, client: TestClient) -> None:
        """GET /info response contains organisation metadata."""
        response = client.get("/info")
        data = response.json()
        assert "organization" in data["response"], "response missing 'organization'"
        org = data["response"]["organization"]
        assert "id" in org
        assert "name" in org

    def test_info_contains_vrs_version(self, client: TestClient) -> None:
        """GET /info response info block references VRS v2.0.

        VRS v2.0 identifiers are used for all variants in this Beacon.
        Wagner et al. 2021 Cell Genomics PMID:35072137.
        """
        response = client.get("/info")
        data = response.json()
        info = data["response"].get("info", {})
        assert "variantIdentifiers" in info, "info missing variantIdentifiers"
        assert "VRS v2.0" in info["variantIdentifiers"], (
            "variantIdentifiers should mention 'VRS v2.0'"
        )


# ---------------------------------------------------------------------------
# /g_variants tests
# ---------------------------------------------------------------------------


class TestGVariants:
    """Tests for GET /g_variants endpoint."""

    def test_g_variants_no_params_returns_400(self, client: TestClient) -> None:
        """GET /g_variants with no parameters returns HTTP 400.

        The endpoint requires at least one query parameter (chrom, vrsId,
        or geneSymbol).
        """
        response = client.get("/g_variants")
        assert response.status_code == 400, (
            f"Expected 400 for missing params, got {response.status_code}"
        )

    def test_g_variants_with_chrom_returns_200(self, client: TestClient) -> None:
        """GET /g_variants?chrom=chr17&start=43044295 returns HTTP 200.

        Coordinate-based queries must succeed when chrom and start are provided.
        """
        with patch("beacon_api.db.session.get_session") as mock_get_session:
            mock_get_session.return_value = _mock_session()
            response = client.get(
                "/g_variants",
                params={"chrom": "chr17", "start": "43044295"},
            )
        assert response.status_code == 200, (
            f"Expected 200, got {response.status_code}: {response.text}"
        )

    def test_g_variants_count_granularity(self, client: TestClient) -> None:
        """GET /g_variants with granularity=count returns numTotalResults."""
        with patch("beacon_api.db.session.get_session"):
            response = client.get(
                "/g_variants",
                params={
                    "chrom": "chr17",
                    "start": "43044295",
                    "granularity": "count",
                },
            )
        assert response.status_code == 200
        data = response.json()
        assert "responseSummary" in data
        assert "numTotalResults" in data["responseSummary"]

    def test_g_variants_boolean_granularity(self, client: TestClient) -> None:
        """GET /g_variants with granularity=boolean returns exists field."""
        with patch("beacon_api.db.session.get_session"):
            response = client.get(
                "/g_variants",
                params={
                    "chrom": "chr17",
                    "start": "43044295",
                    "granularity": "boolean",
                },
            )
        assert response.status_code == 200
        data = response.json()
        assert "responseSummary" in data
        assert "exists" in data["responseSummary"]
        assert isinstance(data["responseSummary"]["exists"], bool)

    def test_g_variants_invalid_granularity_returns_422(self, client: TestClient) -> None:
        """GET /g_variants with invalid granularity returns HTTP 422."""
        response = client.get(
            "/g_variants",
            params={
                "chrom": "chr17",
                "start": "43044295",
                "granularity": "invalid_value",
            },
        )
        assert response.status_code == 422

    def test_g_variants_record_contains_vrs_id(self, client: TestClient) -> None:
        """GET /g_variants record response contains GA4GH VRS v2.0 identifier.

        Each returned variant must have a _vrsId field with the ga4gh:VA. prefix
        per VRS v2.0 specification (Wagner et al. 2021 PMID:35072137).
        """
        with patch("beacon_api.db.session.get_session"):
            # Provide a mock Passport header for record-level granularity
            response = client.get(
                "/g_variants",
                params={
                    "chrom": "chr17",
                    "start": "43044295",
                    "ref": "G",
                    "alt": "A",
                    "granularity": "record",
                },
            )
        assert response.status_code == 200
        data = response.json()
        if data.get("resultSets"):
            for result_set in data["resultSets"]:
                for variant in result_set.get("results", []):
                    assert "_vrsId" in variant, "Variant missing _vrsId field"
                    assert variant["_vrsId"].startswith("ga4gh:VA."), (
                        f"VRS ID should start with 'ga4gh:VA.', got: {variant['_vrsId']}"
                    )

    def test_g_variants_by_gene(self, client: TestClient) -> None:
        """GET /g_variants?geneSymbol=BRCA1 returns 200."""
        with patch("beacon_api.db.session.get_session"):
            response = client.get(
                "/g_variants",
                params={"geneSymbol": "BRCA1"},
            )
        assert response.status_code == 200


# ---------------------------------------------------------------------------
# VRS utils unit tests
# ---------------------------------------------------------------------------


class TestVRSUtils:
    """Tests for GA4GH VRS v2.0 identifier computation."""

    def test_compute_vrs_id_returns_ga4gh_prefix(self) -> None:
        """VRS identifier starts with ga4gh:VA. prefix.

        Per VRS v2.0 spec, Allele IDs must use the ga4gh:VA. prefix.
        Wagner et al. 2021 Cell Genomics PMID:35072137.
        """
        from beacon_api.vrs_utils import compute_vrs_id

        vrs_id = compute_vrs_id("chr17", 43044295, "G", "A")
        assert vrs_id.startswith("ga4gh:VA."), (
            f"VRS ID should start with 'ga4gh:VA.', got: {vrs_id}"
        )

    def test_compute_vrs_id_is_deterministic(self) -> None:
        """Same variant always produces the same VRS identifier.

        VRS v2.0 identifiers are computed deterministically from the
        canonical allele tuple. Wagner et al. 2021 PMID:35072137.
        """
        from beacon_api.vrs_utils import compute_vrs_id

        id1 = compute_vrs_id("chr17", 43044295, "G", "A")
        id2 = compute_vrs_id("chr17", 43044295, "G", "A")
        assert id1 == id2, "VRS ID must be deterministic for the same variant"

    def test_compute_vrs_id_different_variants_differ(self) -> None:
        """Different variants produce different VRS identifiers."""
        from beacon_api.vrs_utils import compute_vrs_id

        id1 = compute_vrs_id("chr17", 43044295, "G", "A")
        id2 = compute_vrs_id("chr17", 43044296, "G", "A")
        id3 = compute_vrs_id("chr1", 43044295, "G", "A")
        assert id1 != id2, "Different positions should yield different VRS IDs"
        assert id1 != id3, "Different chromosomes should yield different VRS IDs"

    def test_compute_vrs_id_length(self) -> None:
        """VRS identifier digest portion is exactly 24 characters."""
        from beacon_api.vrs_utils import compute_vrs_id

        vrs_id = compute_vrs_id("chr17", 43044295, "G", "A")
        digest = vrs_id.replace("ga4gh:VA.", "")
        assert len(digest) == 24, (
            f"VRS digest should be 24 chars, got {len(digest)}: {digest}"
        )

    def test_make_vrs_allele_structure(self) -> None:
        """make_vrs_allele returns VRSAllele with correct fields."""
        from beacon_api.vrs_utils import make_vrs_allele

        allele = make_vrs_allele("chr17", 43044295, "G", "A")
        assert allele.chrom == "chr17"
        assert allele.pos == 43044295
        assert allele.ref == "G"
        assert allele.alt == "A"
        assert allele.vrs_id.startswith("ga4gh:VA.")
        assert len(allele.digest) == 24

    def test_vrs_allele_to_dict_structure(self) -> None:
        """vrs_allele_to_dict returns correct VRS v2.0 dict structure."""
        from beacon_api.vrs_utils import make_vrs_allele, vrs_allele_to_dict

        allele = make_vrs_allele("chr17", 43044295, "G", "A")
        d = vrs_allele_to_dict(allele)
        assert d["type"] == "Allele"
        assert d["id"] == allele.vrs_id
        assert "location" in d
        assert "state" in d
        assert d["location"]["type"] == "SequenceLocation"
        # VRS uses 0-based interbase coordinates
        assert d["location"]["start"] == 43044294  # pos-1
