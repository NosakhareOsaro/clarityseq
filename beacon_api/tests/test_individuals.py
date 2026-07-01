"""
beacon_api.tests.test_individuals
==================================
pytest tests for GET /individuals — the GA4GH Passport-gated individual-level
query endpoint in beacon_api.routers.individuals.

Tests cover:
    - _build_individual_response(): field mapping for sex/ethnicity/diseases/
      phenotypic features.
    - The endpoint requires a valid GA4GH Passport (verify_passport
      dependency) — unauthenticated requests get 401.
    - Missing filter params (phenotypeId/diseaseId/sex) -> 400.
    - Visa-based access check branches (with/without ResearcherStatus or
      ControlledAccessGrants visa).
    - Pagination (skip/limit) and result shape.

Follows the TestClient + dependency_overrides pattern established in
test_beacon_api.py.
"""

from __future__ import annotations

from unittest.mock import AsyncMock

import pytest
from fastapi.testclient import TestClient

from beacon_api.routers.individuals import _build_individual_response


# ---------------------------------------------------------------------------
# App fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def client() -> TestClient:
    """TestClient with DB session mocked out; Passport auth left real."""
    from beacon_api.main import app
    from beacon_api.db.session import get_session

    async def _mock_get_session():  # type: ignore[return]
        yield AsyncMock()

    app.dependency_overrides[get_session] = _mock_get_session
    with TestClient(app, raise_server_exceptions=False) as c:
        yield c
    app.dependency_overrides.clear()


@pytest.fixture()
def authed_client(client: TestClient):
    """TestClient with verify_passport overridden to a valid claims dict.

    Yields a tuple of (client, claims_dict) so tests can mutate the visa
    list per-test before making requests.
    """
    from beacon_api.main import app
    from beacon_api.auth.passports import verify_passport

    claims: dict = {
        "sub": "test-subject",
        "iss": "https://login.elixir-czech.org/oidc/",
        "ga4gh_passport_v1": [{"type": "ResearcherStatus"}],
    }

    async def _mock_verify_passport():
        return claims

    app.dependency_overrides[verify_passport] = _mock_verify_passport
    yield client, claims
    del app.dependency_overrides[verify_passport]


# ---------------------------------------------------------------------------
# _build_individual_response unit tests
# ---------------------------------------------------------------------------


class TestBuildIndividualResponse:
    def test_female_sex_maps_to_gsso_009523(self) -> None:
        """FEMALE sex maps to GSSO:009523."""
        row = {"individual_id": "IND-1", "sex": "FEMALE"}
        result = _build_individual_response(row)
        assert result["sex"]["id"] == "GSSO:009523"
        assert result["sex"]["label"] == "FEMALE"

    def test_non_female_sex_maps_to_gsso_009521(self) -> None:
        """MALE (or any non-FEMALE) sex maps to GSSO:009521."""
        row = {"individual_id": "IND-2", "sex": "MALE"}
        result = _build_individual_response(row)
        assert result["sex"]["id"] == "GSSO:009521"

    def test_missing_sex_defaults_to_unknown(self) -> None:
        """Missing sex field defaults label to UNKNOWN_SEX."""
        row = {"individual_id": "IND-3"}
        result = _build_individual_response(row)
        assert result["sex"]["label"] == "UNKNOWN_SEX"

    def test_diseases_and_phenotypic_features_mapped(self) -> None:
        """Diseases and phenotypic features are mapped into Beacon schema shape."""
        row = {
            "individual_id": "IND-4",
            "sex": "FEMALE",
            "ethnicity": "Not specified",
            "ethnicity_id": "HANCESTRO:0004",
            "diseases": [{"omim_id": "OMIM:114480", "label": "Breast cancer", "stage": "IV"}],
            "phenotypic_features": [{"hpo_id": "HP:0001250", "label": "Seizure", "excluded": False}],
        }
        result = _build_individual_response(row)
        assert result["diseases"][0]["diseaseCode"]["id"] == "OMIM:114480"
        assert result["diseases"][0]["stage"] == "IV"
        assert result["phenotypicFeatures"][0]["featureType"]["id"] == "HP:0001250"
        assert result["phenotypicFeatures"][0]["excluded"] is False
        assert result["info"]["dataAccessLevel"] == "controlled"
        assert result["info"]["requiresPassport"] is True

    def test_empty_diseases_and_features_default_to_empty_lists(self) -> None:
        """Absent diseases/phenotypic_features keys produce empty lists, not errors."""
        row = {"individual_id": "IND-5"}
        result = _build_individual_response(row)
        assert result["diseases"] == []
        assert result["phenotypicFeatures"] == []


# ---------------------------------------------------------------------------
# Endpoint authentication tests
# ---------------------------------------------------------------------------


class TestIndividualsAuth:
    def test_no_authorization_header_returns_401(self, client: TestClient) -> None:
        """Without a Bearer token, verify_passport rejects with 401."""
        response = client.get("/individuals", params={"sex": "FEMALE"})
        assert response.status_code == 401

    def test_invalid_bearer_token_returns_401(self, client: TestClient) -> None:
        """A malformed bearer token is rejected with 401."""
        response = client.get(
            "/individuals",
            params={"sex": "FEMALE"},
            headers={"Authorization": "Bearer not-a-real-jwt"},
        )
        assert response.status_code == 401


# ---------------------------------------------------------------------------
# Endpoint behaviour tests (authenticated)
# ---------------------------------------------------------------------------


class TestIndividualsQuery:
    def test_missing_filters_returns_400(self, authed_client) -> None:
        """No phenotypeId/diseaseId/sex filter provided -> 400."""
        client, _claims = authed_client
        response = client.get("/individuals")
        assert response.status_code == 400

    def test_phenotype_filter_returns_200_with_results(self, authed_client) -> None:
        """phenotypeId filter returns a mock individual with that HPO term."""
        client, _claims = authed_client
        response = client.get("/individuals", params={"phenotypeId": "HP:0001250"})
        assert response.status_code == 200
        data = response.json()
        assert data["responseSummary"]["numTotalResults"] == 1
        assert data["responseSummary"]["exists"] is True
        results = data["resultSets"][0]["results"]
        assert len(results) == 1
        assert results[0]["phenotypicFeatures"][0]["featureType"]["id"] == "HP:0001250"
        assert data["_passportSubject"] == "test-subject"

    def test_disease_only_filter_returns_no_mock_results(self, authed_client) -> None:
        """diseaseId-only filter passes the >=1 filter check but yields zero mock rows."""
        client, _claims = authed_client
        response = client.get("/individuals", params={"diseaseId": "OMIM:114480"})
        assert response.status_code == 200
        data = response.json()
        assert data["responseSummary"]["numTotalResults"] == 0
        assert data["responseSummary"]["exists"] is False

    def test_sex_only_filter_returns_200(self, authed_client) -> None:
        """sex-only filter is accepted as a valid filter."""
        client, _claims = authed_client
        response = client.get("/individuals", params={"sex": "MALE"})
        assert response.status_code == 200

    def test_pagination_skip_beyond_results(self, authed_client) -> None:
        """skip beyond the single mock result yields an empty page."""
        client, _claims = authed_client
        response = client.get(
            "/individuals", params={"phenotypeId": "HP:0001250", "skip": 5}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["resultSets"][0]["results"] == []
        assert data["resultSets"][0]["resultsCount"] == 0

    def test_response_meta_schema(self, authed_client) -> None:
        """Response meta block conforms to Beacon v2.1.1 individuals schema."""
        client, _claims = authed_client
        response = client.get("/individuals", params={"phenotypeId": "HP:0001250"})
        data = response.json()
        assert data["meta"]["beaconId"] == "org.clarityseq.beacon"
        assert data["meta"]["apiVersion"] == "v2.1.1"
        assert data["meta"]["returnedSchemas"][0]["entityType"] == "individuals"

    def test_no_visas_denied_403(self, authed_client) -> None:
        """Passport with empty visa list is rejected with 403 Forbidden.

        beacon_api/routers/individuals.py enforces that the Passport must
        carry a ResearcherStatus or ControlledAccessGrants visa before
        individual-level (PHI) records are returned. A validly-signed
        Passport with no visas is authenticated but not authorized.
        """
        client, claims = authed_client
        claims["ga4gh_passport_v1"] = []  # no visas at all
        response = client.get("/individuals", params={"phenotypeId": "HP:0001250"})
        assert response.status_code == 403

    def test_visa_without_recognized_type_denied_403(self, authed_client) -> None:
        """A visa present but of an unrecognized type also fails the has_access check."""
        client, claims = authed_client
        claims["ga4gh_passport_v1"] = [{"type": "AffiliationAndRole"}]
        response = client.get("/individuals", params={"phenotypeId": "HP:0001250"})
        assert response.status_code == 403

    def test_controlled_access_grants_visa_allows_access(self, authed_client) -> None:
        """A ControlledAccessGrants visa satisfies has_access=True."""
        client, claims = authed_client
        claims["ga4gh_passport_v1"] = [{"type": "ControlledAccessGrants"}]
        response = client.get("/individuals", params={"phenotypeId": "HP:0001250"})
        assert response.status_code == 200
