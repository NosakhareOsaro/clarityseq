"""
beacon_api.tests.test_passports
================================
pytest tests for GA4GH Passport JWT validation in beacon_api.auth.passports.

Tests cover:
    - JWKS fetching/caching (_fetch_jwks).
    - Public key resolution (_get_public_key): issuer allowlist, JWKS fetch
      failures, missing kid match, PyJWT-not-installed fallback.
    - verify_passport(): missing header, expired/malformed/wrong-audience
      JWTs, valid Passport acceptance.
    - optional_passport(): None on missing/invalid token, claims on valid.

Real RSA keypairs (via `cryptography`) and real PyJWT encode/decode are used
so that signature verification genuinely exercises the RS256 path; only the
network call to the JWKS endpoint (httpx.get) is mocked with unittest.mock.
"""

from __future__ import annotations

import time
from unittest.mock import MagicMock, patch

import httpx
import jwt as pyjwt
import pytest
from cryptography.hazmat.primitives.asymmetric import rsa
from fastapi.security import HTTPAuthorizationCredentials
from jwt.algorithms import RSAAlgorithm

from beacon_api.auth import passports

ALLOWED_ISSUER = passports._ALLOWED_ISSUERS[0]
AUDIENCE = passports._EXPECTED_AUDIENCE
KID = "test-key-1"


@pytest.fixture(autouse=True)
def _clear_jwks_cache():
    """Ensure the lru_cache on _fetch_jwks doesn't leak across tests."""
    passports._fetch_jwks.cache_clear()
    yield
    passports._fetch_jwks.cache_clear()


@pytest.fixture(scope="module")
def rsa_private_key():
    """Generate a throwaway RSA keypair for signing test Passport JWTs."""
    return rsa.generate_private_key(public_exponent=65537, key_size=2048)


@pytest.fixture(scope="module")
def jwks_dict(rsa_private_key):
    """Build a JWKS document containing the public half of the test key."""
    jwk_json = RSAAlgorithm.to_jwk(rsa_private_key.public_key())
    import json

    jwk = json.loads(jwk_json)
    jwk["kid"] = KID
    jwk["use"] = "sig"
    jwk["alg"] = "RS256"
    return {"keys": [jwk]}


def _make_token(private_key, *, kid=KID, issuer=ALLOWED_ISSUER, audience=AUDIENCE,
                 exp_delta=3600, extra_claims=None):
    """Sign an RS256 JWT with the given claims for use as a test Passport."""
    now = int(time.time())
    payload = {
        "sub": "test-user-001",
        "iss": issuer,
        "aud": audience,
        "iat": now,
        "exp": now + exp_delta,
    }
    if extra_claims:
        payload.update(extra_claims)
    headers = {"kid": kid} if kid else {}
    return pyjwt.encode(payload, private_key, algorithm="RS256", headers=headers)


def _mock_httpx_response(json_data):
    """Return a MagicMock standing in for an httpx.Response."""
    resp = MagicMock()
    resp.json.return_value = json_data
    resp.raise_for_status.return_value = None
    return resp


# ---------------------------------------------------------------------------
# _fetch_jwks
# ---------------------------------------------------------------------------


class TestFetchJwks:
    def test_fetches_and_returns_jwks(self, jwks_dict) -> None:
        """_fetch_jwks calls httpx.get and returns the parsed JSON body."""
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)) as mock_get:
            result = passports._fetch_jwks("https://issuer.example.com/.well-known/jwks.json")
        mock_get.assert_called_once()
        assert result == jwks_dict

    def test_caches_repeated_calls(self, jwks_dict) -> None:
        """A second call with the same URI does not re-hit the network (lru_cache)."""
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)) as mock_get:
            passports._fetch_jwks("https://issuer.example.com/cached/jwks.json")
            passports._fetch_jwks("https://issuer.example.com/cached/jwks.json")
        assert mock_get.call_count == 1


# ---------------------------------------------------------------------------
# _get_public_key
# ---------------------------------------------------------------------------


class TestGetPublicKey:
    def test_raises_503_when_pyjwt_unavailable(self) -> None:
        """When PyJWT is not installed, a 503 is raised."""
        with patch.object(passports, "_JWT_AVAILABLE", False):
            with pytest.raises(Exception) as exc_info:
                passports._get_public_key("irrelevant")
        assert exc_info.value.status_code == 503

    def test_raises_401_for_malformed_token(self) -> None:
        """A non-JWT string raises 401 'Invalid Passport JWT format'."""
        with pytest.raises(Exception) as exc_info:
            passports._get_public_key("not-a-valid-jwt-string")
        assert exc_info.value.status_code == 401
        assert "Invalid Passport JWT format" in exc_info.value.detail

    def test_raises_401_for_disallowed_issuer(self, rsa_private_key) -> None:
        """A JWT with an issuer outside the allowlist is rejected with 401."""
        token = _make_token(rsa_private_key, issuer="https://evil.example.com/oidc/")
        with pytest.raises(Exception) as exc_info:
            passports._get_public_key(token)
        assert exc_info.value.status_code == 401
        assert "not in the allowed issuers list" in exc_info.value.detail

    def test_raises_401_on_jwks_fetch_network_error(self, rsa_private_key) -> None:
        """If the JWKS endpoint is unreachable, 401 is raised."""
        token = _make_token(rsa_private_key)
        with patch.object(passports.httpx, "get", side_effect=httpx.RequestError("boom")):
            with pytest.raises(Exception) as exc_info:
                passports._get_public_key(token)
        assert exc_info.value.status_code == 401
        assert "Failed to fetch JWKS" in exc_info.value.detail

    def test_raises_401_on_jwks_http_status_error(self, rsa_private_key) -> None:
        """If the JWKS endpoint returns a non-2xx status, 401 is raised."""
        token = _make_token(rsa_private_key)
        resp = MagicMock()
        resp.raise_for_status.side_effect = httpx.HTTPStatusError(
            "500", request=MagicMock(), response=MagicMock()
        )
        with patch.object(passports.httpx, "get", return_value=resp):
            with pytest.raises(Exception) as exc_info:
                passports._get_public_key(token)
        assert exc_info.value.status_code == 401
        assert "Failed to fetch JWKS" in exc_info.value.detail

    def test_raises_401_when_no_matching_kid(self, rsa_private_key, jwks_dict) -> None:
        """If no JWKS key matches the token's kid, 401 is raised."""
        token = _make_token(rsa_private_key, kid="does-not-exist")
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            with pytest.raises(Exception) as exc_info:
                passports._get_public_key(token)
        assert exc_info.value.status_code == 401
        assert "No matching public key" in exc_info.value.detail

    def test_returns_key_on_kid_match(self, rsa_private_key, jwks_dict) -> None:
        """A token whose kid matches a JWKS entry resolves to a usable public key."""
        token = _make_token(rsa_private_key)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            key = passports._get_public_key(token)
        assert key is not None

    def test_returns_key_when_token_has_no_kid(self, rsa_private_key, jwks_dict) -> None:
        """When the token omits kid, the first JWKS key is used."""
        token = _make_token(rsa_private_key, kid=None)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            key = passports._get_public_key(token)
        assert key is not None


# ---------------------------------------------------------------------------
# verify_passport
# ---------------------------------------------------------------------------


class TestVerifyPassport:
    @pytest.mark.asyncio
    async def test_missing_credentials_raises_401(self) -> None:
        """No Authorization header -> 401 with WWW-Authenticate header."""
        with pytest.raises(Exception) as exc_info:
            await passports.verify_passport(credentials=None)
        assert exc_info.value.status_code == 401
        assert exc_info.value.headers["WWW-Authenticate"] == "Bearer"

    @pytest.mark.asyncio
    async def test_pyjwt_unavailable_returns_mock_claims(self) -> None:
        """When PyJWT isn't installed, dev-mode mock claims are returned."""
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials="whatever")
        with patch.object(passports, "_JWT_AVAILABLE", False):
            claims = await passports.verify_passport(credentials=creds)
        assert claims["sub"] == "dev-user"
        assert claims["ga4gh_passport_v1"] == []

    @pytest.mark.asyncio
    async def test_expired_token_raises_401(self, rsa_private_key, jwks_dict) -> None:
        """An expired JWT raises 401 'has expired'."""
        token = _make_token(rsa_private_key, exp_delta=-3600)
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials=token)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            with pytest.raises(Exception) as exc_info:
                await passports.verify_passport(credentials=creds)
        assert exc_info.value.status_code == 401
        assert "expired" in exc_info.value.detail

    @pytest.mark.asyncio
    async def test_wrong_audience_raises_401(self, rsa_private_key, jwks_dict) -> None:
        """A JWT with the wrong audience raises 401 'audience mismatch'."""
        token = _make_token(rsa_private_key, audience="some.other.beacon")
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials=token)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            with pytest.raises(Exception) as exc_info:
                await passports.verify_passport(credentials=creds)
        assert exc_info.value.status_code == 401
        assert "audience mismatch" in exc_info.value.detail

    @pytest.mark.asyncio
    async def test_bad_signature_raises_generic_401(self, rsa_private_key, jwks_dict) -> None:
        """A JWT signed with a different key fails signature verification."""
        other_key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
        token = _make_token(other_key)  # signed with a key NOT in the JWKS
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials=token)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            with pytest.raises(Exception) as exc_info:
                await passports.verify_passport(credentials=creds)
        assert exc_info.value.status_code == 401
        assert "validation failed" in exc_info.value.detail

    @pytest.mark.asyncio
    async def test_valid_passport_returns_claims(self, rsa_private_key, jwks_dict) -> None:
        """A valid, signed, non-expired, correct-audience Passport is accepted."""
        token = _make_token(
            rsa_private_key,
            extra_claims={
                "ga4gh_passport_v1": [{"type": "ResearcherStatus"}],
            },
        )
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials=token)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            claims = await passports.verify_passport(credentials=creds)
        assert claims["sub"] == "test-user-001"
        assert claims["iss"] == ALLOWED_ISSUER
        assert claims["ga4gh_passport_v1"][0]["type"] == "ResearcherStatus"


# ---------------------------------------------------------------------------
# optional_passport
# ---------------------------------------------------------------------------


class TestOptionalPassport:
    @pytest.mark.asyncio
    async def test_none_credentials_returns_none(self) -> None:
        """No Authorization header -> None (no exception raised)."""
        result = await passports.optional_passport(credentials=None)
        assert result is None

    @pytest.mark.asyncio
    async def test_valid_credentials_returns_claims(self, rsa_private_key, jwks_dict) -> None:
        """A valid Passport returns the decoded claims dict."""
        token = _make_token(rsa_private_key)
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials=token)
        with patch.object(passports.httpx, "get", return_value=_mock_httpx_response(jwks_dict)):
            result = await passports.optional_passport(credentials=creds)
        assert result is not None
        assert result["sub"] == "test-user-001"

    @pytest.mark.asyncio
    async def test_invalid_credentials_returns_none(self) -> None:
        """An invalid/malformed token results in None rather than raising."""
        creds = HTTPAuthorizationCredentials(scheme="Bearer", credentials="not-a-jwt")
        result = await passports.optional_passport(credentials=creds)
        assert result is None
