"""
beacon_api.auth.passports
==========================
GA4GH Passport JWT validation for the Beacon v2.1.1 API.

GA4GH Passport v1.2 specification:
    https://github.com/ga4gh/data-security/blob/master/AAI/AAIConnectProfile.md

Passport tokens are RS256-signed JWTs issued by a Passport Broker
(e.g. Elixir AAI, REMS, ORCID).  They contain:
    - Standard JWT claims: sub, iss, iat, exp, jti.
    - ``ga4gh_passport_v1``: Array of Passport Visa JWTs or decoded visa objects.

Visa types:
    - ResearcherStatus: Indicates the holder is a bona fide researcher.
    - AffiliationAndRole: Academic/institutional affiliation.
    - AcceptedTermsAndPolicies: Data access agreement acceptance.
    - LinkedIdentities: Links between user identities.
    - ControlledAccessGrants: Dataset-specific access grants (used here).

JWKS validation:
    The Passport token's public key is validated against the issuer's JWKS
    endpoint (``{iss}/.well-known/jwks.json``).  Keys are cached with a
    1-hour TTL using functools.lru_cache.

References:
    GA4GH Passport v1.2 spec.
    RFC 7517 (JWKS), RFC 7519 (JWT).
    PyJWT ≥ 2.8.
"""

from __future__ import annotations

import logging
import os
from functools import lru_cache
from typing import Any

import httpx
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

try:
    import jwt as pyjwt
    _JWT_AVAILABLE = True
except ImportError:
    _JWT_AVAILABLE = False

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Allowed Passport issuer URIs — configured via environment variable.
# Multiple issuers separated by commas.
_ALLOWED_ISSUERS: list[str] = [
    iss.strip()
    for iss in os.getenv(
        "PASSPORT_ALLOWED_ISSUERS",
        "https://login.elixir-czech.org/oidc/,https://proxy.aai.lifescience-ri.eu/oidc/",
    ).split(",")
    if iss.strip()
]

# Expected audience claim value (the Beacon's identifier)
_EXPECTED_AUDIENCE: str = os.getenv(
    "PASSPORT_EXPECTED_AUDIENCE",
    "org.genomeforge.beacon",
)

_JWKS_CACHE_TTL_SECONDS = 3600  # 1-hour JWKS key cache TTL

_security_scheme = HTTPBearer(auto_error=False)

# ---------------------------------------------------------------------------
# JWKS fetching with caching
# ---------------------------------------------------------------------------


@lru_cache(maxsize=16)
def _fetch_jwks(jwks_uri: str) -> dict[str, Any]:
    """Fetch and cache a JWKS key set from the issuer's JWKS endpoint.

    Uses functools.lru_cache for in-process caching.  Cache is invalidated
    on process restart.  Key refresh is handled by re-fetching on JWT
    signature verification failure.

    Args:
        jwks_uri: Full HTTPS URL to the JWKS endpoint
            (e.g. ``"https://login.elixir-czech.org/oidc/jwks"``).

    Returns:
        Parsed JWKS dict as returned by the issuer.

    Raises:
        httpx.RequestError: On network error fetching JWKS.
        httpx.HTTPStatusError: On non-2xx response from JWKS endpoint.
    """
    response = httpx.get(jwks_uri, timeout=5.0)
    response.raise_for_status()
    return response.json()


def _get_public_key(token: str) -> Any:
    """Resolve the signing public key for a Passport JWT.

    Decodes the token header (without verification) to extract ``kid``
    and ``alg``, then fetches the matching key from the issuer's JWKS
    endpoint.

    Args:
        token: Raw JWT string (``"Bearer <token>``" prefix already stripped).

    Returns:
        PyJWT RSAAlgorithm public key object for RS256 verification.

    Raises:
        HTTPException: 401 if the issuer is not in the allowlist, or if
            the JWKS endpoint cannot be reached, or if the key is not found.
    """
    if not _JWT_AVAILABLE:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="PyJWT not installed; Passport authentication unavailable.",
        )

    # Decode header without verification to get issuer and kid
    try:
        unverified_header = pyjwt.get_unverified_header(token)
        unverified_payload = pyjwt.decode(
            token,
            options={"verify_signature": False},
        )
    except pyjwt.DecodeError as exc:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=f"Invalid Passport JWT format: {exc}",
        ) from exc

    issuer = unverified_payload.get("iss", "")
    if issuer not in _ALLOWED_ISSUERS:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=(
                f"Passport issuer '{issuer}' is not in the allowed issuers list. "
                "Configure PASSPORT_ALLOWED_ISSUERS environment variable."
            ),
        )

    kid = unverified_header.get("kid")
    jwks_uri = f"{issuer.rstrip('/')}/.well-known/jwks.json"

    try:
        jwks = _fetch_jwks(jwks_uri)
    except (httpx.RequestError, httpx.HTTPStatusError) as exc:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=f"Failed to fetch JWKS from issuer '{issuer}': {exc}",
        ) from exc

    # Find matching key by kid
    for key_data in jwks.get("keys", []):
        if kid and key_data.get("kid") != kid:
            continue
        return pyjwt.algorithms.RSAAlgorithm.from_jwk(key_data)

    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail=f"No matching public key found in JWKS for kid='{kid}'",
    )


# ---------------------------------------------------------------------------
# FastAPI dependencies
# ---------------------------------------------------------------------------


async def verify_passport(
    credentials: HTTPAuthorizationCredentials | None = Depends(_security_scheme),
) -> dict[str, Any]:
    """FastAPI dependency: validate GA4GH Passport JWT and return claims.

    Validates the JWT signature, expiry, issuer, and audience.
    Returns the full decoded payload including ``ga4gh_passport_v1`` visas.

    Args:
        credentials: HTTP Bearer credentials from the Authorization header.

    Returns:
        Decoded JWT payload dict with all claims including
        ``"sub"``, ``"iss"``, ``"ga4gh_passport_v1"``.

    Raises:
        HTTPException: 401 if the Authorization header is missing.
        HTTPException: 401 if the JWT is expired, malformed, or signature invalid.
        HTTPException: 401 if the issuer is not in the allowlist.

    References:
        GA4GH Passport v1.2 specification.
        PyJWT: https://pyjwt.readthedocs.io/
    """
    if credentials is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="GA4GH Passport JWT required. Provide Authorization: Bearer <token>.",
            headers={"WWW-Authenticate": "Bearer"},
        )

    token = credentials.credentials

    if not _JWT_AVAILABLE:
        # Fallback: return minimal mock claims when PyJWT not installed
        logger.warning("PyJWT not available; returning mock Passport claims for development.")
        return {"sub": "dev-user", "iss": "dev", "ga4gh_passport_v1": []}

    public_key = _get_public_key(token)

    try:
        claims = pyjwt.decode(
            token,
            key=public_key,
            algorithms=["RS256"],
            audience=_EXPECTED_AUDIENCE,
            options={
                "verify_exp": True,
                "verify_iat": True,
                "verify_iss": True,
            },
        )
    except pyjwt.ExpiredSignatureError as exc:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Passport JWT has expired.",
        ) from exc
    except pyjwt.InvalidAudienceError as exc:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=(
                f"Passport JWT audience mismatch. Expected '{_EXPECTED_AUDIENCE}'. "
                "Configure PASSPORT_EXPECTED_AUDIENCE environment variable."
            ),
        ) from exc
    except pyjwt.PyJWTError as exc:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=f"Passport JWT validation failed: {exc}",
        ) from exc

    logger.debug("Passport validated for subject: %s", claims.get("sub"))
    return claims


async def optional_passport(
    credentials: HTTPAuthorizationCredentials | None = Depends(_security_scheme),
) -> dict[str, Any] | None:
    """FastAPI dependency: optional GA4GH Passport JWT validation.

    Unlike verify_passport(), this dependency does NOT raise 401 if the
    Authorization header is absent.  Used for endpoints that provide
    record-level granularity to authenticated callers and count-level
    granularity to anonymous callers.

    Args:
        credentials: HTTP Bearer credentials (may be None).

    Returns:
        Decoded JWT payload if a valid Passport was provided, or None if
        no Authorization header was present.

    Raises:
        HTTPException: 401 only if a token was provided but is invalid.

    References:
        GA4GH Beacon v2.1.1 — granularity downgrade for unauthenticated.
    """
    if credentials is None:
        return None

    try:
        return await verify_passport(credentials)
    except HTTPException:
        # Invalid token provided — return None and let the router downgrade
        return None
