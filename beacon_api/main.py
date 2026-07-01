"""
beacon_api.main
===============
GA4GH Beacon v2.1.1 API (released December 13, 2024).
VRS v2.0 variant identifiers. GA4GH Passport authentication.

Beacon v2.1.1 specification:
    https://github.com/ga4gh-beacon/beacon-v2
    Released December 13, 2024.

GA4GH Passport:
    https://github.com/ga4gh/data-security/blob/master/AAI/AAIConnectProfile.md
    JWT-based access control for sensitive genomic data.

VRS v2.0:
    GA4GH Variant Representation Specification v2.0.
    Wagner et al. 2021 Cell Genomics PMID:35072137.
    https://vrs.ga4gh.org/
"""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from beacon_api.routers import info, g_variants, individuals

app = FastAPI(
    title="ClaritySeq Beacon API",
    description=(
        "GA4GH Beacon v2.1.1 with VRS v2.0 identifiers. "
        "Implements the Beacon v2 specification (December 13, 2024)."
    ),
    version="2.1.1",
    contact={
        "name": "ClaritySeq",
        "url": "https://github.com/clarityseq/clarityseq",
    },
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0",
    },
)

# ---------------------------------------------------------------------------
# CORS middleware
# ---------------------------------------------------------------------------
# Allow cross-origin requests from trusted GA4GH portal origins.
# In production, replace ["*"] with specific allowed origins.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["GET", "POST"],
    allow_headers=["Authorization", "Content-Type"],
)

# ---------------------------------------------------------------------------
# Routers
# ---------------------------------------------------------------------------
app.include_router(info.router, tags=["Info"])
app.include_router(g_variants.router, tags=["Genomic Variants"])
app.include_router(individuals.router, tags=["Individuals"])


@app.get("/", include_in_schema=False)
async def root() -> dict[str, str]:
    """Redirect hint for API root.

    Returns:
        Dict with link to /info endpoint.
    """
    return {"message": "GA4GH Beacon v2.1.1. See /info for metadata or /docs for OpenAPI."}
