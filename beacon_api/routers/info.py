"""
beacon_api.routers.info
=======================
GET /info — GA4GH Beacon v2.1.1 metadata endpoint.

The /info endpoint returns the Beacon's metadata including its identifier,
name, API version, supported granularities, and dataset information.

Beacon v2.1.1 specification:
    https://github.com/ga4gh-beacon/beacon-v2/blob/main/framework/src/responses/beaconInfoResponse.yaml

References:
    Rambla et al. 2022 Human Mutation PMID:35297560 (Beacon v2 overview).
    GA4GH Beacon v2.1.1 released December 13, 2024.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from fastapi import APIRouter

router = APIRouter(prefix="/info")

# ---------------------------------------------------------------------------
# Static Beacon metadata — per GA4GH Beacon v2.1.1 schema
# ---------------------------------------------------------------------------

_BEACON_INFO: dict[str, Any] = {
    "id": "org.genomeforge.beacon",
    "name": "GenomeForge Genomic Beacon",
    "apiVersion": "v2.1.1",
    "environment": "prod",
    "organization": {
        "id": "org.genomeforge",
        "name": "GenomeForge",
        "description": "NHS GMS-compliant clinical whole-genome sequencing pipeline.",
        "address": "United Kingdom",
        "welcomeUrl": "https://github.com/genomeforge/genomeforge",
        "contactUrl": "https://github.com/genomeforge/genomeforge/issues",
        "logoUrl": "https://github.com/genomeforge/genomeforge/logo.png",
    },
    "description": (
        "GenomeForge GA4GH Beacon v2.1.1 providing access to variant-level "
        "genomic data from NHS GMS whole-genome sequencing. "
        "Variant identifiers use VRS v2.0 (Wagner et al. 2021 PMID:35072137). "
        "Population frequency from gnomAD v4.1 (April 2024, 807,162 individuals)."
    ),
    "version": "2.1.1",
    "welcomeUrl": "https://github.com/genomeforge/genomeforge",
    "alternativeUrl": None,
    "createDateTime": "2024-12-13T00:00:00Z",
    "updateDateTime": datetime.now(timezone.utc).isoformat(),
    "datasets": [
        {
            "id": "genomeforge.wgs.grch38",
            "name": "GenomeForge WGS GRCh38",
            "description": (
                "Whole-genome sequencing variants on GRCh38 reference assembly. "
                "Called with DRAGEN-GATK 4.6.0.0. "
                "Annotated with VEP 111, gnomAD v4.1, AlphaMissense."
            ),
            "assemblyId": "GRCh38",
            "createDateTime": "2024-12-13T00:00:00Z",
            "updateDateTime": datetime.now(timezone.utc).isoformat(),
            "version": "1.0",
            "variantCount": None,  # populated at runtime from DB
            "callCount": None,
            "sampleCount": None,
            "externalUrl": None,
            "info": {
                "caller": "DRAGEN-GATK 4.6.0.0",
                "annotator": "VEP 111",
                "gnomadVersion": "4.1",
                "variantIdScheme": "GA4GH VRS v2.0",
            },
        }
    ],
    "info": {
        "variantIdentifiers": "GA4GH VRS v2.0 (Wagner et al. 2021 PMID:35072137)",
        "authentication": "GA4GH Passport JWT (AAI Connect Profile)",
        "complianceStandards": [
            "NHS GMS",
            "ACGS 2024 v1.2",
            "GA4GH Beacon v2.1.1",
        ],
    },
}

# ---------------------------------------------------------------------------
# Beacon v2.1.1 granularity levels
# ---------------------------------------------------------------------------

_GRANULARITY = {
    "boolean": {
        "level": "boolean",
        "description": "YES/NO response indicating variant presence.",
    },
    "count": {
        "level": "count",
        "description": "Number of datasets/individuals with the queried variant.",
    },
    "record": {
        "level": "record",
        "description": "Full variant records with VRS v2.0 identifiers.",
        "authentication": "GA4GH Passport required for individual-level records.",
    },
}


@router.get(
    "",
    summary="Beacon metadata",
    description=(
        "Returns metadata about this Beacon: identifier, name, API version, "
        "supported granularities, datasets, and organisation information. "
        "Per GA4GH Beacon v2.1.1 specification (December 13, 2024)."
    ),
    response_model=None,
)
async def get_info() -> dict[str, Any]:
    """Return GA4GH Beacon v2.1.1 metadata.

    Returns:
        Dict conforming to the BeaconInfoResponse schema.  Contains Beacon
        identifier, API version, datasets, organisation metadata, and
        supported granularities.

    References:
        GA4GH Beacon v2.1.1 spec (December 13, 2024).
        Rambla et al. 2022 Human Mutation PMID:35297560.
    """
    return {
        "meta": {
            "beaconId": _BEACON_INFO["id"],
            "apiVersion": _BEACON_INFO["apiVersion"],
            "returnedSchemas": [
                {
                    "entityType": "info",
                    "schema": "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/main/framework/src/responses/beaconInfoResponse.yaml",
                }
            ],
        },
        "response": {
            **_BEACON_INFO,
            "granularity": _GRANULARITY,
        },
    }
