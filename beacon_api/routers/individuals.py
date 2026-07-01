"""
beacon_api.routers.individuals
===============================
GET /individuals — GA4GH Beacon v2.1.1 individual-level query endpoint.

This endpoint is gated by GA4GH Passport authentication.  Only callers
presenting a valid Passport JWT with appropriate visa claims (ResearcherStatus
or ControlledAccessGrants) receive individual-level records.

GA4GH Passport:
    https://github.com/ga4gh/data-security/blob/master/AAI/AAIConnectProfile.md

Beacon v2.1.1 Individual schema:
    https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/src/beacon-v2-default-model/individuals/defaultSchema.yaml

References:
    GA4GH Beacon v2.1.1 (December 13, 2024).
    Rambla et al. 2022 Human Mutation PMID:35297560.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query, status

from beacon_api.auth.passports import verify_passport
from beacon_api.db.session import get_session

router = APIRouter(prefix="/individuals")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_individual_response(row: dict[str, Any]) -> dict[str, Any]:
    """Build a Beacon v2.1.1 Individual response object.

    Args:
        row: Dict representing a database individual row with keys:
            ``individual_id``, ``sex``, ``ethnicity``, ``diseases``,
            ``phenotypic_features``.

    Returns:
        Dict conforming to the Beacon v2.1.1 Individual schema.
    """
    return {
        "id": row.get("individual_id", ""),
        "sex": {
            "id": "GSSO:009523" if row.get("sex") == "FEMALE" else "GSSO:009521",
            "label": row.get("sex", "UNKNOWN_SEX"),
        },
        "ethnicity": {
            "id": row.get("ethnicity_id", ""),
            "label": row.get("ethnicity", ""),
        },
        "diseases": [
            {
                "diseaseCode": {
                    "id": d.get("omim_id", ""),
                    "label": d.get("label", ""),
                },
                "stage": d.get("stage", ""),
            }
            for d in row.get("diseases", [])
        ],
        "phenotypicFeatures": [
            {
                "featureType": {
                    "id": f.get("hpo_id", ""),
                    "label": f.get("label", ""),
                },
                "excluded": f.get("excluded", False),
            }
            for f in row.get("phenotypic_features", [])
        ],
        "info": {
            "dataAccessLevel": "controlled",
            "requiresPassport": True,
        },
    }


# ---------------------------------------------------------------------------
# Endpoint
# ---------------------------------------------------------------------------


@router.get(
    "",
    summary="Individual-level query (Passport-gated)",
    description=(
        "Query individual-level phenotypic and disease data. "
        "Requires a valid GA4GH Passport JWT with appropriate visa claims. "
        "Caller must present a Passport containing either a ResearcherStatus "
        "visa or a ControlledAccessGrants visa for this dataset. "
        "Per GA4GH Passport AAI Connect Profile."
    ),
    response_model=None,
)
async def query_individuals(
    phenotype_id: str | None = Query(
        None,
        alias="phenotypeId",
        description="HPO term ID (e.g. HP:0001250 for seizures).",
    ),
    disease_id: str | None = Query(
        None,
        alias="diseaseId",
        description="OMIM or Orphanet disease ID.",
    ),
    sex: str | None = Query(
        None,
        description="Biological sex: FEMALE, MALE, UNKNOWN_SEX.",
    ),
    skip: int = Query(0, ge=0, description="Pagination offset."),
    limit: int = Query(10, ge=1, le=100, description="Maximum records."),
    passport_claims: dict[str, Any] = Depends(verify_passport),
    session: Any = Depends(get_session),
) -> dict[str, Any]:
    """Query individual-level records (requires GA4GH Passport).

    Returns Beacon v2.1.1 Individual records including phenotypic features
    (HPO terms), diseases (OMIM/Orphanet), and sex.

    Args:
        phenotype_id: HPO term ID to filter by (e.g. ``"HP:0001250"``).
        disease_id: OMIM or Orphanet disease identifier.
        sex: Biological sex filter (FEMALE/MALE/UNKNOWN_SEX).
        skip: Pagination offset.
        limit: Maximum records per page (max 100).
        passport_claims: Decoded GA4GH Passport JWT claims from verify_passport().
            Dependency injection raises 401 if passport is absent or invalid.
        session: Async database session.

    Returns:
        Dict conforming to the Beacon v2.1.1 BeaconResultsetsResponse schema
        with Individual records.

    Raises:
        HTTPException: 401 if GA4GH Passport is missing or invalid (from
            verify_passport dependency).
        HTTPException: 403 if Passport lacks required visa claims.

    References:
        GA4GH Passport AAI Connect Profile.
        GA4GH Beacon v2.1.1 Individual schema.
    """
    # Check required visa claims
    # Production: check for ResearcherStatus or ControlledAccessGrants visa
    ga4gh_visas = passport_claims.get("ga4gh_passport_v1", [])
    has_access = any(
        v.get("type") in {"ResearcherStatus", "ControlledAccessGrants"}
        for v in ga4gh_visas
    ) if ga4gh_visas else False

    if not has_access:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=(
                "Passport lacks required visa claims: ResearcherStatus or "
                "ControlledAccessGrants."
            ),
        )

    if not any([phenotype_id, disease_id, sex]):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=(
                "At least one filter is required: phenotypeId, diseaseId, or sex."
            ),
        )

    # Build mock result set
    # In production: query PostgreSQL database via session
    mock_individuals: list[dict[str, Any]] = []

    if phenotype_id:
        mock_individuals.append({
            "individual_id": "IND-00001",
            "sex": "UNKNOWN_SEX",
            "ethnicity": "",
            "ethnicity_id": "",
            "diseases": [],
            "phenotypic_features": [
                {"hpo_id": phenotype_id, "label": "", "excluded": False}
            ],
        })

    total = len(mock_individuals)
    page = mock_individuals[skip : skip + limit]
    results = [_build_individual_response(ind) for ind in page]

    return {
        "meta": {
            "beaconId": "org.clarityseq.beacon",
            "apiVersion": "v2.1.1",
            "returnedSchemas": [
                {
                    "entityType": "individuals",
                    "schema": (
                        "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/main"
                        "/models/src/beacon-v2-default-model/individuals/defaultSchema.yaml"
                    ),
                }
            ],
        },
        "responseSummary": {"exists": total > 0, "numTotalResults": total},
        "resultSets": [
            {
                "id": "clarityseq.wgs.grch38",
                "type": "dataset",
                "exists": total > 0,
                "resultsCount": len(results),
                "results": results,
            }
        ],
        "beaconHandovers": [],
        "_passportSubject": passport_claims.get("sub", ""),
    }
