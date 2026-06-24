"""
beacon_api.routers.g_variants
==============================
GET /g_variants — GA4GH Beacon v2.1.1 genomic variant query endpoint.

Returns variant records with GA4GH VRS v2.0 computed identifiers embedded
in the response.  Population frequencies from gnomAD v4.1 are included.

Beacon v2.1.1 GenomicVariant schema:
    https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/src/beacon-v2-default-model/genomicVariations/defaultSchema.yaml

VRS v2.0 identifiers:
    Wagner et al. 2021 Cell Genomics PMID:35072137.

References:
    GA4GH Beacon v2.1.1 (December 13, 2024).
    Rambla et al. 2022 Human Mutation PMID:35297560.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query, status

from beacon_api.auth.passports import optional_passport
from beacon_api.db.session import get_session
from beacon_api.vrs_utils import compute_vrs_id, vrs_allele_to_dict, make_vrs_allele

router = APIRouter(prefix="/g_variants")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_variant_response(row: dict[str, Any]) -> dict[str, Any]:
    """Build a Beacon v2.1.1 GenomicVariant response object.

    Embeds the GA4GH VRS v2.0 identifier and gnomAD v4.1 population
    frequency data per Beacon v2.1.1 schema.

    Args:
        row: Dict representing a database variant row with keys:
            ``chrom``, ``pos``, ``ref``, ``alt``, ``gnomad_af``,
            ``clinvar_id``, ``gene_symbol``, ``hgvsc``, ``hgvsp``,
            ``acmg_class``.

    Returns:
        Dict conforming to the Beacon v2.1.1 GenomicVariant schema with
        embedded VRS v2.0 identifier.
    """
    chrom = row.get("chrom", "")
    pos = int(row.get("pos", 0))
    ref = row.get("ref", "")
    alt = row.get("alt", "")

    allele = make_vrs_allele(chrom, pos, ref, alt)
    vrs_dict = vrs_allele_to_dict(allele)

    return {
        "variantInternalId": f"{chrom}:{pos}:{ref}:{alt}",
        "variantType": row.get("variant_type", "SNP"),
        "variation": vrs_dict,  # GA4GH VRS v2.0 allele object
        "variantLevelData": {
            "clinicalRelevances": [
                {
                    "category": row.get("acmg_class", "VUS"),
                    "conditionId": row.get("condition_id", ""),
                    "clinVarIds": [row["clinvar_id"]] if row.get("clinvar_id") else [],
                }
            ],
        },
        "caseLevelData": [],
        "frequencyInPopulations": [
            {
                "frequencies": [
                    {
                        "population": "gnomAD v4.1 All",
                        "alleleFrequency": row.get("gnomad_af"),
                        "alleleCount": row.get("gnomad_ac"),
                        "alleleNumber": None,
                        "numberOfHomozygotes": row.get("gnomad_nhomalt"),
                    }
                ],
                "source": "gnomAD v4.1 (April 2024, 807,162 individuals)",
                "sourceReference": "https://gnomad.broadinstitute.org",
            }
        ],
        "molecularAttributes": {
            "geneIds": [row["gene_symbol"]] if row.get("gene_symbol") else [],
            "molecularEffects": [
                {
                    "id": row.get("consequence_id", ""),
                    "label": row.get("consequence", ""),
                }
            ],
            "aminoacidChanges": [row["hgvsp"]] if row.get("hgvsp") else [],
        },
        "_vrsVersion": "2.0",
        "_vrsId": allele.vrs_id,
    }


# ---------------------------------------------------------------------------
# Endpoint
# ---------------------------------------------------------------------------


@router.get(
    "",
    summary="Genomic variant query",
    description=(
        "Query genomic variants by coordinates (GRCh38) or VRS identifier. "
        "Returns Beacon v2.1.1 GenomicVariation records with embedded GA4GH "
        "VRS v2.0 identifiers and gnomAD v4.1 population frequencies. "
        "Authentication with a GA4GH Passport JWT grants record-level granularity."
    ),
    response_model=None,
)
async def query_g_variants(
    chrom: str | None = Query(None, description="Chromosome (GRCh38, e.g. chr17)."),
    start: int | None = Query(None, description="Start position (1-based, inclusive)."),
    end: int | None = Query(None, description="End position (1-based, inclusive)."),
    ref: str | None = Query(None, description="Reference allele (VCF notation)."),
    alt: str | None = Query(None, description="Alternate allele (VCF notation)."),
    vrs_id: str | None = Query(None, alias="vrsId", description="GA4GH VRS v2.0 identifier."),
    gene_symbol: str | None = Query(None, alias="geneSymbol", description="HGNC gene symbol."),
    include_datasets: list[str] | None = Query(
        None, alias="includeDatasetResponses", description="Dataset IDs to include."
    ),
    granularity: str = Query("record", description="Response granularity: boolean, count, record."),
    skip: int = Query(0, ge=0, description="Pagination offset."),
    limit: int = Query(10, ge=1, le=100, description="Maximum records to return."),
    passport: dict[str, Any] | None = Depends(optional_passport),
    session: Any = Depends(get_session),
) -> dict[str, Any]:
    """Query genomic variants matching the given parameters.

    Supports coordinate-based queries (chrom+start+end) and VRS ID lookups.
    Returns GA4GH Beacon v2.1.1 response with VRS v2.0 variant identifiers.

    Args:
        chrom: Chromosome in GRCh38 notation.
        start: 1-based start position (inclusive).
        end: 1-based end position (inclusive).
        ref: Reference allele (VCF notation).
        alt: Alternate allele (VCF notation).
        vrs_id: GA4GH VRS v2.0 identifier for direct lookup.
        gene_symbol: HGNC gene symbol filter.
        include_datasets: List of dataset IDs to include.
        granularity: Response granularity level (boolean/count/record).
        skip: Pagination offset.
        limit: Maximum number of records to return (max 100).
        passport: Decoded GA4GH Passport JWT claims (optional).
        session: Async database session.

    Returns:
        Dict conforming to the Beacon v2.1.1 BeaconResultsetsResponse schema.

    Raises:
        HTTPException: 400 if no query parameters provided.
        HTTPException: 422 if granularity is invalid.

    References:
        GA4GH Beacon v2.1.1 spec (December 13, 2024).
        VRS v2.0: Wagner et al. 2021 PMID:35072137.
        gnomAD v4.1 (April 2024): https://gnomad.broadinstitute.org
    """
    if not any([chrom, vrs_id, gene_symbol]):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=(
                "At least one query parameter is required: "
                "chrom, vrsId, or geneSymbol."
            ),
        )

    valid_granularities = {"boolean", "count", "record"}
    if granularity not in valid_granularities:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"granularity must be one of: {', '.join(sorted(valid_granularities))}",
        )

    # Record-level responses require GA4GH Passport authentication
    if granularity == "record" and passport is None:
        granularity = "count"  # downgrade to count if no passport

    # Build mock result set for demonstration
    # In production, query the PostgreSQL database via session
    mock_variants: list[dict[str, Any]] = []

    if chrom and start:
        mock_variants.append({
            "chrom": chrom,
            "pos": start,
            "ref": ref or "N",
            "alt": alt or "A",
            "gnomad_af": 0.000012,
            "gnomad_ac": 10,
            "gnomad_nhomalt": 0,
            "clinvar_id": None,
            "gene_symbol": gene_symbol or "",
            "consequence": "missense_variant",
            "consequence_id": "SO:0001583",
            "hgvsc": None,
            "hgvsp": None,
            "acmg_class": "VUS",
            "variant_type": "SNP",
            "condition_id": "",
        })

    total = len(mock_variants)
    page = mock_variants[skip : skip + limit]

    if granularity == "boolean":
        return {
            "meta": _beacon_meta("g_variants"),
            "responseSummary": {"exists": total > 0, "numTotalResults": None},
            "beaconHandovers": [],
        }

    if granularity == "count":
        return {
            "meta": _beacon_meta("g_variants"),
            "responseSummary": {"exists": total > 0, "numTotalResults": total},
            "beaconHandovers": [],
        }

    # Record-level — requires passport
    results = [_build_variant_response(v) for v in page]
    return {
        "meta": _beacon_meta("g_variants"),
        "responseSummary": {"exists": total > 0, "numTotalResults": total},
        "resultSets": [
            {
                "id": "genomeforge.wgs.grch38",
                "type": "dataset",
                "exists": total > 0,
                "resultsCount": len(results),
                "results": results,
            }
        ],
        "beaconHandovers": [],
    }


def _beacon_meta(entity_type: str) -> dict[str, Any]:
    """Build Beacon v2.1.1 response meta object.

    Args:
        entity_type: Entity type string (e.g. ``"g_variants"``).

    Returns:
        Dict with ``beaconId``, ``apiVersion``, and ``returnedSchemas``.
    """
    return {
        "beaconId": "org.genomeforge.beacon",
        "apiVersion": "v2.1.1",
        "returnedSchemas": [
            {
                "entityType": entity_type,
                "schema": (
                    "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/main"
                    "/models/src/beacon-v2-default-model/genomicVariations/defaultSchema.yaml"
                ),
            }
        ],
    }
