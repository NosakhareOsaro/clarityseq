"""
prioritisation.exomiser_client
================================
Exomiser 14 REST API client for variant prioritisation.

Exomiser 14 accepts Phenopackets v2.0 input and returns ranked candidate
variants with HPO similarity scores, variant pathogenicity scores, and
inheritance mode compatibility.

Exomiser 14 (July 2023):
    - Phenopackets v2.0 native input support.
    - Updated human phenotype network (HPN) from HPO 2023 release.
    - OMIM/Orphanet disease phenotype associations (monthly updates).
    - ClinVar, gnomAD, and CADD integration.
    - HIPHY algorithm for HPO phenotypic similarity.

REST API (Exomiser 14 Spring Boot server):
    POST /api/v1/analysis   — submit Phenopacket + analysis settings.
    GET  /api/v1/results/{id} — poll for analysis results.

References:
    Robinson et al. 2023 Nature Genetics PMID:37604970 (Exomiser).
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    Exomiser docs: https://exomiser.readthedocs.io/
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass, field
from typing import Any

import httpx
from tenacity import retry, stop_after_attempt, wait_exponential

logger = logging.getLogger(__name__)

_EXOMISER_BASE_URL: str = os.getenv(
    "EXOMISER_API_URL",
    "http://localhost:8080/exomiser",
)
_EXOMISER_API_KEY: str | None = os.getenv("EXOMISER_API_KEY")
_POLL_INTERVAL_SECONDS: float = 5.0
_MAX_POLL_ATTEMPTS: int = 120  # max 10 minutes polling


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class ExomiserResult:
    """A ranked candidate variant from Exomiser 14.

    Attributes:
        rank: Overall prioritisation rank (1 = highest priority).
        gene_symbol: HGNC gene symbol.
        combined_score: Exomiser combined score (variant + phenotype).
        phenotype_score: Phenotypic similarity score (HPO BMA Resnik).
        variant_score: Variant pathogenicity score.
        acmg_class: ACMG classification from Exomiser.
        contributing_variants: List of dicts with variant details.
        disease_associations: List of dicts with OMIM/Orphanet disease associations.
        inheritance_mode: Inheritance mode from Exomiser analysis.
    """

    rank: int
    gene_symbol: str
    combined_score: float
    phenotype_score: float
    variant_score: float
    acmg_class: str = "VUS"
    contributing_variants: list[dict[str, Any]] = field(default_factory=list)
    disease_associations: list[dict[str, Any]] = field(default_factory=list)
    inheritance_mode: str = ""


def _build_analysis_request(
    phenopacket: dict[str, Any],
    vcf_path: str | None = None,
    inheritance_modes: list[str] | None = None,
    min_frequency: float = 1.0,
    gene_panel: list[str] | None = None,
) -> dict[str, Any]:
    """Build the Exomiser 14 REST API analysis request body.

    Args:
        phenopacket: Phenopackets v2.0 dict (Exomiser 14 native input).
        vcf_path: Optional path to a VCF file for variant input.
            If None, variants are read from the phenopacket.
        inheritance_modes: List of inheritance modes to test
            (e.g. ``["AUTOSOMAL_DOMINANT", "AUTOSOMAL_RECESSIVE"]``).
            Default: all modes.
        min_frequency: Maximum allele frequency filter (default 1.0%).
        gene_panel: Optional list of gene symbols to restrict analysis.

    Returns:
        Dict conforming to the Exomiser 14 analysis settings API.
    """
    modes = inheritance_modes or [
        "AUTOSOMAL_DOMINANT",
        "AUTOSOMAL_RECESSIVE",
        "X_DOMINANT",
        "X_RECESSIVE",
        "MITOCHONDRIAL",
    ]

    analysis: dict[str, Any] = {
        "phenopacket": phenopacket,
        "analysisMode": "PASS_ONLY",
        "inheritanceModes": {mode: 0.0 for mode in modes},
        "frequencySources": [
            "THOUSAND_GENOMES",
            "TOPMED",
            "UK10K",
            "ESP_EUROPEAN_AMERICAN",
            "ESP_AFRICAN_AMERICAN",
            "EXAC_NON_FINNISH_EUROPEAN",
            "GNOMAD_E_NFE",
            "GNOMAD_G_NFE",
        ],
        "pathogenicitySources": ["CLINVAR", "CADD", "REVEL"],
        "steps": [
            {
                "failedVariantFilter": {},
            },
            {
                "variantEffectFilter": {
                    "remove": [
                        "FIVE_PRIME_UTR_EXON_VARIANT",
                        "THREE_PRIME_UTR_EXON_VARIANT",
                        "SYNONYMOUS_VARIANT",
                        "INTERGENIC_VARIANT",
                        "NON_CODING_TRANSCRIPT_INTRON_VARIANT",
                    ]
                },
            },
            {
                "frequencyFilter": {
                    "maxFrequency": min_frequency,
                },
            },
            {
                "pathogenicityFilter": {
                    "keepNonPathogenic": True,
                },
            },
            {
                "inheritanceFilter": {},
            },
            {
                "omimPrioritiser": {},
            },
            {
                "hiPhivePrioritiser": {
                    "runParams": "human",
                },
            },
        ],
    }

    if vcf_path:
        analysis["vcf"] = vcf_path

    if gene_panel:
        analysis["steps"].insert(
            2,
            {
                "genePanel": {
                    "geneSymbols": gene_panel,
                },
            },
        )

    return analysis


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    reraise=True,
)
def _submit_analysis(analysis_request: dict[str, Any]) -> str:
    """Submit an Exomiser 14 analysis job via REST API.

    Args:
        analysis_request: Analysis settings dict from _build_analysis_request().

    Returns:
        Analysis job ID for polling.

    Raises:
        httpx.HTTPStatusError: On non-2xx response.
        httpx.RequestError: On network error.
    """
    url = f"{_EXOMISER_BASE_URL}/api/v1/analysis"
    headers = {"Content-Type": "application/json"}
    if _EXOMISER_API_KEY:
        headers["Authorization"] = f"Bearer {_EXOMISER_API_KEY}"

    response = httpx.post(url, json=analysis_request, headers=headers, timeout=30.0)
    response.raise_for_status()
    data = response.json()
    job_id = data.get("id", data.get("jobId", ""))
    if not job_id:
        raise RuntimeError(f"Exomiser API returned no job ID: {data}")
    return str(job_id)


def _poll_results(job_id: str) -> dict[str, Any]:
    """Poll for Exomiser analysis results until completion.

    Args:
        job_id: Analysis job ID from _submit_analysis().

    Returns:
        Exomiser results dict with ranked gene list.

    Raises:
        RuntimeError: If analysis fails or times out.
    """
    url = f"{_EXOMISER_BASE_URL}/api/v1/results/{job_id}"
    headers: dict[str, str] = {}
    if _EXOMISER_API_KEY:
        headers["Authorization"] = f"Bearer {_EXOMISER_API_KEY}"

    for attempt in range(_MAX_POLL_ATTEMPTS):
        response = httpx.get(url, headers=headers, timeout=15.0)
        response.raise_for_status()
        data = response.json()

        status = data.get("status", "")
        if status == "COMPLETED":
            return data
        elif status in ("FAILED", "ERROR"):
            raise RuntimeError(f"Exomiser analysis {job_id} failed: {data.get('message', '')}")

        logger.debug(
            "Exomiser job %s status=%s (attempt %d/%d)",
            job_id, status, attempt + 1, _MAX_POLL_ATTEMPTS,
        )
        time.sleep(_POLL_INTERVAL_SECONDS)

    raise RuntimeError(
        f"Exomiser analysis {job_id} timed out after "
        f"{_MAX_POLL_ATTEMPTS * _POLL_INTERVAL_SECONDS}s."
    )


def _parse_exomiser_results(data: dict[str, Any]) -> list[ExomiserResult]:
    """Parse Exomiser 14 results into ExomiserResult objects.

    Args:
        data: Completed analysis result dict from Exomiser API.

    Returns:
        List of ExomiserResult sorted by rank.
    """
    results: list[ExomiserResult] = []
    genes = data.get("results", data.get("genes", []))

    for i, gene in enumerate(genes):
        results.append(
            ExomiserResult(
                rank=gene.get("rank", i + 1),
                gene_symbol=gene.get("geneSymbol", gene.get("geneName", "")),
                combined_score=float(gene.get("combinedScore", 0.0)),
                phenotype_score=float(gene.get("phenotypeScore", 0.0)),
                variant_score=float(gene.get("variantScore", 0.0)),
                acmg_class=gene.get("acmgClassification", "VUS"),
                contributing_variants=gene.get("contributingVariants", []),
                disease_associations=gene.get("omimResults", []),
                inheritance_mode=gene.get("modeOfInheritance", ""),
            )
        )

    results.sort(key=lambda r: r.rank)
    return results


def run_exomiser(
    phenopacket: dict[str, Any],
    vcf_path: str | None = None,
    inheritance_modes: list[str] | None = None,
    gene_panel: list[str] | None = None,
    min_frequency: float = 1.0,
) -> list[ExomiserResult]:
    """Run Exomiser 14 variant prioritisation for a patient phenopacket.

    Submits the phenopacket to the Exomiser 14 REST API, polls for results,
    and returns ranked candidate genes/variants.

    Args:
        phenopacket: Phenopackets v2.0 dict (native Exomiser 14 input).
            Must include phenotypicFeatures with HPO terms.
        vcf_path: Path to proband VCF (optional; variants from phenopacket
            if not provided).
        inheritance_modes: Inheritance modes to analyse.
            Default: all (AD, AR, XL, Mito).
        gene_panel: Optional gene panel to restrict analysis
            (e.g. NHS GMS panel gene list).
        min_frequency: Maximum allele frequency for rare variant filter
            (default 1.0% = rare variant threshold).

    Returns:
        List of ExomiserResult objects ranked by combined score.

    Raises:
        httpx.HTTPStatusError: On Exomiser API errors.
        RuntimeError: On analysis failure or timeout.

    References:
        Robinson et al. 2023 Nature Genetics PMID:37604970.
        Exomiser 14 REST API: https://exomiser.readthedocs.io/en/latest/api.html
    """
    analysis_request = _build_analysis_request(
        phenopacket=phenopacket,
        vcf_path=vcf_path,
        inheritance_modes=inheritance_modes,
        min_frequency=min_frequency,
        gene_panel=gene_panel,
    )

    logger.info(
        "Submitting Exomiser 14 analysis for phenopacket %s",
        phenopacket.get("id", "unknown"),
    )

    job_id = _submit_analysis(analysis_request)
    logger.info("Exomiser job submitted: %s", job_id)

    raw_results = _poll_results(job_id)
    results = _parse_exomiser_results(raw_results)

    logger.info(
        "Exomiser analysis %s complete: %d candidate genes ranked.",
        job_id,
        len(results),
    )
    return results
