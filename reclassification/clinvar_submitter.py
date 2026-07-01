"""NHS-mandated ClinVar submission client for ClaritySeq.

This module implements the NCBI ClinVar API submission pipeline as required
by the NHS Genomic Medicine Service participation agreement.

Mandate reference:
    ACGS 2024 Introduction: "All NHS WGS laboratories are required to submit
    clinically interpreted variants to ClinVar as a condition of participation
    in the NHS Genomic Medicine Service. Submission of pathogenic (P) and
    likely pathogenic (LP) variants must occur within 3 months of the date
    of the clinical report. Variants of uncertain significance (VUS) must be
    submitted within 6 months."

NCBI ClinVar Submission API:
    - REST endpoint: https://submit.ncbi.nlm.nih.gov/api/2.0/files/
    - API documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/
    - Authentication: NCBI API key in X-API-KEY request header
    - Submission format: JSON (preferred) or XML; this module uses JSON
    - MANE Select HGVSc: Required by NCBI for unambiguous transcript
      identification (Morales et al. 2022, Nat Methods PMID:35379937).
      Must use RefSeq NM_ accession matching MANE Select v1.3 or later.

Rate limits:
    - NCBI API: 3 requests/second without API key; 10 requests/second with.
    - Batch submissions up to 10,000 variants per submission file.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Optional
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import Element, SubElement

import requests

from reclassification.models import ClinVarSubmissionQueue, SubmissionStatus

logger = logging.getLogger(__name__)

# NCBI ClinVar submission API endpoints
NCBI_SUBMISSION_API_BASE = "https://submit.ncbi.nlm.nih.gov/api/2.0"
NCBI_SUBMISSION_ENDPOINT = f"{NCBI_SUBMISSION_API_BASE}/files/"
NCBI_STATUS_ENDPOINT = f"{NCBI_SUBMISSION_API_BASE}/submissions/{{submission_id}}/actions/"

# NCBI rate limiting: 10 req/s with API key, 3 req/s without
NCBI_REQUEST_INTERVAL_WITH_KEY = 0.12   # seconds between requests
NCBI_REQUEST_INTERVAL_WITHOUT_KEY = 0.35

# ClinVar submission XML namespace
CLINVAR_XML_NS = "http://www.ncbi.nlm.nih.gov/clinvar"

# ClaritySeq organisation identifiers for ClinVar
CLARITYSEQ_ORG_ID = "507000"  # NCBI-assigned organisation ID (placeholder)
CLARITYSEQ_ORG_NAME = "ClaritySeq NHS WGS Laboratory"

# Mapping from internal significance to ClinVar controlled vocabulary
SIGNIFICANCE_TO_CLINVAR: dict[str, str] = {
    "Pathogenic": "Pathogenic",
    "Likely pathogenic": "Likely pathogenic",
    "Uncertain significance": "Uncertain significance",
    "Likely benign": "Likely benign",
    "Benign": "Benign",
    "Conflicting interpretations": "Conflicting interpretations of pathogenicity",
    "risk factor": "risk factor",
}

# Mapping from internal significance to ClinVar assertion type
ASSERTION_TYPE_MAP: dict[str, str] = {
    "Pathogenic": "variation to disease",
    "Likely pathogenic": "variation to disease",
    "Uncertain significance": "variation to disease",
    "Likely benign": "variation to disease",
    "Benign": "variation to disease",
}


@dataclass
class SubmissionResult:
    """Result of a ClinVar API submission attempt.

    Attributes:
        success: True if submission was accepted by NCBI.
        submission_id: NCBI-assigned submission batch ID (on success).
        status: Current SubmissionStatus enum value.
        ncbi_response_raw: Raw JSON response from NCBI API.
        error_message: Human-readable error description (on failure).
        submitted_at: Timestamp of the API call.
    """

    success: bool
    submission_id: Optional[str]
    status: SubmissionStatus
    ncbi_response_raw: Optional[dict[str, Any]]
    error_message: Optional[str]
    submitted_at: datetime


def _get_request_session(api_key: Optional[str] = None) -> requests.Session:
    """Create a configured requests.Session for NCBI API calls.

    Sets standard headers including optional X-API-KEY for increased
    rate limits.

    Args:
        api_key: Optional NCBI API key. If provided, enables 10 req/s rate
            limit instead of the default 3 req/s.

    Returns:
        Configured requests.Session object.
    """
    session = requests.Session()
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "ClaritySeq-ClinVarSubmitter/1.0 (NHS WGS; contact@clarityseq.nhs.uk)",
    }
    if api_key:
        headers["X-API-KEY"] = api_key
    session.headers.update(headers)
    return session


def build_submission_json(submission: ClinVarSubmissionQueue) -> dict[str, Any]:
    """Build the NCBI ClinVar API JSON submission payload.

    Constructs the JSON structure required by the NCBI ClinVar Submission
    API v2.0. Uses MANE Select HGVSc for transcript identification when
    available.

    Args:
        submission: ClinVarSubmissionQueue ORM object containing all
            variant and classification details.

    Returns:
        JSON-serialisable dictionary conforming to the NCBI ClinVar
        submission schema. Contains:
        - clinvarSubmission[]: array with one variant submission
        - submissionName: unique identifier for this batch
        - assertionCriteria: reference to ACMG/AMP 2015 guidelines

    Raises:
        ValueError: If required fields (gene_symbol, clinical_significance,
            condition_name) are missing from the submission.

    Example:
        >>> payload = build_submission_json(queue_item)
        >>> json.dumps(payload, indent=2)
    """
    if not submission.gene_symbol:
        raise ValueError(f"Missing gene_symbol for submission {submission.id}")
    if not submission.clinical_significance:
        raise ValueError(
            f"Missing clinical_significance for submission {submission.id}"
        )
    if not submission.condition_name:
        raise ValueError(f"Missing condition_name for submission {submission.id}")

    clinvar_significance = SIGNIFICANCE_TO_CLINVAR.get(
        submission.clinical_significance,
        submission.clinical_significance,
    )

    # Build variant specification — prefer MANE Select HGVSc if available
    variant_spec: dict[str, Any]
    if submission.mane_select_hgvsc:
        variant_spec = {
            "hgvs": submission.mane_select_hgvsc,
            # MANE Select transcript provides unambiguous gene context
        }
    else:
        # Fall back to chromosomal coordinates
        variant_spec = {
            "chromosomeCoordinates": {
                "assembly": "GRCh38",
                "chromosome": submission.chromosome,
                "start": submission.position_grch38,
                "referenceAllele": submission.ref_allele,
                "alternateAllele": submission.alt_allele,
            }
        }

    # Build condition specification
    condition_spec: dict[str, Any] = {
        "db": "MedGen" if submission.condition_id else "General",
        "name": submission.condition_name,
    }
    if submission.condition_id:
        condition_spec["id"] = submission.condition_id

    # Build evidence/citation (BayesACMG posterior probability as comment)
    observation: dict[str, Any] = {
        "affectedStatus": "yes",
        "alleleOrigin": "germline",
        "collectionMethod": "clinical testing",
    }
    if submission.bayesacmg_probability is not None:
        observation["comment"] = (
            f"BayesACMG posterior probability of pathogenicity: "
            f"{submission.bayesacmg_probability:.4f} "
            f"(Tavtigian et al. 2020, PMID:31479589). "
            f"ACMG criteria applied: {submission.evidence_codes or 'not specified'}."
        )

    payload: dict[str, Any] = {
        "actions": [
            {
                "type": "AddData",
                "targetDb": "clinvar",
                "data": {
                    "content": {
                        "clinvarSubmission": [
                            {
                                "clinicalSignificance": {
                                    "clinicalSignificanceDescription": clinvar_significance,
                                    "comment": (
                                        f"Classified per ACMG/AMP 2015 criteria "
                                        f"(Richards et al. 2015, PMID:25741868) "
                                        f"and ACGS 2024 Best Practice Guidelines. "
                                        f"Evidence codes: "
                                        f"{submission.evidence_codes or 'not specified'}."
                                    ),
                                    "dateLastEvaluated": datetime.now(timezone.utc)
                                    .date()
                                    .isoformat(),
                                },
                                "conditionSet": {
                                    "condition": [condition_spec],
                                },
                                "observedIn": [observation],
                                "variantSet": {
                                    "variant": [
                                        {
                                            **variant_spec,
                                            "gene": [
                                                {
                                                    "symbol": submission.gene_symbol
                                                }
                                            ],
                                        }
                                    ]
                                },
                                "assertionCriteria": {
                                    "db": "PubMed",
                                    "id": "25741868",  # Richards et al. 2015
                                },
                                "localID": submission.variant_id,
                                "localKey": f"clarityseq-{submission.id}",
                            }
                        ]
                    }
                },
            }
        ]
    }

    return payload


def build_submission_xml(submission: ClinVarSubmissionQueue) -> str:
    """Generate ClinVar XML submission format (legacy API).

    Produces a ClinVarSubmission XML document compatible with the NCBI
    ClinVar legacy XML upload format. Prefer build_submission_json() for
    new submissions (JSON API v2.0 is the current standard), but this
    XML format is retained for compatibility with laboratory LIMS systems
    that use the older format.

    Args:
        submission: ClinVarSubmissionQueue ORM object.

    Returns:
        UTF-8 encoded XML string conforming to the ClinVar submission
        XML schema (ClinVarSubmission.xsd v1.7).

    Raises:
        ValueError: If required submission fields are missing.

    Example:
        >>> xml_str = build_submission_xml(queue_item)
        >>> print(xml_str[:200])
        <?xml version='1.0' encoding='utf-8'?>
        <ClinVarSubmission ...>
    """
    if not submission.gene_symbol or not submission.clinical_significance:
        raise ValueError(
            f"Missing required fields for submission {submission.id}"
        )

    clinvar_significance = SIGNIFICANCE_TO_CLINVAR.get(
        submission.clinical_significance, submission.clinical_significance
    )

    # Root element
    root = Element("ClinVarSubmission")
    root.set("xmlns", CLINVAR_XML_NS)
    root.set("SubmissionDate", datetime.now(timezone.utc).date().isoformat())

    # Submission header
    header = SubElement(root, "SubmissionHeader")
    SubElement(header, "OrganizationID").text = CLARITYSEQ_ORG_ID
    SubElement(header, "OrganizationName").text = CLARITYSEQ_ORG_NAME
    SubElement(header, "SubmissionName").text = f"clarityseq-{submission.id}"

    # Clinical assertion
    assertion = SubElement(root, "ClinicalAssertion")
    assertion.set("LocalID", submission.variant_id)

    # Variant description
    variant_elem = SubElement(assertion, "SimpleAllele")
    gene_list = SubElement(variant_elem, "GeneList")
    gene_elem = SubElement(gene_list, "Gene")
    SubElement(gene_elem, "Symbol").text = submission.gene_symbol

    # Use MANE Select HGVSc if available
    if submission.mane_select_hgvsc:
        attribute_list = SubElement(variant_elem, "AttributeList")
        attr = SubElement(attribute_list, "Attribute")
        attr.set("Type", "HGVS, coding, RefSeq")
        attr.text = submission.mane_select_hgvsc
    else:
        location = SubElement(variant_elem, "Location")
        seq_loc = SubElement(location, "SequenceLocation")
        seq_loc.set("Assembly", "GRCh38")
        seq_loc.set("Chr", submission.chromosome)
        seq_loc.set("start", str(submission.position_grch38))
        seq_loc.set("referenceAllele", submission.ref_allele)
        seq_loc.set("alternateAllele", submission.alt_allele)

    # Clinical significance
    clinical_sig = SubElement(assertion, "ClinicalSignificance")
    sig_desc = SubElement(clinical_sig, "Description")
    sig_desc.text = clinvar_significance
    date_evaluated = SubElement(clinical_sig, "DateLastEvaluated")
    date_evaluated.text = datetime.now(timezone.utc).date().isoformat()

    comment = SubElement(clinical_sig, "Comment")
    comment.text = (
        f"Classified per ACMG/AMP 2015 criteria. "
        f"BayesACMG probability: "
        f"{submission.bayesacmg_probability:.4f}"
        if submission.bayesacmg_probability is not None
        else "Classified per ACMG/AMP 2015 criteria."
    )

    # Condition
    trait_set = SubElement(assertion, "TraitSet")
    trait_set.set("Type", "Disease")
    trait = SubElement(trait_set, "Trait")
    trait.set("Type", "Disease")
    trait_name = SubElement(trait, "Name")
    value_elem = SubElement(trait_name, "ElementValue")
    value_elem.set("Type", "Preferred")
    value_elem.text = submission.condition_name

    if submission.condition_id:
        xref = SubElement(trait, "XRef")
        xref.set("DB", "MedGen")
        xref.set("ID", submission.condition_id)

    # Observed in
    observed = SubElement(assertion, "ObservedInList")
    observed_in = SubElement(observed, "ObservedIn")
    sample = SubElement(observed_in, "Sample")
    SubElement(sample, "Origin").text = "germline"
    SubElement(sample, "Tissue").text = "blood"
    SubElement(sample, "AffectedStatus").text = "yes"
    SubElement(observed_in, "Method").text = "clinical testing"

    # Assertion criteria reference (ACMG/AMP 2015)
    criteria = SubElement(assertion, "AssertionCriteria")
    SubElement(criteria, "Database").text = "PubMed"
    SubElement(criteria, "ID").text = "25741868"

    # Serialise to string with XML declaration
    ET.indent(root, space="  ")
    return ET.tostring(root, encoding="unicode", xml_declaration=True)


def submit_variant(
    submission: ClinVarSubmissionQueue,
    api_key: Optional[str] = None,
    dry_run: bool = False,
) -> SubmissionResult:
    """Submit a variant to the NCBI ClinVar API.

    Submits a single variant classification to the NCBI ClinVar Submission
    API v2.0. Handles authentication, rate limiting, and response parsing.
    On success, populates submission.ncbi_submission_id for status tracking.

    Args:
        submission: ClinVarSubmissionQueue ORM object with complete variant
            and classification information. Must have mane_select_hgvsc or
            chromosome/position/ref/alt fields populated.
        api_key: NCBI API key for increased rate limits. Should be stored
            in environment variable NCBI_API_KEY and passed here. Allows
            10 req/s vs 3 req/s for anonymous access.
        dry_run: If True, build the submission payload and log it but do
            not make the actual API call. Useful for testing.

    Returns:
        SubmissionResult with success status, NCBI submission ID (on success),
        and the raw API response.

    Raises:
        requests.RequestException: On network-level errors.
        ValueError: If submission payload fails validation.

    Example:
        >>> result = submit_variant(queue_item, api_key=os.environ["NCBI_API_KEY"])
        >>> if result.success:
        ...     print(f"Submitted as {result.submission_id}")
    """
    now = datetime.now(timezone.utc)

    # Build submission payload
    try:
        payload = build_submission_json(submission)
    except ValueError as exc:
        logger.error("Submission payload error for %s: %s", submission.id, exc)
        return SubmissionResult(
            success=False,
            submission_id=None,
            status=SubmissionStatus.ERROR,
            ncbi_response_raw=None,
            error_message=str(exc),
            submitted_at=now,
        )

    if dry_run:
        logger.info(
            "DRY RUN: Would submit variant %s to ClinVar:\n%s",
            submission.variant_id,
            json.dumps(payload, indent=2),
        )
        return SubmissionResult(
            success=True,
            submission_id="dry-run-id",
            status=SubmissionStatus.SUBMITTED,
            ncbi_response_raw={"dryRun": True},
            error_message=None,
            submitted_at=now,
        )

    session = _get_request_session(api_key)

    # Respect NCBI rate limits
    interval = (
        NCBI_REQUEST_INTERVAL_WITH_KEY
        if api_key
        else NCBI_REQUEST_INTERVAL_WITHOUT_KEY
    )
    time.sleep(interval)

    try:
        response = session.post(
            NCBI_SUBMISSION_ENDPOINT,
            json=payload,
            timeout=30,
        )
        response.raise_for_status()
    except requests.HTTPError as exc:
        error_body = ""
        try:
            error_body = response.json()
        except Exception:
            error_body = response.text

        logger.error(
            "NCBI ClinVar HTTP error for submission %s: %s %s",
            submission.id, response.status_code, error_body,
        )
        return SubmissionResult(
            success=False,
            submission_id=None,
            status=SubmissionStatus.ERROR,
            ncbi_response_raw=error_body if isinstance(error_body, dict) else None,
            error_message=(
                f"HTTP {response.status_code}: {exc}"
            ),
            submitted_at=now,
        )
    except requests.RequestException as exc:
        logger.error(
            "Network error submitting to ClinVar: %s", exc
        )
        return SubmissionResult(
            success=False,
            submission_id=None,
            status=SubmissionStatus.ERROR,
            ncbi_response_raw=None,
            error_message=str(exc),
            submitted_at=now,
        )

    response_data = response.json()
    submission_id = response_data.get("id")

    logger.info(
        "ClinVar submission accepted for variant %s: submission_id=%s",
        submission.variant_id, submission_id,
    )

    return SubmissionResult(
        success=True,
        submission_id=submission_id,
        status=SubmissionStatus.SUBMITTED,
        ncbi_response_raw=response_data,
        error_message=None,
        submitted_at=now,
    )


def check_submission_status(
    submission_id: str,
    api_key: Optional[str] = None,
) -> str:
    """Poll NCBI ClinVar API for the status of a previously submitted batch.

    Queries the NCBI ClinVar Submission API to determine the current
    processing status of a submission batch. NCBI typically processes
    submissions within 2–5 business days.

    Args:
        submission_id: NCBI-assigned submission batch ID returned by
            submit_variant() or the NCBI web interface.
        api_key: Optional NCBI API key for authentication.

    Returns:
        One of the SubmissionStatus values as a string:
        - 'processing': NCBI is still processing the submission
        - 'accepted': Submission accepted and accession numbers assigned
        - 'rejected': Submission rejected (check ncbi_response for details)
        - 'error': Could not determine status

    Raises:
        requests.RequestException: On network-level errors.

    Example:
        >>> status = check_submission_status("SUB123456")
        >>> print(status)
        'accepted'
    """
    url = NCBI_STATUS_ENDPOINT.format(submission_id=submission_id)
    session = _get_request_session(api_key)

    interval = (
        NCBI_REQUEST_INTERVAL_WITH_KEY
        if api_key
        else NCBI_REQUEST_INTERVAL_WITHOUT_KEY
    )
    time.sleep(interval)

    try:
        response = session.get(url, timeout=30)
        response.raise_for_status()
    except requests.HTTPError as exc:
        logger.error(
            "Error checking ClinVar submission status %s: %s",
            submission_id, exc,
        )
        return SubmissionStatus.ERROR.value
    except requests.RequestException as exc:
        logger.error("Network error checking submission %s: %s", submission_id, exc)
        return SubmissionStatus.ERROR.value

    data = response.json()

    # Parse NCBI response structure
    actions = data.get("actions", [])
    if not actions:
        logger.warning(
            "No actions in status response for submission %s", submission_id
        )
        return SubmissionStatus.PROCESSING.value

    action = actions[0]
    status_str = action.get("status", "").lower()

    status_map = {
        "submitted": SubmissionStatus.SUBMITTED.value,
        "processing": SubmissionStatus.PROCESSING.value,
        "processed": SubmissionStatus.ACCEPTED.value,
        "error": SubmissionStatus.REJECTED.value,
        "failed": SubmissionStatus.ERROR.value,
    }

    mapped_status = status_map.get(status_str, SubmissionStatus.PROCESSING.value)

    logger.info(
        "Submission %s status: %s -> %s",
        submission_id, status_str, mapped_status,
    )

    return mapped_status
