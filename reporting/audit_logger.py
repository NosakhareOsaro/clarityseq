"""
reporting.audit_logger
========================
JSON-LD audit trail writer for GenomeForge clinical reports.

JSON-LD provides a linked data representation of the audit record,
suitable for FHIR Provenance resources and NHS GMS audit requirements.

Audit trail includes:
    - Report generation date and operator.
    - Pipeline version and component versions.
    - Data sources: gnomAD v4.1, VEP 111, AlphaMissense, ClinVar.
    - Classification rules applied (ACMG/AMP rule IDs and strengths).
    - PM2 evidence weight: Supporting (1 pt) per ClinGen SVI 2024.
    - VUS review dates (ACGS 2024 §9).
    - Novel P/LP pending ClinVar submission flags.
    - Bayesian model posterior probabilities and 95% HDI bounds.

JSON-LD context:
    Uses schema.org and prov-o terms for interoperability.
    FHIR Provenance mapping available via fhir_mapper.py.

References:
    JSON-LD 1.1: https://www.w3.org/TR/json-ld11/
    PROV-O: https://www.w3.org/TR/prov-o/
    NHS GMS audit requirements.
    ACGS 2024 v1.2 §5 — classification audit trail.
"""

from __future__ import annotations

import json
import logging
import os
import socket
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# JSON-LD context
# ---------------------------------------------------------------------------

_JSON_LD_CONTEXT: dict[str, Any] = {
    "@vocab": "https://schema.org/",
    "prov": "http://www.w3.org/ns/prov#",
    "xsd": "http://www.w3.org/2001/XMLSchema#",
    "ga4gh": "https://ga4gh.org/terms/",
    "genomeforge": "https://genomeforge.github.io/terms/",
    "report_date": {"@type": "xsd:date"},
    "generated_at": {"@type": "xsd:dateTime"},
    "pipeline_version": "genomeforge:pipelineVersion",
    "classification_scheme": "genomeforge:classificationScheme",
    "pm2_weight": "genomeforge:pm2Weight",
    "vus_review_date": "genomeforge:vusReviewDate",
    "pending_clinvar_submission": "genomeforge:pendingClinVarSubmission",
    "posterior_p": "genomeforge:bayesianPosteriorP",
}


def write_audit_log(
    audit_data: dict[str, Any],
    output_path: Path,
    pretty: bool = True,
) -> None:
    """Write a JSON-LD audit trail file for a clinical report.

    Args:
        audit_data: Audit data dict from ReportGenerator._write_json_ld_audit().
            Expected keys: ``report_date``, ``patient_id``, ``sample_id``,
            ``pipeline_version``, ``tools``, ``classification_scheme``,
            ``variants``.
        output_path: Path to write the JSON-LD audit file (``*.jsonld``).
        pretty: If True, indent JSON output (default True).

    Returns:
        None.

    References:
        JSON-LD 1.1: https://www.w3.org/TR/json-ld11/
        PROV-O: https://www.w3.org/TR/prov-o/
    """
    now_iso = datetime.now(timezone.utc).isoformat()

    json_ld: dict[str, Any] = {
        "@context": _JSON_LD_CONTEXT,
        "@type": "prov:Activity",
        "@id": f"urn:genomeforge:audit:{audit_data.get('sample_id', 'unknown')}:{audit_data.get('report_date', '')}",
        "generated_at": now_iso,
        "report_date": audit_data.get("report_date", ""),
        "prov:wasAssociatedWith": {
            "@type": "prov:Agent",
            "name": "GenomeForge",
            "version": audit_data.get("pipeline_version", "unknown"),
            "software": "https://github.com/genomeforge/genomeforge",
        },
        "sample": {
            "@type": "genomeforge:Sample",
            "patient_id": audit_data.get("patient_id", ""),
            "sample_id": audit_data.get("sample_id", ""),
        },
        "data_sources": {
            "genome_assembly": audit_data.get("assembly", "GRCh38"),
            "vep_version": audit_data.get("tools", {}).get("vep", "111"),
            "gnomad_version": audit_data.get("tools", {}).get("gnomad", "4.1"),
            "gnomad_date": "April 2024",
            "gnomad_individuals": 807162,
            "alphamissense": "Cheng et al. 2023 Science PMID:37703350",
            "clinvar_access_date": audit_data.get("report_date", ""),
            "mane_select": "Morales et al. 2022 Nature Methods PMID:35356062",
        },
        "classification_scheme": audit_data.get(
            "classification_scheme",
            {
                "framework": "ACGS 2024 v1.2",
                "point_system": "Tavtigian et al. 2020 PMID:32645316",
                "pm2_weight": "Supporting (1 pt) — ClinGen SVI 2024",
                "pp3_bp4_primary": "AlphaMissense ≥0.564→PP3, ≤0.340→BP4 (Cheng 2023)",
            },
        ),
        "variants": [
            {
                "@type": "genomeforge:VariantClassification",
                "gene": v.get("gene", ""),
                "transcript": v.get("transcript", ""),
                "hgvsc": v.get("hgvsc", ""),
                "hgvsp": v.get("hgvsp", ""),
                "acmg_class": v.get("acmg_class", ""),
                "rules_applied": v.get("rules_applied", []),
                "gnomad_af": v.get("gnomad_af"),
                "alphamissense_score": v.get("alphamissense_score"),
                "posterior_p": v.get("posterior_p"),
                "vus_review_date": v.get("vus_review_date"),
                "pending_clinvar_submission": v.get("pending_clinvar_submission", False),
            }
            for v in audit_data.get("variants", [])
        ],
        "system_info": {
            "hostname": _get_hostname(),
            "user": os.getenv("USER", os.getenv("USERNAME", "unknown")),
        },
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    indent = 2 if pretty else None
    output_path.write_text(
        json.dumps(json_ld, indent=indent, default=str, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info("JSON-LD audit trail written: %s", output_path)


def _get_hostname() -> str:
    """Return the current hostname for audit trail system_info.

    Returns:
        Hostname string, or ``"unknown"`` on failure.
    """
    try:
        return socket.gethostname()
    except Exception:
        return "unknown"


def append_audit_event(
    audit_path: Path,
    event_type: str,
    event_data: dict[str, Any],
) -> None:
    """Append an event to an existing JSON-LD audit trail file.

    Used to record post-report events such as:
    - ClinVar submission confirmation.
    - VUS reclassification after review.
    - Report amendment.

    Args:
        audit_path: Path to the existing JSON-LD audit file.
        event_type: Event type string (e.g. ``"clinvar_submission"``,
            ``"vus_reclassification"``).
        event_data: Dict with event-specific data.

    Returns:
        None.

    Raises:
        FileNotFoundError: If audit_path does not exist.
    """
    if not audit_path.exists():
        raise FileNotFoundError(f"Audit file not found: {audit_path}")

    with audit_path.open("r", encoding="utf-8") as fh:
        existing = json.load(fh)

    events: list[dict[str, Any]] = existing.setdefault("events", [])
    events.append({
        "@type": f"genomeforge:{event_type}",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "data": event_data,
    })

    audit_path.write_text(
        json.dumps(existing, indent=2, default=str, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info("Audit event '%s' appended to %s", event_type, audit_path)
