"""FHIR R4 Task resource generation for variant reclassification recontact.

This module generates HL7 FHIR R4 Task resources to orchestrate the clinical
recontact workflow when a variant reclassification event is detected.

Specification references:
    - HL7 FHIR Genomics Reporting Implementation Guide v3.0.0 (2024):
      Defines the 'recontact' Task profile (StructureDefinition
      genomics-report-task-rec-followup) for communicating reclassification
      events to ordering clinicians.
      URL: https://hl7.org/fhir/uv/genomics-reporting/STU3/
    - FHIR R4 Task resource: https://hl7.org/fhir/R4/task.html
    - GA4GH VRS v2.0 (Wagner et al. 2025): Variant representation identifiers
      are embedded in Task.input to provide an unambiguous, computable
      reference to the reclassified variant. VRS identifiers use a 24-character
      truncated SHA512 digest of the canonical variant representation.
    - ACGS 2024 §9: The recontact process must be documented and traceable;
      FHIR Task resources provide the audit trail required by the guideline.

Task workflow:
    1. clinvar_diff.py detects a reclassification event.
    2. create_recontact_task() generates a FHIR R4 Task resource.
    3. Task is POSTed to the clinical FHIR server (EHR integration).
    4. Task.id is stored back in ReclassificationEvent.fhir_task_id.
    5. Clinical scientists action the Task to generate a recontact letter.
"""

from __future__ import annotations

import json
import logging
import uuid
from datetime import datetime, timezone
from typing import Any, Optional

from reclassification.models import ClinicalSignificance, ReclassificationEvent

logger = logging.getLogger(__name__)

# HL7 FHIR R4 system URIs
FHIR_TASK_PROFILE = (
    "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/"
    "genomics-report-task-rec-followup"
)
FHIR_GENOMICS_REPORTING_CS = (
    "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs"
)
FHIR_LOINC_SYSTEM = "http://loinc.org"
FHIR_SNOMED_SYSTEM = "http://snomed.info/sct"

# VRS v2.0 namespace URI for Task.input identifiers
VRS_IDENTIFIER_SYSTEM = "https://identifiers.org/ga4gh.vrs"

# LOINC codes for genomics Task inputs
LOINC_VARIANT_ANALYZED = "69548-6"  # Genetic variant assessment
LOINC_RECLASSIFICATION_REASON = "LA26333-7"  # Reclassification reason

# SNOMED code for clinical recontact activity
SNOMED_RECONTACT = "185087000"  # Telephone follow-up (generic recontact proxy)

# Mapping from ClinicalSignificance to LOINC answer codes for Task output
SIGNIFICANCE_LOINC_MAP: dict[str, str] = {
    ClinicalSignificance.PATHOGENIC.value: "LA6668-3",
    ClinicalSignificance.LIKELY_PATHOGENIC.value: "LA26332-9",
    ClinicalSignificance.VUS.value: "LA26333-7",
    ClinicalSignificance.LIKELY_BENIGN.value: "LA26334-5",
    ClinicalSignificance.BENIGN.value: "LA6675-8",
    ClinicalSignificance.CONFLICTING.value: "LA26333-7",  # Use VUS code for conflicting
}


def _vrs_identifier_from_variant_id(variant_id: str) -> str:
    """Construct a GA4GH VRS v2.0-style identifier from a variant key.

    In production, this would call compute_vrs_id() from beacon_api.vrs_utils
    after parsing the variant key into (chrom, pos, ref, alt) components.
    Here we produce a deterministic placeholder using the variant_id.

    The VRS v2.0 identifier format is:
        ga4gh:VA.<24-char-base64url-truncated-SHA512>

    Args:
        variant_id: Internal variant identifier, typically in the form
            'chrN:pos:ref:alt'.

    Returns:
        VRS v2.0-format identifier string.
    """
    import hashlib
    import base64

    digest = hashlib.sha512(variant_id.encode()).digest()
    b64 = base64.urlsafe_b64encode(digest).decode().rstrip("=")
    return f"ga4gh:VA.{b64[:24]}"


def _build_task_input_variant(
    variant_id: str,
    vrs_id: str,
) -> dict[str, Any]:
    """Build the Task.input component for the reclassified variant.

    Per FHIR Genomics Reporting IG v3.0.0, Task.input[type=variant]
    provides the machine-readable variant reference using VRS identifier.

    Args:
        variant_id: Internal GenomeForge variant identifier.
        vrs_id: GA4GH VRS v2.0 identifier for the variant.

    Returns:
        FHIR Task.input component dictionary.
    """
    return {
        "type": {
            "coding": [
                {
                    "system": FHIR_LOINC_SYSTEM,
                    "code": LOINC_VARIANT_ANALYZED,
                    "display": "Genetic variant assessment",
                }
            ]
        },
        "valueReference": {
            "identifier": {
                "system": VRS_IDENTIFIER_SYSTEM,
                "value": vrs_id,
            },
            "display": f"Variant: {variant_id}",
        },
        "extension": [
            {
                "url": (
                    "http://hl7.org/fhir/uv/genomics-reporting/"
                    "StructureDefinition/vrs-identifier"
                ),
                "valueString": vrs_id,
            }
        ],
    }


def _build_task_input_reclassification(
    old_class: str,
    new_class: str,
) -> dict[str, Any]:
    """Build the Task.input component describing the reclassification change.

    Args:
        old_class: Previous ClinVar clinical significance value.
        new_class: New ClinVar clinical significance value.

    Returns:
        FHIR Task.input component dictionary for the reclassification reason.
    """
    old_loinc = SIGNIFICANCE_LOINC_MAP.get(old_class, "LA26333-7")
    new_loinc = SIGNIFICANCE_LOINC_MAP.get(new_class, "LA26333-7")

    return {
        "type": {
            "coding": [
                {
                    "system": FHIR_GENOMICS_REPORTING_CS,
                    "code": "rec-followup-reason",
                    "display": "Recommendation follow-up reason",
                }
            ]
        },
        "valueCodeableConcept": {
            "coding": [
                {
                    "system": FHIR_GENOMICS_REPORTING_CS,
                    "code": "variant-reclassified",
                    "display": "Variant reclassified",
                }
            ],
            "text": (
                f"Variant reclassified from '{old_class}' to '{new_class}' "
                f"in ClinVar. Previous LOINC code: {old_loinc}; "
                f"New LOINC code: {new_loinc}."
            ),
        },
    }


def create_recontact_task(
    event: ReclassificationEvent,
    requester_reference: Optional[str] = None,
    patient_reference: Optional[str] = None,
    due_date: Optional[str] = None,
) -> dict[str, Any]:
    """Generate a FHIR R4 Task resource for clinical recontact workflow.

    Creates a structured FHIR R4 Task conforming to the HL7 Genomics
    Reporting IG v3.0.0 'recontact' Task profile. The Task encodes:
    - The reclassified variant using GA4GH VRS v2.0 identifier in Task.input
    - The old and new classifications as coded values
    - The clinical urgency based on the direction of reclassification
    - Requester (laboratory) and owner (ordering clinician) references

    Args:
        event: ReclassificationEvent containing variant_id, old_class,
            new_class, clinvar_date, and clinvar_accession.
        requester_reference: FHIR Reference to the requesting organisation
            (e.g. 'Organization/nhs-genomics-lab-123'). Defaults to a
            placeholder if not provided.
        patient_reference: FHIR Reference to the patient resource
            (e.g. 'Patient/gms-patient-456'). Defaults to a placeholder.
        due_date: ISO 8601 datetime string for Task completion deadline.
            Defaults to 30 days from now for P/LP upgrades, 90 days otherwise.

    Returns:
        FHIR R4 Task resource as a Python dictionary (JSON-serialisable).
        The returned dict has the full Task structure including:
        - resourceType, id, meta, status, intent, priority
        - code (recontact action code)
        - description
        - input[]: variant VRS identifier + reclassification details
        - output[]: placeholder for recording recontact completion

    Example:
        >>> task = create_recontact_task(event)
        >>> json.dumps(task, indent=2)
        '{ "resourceType": "Task", ... }'
    """
    task_id = str(uuid.uuid4())
    now_utc = datetime.now(timezone.utc).isoformat()

    # Compute VRS v2.0 identifier for the variant
    vrs_id = _vrs_identifier_from_variant_id(event.variant_id)

    # Determine task priority based on reclassification direction
    # Upgrades to P/LP require urgent recontact; downgrades are less urgent.
    actionable = {
        ClinicalSignificance.PATHOGENIC.value,
        ClinicalSignificance.LIKELY_PATHOGENIC.value,
    }
    is_upgrade = event.new_class in actionable and event.old_class not in actionable
    is_downgrade = event.old_class in actionable and event.new_class not in actionable
    priority = "urgent" if is_upgrade else "routine"

    # Determine due date
    if due_date is None:
        from datetime import timedelta
        days_offset = 30 if is_upgrade else 90
        due_dt = datetime.now(timezone.utc) + timedelta(days=days_offset)
        due_date = due_dt.date().isoformat()

    # Default references if not supplied
    if requester_reference is None:
        requester_reference = "Organization/nhs-genomics-laboratory"
    if patient_reference is None:
        patient_reference = "Patient/unknown"

    # Build human-readable description
    direction = "upgraded" if is_upgrade else ("downgraded" if is_downgrade else "reclassified")
    description = (
        f"Variant {event.variant_id} has been {direction} in ClinVar "
        f"from '{event.old_class}' to '{event.new_class}' "
        f"(ClinVar date: {event.clinvar_date}). "
        f"Patient recontact {'is required' if event.recontact_required else 'may be required'} "
        f"per ACGS 2024 §9."
    )
    if event.clinvar_accession:
        description += f" ClinVar accession: {event.clinvar_accession}."

    task: dict[str, Any] = {
        "resourceType": "Task",
        "id": task_id,
        "meta": {
            "profile": [FHIR_TASK_PROFILE],
            "lastUpdated": now_utc,
            "source": "https://genomeforge.nhs.uk/reclassification-daemon",
        },
        # Task is 'requested' — requires action from a clinical scientist
        "status": "requested",
        # This is a proposal (filler will decide whether to action)
        "intent": "proposal",
        "priority": priority,
        "code": {
            "coding": [
                {
                    "system": FHIR_GENOMICS_REPORTING_CS,
                    "code": "recontact",
                    "display": "Recontact for variant reclassification",
                }
            ],
            "text": "Recontact patient regarding variant reclassification",
        },
        "description": description,
        "authoredOn": now_utc,
        "lastModified": now_utc,
        "requester": {
            "reference": requester_reference,
            "display": "NHS Genomics Laboratory",
        },
        "for": {
            # The patient who needs to be recontacted
            "reference": patient_reference,
        },
        "restriction": {
            "repetitions": 1,
            "period": {
                "end": due_date,
            },
        },
        "input": [
            # Input 1: The reclassified variant with VRS identifier
            _build_task_input_variant(event.variant_id, vrs_id),
            # Input 2: Description of the reclassification event
            _build_task_input_reclassification(event.old_class, event.new_class),
            # Input 3: ClinVar accession for traceability
            {
                "type": {
                    "coding": [
                        {
                            "system": FHIR_LOINC_SYSTEM,
                            "code": "81252-9",
                            "display": "Discrete genetic variant",
                        }
                    ]
                },
                "valueString": event.clinvar_accession or "Unknown",
                "_valueString": {
                    "extension": [
                        {
                            "url": "http://hl7.org/fhir/StructureDefinition/data-absent-reason",
                            "valueCode": "unknown" if not event.clinvar_accession else None,
                        }
                    ] if not event.clinvar_accession else []
                },
            },
            # Input 4: ClinVar evaluation date
            {
                "type": {
                    "coding": [
                        {
                            "system": FHIR_LOINC_SYSTEM,
                            "code": "93044-6",
                            "display": "Level of evidence",
                        }
                    ]
                },
                "valueDate": str(event.clinvar_date),
            },
        ],
        "output": [
            # Placeholder output to be populated when recontact is completed
            {
                "type": {
                    "coding": [
                        {
                            "system": FHIR_GENOMICS_REPORTING_CS,
                            "code": "recontact-completed",
                            "display": "Recontact action completed",
                        }
                    ]
                },
                # valueBoolean will be set to true when Task is completed
            }
        ],
        "note": [
            {
                "text": (
                    f"Generated by GenomeForge reclassification daemon. "
                    f"VRS v2.0 identifier: {vrs_id}. "
                    f"Reclassification detected: "
                    f"{event.old_class} -> {event.new_class}."
                )
            }
        ],
    }

    logger.info(
        "Created FHIR R4 Task %s for reclassification event: %s -> %s (variant=%s)",
        task_id, event.old_class, event.new_class, event.variant_id,
    )

    return task


def task_to_json(task: dict[str, Any], indent: int = 2) -> str:
    """Serialise a FHIR Task dictionary to JSON string.

    Args:
        task: FHIR Task dictionary as returned by create_recontact_task().
        indent: JSON indentation level. Defaults to 2.

    Returns:
        JSON string representation of the FHIR Task resource.
    """
    return json.dumps(task, indent=indent, default=str)


def validate_task_structure(task: dict[str, Any]) -> list[str]:
    """Validate a FHIR Task dictionary against basic structural requirements.

    Performs lightweight validation of the generated Task resource to catch
    structural errors before POSTing to the FHIR server. Does not perform
    full HL7 validation (use HAPI FHIR validator for that).

    Args:
        task: FHIR Task dictionary to validate.

    Returns:
        List of validation error messages. Empty list means no errors found.

    Example:
        >>> errors = validate_task_structure(task)
        >>> assert errors == [], f"Task validation failed: {errors}"
    """
    errors: list[str] = []

    # Required top-level fields
    required_fields = [
        "resourceType", "id", "status", "intent", "code", "input"
    ]
    for field in required_fields:
        if field not in task:
            errors.append(f"Missing required field: {field}")

    # resourceType must be 'Task'
    if task.get("resourceType") != "Task":
        errors.append(
            f"resourceType must be 'Task', got: {task.get('resourceType')!r}"
        )

    # Status must be a valid FHIR Task status
    valid_statuses = {
        "draft", "requested", "received", "accepted", "rejected",
        "ready", "cancelled", "in-progress", "on-hold", "failed",
        "completed", "entered-in-error",
    }
    if task.get("status") not in valid_statuses:
        errors.append(f"Invalid Task.status: {task.get('status')!r}")

    # Intent must be valid
    valid_intents = {
        "unknown", "proposal", "plan", "directive", "order",
        "original-order", "reflex-order", "filler-order", "instance-order",
        "option",
    }
    if task.get("intent") not in valid_intents:
        errors.append(f"Invalid Task.intent: {task.get('intent')!r}")

    # Must have at least one input
    if not task.get("input"):
        errors.append("Task must have at least one input component")

    # Check VRS identifier is present in first input
    inputs = task.get("input", [])
    if inputs:
        first_input = inputs[0]
        val_ref = first_input.get("valueReference", {})
        if not val_ref.get("identifier", {}).get("value", "").startswith("ga4gh:VA."):
            errors.append(
                "First Task.input must contain a GA4GH VRS v2.0 identifier "
                "(ga4gh:VA.<digest>)"
            )

    return errors
