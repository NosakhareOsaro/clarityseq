"""
phenopackets_input.fhir_mapper
================================
Bidirectional mapping between Phenopackets v2.0 and FHIR R4.

Implements bidirectional conversion:
    - Phenopackets v2.0 → FHIR R4 resources (Patient, Condition, Observation).
    - FHIR R4 resources → Phenopackets v2.0.

Why both formats:
    Phenopackets: exchange format optimised for rare disease / genomics.
    FHIR R4: HL7 interoperability standard used by NHS/EHR systems.
    GenomeForge supports both for NHS GMS data exchange.

Mapping strategy:
    Phenopackets subject   ↔  FHIR Patient resource.
    Phenopackets diseases  ↔  FHIR Condition resources.
    Phenotypic features    ↔  FHIR Observation resources (HPO-coded).
    metaData.created_by    ↔  FHIR Provenance.agent.

Limitations:
    - Phenopackets family pedigree → FHIR FamilyMemberHistory (partial).
    - Complex structural variants (VRS) require FHIR MolecularSequence.
    - FHIR Genomics Reporting IG terms used where available.

References:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    HL7 FHIR R4: https://www.hl7.org/fhir/R4/
    FHIR Genomics Reporting IG: https://build.fhir.org/ig/HL7/genomics-reporting/
    NHS Digital FHIR profile: https://simplifier.net/nhsdigital
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# FHIR R4 system URIs
# ---------------------------------------------------------------------------

_HPO_SYSTEM = "http://purl.obolibrary.org/obo/hp.owl"
_OMIM_SYSTEM = "https://www.omim.org"
_ORPHANET_SYSTEM = "http://www.orphadata.org/cgi-bin/rare_main.php"
_MONDO_SYSTEM = "http://purl.obolibrary.org/obo/mondo.owl"
_SNOMED_SYSTEM = "http://snomed.info/sct"
_LOINC_SYSTEM = "http://loinc.org"

# HPO to FHIR clinical status mapping
_HPO_CLINICAL_STATUS = "http://terminology.hl7.org/CodeSystem/condition-clinical"

# Sex mapping Phenopackets → FHIR
_SEX_TO_FHIR: dict[str, str] = {
    "FEMALE": "female",
    "MALE": "male",
    "OTHER_SEX": "other",
    "UNKNOWN_SEX": "unknown",
}

# Sex mapping FHIR → Phenopackets
_FHIR_TO_SEX: dict[str, str] = {v: k for k, v in _SEX_TO_FHIR.items()}


# ---------------------------------------------------------------------------
# Phenopackets → FHIR R4
# ---------------------------------------------------------------------------


def phenopacket_to_fhir(
    phenopacket: dict[str, Any],
    fhir_server_base: str = "https://fhir.example.org",
) -> dict[str, list[dict[str, Any]]]:
    """Convert a Phenopackets v2.0 document to a bundle of FHIR R4 resources.

    Converts the phenopacket subject to a FHIR Patient, diseases to
    FHIR Conditions, and phenotypic features to FHIR Observations.

    Args:
        phenopacket: Parsed Phenopackets v2.0 JSON dict.
        fhir_server_base: Base URL for FHIR resource references
            (e.g. ``"https://fhir.example.org"``).

    Returns:
        Dict with keys ``"Patient"``, ``"Condition"``, ``"Observation"``,
        each mapping to a list of FHIR R4 resource dicts.

    References:
        Jacobsen et al. 2022 PMID:35705716 (Phenopackets v2).
        HL7 FHIR R4 Patient: https://www.hl7.org/fhir/R4/patient.html
        HL7 FHIR R4 Condition: https://www.hl7.org/fhir/R4/condition.html
    """
    result: dict[str, list[dict[str, Any]]] = {
        "Patient": [],
        "Condition": [],
        "Observation": [],
    }

    subject = phenopacket.get("subject", {})
    phenopacket_id = phenopacket.get("id", "")
    patient_id = subject.get("id", phenopacket_id)

    # --- Patient resource ---
    patient = _subject_to_fhir_patient(subject, patient_id)
    result["Patient"].append(patient)

    patient_ref = {"reference": f"Patient/{patient_id}"}

    # --- Conditions (diseases) ---
    for disease in phenopacket.get("diseases", []):
        condition = _disease_to_fhir_condition(disease, patient_ref, phenopacket_id)
        if condition:
            result["Condition"].append(condition)

    # --- Observations (phenotypic features) ---
    for i, feature in enumerate(phenopacket.get("phenotypicFeatures", [])):
        obs = _feature_to_fhir_observation(feature, patient_ref, phenopacket_id, i)
        if obs:
            result["Observation"].append(obs)

    return result


def _subject_to_fhir_patient(
    subject: dict[str, Any],
    patient_id: str,
) -> dict[str, Any]:
    """Convert a Phenopackets v2.0 subject to a FHIR R4 Patient resource.

    Args:
        subject: Phenopackets subject dict.
        patient_id: FHIR Patient resource ID.

    Returns:
        FHIR R4 Patient resource dict.
    """
    sex = subject.get("sex", "UNKNOWN_SEX")
    fhir_gender = _SEX_TO_FHIR.get(sex, "unknown")

    patient: dict[str, Any] = {
        "resourceType": "Patient",
        "id": patient_id,
        "meta": {
            "profile": [
                "https://fhir.hl7.org.uk/StructureDefinition/UKCore-Patient"
            ]
        },
        "gender": fhir_gender,
    }

    # Date of birth
    dob = subject.get("dateOfBirth", "")
    if dob:
        patient["birthDate"] = dob[:10]  # Take date part only

    # VitalStatus
    vital_status = subject.get("vitalStatus", {})
    if vital_status.get("status") == "DECEASED":
        patient["deceasedBoolean"] = True

    return patient


def _disease_to_fhir_condition(
    disease: dict[str, Any],
    patient_ref: dict[str, str],
    phenopacket_id: str,
) -> dict[str, Any] | None:
    """Convert a Phenopackets v2.0 disease to a FHIR R4 Condition resource.

    Args:
        disease: Phenopackets disease dict with ``term`` and optional ``onset``.
        patient_ref: FHIR Patient reference dict.
        phenopacket_id: Parent phenopacket ID for resource ID generation.

    Returns:
        FHIR R4 Condition resource dict, or None if disease term is absent.
    """
    term = disease.get("term", {})
    term_id = term.get("id", "")
    if not term_id:
        return None

    # Map disease ID to FHIR coding system
    if term_id.startswith("OMIM:"):
        system = _OMIM_SYSTEM
    elif term_id.startswith("Orphanet:"):
        system = _ORPHANET_SYSTEM
    elif term_id.startswith("MONDO:"):
        system = _MONDO_SYSTEM
    else:
        system = ""

    condition: dict[str, Any] = {
        "resourceType": "Condition",
        "id": f"{phenopacket_id}-condition-{term_id.replace(':', '-')}",
        "meta": {
            "profile": [
                "https://fhir.hl7.org.uk/StructureDefinition/UKCore-Condition"
            ]
        },
        "clinicalStatus": {
            "coding": [
                {
                    "system": _HPO_CLINICAL_STATUS,
                    "code": "active",
                    "display": "Active",
                }
            ]
        },
        "code": {
            "coding": [
                {
                    "system": system,
                    "code": term_id,
                    "display": term.get("label", ""),
                }
            ],
            "text": term.get("label", ""),
        },
        "subject": patient_ref,
    }

    # Onset
    onset = disease.get("onset", {})
    onset_class = onset.get("ontologyClass", {})
    if onset_class.get("id"):
        condition["onsetString"] = onset_class.get("label", onset_class.get("id"))

    return condition


def _feature_to_fhir_observation(
    feature: dict[str, Any],
    patient_ref: dict[str, str],
    phenopacket_id: str,
    index: int,
) -> dict[str, Any] | None:
    """Convert a Phenopackets phenotypicFeature to a FHIR R4 Observation.

    Args:
        feature: Phenopackets phenotypicFeature dict with ``type`` and
            optional ``excluded``.
        patient_ref: FHIR Patient reference dict.
        phenopacket_id: Parent phenopacket ID.
        index: Feature index for unique ID generation.

    Returns:
        FHIR R4 Observation resource dict, or None if term ID is absent.
    """
    term = feature.get("type", {})
    term_id = term.get("id", "")
    if not term_id:
        return None

    excluded = feature.get("excluded", False)

    # FHIR Observation status for HPO features
    # Present → LOINC 8302-2 "Body height" style observation; absent → similar
    obs_status = "final"
    value_boolean = not excluded  # True = present, False = absent

    obs: dict[str, Any] = {
        "resourceType": "Observation",
        "id": f"{phenopacket_id}-obs-{index}",
        "status": obs_status,
        "category": [
            {
                "coding": [
                    {
                        "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                        "code": "exam",
                        "display": "Exam",
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": _HPO_SYSTEM,
                    "code": term_id,
                    "display": term.get("label", ""),
                }
            ],
            "text": term.get("label", ""),
        },
        "subject": patient_ref,
        "valueBoolean": value_boolean,
    }

    if excluded:
        obs["interpretation"] = [
            {
                "coding": [
                    {
                        "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                        "code": "N",
                        "display": "Normal (absent)",
                    }
                ]
            }
        ]

    return obs


# ---------------------------------------------------------------------------
# FHIR R4 → Phenopackets v2.0
# ---------------------------------------------------------------------------


def fhir_to_phenopacket(
    patient: dict[str, Any],
    conditions: list[dict[str, Any]] | None = None,
    observations: list[dict[str, Any]] | None = None,
    phenopacket_id: str | None = None,
    created_by: str = "GenomeForge FHIR mapper",
) -> dict[str, Any]:
    """Convert FHIR R4 resources to a Phenopackets v2.0 document.

    Args:
        patient: FHIR R4 Patient resource dict.
        conditions: List of FHIR R4 Condition resource dicts (optional).
        observations: List of FHIR R4 Observation resource dicts (optional).
        phenopacket_id: ID for the generated phenopacket.  Defaults to
            the Patient.id value.
        created_by: Creator identifier for metaData.createdBy.

    Returns:
        Phenopackets v2.0 JSON dict.

    References:
        Jacobsen et al. 2022 PMID:35705716 (Phenopackets v2).
        HL7 FHIR R4: https://www.hl7.org/fhir/R4/
    """
    conditions = conditions or []
    observations = observations or []
    patient_id = patient.get("id", "unknown")
    pp_id = phenopacket_id or patient_id

    # Subject
    fhir_gender = patient.get("gender", "unknown")
    pp_sex = _FHIR_TO_SEX.get(fhir_gender, "UNKNOWN_SEX")

    subject: dict[str, Any] = {
        "id": patient_id,
        "sex": pp_sex,
    }
    if patient.get("birthDate"):
        subject["dateOfBirth"] = patient["birthDate"] + "T00:00:00Z"

    # Phenotypic features from FHIR Observations
    phenotypic_features: list[dict[str, Any]] = []
    for obs in observations:
        code = obs.get("code", {})
        codings = code.get("coding", [])
        for coding in codings:
            if _HPO_SYSTEM in coding.get("system", ""):
                excluded = not obs.get("valueBoolean", True)
                phenotypic_features.append({
                    "type": {
                        "id": coding.get("code", ""),
                        "label": coding.get("display", code.get("text", "")),
                    },
                    "excluded": excluded,
                })
                break

    # Diseases from FHIR Conditions
    diseases: list[dict[str, Any]] = []
    for cond in conditions:
        code = cond.get("code", {})
        codings = code.get("coding", [])
        for coding in codings:
            diseases.append({
                "term": {
                    "id": coding.get("code", ""),
                    "label": coding.get("display", code.get("text", "")),
                }
            })
            break

    # metaData
    metadata: dict[str, Any] = {
        "created": datetime.now(timezone.utc).isoformat(),
        "createdBy": created_by,
        "phenopacketSchemaVersion": "2.0",
        "resources": [
            {
                "id": "hp",
                "name": "Human Phenotype Ontology",
                "url": "http://purl.obolibrary.org/obo/hp.owl",
                "namespacePrefix": "HP",
                "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
            }
        ],
    }

    return {
        "id": pp_id,
        "subject": subject,
        "phenotypicFeatures": phenotypic_features,
        "diseases": diseases,
        "metaData": metadata,
    }
