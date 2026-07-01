"""
phenopackets_input.tests.test_fhir_mapper
==========================================
pytest tests for Phenopackets v2.0 ↔ FHIR R4 bidirectional mapping.

Tests cover:
    - phenopacket_to_fhir: produces Patient, Condition, Observation resources.
    - fhir_to_phenopacket: reconstructs phenopacket from FHIR resources.
    - _subject_to_fhir_patient: sex mapping, date of birth, vital status.
    - _disease_to_fhir_condition: OMIM/Orphanet/MONDO system URIs, onset.
    - _feature_to_fhir_observation: present/excluded HPO features.

References:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    HL7 FHIR R4: https://www.hl7.org/fhir/R4/
"""

from __future__ import annotations

import pytest

from phenopackets_input.fhir_mapper import (
    _disease_to_fhir_condition,
    _feature_to_fhir_observation,
    _subject_to_fhir_patient,
    fhir_to_phenopacket,
    phenopacket_to_fhir,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_PHENOPACKET = {
    "id": "PP-BRCA1-001",
    "subject": {
        "id": "PATIENT-001",
        "sex": "FEMALE",
        "dateOfBirth": "1980-06-15T00:00:00Z",
    },
    "phenotypicFeatures": [
        {"type": {"id": "HP:0003002", "label": "Breast carcinoma"}, "excluded": False},
        {"type": {"id": "HP:0002664", "label": "Neoplasm"}, "excluded": True},
    ],
    "diseases": [
        {
            "term": {
                "id": "OMIM:604370",
                "label": "Hereditary Breast and Ovarian Cancer",
            }
        }
    ],
    "metaData": {"created": "2024-01-01T00:00:00Z", "createdBy": "test"},
}

SAMPLE_FHIR_PATIENT = {
    "resourceType": "Patient",
    "id": "PATIENT-001",
    "gender": "female",
    "birthDate": "1980-06-15",
}


# ---------------------------------------------------------------------------
# phenopacket_to_fhir tests
# ---------------------------------------------------------------------------


class TestPhenopacketToFhir:
    """Tests for phenopacket_to_fhir()."""

    def test_returns_patient_condition_observation_keys(self) -> None:
        """Result has Patient, Condition, Observation keys."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert "Patient" in result
        assert "Condition" in result
        assert "Observation" in result

    def test_patient_resource_created(self) -> None:
        """One Patient resource is created."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert len(result["Patient"]) == 1

    def test_patient_resource_type(self) -> None:
        """Patient resourceType is 'Patient'."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert result["Patient"][0]["resourceType"] == "Patient"

    def test_patient_gender_female(self) -> None:
        """FEMALE sex maps to 'female' FHIR gender."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert result["Patient"][0]["gender"] == "female"

    def test_patient_id_from_subject(self) -> None:
        """Patient.id matches subject.id."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert result["Patient"][0]["id"] == "PATIENT-001"

    def test_condition_created_for_disease(self) -> None:
        """One Condition resource is created per disease."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert len(result["Condition"]) == 1

    def test_condition_resource_type(self) -> None:
        """Condition resourceType is 'Condition'."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert result["Condition"][0]["resourceType"] == "Condition"

    def test_observations_created_for_features(self) -> None:
        """Two Observation resources created for two phenotypicFeatures."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert len(result["Observation"]) == 2

    def test_observation_resource_type(self) -> None:
        """Observation resourceType is 'Observation'."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        assert result["Observation"][0]["resourceType"] == "Observation"

    def test_excluded_feature_has_value_boolean_false(self) -> None:
        """Excluded feature → valueBoolean=False in Observation."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        excluded_obs = [o for o in result["Observation"] if not o.get("valueBoolean", True)]
        assert len(excluded_obs) == 1

    def test_present_feature_has_value_boolean_true(self) -> None:
        """Present (non-excluded) feature → valueBoolean=True."""
        result = phenopacket_to_fhir(SAMPLE_PHENOPACKET)
        present_obs = [o for o in result["Observation"] if o.get("valueBoolean") is True]
        assert len(present_obs) == 1

    def test_empty_phenopacket_produces_empty_patient(self) -> None:
        """Empty phenopacket still produces a Patient resource."""
        result = phenopacket_to_fhir({})
        assert len(result["Patient"]) == 1

    def test_no_diseases_produces_no_conditions(self) -> None:
        """No diseases → empty Condition list."""
        pp = dict(SAMPLE_PHENOPACKET)
        pp["diseases"] = []
        result = phenopacket_to_fhir(pp)
        assert result["Condition"] == []


# ---------------------------------------------------------------------------
# _subject_to_fhir_patient tests
# ---------------------------------------------------------------------------


class TestSubjectToFhirPatient:
    """Tests for _subject_to_fhir_patient()."""

    def test_male_sex_mapping(self) -> None:
        """MALE sex maps to 'male' FHIR gender."""
        subject = {"id": "P1", "sex": "MALE"}
        patient = _subject_to_fhir_patient(subject, "P1")
        assert patient["gender"] == "male"

    def test_unknown_sex_mapping(self) -> None:
        """UNKNOWN_SEX maps to 'unknown' FHIR gender."""
        subject = {"id": "P1", "sex": "UNKNOWN_SEX"}
        patient = _subject_to_fhir_patient(subject, "P1")
        assert patient["gender"] == "unknown"

    def test_date_of_birth_extracted(self) -> None:
        """dateOfBirth is truncated to date part."""
        subject = {"id": "P1", "sex": "FEMALE", "dateOfBirth": "1990-03-25T00:00:00Z"}
        patient = _subject_to_fhir_patient(subject, "P1")
        assert patient["birthDate"] == "1990-03-25"

    def test_deceased_patient(self) -> None:
        """DECEASED vital status sets deceasedBoolean=True."""
        subject = {"id": "P1", "sex": "MALE", "vitalStatus": {"status": "DECEASED"}}
        patient = _subject_to_fhir_patient(subject, "P1")
        assert patient["deceasedBoolean"] is True

    def test_no_date_of_birth_absent(self) -> None:
        """Missing dateOfBirth → birthDate key absent from Patient."""
        subject = {"id": "P1", "sex": "FEMALE"}
        patient = _subject_to_fhir_patient(subject, "P1")
        assert "birthDate" not in patient

    def test_patient_id_set(self) -> None:
        """Patient.id matches the provided patient_id."""
        subject = {"id": "SUBJ-1", "sex": "MALE"}
        patient = _subject_to_fhir_patient(subject, "PATIENT-1")
        assert patient["id"] == "PATIENT-1"


# ---------------------------------------------------------------------------
# _disease_to_fhir_condition tests
# ---------------------------------------------------------------------------


class TestDiseaseToFhirCondition:
    """Tests for _disease_to_fhir_condition()."""

    def test_omim_system_uri(self) -> None:
        """OMIM: prefix maps to OMIM system URI."""
        disease = {"term": {"id": "OMIM:604370", "label": "HBOC"}}
        condition = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert condition is not None
        codings = condition["code"]["coding"]
        assert any("omim" in c.get("system", "").lower() for c in codings)

    def test_orphanet_system_uri(self) -> None:
        """Orphanet: prefix maps to Orphanet system URI."""
        disease = {"term": {"id": "Orphanet:119", "label": "Familial hypercholesterolaemia"}}
        condition = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert condition is not None
        codings = condition["code"]["coding"]
        assert any("orphadata" in c.get("system", "").lower() for c in codings)

    def test_mondo_system_uri(self) -> None:
        """MONDO: prefix maps to MONDO system URI."""
        disease = {"term": {"id": "MONDO:0010297", "label": "BRCA1 syndrome"}}
        condition = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert condition is not None
        codings = condition["code"]["coding"]
        assert any("mondo" in c.get("system", "").lower() for c in codings)

    def test_returns_none_when_no_term_id(self) -> None:
        """Returns None when term.id is absent."""
        disease = {"term": {}}
        result = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert result is None

    def test_condition_subject_set(self) -> None:
        """Condition.subject references the Patient."""
        disease = {"term": {"id": "OMIM:604370", "label": "HBOC"}}
        patient_ref = {"reference": "Patient/PATIENT-001"}
        condition = _disease_to_fhir_condition(disease, patient_ref, "PP-001")
        assert condition is not None
        assert condition["subject"] == patient_ref

    def test_unknown_prefix_uses_empty_system(self) -> None:
        """A disease term ID with an unrecognised prefix maps to an empty system URI."""
        disease = {"term": {"id": "ICD10:E11.9", "label": "Type 2 diabetes"}}
        condition = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert condition is not None
        codings = condition["code"]["coding"]
        assert codings[0]["system"] == ""
        assert codings[0]["code"] == "ICD10:E11.9"

    def test_onset_string_set(self) -> None:
        """Onset label is set as onsetString when present."""
        disease = {
            "term": {"id": "OMIM:604370", "label": "HBOC"},
            "onset": {"ontologyClass": {"id": "HP:0011462", "label": "Young adult onset"}},
        }
        condition = _disease_to_fhir_condition(disease, {"reference": "Patient/P1"}, "PP-001")
        assert condition is not None
        assert condition.get("onsetString") == "Young adult onset"


# ---------------------------------------------------------------------------
# _feature_to_fhir_observation tests
# ---------------------------------------------------------------------------


class TestFeatureToFhirObservation:
    """Tests for _feature_to_fhir_observation()."""

    def test_returns_none_when_no_term_id(self) -> None:
        """Returns None when feature type.id is absent."""
        feature = {"type": {}, "excluded": False}
        result = _feature_to_fhir_observation(feature, {"reference": "Patient/P1"}, "PP-001", 0)
        assert result is None

    def test_excluded_feature_has_interpretation(self) -> None:
        """Excluded feature (absent) has interpretation coding."""
        feature = {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": True}
        obs = _feature_to_fhir_observation(feature, {"reference": "Patient/P1"}, "PP-001", 0)
        assert obs is not None
        assert "interpretation" in obs

    def test_present_feature_no_interpretation(self) -> None:
        """Present feature does not have interpretation key."""
        feature = {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False}
        obs = _feature_to_fhir_observation(feature, {"reference": "Patient/P1"}, "PP-001", 0)
        assert obs is not None
        assert "interpretation" not in obs

    def test_observation_id_uses_index(self) -> None:
        """Observation ID includes the provided index."""
        feature = {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False}
        obs = _feature_to_fhir_observation(feature, {"reference": "Patient/P1"}, "PP-001", 3)
        assert obs is not None
        assert "3" in obs["id"]

    def test_observation_hpo_system_in_code(self) -> None:
        """HPO system URI appears in Observation.code.coding."""
        feature = {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False}
        obs = _feature_to_fhir_observation(feature, {"reference": "Patient/P1"}, "PP-001", 0)
        assert obs is not None
        codings = obs["code"]["coding"]
        assert any("hp" in c.get("system", "").lower() for c in codings)


# ---------------------------------------------------------------------------
# fhir_to_phenopacket tests
# ---------------------------------------------------------------------------


class TestFhirToPhenopacket:
    """Tests for fhir_to_phenopacket()."""

    def test_produces_required_phenopacket_keys(self) -> None:
        """Output has id, subject, phenotypicFeatures, diseases, metaData."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert "id" in pp
        assert "subject" in pp
        assert "phenotypicFeatures" in pp
        assert "diseases" in pp
        assert "metaData" in pp

    def test_subject_id_from_patient(self) -> None:
        """subject.id comes from Patient.id."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert pp["subject"]["id"] == "PATIENT-001"

    def test_female_gender_maps_to_female_sex(self) -> None:
        """FHIR 'female' gender maps back to Phenopackets 'FEMALE'."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert pp["subject"]["sex"] == "FEMALE"

    def test_male_gender_maps_to_male_sex(self) -> None:
        """FHIR 'male' maps back to Phenopackets 'MALE'."""
        patient = {"id": "P2", "gender": "male"}
        pp = fhir_to_phenopacket(patient)
        assert pp["subject"]["sex"] == "MALE"

    def test_birth_date_transferred(self) -> None:
        """Patient.birthDate is transferred to subject.dateOfBirth."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert "dateOfBirth" in pp["subject"]
        assert pp["subject"]["dateOfBirth"].startswith("1980-06-15")

    def test_phenopacket_id_defaults_to_patient_id(self) -> None:
        """Phenopacket id defaults to Patient.id when not specified."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert pp["id"] == "PATIENT-001"

    def test_custom_phenopacket_id(self) -> None:
        """Custom phenopacket_id overrides Patient.id."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT, phenopacket_id="PP-CUSTOM-001")
        assert pp["id"] == "PP-CUSTOM-001"

    def test_hpo_observations_converted_to_features(self) -> None:
        """FHIR Observations with HPO coding become phenotypicFeatures."""
        observations = [
            {
                "resourceType": "Observation",
                "code": {
                    "coding": [
                        {
                            "system": "http://purl.obolibrary.org/obo/hp.owl",
                            "code": "HP:0001250",
                            "display": "Seizures",
                        }
                    ]
                },
                "valueBoolean": True,
            }
        ]
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT, observations=observations)
        assert len(pp["phenotypicFeatures"]) == 1
        assert pp["phenotypicFeatures"][0]["type"]["id"] == "HP:0001250"
        assert pp["phenotypicFeatures"][0]["excluded"] is False

    def test_conditions_converted_to_diseases(self) -> None:
        """FHIR Conditions become diseases in the phenopacket."""
        conditions = [
            {
                "resourceType": "Condition",
                "code": {
                    "coding": [
                        {"code": "OMIM:604370", "display": "HBOC"}
                    ]
                },
            }
        ]
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT, conditions=conditions)
        assert len(pp["diseases"]) == 1
        assert pp["diseases"][0]["term"]["id"] == "OMIM:604370"

    def test_metadata_schema_version_is_v2(self) -> None:
        """metaData.phenopacketSchemaVersion is '2.0'."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT)
        assert pp["metaData"]["phenopacketSchemaVersion"] == "2.0"

    def test_metadata_created_by_set(self) -> None:
        """metaData.createdBy reflects the created_by argument."""
        pp = fhir_to_phenopacket(SAMPLE_FHIR_PATIENT, created_by="test-system")
        assert pp["metaData"]["createdBy"] == "test-system"
