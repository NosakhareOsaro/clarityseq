"""
phenopackets_input.tests.test_schema_validator
================================================
pytest tests for Phenopackets v2.0 schema validation.

Tests cover:
    - validate_phenopacket_dict: valid and invalid dicts.
    - validate_phenopacket: file-based validation with mocked phenopacket-tools.
    - HPO term format checking.
    - Disease term prefix validation.
    - metaData required fields.

References:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def valid_phenopacket_v2() -> dict[str, Any]:
    """Return a minimal valid Phenopackets v2.0 dict.

    Returns:
        Dict conforming to Phenopackets v2.0 required fields.
    """
    return {
        "id": "test-phenopacket-001",
        "subject": {
            "id": "patient-001",
            "sex": "FEMALE",
        },
        "phenotypicFeatures": [
            {
                "type": {
                    "id": "HP:0001250",
                    "label": "Seizure",
                },
                "excluded": False,
            },
            {
                "type": {
                    "id": "HP:0004322",
                    "label": "Short stature",
                },
                "excluded": False,
            },
        ],
        "diseases": [
            {
                "term": {
                    "id": "OMIM:605711",
                    "label": "DRAVET SYNDROME",
                },
            }
        ],
        "metaData": {
            "created": "2024-12-13T00:00:00Z",
            "createdBy": "GenomeForge",
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
        },
    }


@pytest.fixture
def missing_metadata_phenopacket() -> dict[str, Any]:
    """Return a phenopacket dict with missing metaData keys.

    Returns:
        Dict missing required metaData fields.
    """
    return {
        "id": "bad-pp-001",
        "subject": {"id": "patient-002"},
        "metaData": {
            "created": "2024-01-01T00:00:00Z",
            # Missing: createdBy, phenopacketSchemaVersion, resources
        },
    }


@pytest.fixture
def invalid_hpo_phenopacket() -> dict[str, Any]:
    """Return a phenopacket with invalid HPO term ID format.

    Returns:
        Dict with non-HP: prefixed phenotypic feature type.
    """
    return {
        "id": "bad-hpo-001",
        "subject": {"id": "patient-003"},
        "phenotypicFeatures": [
            {
                "type": {
                    "id": "SNOMED:44054006",  # SNOMED not HPO — should warn/error
                    "label": "Diabetes mellitus",
                },
                "excluded": False,
            }
        ],
        "metaData": {
            "created": "2024-01-01T00:00:00Z",
            "createdBy": "test",
            "phenopacketSchemaVersion": "2.0",
            "resources": [],
        },
    }


# ---------------------------------------------------------------------------
# validate_phenopacket_dict tests
# ---------------------------------------------------------------------------


class TestValidatePhenopacketDict:
    """Tests for validate_phenopacket_dict (in-memory validation)."""

    def test_valid_phenopacket_returns_valid(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """Valid Phenopackets v2.0 dict returns valid=True.

        Reference: Jacobsen et al. 2022 PMID:35705716.
        """
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        result = validate_phenopacket_dict(valid_phenopacket_v2)
        assert result.valid is True, (
            f"Valid phenopacket should return valid=True. Errors: {result.errors}"
        )
        assert result.phenopacket_id == "test-phenopacket-001"

    def test_missing_subject_id_returns_error(self) -> None:
        """Phenopacket missing subject.id returns validation error."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp = {
            "id": "pp-no-subject-id",
            "subject": {},  # Missing id
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": [],
            },
        }
        result = validate_phenopacket_dict(pp)
        assert result.valid is False
        assert any("subject.id" in e for e in result.errors), (
            f"Expected 'subject.id' error, got: {result.errors}"
        )

    def test_missing_top_level_keys_returns_error(self) -> None:
        """Phenopacket missing required top-level keys returns errors."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp: dict[str, Any] = {}  # Missing all required keys
        result = validate_phenopacket_dict(pp)
        assert result.valid is False
        # Should report missing id, subject, metaData
        assert len(result.errors) >= 1

    def test_missing_metadata_keys_returns_errors(
        self, missing_metadata_phenopacket: dict[str, Any]
    ) -> None:
        """Phenopacket with incomplete metaData returns validation errors."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        result = validate_phenopacket_dict(missing_metadata_phenopacket)
        assert result.valid is False
        assert any("metaData" in e or "createdBy" in e or "resources" in e
                   for e in result.errors), (
            f"Expected metaData-related errors, got: {result.errors}"
        )

    def test_invalid_hpo_format_returns_error(
        self, invalid_hpo_phenopacket: dict[str, Any]
    ) -> None:
        """Non-HP: prefixed term ID returns validation error.

        HPO terms must use HP:XXXXXXX format per Phenopackets v2.0.
        """
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        result = validate_phenopacket_dict(invalid_hpo_phenopacket)
        assert result.valid is False
        assert any("HP:" in e or "HPO" in e for e in result.errors), (
            f"Expected HPO format error, got: {result.errors}"
        )

    def test_excluded_hpo_feature_accepted(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """Phenopacket with excluded=True HPO features is valid."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp = dict(valid_phenopacket_v2)
        pp["phenotypicFeatures"] = [
            {
                "type": {"id": "HP:0001250", "label": "Seizure"},
                "excluded": True,  # Patient does NOT have seizures
            }
        ]
        result = validate_phenopacket_dict(pp)
        assert result.valid is True, (
            f"Excluded HPO feature should be valid. Errors: {result.errors}"
        )

    def test_schema_version_v1_raises_error(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """Phenopackets v1 schema version string returns an error.

        This validator targets Phenopackets v2.0 only.
        """
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp = dict(valid_phenopacket_v2)
        pp["metaData"] = dict(pp["metaData"])
        pp["metaData"]["phenopacketSchemaVersion"] = "1.0.0-RC3"

        result = validate_phenopacket_dict(pp)
        assert result.valid is False
        assert any("schema" in e.lower() or "version" in e.lower() for e in result.errors)

    def test_valid_orphanet_disease_id(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """Phenopacket with Orphanet disease ID passes validation."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp = dict(valid_phenopacket_v2)
        pp["diseases"] = [
            {"term": {"id": "Orphanet:98895", "label": "Dravet syndrome"}}
        ]
        result = validate_phenopacket_dict(pp)
        assert result.valid is True, (
            f"Orphanet disease ID should be valid. Errors: {result.errors}"
        )

    def test_valid_mondo_disease_id(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """Phenopacket with MONDO disease ID passes validation."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        pp = dict(valid_phenopacket_v2)
        pp["diseases"] = [
            {"term": {"id": "MONDO:0100135", "label": "Dravet syndrome"}}
        ]
        result = validate_phenopacket_dict(pp)
        assert result.valid is True, (
            f"MONDO disease ID should be valid. Errors: {result.errors}"
        )

    def test_validation_result_str_method(
        self, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """ValidationResult __str__ includes phenopacket ID and status."""
        from phenopackets_input.schema_validator import validate_phenopacket_dict

        result = validate_phenopacket_dict(valid_phenopacket_v2)
        s = str(result)
        assert "test-phenopacket-001" in s
        assert "VALID" in s


# ---------------------------------------------------------------------------
# File-based validation tests
# ---------------------------------------------------------------------------


class TestValidatePhenopacketFile:
    """Tests for validate_phenopacket (file path-based)."""

    def test_file_not_found_raises(self) -> None:
        """validate_phenopacket raises FileNotFoundError for missing file."""
        from phenopackets_input.schema_validator import validate_phenopacket

        with pytest.raises(FileNotFoundError):
            validate_phenopacket(Path("/nonexistent/path/pp.json"))

    def test_invalid_json_returns_error(self, tmp_path: Path) -> None:
        """validate_phenopacket with malformed JSON returns valid=False."""
        from phenopackets_input.schema_validator import validate_phenopacket

        bad_file = tmp_path / "bad.json"
        bad_file.write_text("{not: valid json}", encoding="utf-8")

        result = validate_phenopacket(bad_file)
        assert result.valid is False
        assert any("JSON" in e for e in result.errors)

    def test_valid_file_returns_valid(
        self,
        tmp_path: Path,
        valid_phenopacket_v2: dict[str, Any],
    ) -> None:
        """validate_phenopacket with valid JSON file returns valid=True.

        phenopacket-tools JAR is mocked as unavailable.
        """
        from phenopackets_input.schema_validator import validate_phenopacket

        pp_file = tmp_path / "valid.json"
        pp_file.write_text(json.dumps(valid_phenopacket_v2), encoding="utf-8")

        # Mock phenopacket-tools JAR as not available (common in CI)
        with patch("pathlib.Path.exists", side_effect=lambda: False):
            result = validate_phenopacket(pp_file)

        assert result.valid is True, (
            f"Valid phenopacket file should pass. Errors: {result.errors}"
        )
