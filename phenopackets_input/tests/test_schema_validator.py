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
import subprocess
import tempfile
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

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
            "createdBy": "ClaritySeq",
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

        # phenopacket-tools JAR is absent in CI — validation still passes without it
        result = validate_phenopacket(pp_file)

        assert result.valid is True, (
            f"Valid phenopacket file should pass. Errors: {result.errors}"
        )


# ---------------------------------------------------------------------------
# ValidationResult.__str__ tests
# ---------------------------------------------------------------------------


class TestValidationResultStr:
    """Tests for ValidationResult.__str__ error/warning formatting."""

    def test_str_includes_error_lines(self) -> None:
        """__str__ lists each error message when errors are present."""
        from phenopackets_input.schema_validator import ValidationResult

        result = ValidationResult(
            valid=False,
            errors=["error one", "error two"],
            warnings=[],
            validated_by=["Python structural validator"],
            phenopacket_id="pp-err",
        )
        s = str(result)
        assert "INVALID" in s
        assert "Errors (2):" in s
        assert "error one" in s
        assert "error two" in s

    def test_str_includes_warning_lines(self) -> None:
        """__str__ lists each warning message when warnings are present."""
        from phenopackets_input.schema_validator import ValidationResult

        result = ValidationResult(
            valid=True,
            errors=[],
            warnings=["warn one"],
            validated_by=["Python structural validator"],
            phenopacket_id="pp-warn",
        )
        s = str(result)
        assert "Warnings (1):" in s
        assert "warn one" in s


# ---------------------------------------------------------------------------
# _validate_json_schema edge case tests
# ---------------------------------------------------------------------------


class TestValidateJsonSchemaEdgeCases:
    """Tests for _validate_json_schema() structural edge cases."""

    def test_subject_not_dict_returns_error(self) -> None:
        """'subject' as a non-dict value returns an error."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": "not-a-dict",
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": [],
            },
        }
        errors, warnings = _validate_json_schema(data)
        assert any("subject" in e and "JSON object" in e for e in errors)

    def test_metadata_not_dict_returns_error(self) -> None:
        """'metaData' as a non-dict value returns an error."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": {"id": "p1"},
            "metaData": "not-a-dict",
        }
        errors, warnings = _validate_json_schema(data)
        assert any("metaData" in e and "JSON object" in e for e in errors)

    def test_resources_not_list_returns_error(self) -> None:
        """metaData.resources as a non-list value returns an error."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": {"id": "p1"},
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": "hp",
            },
        }
        errors, warnings = _validate_json_schema(data)
        assert any("resources must be a JSON array" in e for e in errors)

    def test_phenotypic_feature_not_dict_returns_error(self) -> None:
        """A non-dict entry in phenotypicFeatures returns an error and is skipped."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": {"id": "p1"},
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": [{"id": "hp"}],
            },
            "phenotypicFeatures": ["not-a-dict", {"type": {"id": "HP:0001250"}}],
        }
        errors, warnings = _validate_json_schema(data)
        assert any(
            "phenotypicFeatures[0]" in e and "JSON object" in e for e in errors
        )

    def test_disease_not_dict_returns_error(self) -> None:
        """A non-dict entry in diseases returns an error and is skipped."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": {"id": "p1"},
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": [{"id": "hp"}],
            },
            "diseases": ["not-a-dict"],
        }
        errors, warnings = _validate_json_schema(data)
        assert any("diseases[0]" in e and "JSON object" in e for e in errors)

    def test_disease_unexpected_prefix_returns_warning(self) -> None:
        """A disease term ID with an unexpected prefix returns a warning, not an error."""
        from phenopackets_input.schema_validator import _validate_json_schema

        data = {
            "id": "pp-1",
            "subject": {"id": "p1"},
            "metaData": {
                "created": "2024-01-01T00:00:00Z",
                "createdBy": "test",
                "phenopacketSchemaVersion": "2.0",
                "resources": [{"id": "hp"}],
            },
            "diseases": [{"term": {"id": "ICD10:E11.9", "label": "Diabetes"}}],
        }
        errors, warnings = _validate_json_schema(data)
        assert not any("diseases[0]" in e for e in errors)
        assert any("unexpected prefix" in w for w in warnings)


# ---------------------------------------------------------------------------
# _run_phenopacket_tools tests (mocked subprocess)
# ---------------------------------------------------------------------------


class TestRunPhenopacketTools:
    """Tests for _run_phenopacket_tools() with a mocked Java subprocess."""

    def test_jar_missing_returns_warning(self, tmp_path: Path) -> None:
        """When the JAR path does not exist, returns a 'not found' warning."""
        from phenopackets_input import schema_validator

        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        with patch.object(
            schema_validator, "_PHENOPACKET_TOOLS_JAR", str(tmp_path / "missing.jar")
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert errors == []
        assert any("not found" in w for w in warnings)

    def test_successful_run_no_errors_or_warnings(self, tmp_path: Path) -> None:
        """Clean stdout/stderr with returncode 0 yields no errors or warnings."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = ""
        mock_proc.stderr = ""

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert errors == []
        assert warnings == []

    def test_error_line_in_output_is_captured(self, tmp_path: Path) -> None:
        """A stdout line containing 'error' is captured as a validation error."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "ERROR: invalid HPO term HP:9999999\n"
        mock_proc.stderr = ""

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert any("invalid HPO term" in e for e in errors)

    def test_nonzero_returncode_marks_all_lines_as_errors(self, tmp_path: Path) -> None:
        """Non-zero returncode marks every non-empty output line as an error."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stdout = "validation failed for subject\n"
        mock_proc.stderr = ""

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert any("validation failed for subject" in e for e in errors)

    def test_warning_line_in_output_is_captured(self, tmp_path: Path) -> None:
        """A stdout line containing 'warn' (returncode 0) is captured as a warning."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "WARN: resource version outdated\n"
        mock_proc.stderr = ""

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert any("resource version outdated" in w for w in warnings)

    def test_blank_lines_are_skipped(self, tmp_path: Path) -> None:
        """Blank lines in stdout/stderr do not produce error/warning entries."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "\n   \n"
        mock_proc.stderr = "\n"

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert errors == []
        assert warnings == []

    def test_timeout_returns_error(self, tmp_path: Path) -> None:
        """subprocess.TimeoutExpired is converted into a timeout error message."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch(
                "subprocess.run",
                side_effect=subprocess.TimeoutExpired(cmd="java", timeout=60),
            ),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert any("timed out" in e for e in errors)
        assert warnings == []

    def test_java_binary_not_found_returns_warning(self, tmp_path: Path) -> None:
        """FileNotFoundError (missing Java binary) is converted into a warning."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        json_path = tmp_path / "pp.json"
        json_path.write_text("{}")

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", side_effect=FileNotFoundError()),
        ):
            errors, warnings = schema_validator._run_phenopacket_tools(json_path)

        assert errors == []
        assert any("Java binary" in w and "not found" in w for w in warnings)


# ---------------------------------------------------------------------------
# validate_phenopacket with phenopacket-tools available (full integration)
# ---------------------------------------------------------------------------


class TestValidatePhenopacketWithToolsAvailable:
    """Tests for validate_phenopacket() when the phenopacket-tools JAR is present."""

    def test_validated_by_includes_phenopacket_tools_on_success(
        self, tmp_path: Path, valid_phenopacket_v2: dict[str, Any]
    ) -> None:
        """validated_by includes 'phenopacket-tools (Java)' when the JAR runs cleanly."""
        from phenopackets_input import schema_validator

        jar_path = tmp_path / "phenopacket-tools.jar"
        jar_path.write_text("fake jar")
        pp_file = tmp_path / "valid.json"
        pp_file.write_text(json.dumps(valid_phenopacket_v2), encoding="utf-8")

        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = ""
        mock_proc.stderr = ""

        with (
            patch.object(schema_validator, "_PHENOPACKET_TOOLS_JAR", str(jar_path)),
            patch("subprocess.run", return_value=mock_proc),
        ):
            result = schema_validator.validate_phenopacket(pp_file)

        assert result.valid is True
        assert "phenopacket-tools (Java)" in result.validated_by
        assert "Python structural validator" in result.validated_by

    def test_logs_warning_when_invalid(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """validate_phenopacket logs a warning listing errors when invalid."""
        from phenopackets_input.schema_validator import validate_phenopacket

        pp_file = tmp_path / "invalid.json"
        pp_file.write_text(
            json.dumps(
                {
                    "id": "pp-invalid",
                    "subject": {},  # missing subject.id -> structural error
                    "metaData": {
                        "created": "2024-01-01T00:00:00Z",
                        "createdBy": "test",
                        "phenopacketSchemaVersion": "2.0",
                        "resources": [],
                    },
                }
            ),
            encoding="utf-8",
        )

        with caplog.at_level("WARNING", logger="phenopackets_input.schema_validator"):
            result = validate_phenopacket(pp_file)

        assert result.valid is False
        assert any("Phenopacket errors" in rec.message for rec in caplog.records)
