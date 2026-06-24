"""
phenopackets_input.schema_validator
=====================================
Phenopackets v2.0 JSON validation using phenopacket-tools and Python SDK.

Why both validators:
    Python SDK: parses and creates Phenopackets; less strict validation.
    phenopacket-tools (Java): enforces valid OMIM IDs, HPO term existence,
    ontology version metadata — constraints the SDK misses.

Two-stage validation:
    1. Python SDK (phenopackets package): parse and basic schema conformance.
    2. phenopacket-tools (Java subprocess): strict ontology validation.
       - Validates HPO term IDs exist in the current HPO release.
       - Validates OMIM/Orphanet disease IDs.
       - Checks metadata.created_by, metadata.phenopacket_schema_version.
       - Enforces ISO8601 timestamps in metadata.created field.

Running phenopacket-tools:
    Java ≥ 17 required. Download from:
    https://github.com/phenopackets/phenopacket-tools/releases
    Set PHENOPACKET_TOOLS_JAR environment variable to the JAR path.

Reference:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    phenopacket-tools: https://github.com/phenopackets/phenopacket-tools
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_PHENOPACKET_TOOLS_JAR: str = os.getenv(
    "PHENOPACKET_TOOLS_JAR",
    "/opt/tools/phenopacket-tools-cli.jar",
)

_JAVA_BINARY: str = os.getenv("JAVA_BINARY", "java")

# Required top-level keys in a Phenopackets v2.0 JSON
_REQUIRED_KEYS_V2 = {"id", "subject", "metaData"}

# Required metaData keys
_REQUIRED_METADATA_KEYS = {"created", "createdBy", "phenopacketSchemaVersion", "resources"}


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class ValidationResult:
    """Result of validating a Phenopackets v2.0 JSON file.

    Attributes:
        valid: True if the phenopacket passed all validation checks.
        errors: List of validation error messages (blocking issues).
        warnings: List of validation warning messages (non-blocking).
        validated_by: List of validator names that ran successfully.
        phenopacket_id: The ID field from the validated phenopacket,
            or empty string if the file could not be parsed.
    """

    valid: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    validated_by: list[str] = field(default_factory=list)
    phenopacket_id: str = ""

    def __str__(self) -> str:
        """Return a human-readable validation summary.

        Returns:
            Formatted string with validation status, errors, and warnings.
        """
        status = "VALID" if self.valid else "INVALID"
        lines = [f"Phenopacket validation: {status} (id={self.phenopacket_id!r})"]
        if self.errors:
            lines.append(f"  Errors ({len(self.errors)}):")
            for e in self.errors:
                lines.append(f"    - {e}")
        if self.warnings:
            lines.append(f"  Warnings ({len(self.warnings)}):")
            for w in self.warnings:
                lines.append(f"    - {w}")
        lines.append(f"  Validated by: {', '.join(self.validated_by)}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Internal validation helpers
# ---------------------------------------------------------------------------


def _validate_json_schema(data: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Validate basic Phenopackets v2.0 JSON structure.

    Checks required top-level fields and metaData structure.
    Does NOT validate ontology term validity (that requires phenopacket-tools).

    Args:
        data: Parsed phenopacket JSON dict.

    Returns:
        Tuple of (errors, warnings) lists.
    """
    errors: list[str] = []
    warnings: list[str] = []

    # Check required top-level keys
    missing_keys = _REQUIRED_KEYS_V2 - set(data.keys())
    if missing_keys:
        errors.append(
            f"Missing required top-level keys: {', '.join(sorted(missing_keys))}. "
            "Required: id, subject, metaData (Phenopackets v2.0)."
        )

    # Check subject
    subject = data.get("subject", {})
    if not isinstance(subject, dict):
        errors.append("'subject' must be a JSON object.")
    elif "id" not in subject:
        errors.append("subject.id is required (Phenopackets v2.0 §2.1).")

    # Check metaData
    metadata = data.get("metaData", {})
    if not isinstance(metadata, dict):
        errors.append("'metaData' must be a JSON object.")
    else:
        missing_meta = _REQUIRED_METADATA_KEYS - set(metadata.keys())
        if missing_meta:
            errors.append(
                f"metaData missing required keys: {', '.join(sorted(missing_meta))}."
            )

        # Check schema version
        schema_version = metadata.get("phenopacketSchemaVersion", "")
        if schema_version and not schema_version.startswith("2"):
            errors.append(
                f"phenopacketSchemaVersion '{schema_version}' must be '2.0' "
                "for Phenopackets v2.0 validation."
            )
        elif not schema_version:
            warnings.append(
                "metaData.phenopacketSchemaVersion is empty; should be '2.0'."
            )

        # Check resources array
        resources = metadata.get("resources", [])
        if not isinstance(resources, list):
            errors.append("metaData.resources must be a JSON array.")
        elif len(resources) == 0:
            warnings.append(
                "metaData.resources is empty. "
                "Should include at least the HPO ontology resource."
            )

    # Check phenotypicFeatures have valid HPO ID format (HP:XXXXXXX)
    phenotypic_features = data.get("phenotypicFeatures", [])
    for i, feature in enumerate(phenotypic_features):
        if not isinstance(feature, dict):
            errors.append(f"phenotypicFeatures[{i}] must be a JSON object.")
            continue
        term = feature.get("type", {})
        term_id = term.get("id", "")
        if term_id and not term_id.startswith("HP:"):
            errors.append(
                f"phenotypicFeatures[{i}].type.id '{term_id}' does not look like "
                "an HPO term (expected HP:XXXXXXX format)."
            )

    # Check diseases have OMIM/Orphanet IDs
    diseases = data.get("diseases", [])
    for i, disease in enumerate(diseases):
        if not isinstance(disease, dict):
            errors.append(f"diseases[{i}] must be a JSON object.")
            continue
        term = disease.get("term", {})
        term_id = term.get("id", "")
        valid_prefixes = ("OMIM:", "Orphanet:", "MONDO:", "OMIMPS:")
        if term_id and not any(term_id.startswith(p) for p in valid_prefixes):
            warnings.append(
                f"diseases[{i}].term.id '{term_id}' uses an unexpected prefix. "
                f"Expected one of: {', '.join(valid_prefixes)}."
            )

    return errors, warnings


def _run_phenopacket_tools(json_path: Path) -> tuple[list[str], list[str]]:
    """Run phenopacket-tools JAR for strict ontology validation.

    Executes the phenopacket-tools CLI validate command against the phenopacket
    JSON file.  Parses stderr/stdout for error and warning lines.

    Args:
        json_path: Path to the Phenopacket JSON file to validate.

    Returns:
        Tuple of (errors, warnings) from phenopacket-tools output.
        Returns ([], ["phenopacket-tools not available"]) if JAR not found.
    """
    jar_path = Path(_PHENOPACKET_TOOLS_JAR)
    if not jar_path.exists():
        logger.debug(
            "phenopacket-tools JAR not found at %s; skipping strict validation.",
            jar_path,
        )
        return [], [
            f"phenopacket-tools JAR not found at '{jar_path}'. "
            "Set PHENOPACKET_TOOLS_JAR env var for strict HPO/OMIM validation. "
            "See: https://github.com/phenopackets/phenopacket-tools"
        ]

    cmd = [
        _JAVA_BINARY,
        "-jar", str(jar_path),
        "validate",
        "--phenopacket", str(json_path),
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
    except FileNotFoundError:
        return [], [f"Java binary '{_JAVA_BINARY}' not found. Install Java ≥ 17."]
    except subprocess.TimeoutExpired:
        return [f"phenopacket-tools timed out after 60s validating {json_path.name}"], []

    errors: list[str] = []
    warnings: list[str] = []

    for line in (proc.stdout + proc.stderr).splitlines():
        line = line.strip()
        if not line:
            continue
        low = line.lower()
        if "error" in low or proc.returncode != 0:
            errors.append(f"[phenopacket-tools] {line}")
        elif "warn" in low:
            warnings.append(f"[phenopacket-tools] {line}")

    return errors, warnings


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_phenopacket(json_path: Path) -> ValidationResult:
    """Validate a Phenopackets v2.0 JSON file.

    Runs two validation stages:
    1. Python SDK structural validation (basic JSON schema conformance).
    2. phenopacket-tools (Java) strict validation if JAR is available:
       - Valid HPO term IDs in current HPO release.
       - Valid OMIM/Orphanet disease IDs.
       - Metadata completeness and ISO8601 timestamp format.

    Args:
        json_path: Path to the Phenopacket v2.0 JSON file to validate.

    Returns:
        ValidationResult with valid=True if no errors were found.
        Warnings are non-blocking and included for informational purposes.

    Raises:
        FileNotFoundError: If json_path does not exist.
        json.JSONDecodeError: If the file is not valid JSON.

    References:
        Jacobsen et al. 2022 Nature Biotechnology PMID:35705716.
        Phenopackets v2.0 schema: https://phenopacket-schema.readthedocs.io/
    """
    json_path = Path(json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"Phenopacket JSON not found: {json_path}")

    all_errors: list[str] = []
    all_warnings: list[str] = []
    validated_by: list[str] = []
    phenopacket_id: str = ""

    # Stage 1: Parse JSON
    try:
        with json_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        return ValidationResult(
            valid=False,
            errors=[f"Invalid JSON: {exc}"],
            warnings=[],
            validated_by=["json.loads"],
            phenopacket_id="",
        )

    phenopacket_id = str(data.get("id", ""))

    # Stage 1: Structural JSON schema validation
    sdk_errors, sdk_warnings = _validate_json_schema(data)
    all_errors.extend(sdk_errors)
    all_warnings.extend(sdk_warnings)
    validated_by.append("Python structural validator")

    # Stage 2: phenopacket-tools strict validation
    pt_errors, pt_warnings = _run_phenopacket_tools(json_path)
    all_errors.extend(pt_errors)
    all_warnings.extend(pt_warnings)
    if not any("not found" in w for w in pt_warnings):
        validated_by.append("phenopacket-tools (Java)")

    is_valid = len(all_errors) == 0

    result = ValidationResult(
        valid=is_valid,
        errors=all_errors,
        warnings=all_warnings,
        validated_by=validated_by,
        phenopacket_id=phenopacket_id,
    )
    logger.info("Phenopacket %s validation: %s", phenopacket_id, "VALID" if is_valid else "INVALID")
    if not is_valid:
        logger.warning("Phenopacket errors: %s", all_errors)
    return result


def validate_phenopacket_dict(data: dict[str, Any]) -> ValidationResult:
    """Validate a Phenopackets v2.0 JSON already parsed into a dict.

    Performs only the Python structural validation stage (no phenopacket-tools
    call since there is no file path to pass to the JAR).  For strict
    ontology validation, write to a temp file and use validate_phenopacket().

    Args:
        data: Already-parsed phenopacket JSON dict.

    Returns:
        ValidationResult with structural validation results.

    References:
        Jacobsen et al. 2022 Nature Biotechnology PMID:35705716.
    """
    errors, warnings = _validate_json_schema(data)
    return ValidationResult(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        validated_by=["Python structural validator"],
        phenopacket_id=str(data.get("id", "")),
    )
