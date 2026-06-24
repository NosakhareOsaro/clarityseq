"""
phenopackets_input — Phenopackets v2 ingestion, validation, and conversion.

Phenopackets v2 (Jacobsen et al. 2022 Nature Biotechnology) provides a
GA4GH-standardised representation of clinical and genomic data.

Modules
-------
- schema_validator: Validate phenopacket JSON against the v2 schema using
  both the Python SDK and phenopacket-tools for strict compliance.
- converter: Convert HPO terms and clinical data to variant filter parameters.
- fhir_mapper: Bidirectional mapping between Phenopackets v2 and FHIR R4.

References
----------
Jacobsen JOB, et al. "The GA4GH Phenopacket Schema: A computable representation
of clinical data for precision medicine." Nature Biotechnology. 2022;40:817–820.
PMID:35705715. DOI:10.1038/s41587-022-01357-4

Haendel MA, et al. "How many rare diseases are there?" Nature Reviews Drug
Discovery. 2020;19:77–78. PMID:32020066.
"""

from phenopackets_input.schema_validator import ValidationResult, validate_phenopacket

__all__ = [
    "validate_phenopacket",
    "ValidationResult",
]
