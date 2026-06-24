# phenopackets_input/

Phenopackets v2.0 clinical input handling for GenomeForge.

## What are Phenopackets?

Phenopackets v2.0 is the GA4GH standard for sharing disease and phenotype data in a computable format. GenomeForge uses them to capture clinical context (HPO terms, family history, pedigree) and pass it to Exomiser 14 for phenotype-driven variant prioritisation.

Reference: Jacobsen et al. 2022 Nature Biotechnology PMID:35705716

## Components

| File | Description |
|------|-------------|
| `schema_validator.py` | Validates Phenopackets v2 JSON (phenopacket-tools + Python SDK) |
| `converter.py` | Converts HPO terms → variant filter parameters |
| `fhir_mapper.py` | Bidirectional Phenopackets v2 ↔ FHIR R4 mapper |
| `examples/` | Example Phenopackets JSON files |

## Why two validators?

- **Python SDK**: Parses and creates Phenopackets; less strict validation
- **phenopacket-tools** (Java): Enforces additional constraints — valid OMIM IDs, HPO term existence, ontology version metadata. The SDK misses these.

Both are run in sequence; both must pass before a phenopacket is accepted.

## Usage

```python
from phenopackets_input.schema_validator import validate_phenopacket
result = validate_phenopacket("patient.json")
if not result.valid:
    print(result.errors)
```
