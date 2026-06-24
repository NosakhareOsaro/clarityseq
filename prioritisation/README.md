# prioritisation/

Phenotype-driven variant prioritisation using HPO terms and Exomiser 14.

## Ranking strategy

Variants are ranked by a composite score combining:

1. **ACMG classification** (BayesACMG posterior probability)
2. **HPO phenotype match** (Exomiser 14 phenotype score via Phenopackets v2)
3. **Inheritance model** (AD/AR/XL/de_novo filter applied first)
4. **Panel membership** (PanelApp gene panel status)

## Exomiser 14

Exomiser is run via its REST API accepting Phenopackets v2 format.
It scores variants by combining phenotype similarity (HPO ontology traversal) with gene-level constraint scores.

Reference: Smedley et al. 2015 PMID:26562621

## Components

| File | Description |
|------|-------------|
| `hpo_scorer.py` | HPO term scoring against variant consequences |
| `inheritance_filter.py` | Filter by inheritance mode (AD/AR/XL/de_novo) |
| `exomiser_client.py` | Exomiser 14 REST client (Phenopackets v2 input) |
| `ranking.py` | Final composite variant ranking |
