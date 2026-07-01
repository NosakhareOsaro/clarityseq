# multi_ancestry/

Multi-ancestry VQSR model selection via somalier ancestry inference.

## Why per-ancestry VQSR?

Variant quality score recalibration (VQSR) uses training variants from population databases. Using the wrong ancestry training set inflates false positive rates for underrepresented populations. ClaritySeq uses somalier to infer ancestry, then selects the appropriate gnomAD v4.1 ancestry-stratified subset for VQSR training.

## Ancestry-to-VQSR mapping

| Inferred ancestry | gnomAD v4.1 subset | Coverage |
|-------------------|--------------------|---------|
| AFR | gnomAD v4.1 AFR | 43,538 genomes |
| AMR | gnomAD v4.1 AMR | 6,972 genomes |
| EAS | gnomAD v4.1 EAS | 5,586 genomes |
| EUR | gnomAD v4.1 EUR | 33,082 genomes |
| SAS | gnomAD v4.1 SAS | 5,024 genomes |
| Admixed/unknown | gnomAD v4.1 combined | All populations |

## Components

| File | Description |
|------|-------------|
| `somalier_runner.py` | Run somalier relate + ancestry (uses 1000G PCA sites) |
| `ancestry_assigner.py` | Assign population label from somalier PCA output |
| `vqsr_selector.py` | Select gnomAD v4.1 ancestry-stratified training resources |

## Reference

somalier: Pedersen et al. 2020 Genome Biology PMID:32664994
gnomAD v4.1 ancestry groups: Chen et al. 2024 bioRxiv
