# pipelines/

Nextflow DSL2 pipeline entry points for GenomeForge.

| File | Description |
|------|-------------|
| `wgs_grch38.nf` | Primary WGS pipeline: DRAGMAP → DRAGEN-GATK → DeepVariant ensemble → VEP v111 → BayesACMG |
| `wgs_t2t.nf` | T2T-CHM13 v2.0 arm with CrossMap liftover to GRCh38 |
| `wgs_pangenome.nf` | HPRC pangenome arm (vg giraffe v1.1 graph) |
| `mito.nf` | Mitochondrial sub-workflow (Mutect2 mito mode + Haplogrep3; ACGS 2024 §6) |
| `expansions.nf` | Repeat expansion sub-workflow (ExpansionHunter v5.0, 60 loci) |

## Usage

```bash
# Standard GRCh38 (primary)
nextflow run pipelines/wgs_grch38.nf -profile local --sample_sheet samples.csv

# Test profile (chr22, GIAB HG001)
nextflow run pipelines/wgs_grch38.nf -profile test

# BWA-MEM2 fallback (when DRAGMAP hash table unavailable)
nextflow run pipelines/wgs_grch38.nf -profile bwa_mem2

# Pangenome arm (nightly CI only; requires 32+ CPUs, 256 GB RAM)
nextflow run pipelines/wgs_pangenome.nf -profile pangenome
```

## Guidelines

- DRAGMAP (primary) + DRAGEN-GATK mode: BQSR is **not** run (see conf/dragen_gatk.config)
- AlphaMissense used as PRIMARY PP3/BP4 predictor (ClinGen SVI 2024)
- PM2 applied at SUPPORTING weight (ClinGen SVI 2024; see bayesacmg/src/bayesacmg/rules/pathogenic.py)
