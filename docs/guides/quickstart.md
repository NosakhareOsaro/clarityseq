# Quick-start Guide

Run your first GenomeForge WGS analysis in 5 commands.

## Prerequisites

- Docker ≥ 25.x (for `-profile local`)
- Nextflow ≥ 24.x
- 64 GB RAM; 16+ CPUs recommended
- Reference data (see `docs/guides/data_setup.md`)

## 5-command quick-start

```bash
# 1. Clone repository
git clone https://github.com/genome-forge/genome-forge.git && cd genome-forge

# 2. Configure environment
cp .env.example .env
# Edit .env: set DRAGMAP_REFERENCE, GNOMAD_VCF, ALPHAMISSENSE_SCORES_PATH

# 3. Run test (chr22, GIAB HG001, ~30 min)
nextflow run pipelines/wgs_grch38.nf -profile test,docker

# 4. Verify: SNP sensitivity ≥ 99.0%
cat results/test/benchmark/hap_py_summary.txt

# 5. Run with your samples
nextflow run pipelines/wgs_grch38.nf -profile local \
  --input my_samples.csv \
  --dragmap_reference /data/dragmap \
  --gnomad_vcf /data/gnomad/gnomad.genomes.v4.1.sites.chr*.vcf.bgz \
  --alphamissense_tsv /data/AlphaMissense_hg38.tsv.gz \
  --vep_cache_dir /data/vep_cache \
  --outdir results/
```

## Sample sheet format

```csv
sample,fastq_1,fastq_2,sex,affected,ped_file
Patient1,/data/Patient1_R1.fastq.gz,/data/Patient1_R2.fastq.gz,female,true,
Trio_Proband,/data/Proband_R1.fastq.gz,/data/Proband_R2.fastq.gz,male,true,family.ped
```

- `sex`: `male` or `female` (required for X-linked inheritance filter)
- `affected`: `true` or `false`
- `ped_file`: optional; provide for trio analysis (enables DeepTrio)

## Outputs

```
results/
└── Sample1/
    ├── alignment/             BAM + BAI (DRAGMAP-aligned)
    ├── variants/
    │   ├── gatk4/             GATK4 GVCF and genotyped VCF
    │   ├── deepvariant/       DeepVariant VCF
    │   └── ensemble/          Ensemble-merged VCF (INTERSECTION mode)
    ├── annotation/            VEP v111 annotated VCF
    ├── classification/        BayesACMG JSON + TSV
    ├── report/
    │   ├── Sample1_clinical_report.html   NHS GMS-style report
    │   └── Sample1_clinical_report.pdf
    ├── mito/                  Mitochondrial variants + haplogroup
    ├── expansions/            ExpansionHunter repeat expansion calls
    ├── pgx/                   CYP2D6 star alleles + CPIC dosing
    └── qc/                    fastp + mosdepth + multiqc reports
```

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dragmap_reference` | null | DRAGMAP hash table directory (PRIMARY) |
| `--gnomad_vcf` | null | gnomAD v4.1 sites VCF (use v4.1 NOT v4.0) |
| `--alphamissense_tsv` | null | AlphaMissense scores TSV.gz (tabix-indexed) |
| `--vep_cache_dir` | null | VEP v111 cache directory |
| `--run_deepvariant` | true | Run DeepVariant v1.8.0 parallel arm |
| `--run_mito` | true | GATK Mutect2 mito + Haplogrep3 |
| `--run_expansions` | true | ExpansionHunter v5.0 (60 loci) |
| `--run_pgx` | true | Cyrius CYP2D6 pharmacogenomics |
| `--run_pangenome` | false | vg giraffe HPRC arm (requires 256 GB RAM) |

## Profiles

| Profile | Description |
|---------|-------------|
| `test` | chr22 only; GIAB HG001; max 8 CPUs; no pangenome |
| `local` | Docker; DRAGMAP primary; default resources |
| `aws` | AWS Batch + Spot fleet; ECR images; S3 staging |
| `hpc` | Singularity; SLURM executor |
| `bwa_mem2` | **FALLBACK**: BWA-MEM2 + BQSR (when DRAGMAP hash unavailable) |
| `pangenome` | vg giraffe arm; requires 32+ CPUs, 256 GB RAM |
