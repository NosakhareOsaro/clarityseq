# conf/

Nextflow configuration files.

| File | Description |
|------|-------------|
| `base.config` | Default resource limits per process label |
| `dragen_gatk.config` | DRAGEN-GATK mode (BQSR=false); critical comment block |
| `test.config` | chr22 test profile; max 8 CPUs; no pangenome |
| `aws.config` | AWS Batch with Spot fleet; ECR images; S3 staging |
| `hpc.config` | SLURM executor; Singularity containers |
| `resources.config` | Per-process resource overrides |

## Critical: DRAGEN-GATK and BQSR

**`dragen_gatk.config` sets `run_bqsr = false`.** This is intentional and mandatory.

When using DRAGMAP + HaplotypeCaller in DRAGEN-GATK mode, BQSR must NOT be run because:
1. HaplotypeCaller's **BQD (Base Quality Dropoff)** genotyping model replaces BQSR by modelling systematic errors internally
2. Applying BQSR to DRAGMAP-aligned reads produces quality scores incompatible with the BQD model, **reducing accuracy**

Reference: Broad Institute DRAGEN-GATK documentation
https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939

If you need BQSR (e.g., BWA-MEM2 fallback), use `-profile bwa_mem2`.
