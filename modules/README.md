# modules/

Nextflow DSL2 process modules following nf-core module style.

Every module has:
- An 8-line header comment block (see §0.3 of PROJECT_GUIDE.MD)
- Pinned container versions (never `latest`)
- Inline comments on every non-default parameter
- Version note explaining why this specific version was chosen

## Aligner modules

| Module | Container | Role |
|--------|-----------|------|
| `dragmap/` | `nfcore/dragmap:1.3.0` | **PRIMARY aligner** (DRAGEN-GATK Best Practices); BQSR NOT run |
| `bwa_mem2/` | `nfcore/bwa-mem2:2.2.1` | **FALLBACK only** (`-profile bwa_mem2`); BQSR IS run |

## Variant calling modules

| Module | Container | Notes |
|--------|-----------|-------|
| `gatk4/haplotypecaller/` | `broadinstitute/gatk:4.6.0.0` | `--dragen-mode true`; BQD model |
| `deepvariant/` | `google/deepvariant:1.8.0` | CNN pileup; parallel with GATK4 |
| `deeptrio/` | `google/deepvariant:1.8.0` | Trio samples; +15% de novo sensitivity |
| `gatk4/mutect2_mito/` | `broadinstitute/gatk:4.6.0.0` | `--mitochondria-mode` |
| `expansionhunter/` | `clinicalgenomics/expansionhunter:5.0.0` | 60 loci; TRGT NOT used |

## Annotation modules

| Module | Notes |
|--------|-------|
| `vep/` | VEP v111; MANE Select pick order |
| `alphamissense/` | ClinGen SVI 2024: ≥0.564→PP3, ≤0.340→BP4 |

## QC modules

| Module | Notes |
|--------|-------|
| `mosdepth/` | Coverage QC; 30× minimum gate |
| `somalier/` | Ancestry inference → per-ancestry VQSR |
| `fastp/` | Adapter trimming; poly-G for NovaSeq |
