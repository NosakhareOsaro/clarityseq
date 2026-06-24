# singularity/

Singularity/Apptainer container definition for HPC deployment.

## Usage on HPC

```bash
# Build from definition file (requires root or --fakeroot)
singularity build genome-forge.sif singularity/genome-forge.def

# Run pipeline with Singularity
nextflow run pipelines/wgs_grch38.nf -profile hpc --singularity_image genome-forge.sif

# Or with Apptainer (rootless)
apptainer build genome-forge.sif singularity/genome-forge.def
```

## Note on DeepVariant

DeepVariant v1.8.0 uses a separate GPU-optimised Singularity image pulled at runtime. The base image contains GATK, DRAGMAP, samtools, and other tools but NOT DeepVariant (too large for a single image).
