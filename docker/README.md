# docker/

Multi-stage Docker images for GenomeForge components. All based on Ubuntu 24.04 LTS.

## Images

| Dockerfile | Purpose | Key tools |
|------------|---------|-----------|
| `Dockerfile.pipeline` | Nextflow pipeline runner | GATK 4.6.0.0, DRAGMAP 1.3.0, samtools, fastp, Python 3.12 |
| `Dockerfile.beacon` | GA4GH Beacon API | Python 3.12, FastAPI, uvicorn |
| `Dockerfile.daemon` | ClinVar reclassification daemon | Python 3.12, Celery, Redis client |

## Base image rationale

**Ubuntu 24.04** (Noble Numbat, released April 2024):
- LTS support until 2029 (vs 22.04 → 2027)
- GLIBC 2.39 (required for GATK 4.6.0.0)
- OpenSSL 3.x (security requirement)
- System Python 3.12 compatible

## Build

```bash
# Build all images
make docker-build

# Build individually
docker build -f docker/Dockerfile.pipeline -t genomeforge/pipeline:0.1.0 .
docker build -f docker/Dockerfile.beacon -t genomeforge/beacon:0.1.0 .
docker build -f docker/Dockerfile.daemon -t genomeforge/daemon:0.1.0 .
```

## DeepVariant note

DeepVariant v1.8.0 is NOT installed in Dockerfile.pipeline — it is pulled as a
separate Docker image (`google/deepvariant:1.8.0`) in the Nextflow module.
Reason: DeepVariant image is ~10 GB; embedding it would make pipeline image unmanageable.
