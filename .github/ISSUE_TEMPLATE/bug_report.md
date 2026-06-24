---
name: Bug report
about: Report a reproducible bug in GenomeForge
title: "[BUG] "
labels: bug
assignees: ""
---

## Bug description

<!-- A clear, concise description of the bug. What is broken? -->

## Severity

- [ ] P1 — Clinical safety impact (incorrect variant classification or missed pathogenic call)
- [ ] P2 — Pipeline failure (run crashes or produces no output)
- [ ] P3 — Incorrect output (wrong result but pipeline completes)
- [ ] P4 — Minor (cosmetic, documentation, or performance issue)

## Component affected

<!-- Check all that apply -->
- [ ] Pipeline (Nextflow) — module name: `modules/____________/main.nf`
- [ ] BayesACMG variant classifier (`bayesacmg/`)
- [ ] Ensemble caller (`modules/ensemble/ensemble.py`)
- [ ] Annotation stack (`annotation/`, VEP, AlphaMissense)
- [ ] Clinical reporting (`reporting/`)
- [ ] ClinVar reclassification daemon (`reclassification/`)
- [ ] GA4GH Beacon API (`beacon_api/`)
- [ ] Pharmacogenomics — Cyrius CYP2D6 (`modules/cyrius/`)
- [ ] STR expansion calling — ExpansionHunter (`modules/expansionhunter/`)
- [ ] Mitochondrial calling (`modules/gatk4/mutect2_mito/`, `modules/haplogrep3/`)
- [ ] Mosdepth coverage QC (`modules/mosdepth/`)
- [ ] Somalier ancestry/relatedness (`modules/somalier/`)
- [ ] CI / GitHub Actions (`.github/workflows/`)
- [ ] Terraform / infrastructure (`terraform/`)
- [ ] Other: _______________

## Steps to reproduce

<!-- Provide the exact commands needed to reproduce the issue -->

```bash
# 1. Environment setup
nextflow -version
python --version

# 2. Exact command that triggers the bug
nextflow run pipelines/wgs_grch38.nf \
  --input samplesheet.csv \
  --outdir results/ \
  ...

# OR for Python bugs:
python -m pytest bayesacmg/tests/test_<filename>.py -v
```

**Samplesheet (if applicable):**
```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
```

## Expected behaviour

<!-- What should happen according to the documentation, ACGS guidelines, or prior behaviour? -->

## Actual behaviour

<!-- What actually happens? Paste the complete error output. -->

```
ERROR: paste full error output here, including stack trace if Python
```

**Relevant log files:**
```
# .nextflow.log excerpt, or pytest output
```

## ACMG/Clinical impact assessment

<!-- REQUIRED for P1 bugs. Complete for any bug affecting variant classification. -->

- Variant affected (if known): `CHROM:POS:REF:ALT`
- Expected ACMG classification: `P / LP / VUS / LB / B`
- Actual ACMG classification produced: `P / LP / VUS / LB / B`
- ACMG rule affected (if known): `PVS1 / PM2 / PP3 / BA1 / ...`
- PM2 weight returned: `Supporting / Moderate`
  (MUST be Supporting per ClinGen SVI 2024; flag immediately if Moderate is returned)
- Does this affect haplogroup-dependent mito classification? `Y / N`
- Does this affect CYP2D6 metaboliser phenotype assignment? `Y / N`

## Environment

```
GenomeForge version: v0.x.x  (from `git describe --tags`)
Python version:      3.12.x
Nextflow version:    24.x.x  (from `nextflow -version`)
Docker version:      xx.x.x  (from `docker --version`)
OS / Platform:       Ubuntu 24.04 / AWS EC2 r6i.8xlarge / macOS 15
Container runtime:   Docker / Singularity / Apptainer
```

## Additional context

<!-- Add any other context, screenshots, or relevant VCF/log snippets -->

<!-- For Nextflow bugs: attach the .nextflow.log file or relevant excerpt -->
<!-- For variant classification bugs: attach the relevant portion of the annotated VCF -->
