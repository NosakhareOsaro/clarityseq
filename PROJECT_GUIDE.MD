# ClaritySeq — Claude Code Build Prompt (FINAL, UPDATED)

## Complete Start-to-Finish Instructions · June 2026 Edition

### Incorporates: ACGS 2024 v1.2 · ClinGen SVI 2024 · DRAGEN-GATK · DeepVariant 1.8 · gnomAD v4.1 · AlphaMissense · Beacon v2.1.1 · VRS v2.0

---

> **How to use this prompt:**
> Open Claude Code in your terminal (`claude`) and paste this entire document.
> Claude Code will read it top-to-bottom and execute every step in order.
> Each section ends with a `git commit` checkpoint.
> Do not skip sections — later steps depend on earlier ones.
> This is the ONLY prompt you need. Do not use any previous version.

---

## PART 0 — CONTEXT, STACK & NON-NEGOTIABLE STANDARDS

### 0.1 What you are building

You are building **ClaritySeq** — a production-grade, research-novel
whole-genome sequencing (WGS) clinical variant interpretation platform
with seven distinct novel contributions not present in any existing
public pipeline.

**Pipeline inputs:** paired-end FASTQ (Illumina WGS, 30× minimum)
**Pipeline outputs:**

- DRAGMAP + GATK4 DRAGEN-GATK VCF, benchmarked vs GIAB truth sets
- DeepVariant v1.8.0 parallel VCF + ensemble caller output
- HPRC pangenome (vg giraffe) VCF, surjected to GRCh38
- T2T-CHM13 v2.0 VCF with CrossMap liftover to GRCh38
- Three-way benchmark report (GRCh38 vs T2T vs pangenome) on HG001/HG002/HG005
- VEP v111 + AlphaMissense + Pangolin annotated VCF
- Bayesian ACMG/AMP classification (BayesACMG PyPI package)
  implementing all 28 rules + ClinGen SVI 2024 updates + VCEP lookup
- NHS GMS-style HTML + PDF clinical report with MANE Select transcripts,
  Bayesian credible intervals, reclassification review date
- ClinVar reclassification daemon (Celery + FHIR R4 Task IG v3.0.0)
- GA4GH Beacon v2.1.1 API with VRS identifiers + GA4GH Passports
- Phenopackets v2 clinical input + Exomiser 14 integration
- Cyrius CYP2D6 SV calling + PharmVar + CPIC dosing
- GATK Mutect2 mito mode + Haplogrep3 + ACGS 2024 §6 mito criteria
- ExpansionHunter v5.0 catalog (60 loci) + STRipy 2024 thresholds
- ClinVar submission workflow (mandatory for NHS England labs, ACGS 2024)

### 0.2 Complete technology stack

| Layer                      | Technology               | Version                      | Notes                                           |
| -------------------------- | ------------------------ | ---------------------------- | ----------------------------------------------- |
| Pipeline orchestration     | Nextflow DSL2            | Latest stable                | nf-core module style; sarek v3.8.1 as reference |
| Primary aligner            | DRAGMAP                  | v1.3.0+                      | DRAGEN-GATK Best Practices; replaces BWA-MEM2   |
| Fallback aligner           | BWA-MEM2                 | v2.2.1                       | `-profile bwa_mem2` only                        |
| Variant calling (germline) | GATK4 HaplotypeCaller    | **4.6.0.0**                  | DRAGEN-GATK mode; BQSR SKIPPED                  |
| Variant calling (DL)       | DeepVariant              | **1.8.0**                    | WGS model; pangenome-aware                      |
| Trio variant calling       | DeepTrio                 | **1.8.0**                    | For proband+parents samples                     |
| Variant calling (mito)     | GATK Mutect2             | 4.6.0.0                      | --mitochondria-mode                             |
| Repeat expansion           | ExpansionHunter          | v5.0 catalog                 | 60 loci; TRGT not used (short reads)            |
| Pangenome alignment        | vg giraffe               | Latest stable                | HPRC v1.1 graph                                 |
| Annotation                 | VEP                      | **v111**                     | MANE Select priority                            |
| Missense predictor         | AlphaMissense            | 2023 scores                  | PRIMARY PP3/BP4 (ClinGen SVI 2024)              |
| Splicing predictor         | SpliceAI + Pangolin      | Latest                       | ClinGen SVI splicing framework                  |
| Variant scores             | dbNSFP                   | **v4.7**                     | Includes AlphaMissense + ESM1b                  |
| Population frequencies     | gnomAD                   | **v4.1**                     | 807,162 individuals; April 2024                 |
| ACMG classification        | BayesACMG (our package)  | 0.1.0                        | ACGS 2024 + ClinGen SVI 2024                    |
| Clinical guidelines        | ACGS 2024 v1.2           | Feb 2024                     | Durkie et al.; replaces ACGS 2020               |
| ClinGen SVI                | 2024 recommendations     | Current                      | PM2→Supporting; AlphaMissense approved          |
| Clinical input             | Phenopackets v2          | v2.0                         | + phenopacket-tools validation                  |
| Clinical output            | FHIR R4                  | Genomics Reporting IG v3.0.0 | VRS identifiers; recontact Task                 |
| Federation API             | GA4GH Beacon             | **v2.1.1**                   | VRS aligned; Dec 2024 release                   |
| Variant IDs                | GA4GH VRS                | v2.0                         | Canonical variant representation                |
| PGx                        | Cyrius + PharmVar + CPIC | Latest                       | CYP2D6 SV                                       |
| Haplogroup                 | Haplogrep3               | Latest                       | Mito haplogroup                                 |
| Reclassification           | Celery + Redis 7         | —                            | Weekly ClinVar diff daemon                      |
| Database                   | PostgreSQL               | **16**                       | db.t4g.micro on AWS                             |
| Cache/broker               | Redis                    | **7.x**                      | ElastiCache                                     |
| Reporting                  | Jinja2 + weasyprint      | Latest                       | HTML + PDF                                      |
| Infrastructure             | Terraform                | —                            | AWS provider **v5**                             |
| Container (pipeline)       | Docker                   | Ubuntu **24.04**             | Multi-stage                                     |
| Container (HPC)            | Singularity/Apptainer    | —                            | .def file                                       |
| CI/CD                      | GitHub Actions           | —                            | chr22 CI + nightly pangenome                    |
| Testing                    | pytest                   | —                            | ≥ 90% coverage                                  |
| Python                     | —                        | **3.12**                     | All packages                                    |

### 0.3 Non-negotiable coding standards

Apply ALL of the following to **every single file**:

1. **Module-level docstring** (≥ 3 sentences): what the module does, inputs,
   outputs, external dependencies, and which guidelines/papers it implements
   (with PMID or DOI where applicable).

2. **Function docstrings** (Google style): description, Args, Returns, Raises,
   Example. ACMG rule functions must cite the specific guideline document,
   section, and PMID.

3. **Inline comments on every non-obvious line.** Prefer `# Why:` over
   `# What:`. Every flag passed to a bioinformatics tool must have an inline
   comment explaining why that value was chosen.

4. **Type annotations** on every Python function signature.

5. **Nextflow module header comment block** on every process:

   ```
   // ============================================================================
   // Module: <NAME>
   // Description: <what this does and why it is in this pipeline>
   // Guidelines: <which ACGS/ClinGen/GA4GH spec this implements>
   // Inputs:  <channel names and content>
   // Outputs: <channel names and content>
   // Container: <image:tag — must be pinned, never 'latest'>
   // Docs: <URL>
   // Parameter rationale: <explain each non-default param choice>
   // Version note: <why this specific version was chosen>
   // ============================================================================
   ```

6. **Terraform resource comments**: every resource block explains what it
   creates, why those settings, and the cost implication.

7. **Dockerfile comments**: every RUN instruction is commented.

8. No file longer than 500 lines — split into logical sub-modules.

9. Every TODO marked `# TODO(clarityseq):` with a description.

10. Every ACMG rule function cites: Richards 2015 PMID:25741868, the ClinGen
    SVI recommendation, and the ACGS 2024 v1.2 section number.

11. All tool versions **pinned** in Docker and Conda/pip. Never use `latest`.

12. Every numeric threshold (PP3 cut-off, BA1 frequency, coverage minimum)
    must have an inline comment citing the paper or guideline it comes from.

13. Test fixtures for ACMG rules must use **real ClinVar RCV accession numbers**
    and document which variant the fixture represents.

---

## PART 1 — REPOSITORY SETUP

### Step 1.1 — Check prerequisites

Run these checks. Print each version. Document missing tools in
`docs/setup/prerequisites.md`:

```
Required (fail loudly if absent):
- git ≥ 2.40
- gh (GitHub CLI) — must be authenticated (gh auth status)
- nextflow ≥ 24.x (DSL2)
- docker ≥ 25.x
- python3 = 3.12.x exactly (not 3.11, not 3.13)
- conda or mamba
- terraform ≥ 1.7
- java ≥ 17 (required by GATK 4.6.0.0)

Optional (note if absent but do not fail):
- singularity or apptainer ≥ 1.2
- node ≥ 20 (LTS)
```

### Step 1.2 — Create the GitHub repository

Using `gh repo create`:

- Name: `clarityseq`
- Public visibility
- Description: "Production WGS clinical interpretation: DRAGEN-GATK +
  pangenome + Bayesian ACMG (ACGS 2024 v1.2) + ClinVar reclassification
  daemon + GA4GH Beacon v2.1.1 + VRS v2.0 · ClaritySeq 2026"

Clone locally and `cd` into it.

### Step 1.3 — Create directory structure

Create every directory listed below. Every directory gets a `README.md`
explaining what belongs there. Empty directories get `.gitkeep`.

```
clarityseq/
├── .github/
│   ├── workflows/
│   │   ├── ci.yml
│   │   ├── ci_pangenome.yml
│   │   └── release.yml
│   ├── PULL_REQUEST_TEMPLATE.md
│   └── ISSUE_TEMPLATE/
│       ├── bug_report.md
│       └── feature_request.md
├── pipelines/
│   ├── wgs_grch38.nf
│   ├── wgs_t2t.nf
│   ├── wgs_pangenome.nf
│   ├── mito.nf
│   └── expansions.nf
├── modules/
│   ├── dragmap/               # PRIMARY aligner (DRAGEN-GATK)
│   ├── bwa_mem2/              # FALLBACK aligner only
│   ├── fastp/
│   ├── gatk4/
│   │   ├── markduplicates/
│   │   ├── haplotypecaller/   # DRAGEN-GATK mode
│   │   ├── genotypegvcfs/
│   │   ├── vqsr/
│   │   ├── mutect2_mito/
│   │   └── applybqsr/         # bwa_mem2 profile ONLY
│   ├── deepvariant/           # v1.8.0
│   ├── deeptrio/              # v1.8.0 for trio samples
│   ├── ensemble/
│   ├── vep/                   # v111; MANE Select
│   ├── alphamissense/         # ClinGen SVI PP3/BP4
│   ├── vg_giraffe/
│   ├── expansionhunter/       # v5.0 catalog
│   ├── cyrius/
│   ├── haplogrep3/
│   ├── samtools/
│   ├── mosdepth/
│   └── somalier/
├── bayesacmg/
│   ├── src/bayesacmg/
│   │   ├── __init__.py
│   │   ├── models.py
│   │   ├── rules/
│   │   │   ├── __init__.py
│   │   │   ├── pathogenic.py  # PVS1, PS1-4, PM1-6, PP1-5
│   │   │   ├── benign.py      # BA1, BS1-4, BP1-7
│   │   │   ├── splicing.py    # ClinGen SVI splicing framework
│   │   │   └── mito.py        # ACGS 2024 §6 mito rules
│   │   ├── combinations.py    # Novel ClinGen SVI 2024 combinations
│   │   ├── vcep_client.py     # ClinGen CSpec registry lookup
│   │   ├── model.py           # PyMC Bayesian model
│   │   ├── calibration.py
│   │   └── cli.py
│   ├── tests/
│   ├── pyproject.toml
│   └── README.md
├── reclassification/
│   ├── daemon.py
│   ├── clinvar_diff.py
│   ├── fhir_task.py           # FHIR Genomics Reporting IG v3.0.0
│   ├── clinvar_submitter.py   # NHS mandatory ClinVar submission
│   ├── notifier.py
│   ├── models.py
│   └── tests/
├── beacon_api/
│   ├── main.py
│   ├── routers/
│   │   ├── info.py
│   │   ├── g_variants.py      # VRS identifiers in response
│   │   └── individuals.py
│   ├── auth/passports.py
│   ├── vrs_utils.py           # GA4GH VRS v2.0
│   ├── db/
│   │   ├── models.py
│   │   └── session.py
│   └── tests/
├── annotation/
│   ├── vep_runner.py          # VEP v111; MANE Select; dbNSFP v4.7
│   ├── alphamissense_client.py
│   ├── dbnsfp_client.py       # dbNSFP v4.7
│   ├── mane_select.py         # MANE Select transcript utilities
│   ├── panelapp_client.py
│   ├── clinvar_client.py
│   ├── gnomad_client.py       # gnomAD v4.1
│   ├── clingen_client.py
│   └── tests/
├── prioritisation/
│   ├── hpo_scorer.py
│   ├── inheritance_filter.py
│   ├── exomiser_client.py     # Exomiser 14 via Phenopackets v2
│   ├── ranking.py
│   └── tests/
├── reporting/
│   ├── report_generator.py    # ACGS 2024 compliant
│   ├── pdf_renderer.py
│   ├── audit_logger.py
│   ├── templates/
│   │   ├── base.html.j2
│   │   ├── variants.html.j2   # MANE Select + Bayesian CI + review date
│   │   ├── mito.html.j2       # ACGS 2024 §6
│   │   ├── expansions.html.j2
│   │   └── pgx.html.j2
│   └── tests/
├── pgx/
│   ├── cyrius_runner.py
│   ├── pharmvar_client.py
│   ├── cpic_client.py
│   └── tests/
├── phenopackets_input/
│   ├── schema_validator.py    # phenopacket-tools + Python SDK
│   ├── converter.py
│   ├── fhir_mapper.py         # Phenopackets v2 ↔ FHIR R4
│   ├── examples/
│   └── tests/
├── multi_ancestry/
│   ├── somalier_runner.py
│   ├── ancestry_assigner.py
│   ├── vqsr_selector.py
│   └── tests/
├── benchmark/
│   ├── run_hap_py.py
│   ├── parse_results.py
│   └── results/
├── conf/
│   ├── base.config
│   ├── test.config
│   ├── aws.config
│   ├── hpc.config
│   ├── resources.config
│   └── dragen_gatk.config
├── docker/
│   ├── Dockerfile.pipeline    # Ubuntu 24.04; GATK 4.6.0.0
│   ├── Dockerfile.beacon
│   └── Dockerfile.daemon
├── singularity/clarityseq.def
├── terraform/
│   ├── main.tf                # AWS provider v5
│   ├── variables.tf
│   ├── outputs.tf
│   ├── ecr.tf
│   ├── batch.tf
│   ├── s3.tf
│   ├── iam.tf
│   └── rds.tf                 # PostgreSQL 16; db.t4g.micro
├── docs/
│   ├── conf.py
│   ├── index.rst
│   ├── setup/prerequisites.md
│   ├── adr/
│   │   ├── 001-why-dragmap-over-bwa-mem2.md
│   │   ├── 002-why-pymc-over-stan.md
│   │   ├── 003-why-celery-over-apscheduler.md
│   │   ├── 004-why-vg-giraffe-over-vg-map.md
│   │   ├── 005-why-alphamissense-primary-pp3-bp4.md
│   │   └── 006-why-pm2-supporting-not-moderate.md
│   ├── pipeline/
│   ├── api/
│   └── guides/
│       ├── quickstart.md
│       ├── data_setup.md
│       ├── aws_deployment.md
│       ├── acmg_classification.md
│       └── clinvar_submission.md
├── calibration/
│   ├── download_clingen_curations.py
│   ├── run_calibration.py
│   └── results/
├── joss_paper/paper.md
├── nextflow.config
├── CITATION.cff
├── LICENSE
├── README.md
├── Makefile
├── pyproject.toml
├── .pre-commit-config.yaml
├── .gitignore
└── .env.example
```

### Step 1.4 — Root configuration files

**`.gitignore`** — include all of:

- `*.fastq.gz`, `*.bam`, `*.vcf.gz`, `*.cram`, `*.gvcf.gz`
- `.nextflow/`, `work/`, `results/`
- `.terraform/`, `*.tfstate`, `*.tfstate.backup`
- `__pycache__/`, `*.pyc`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`
- `.env`, `*.pem`, `*.key`, `*.p12`
- `venv/`, `.venv/`, `conda-env/`
- `AlphaMissense_hg38.tsv.gz*`
- `*.dragen.hash` (DRAGMAP reference files)

**`.env.example`** — with an explanatory comment for every variable:

```
# AWS
AWS_REGION=eu-west-2
AWS_ACCOUNT_ID=

# Database
POSTGRES_DSN=postgresql+asyncpg://user:pass@localhost:5432/clarityseq
REDIS_URL=redis://localhost:6379/0

# External APIs
CLINVAR_FTP_URL=ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
PANELAPP_BASE_URL=https://panelapp.genomicsengland.co.uk/api/v1
CPIC_API_URL=https://api.cpicpgx.org/v1
PHARMVAR_API_URL=https://www.pharmvar.org/api-service
CLINGEN_CSPEC_URL=https://cspec.genome.network/cspec/api/svi/

# GA4GH
BEACON_JWT_SECRET=
BEACON_EGA_REGISTRATION_URL=https://ega-archive.org/beacon-registration
FHIR_SERVER_URL=
FHIR_GENOMICS_IG_VERSION=3.0.0
VRS_SEQREPO_URL=https://services.genomicmedlab.org/seqrepo

# Notifications
SLACK_WEBHOOK_URL=
SENDGRID_API_KEY=

# ClinVar submission (mandatory for NHS England labs — ACGS 2024)
NCBI_API_KEY=
CLINVAR_SUBMISSION_ORG_ID=

# Annotation data paths
ALPHAMISSENSE_SCORES_PATH=/data/AlphaMissense_hg38.tsv.gz
DBNSFP_DB_PATH=/data/dbNSFP4.7_variant.chr*.gz
```

**`pyproject.toml`** (root) — configure:

- black line length 88
- ruff rules: E, F, I, B, UP, ANN (enforce type annotations)
- mypy strict mode, python_version = "3.12"
- pytest: testpaths, cov-fail-under = 90
- isort profile black

**`.pre-commit-config.yaml`** — hooks:

- pre-commit-hooks (trailing-whitespace, end-of-file-fixer, check-yaml,
  check-added-large-files max_kB=500, check-case-conflict)
- black, ruff, mypy, hadolint, detect-secrets

**`Makefile`** — targets:
`install`, `test`, `lint`, `docs`, `figures`, `docker-build`,
`ci-test`, `beacon-start`, `daemon-start`, `clean`

**`CITATION.cff`** — CFF v1.2.0:

- title: "ClaritySeq: a DRAGEN-GATK pangenome-aware clinical WGS platform
  with Bayesian ACMG classification (ACGS 2024 v1.2), automated ClinVar
  reclassification monitoring, and GA4GH Beacon v2.1.1 federation"
- version: "0.1.0"
- keywords: wgs, dragen-gatk, pangenome, bayesian-acmg, acgs-2024,
  clingen-svi, alphamissense, fhir, beacon, vrs, pharmacogenomics

**`README.md`** — include:

- Badges: CI, coverage, DOI, licence, Python 3.12, nf-core
- "What makes this different from nf-core/sarek" — be specific
- ACGS 2024 v1.2 compliance statement
- Quick-start (≤ 5 commands)
- Mermaid architecture diagram

### Step 1.5 — First commit

```
git add -A
git commit -m "chore: initial repository scaffold

- Full directory structure with README.md in each major directory
- .gitignore: excludes genomic data, DRAGMAP hash, AlphaMissense TSV
- .env.example: all variables documented including NCBI ClinVar submission
- LICENSE (MIT), CITATION.cff (CFF v1.2.0)
- Makefile (10 targets), pyproject.toml (black/ruff-ANN/mypy-3.12/pytest-90%)
- .pre-commit-config.yaml (black+ruff+mypy+hadolint+detect-secrets)
- README.md: 7 novel contributions, ACGS 2024 compliance, nf-core comparison

Stack: DRAGEN-GATK · DeepVariant 1.8.0 · gnomAD v4.1 · ACGS 2024 v1.2 ·
AlphaMissense · VEP v111 (MANE Select) · Beacon v2.1.1 · VRS v2.0 · Python 3.12"

git push origin main
```

---

## PART 2 — NEXTFLOW PIPELINE CORE (DRAGEN-GATK + GRCh38 ARM)

### Step 2.1 — `nextflow.config`

Write a fully commented `nextflow.config`. Comment every parameter.

Key params to include with their comments:

```groovy
params {
    // ── Aligners ──────────────────────────────────────────────────────
    dragmap_reference = null
    // DRAGMAP hash table directory (built with dragen-os --build-hash-table)
    // This is the PRIMARY aligner reference. Different format from FASTA.
    // Pre-built hg38 hash tables: s3://broad-references/hg38/v0/dragmap/
    // BQSR is NOT run when using DRAGMAP (BQD model in HaplotypeCaller
    // replaces BQSR in DRAGEN-GATK mode — see conf/dragen_gatk.config)

    reference_fasta   = null
    // GRCh38 FASTA — used only with -profile bwa_mem2 (fallback)
    // and for T2T liftover step

    // ── gnomAD v4.1 ───────────────────────────────────────────────────
    gnomad_vcf        = null
    // gnomAD v4.1 sites VCF (released April 19, 2024)
    // IMPORTANT: Use v4.1, NOT v4.0. v4.0 had an allele number
    // calculation bug fixed in v4.1. v4.1 contains 807,162 individuals
    // (730,947 exomes + 76,215 genomes including 416,555 UK Biobank).
    // Reference: https://gnomad.broadinstitute.org/news/2024-05-gnomad-v4-1-updates/

    // ── AlphaMissense ─────────────────────────────────────────────────
    alphamissense_tsv = null
    // Pre-computed AlphaMissense scores (Cheng et al. 2023, Science)
    // ClinGen SVI approved for PP3/BP4 evidence in 2024
    // Thresholds: ≥0.564→PP3; ≤0.340→BP4; ambiguous otherwise
    // Download: gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz (~2.7 GB)
    // Must be tabix-indexed before use

    // ── VEP v111 ──────────────────────────────────────────────────────
    vep_cache_dir     = null
    vep_cache_version = 111
    // VEP cache must be v111. Do NOT use older versions.
    // v111 includes improved MANE Select annotation and dbNSFP v4.7 support.
    // MANE Select is used as the primary clinical transcript per
    // ACGS 2024 v1.2 and Morales et al. 2022 Nature PMID:35356062

    // ── Pipeline arms ─────────────────────────────────────────────────
    run_deepvariant  = true
    // Run DeepVariant v1.8.0 in parallel with GATK4; merge via ensemble caller
    run_deeptrio     = false
    // Enable for trio samples (requires PED file); uses DeepTrio v1.8.0
    run_mito         = true
    // GATK Mutect2 mito mode + Haplogrep3; follows ACGS 2024 §6
    run_expansions   = true
    // ExpansionHunter v5.0 catalog (60 loci)
    run_pgx          = true
    // Cyrius CYP2D6 SV + PharmVar + CPIC dosing
    run_pangenome    = false
    // vg giraffe HPRC alignment arm (resource-intensive; nightly CI only)
    run_t2t          = false
    // T2T-CHM13 v2.0 arm with GRCh38 liftover
    ancestry_vqsr    = true
    // Per-ancestry VQSR model selection via somalier ancestry inference
}
```

**`conf/dragen_gatk.config`** — dedicated DRAGEN-GATK settings:

```groovy
// conf/dragen_gatk.config
// ========================
// DRAGEN-GATK mode configuration.
//
// CRITICAL: BQSR (Base Quality Score Recalibration) is DISABLED
// in this config. When using DRAGMAP + HaplotypeCaller DRAGEN-GATK
// mode, BQSR must NOT be run because:
// 1. HaplotypeCaller's BQD (Base Quality Dropoff) genotyping model
//    replaces BQSR by modelling systematic errors internally
// 2. Applying BQSR to DRAGMAP-aligned reads produces quality scores
//    incompatible with the BQD model, REDUCING accuracy
// Reference: Broad Institute DRAGEN-GATK documentation
// https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939
//
// If you need BQSR (e.g., for BWA-MEM2 alignment), use -profile bwa_mem2

params {
    run_bqsr = false  // Disabled in DRAGEN-GATK mode
}
```

**`profiles` block** — with full comments explaining each:

- `test`: chr22 only; DRAGMAP; DeepVariant; no pangenome; max 8 CPUs
- `local`: Docker; DRAGMAP primary
- `aws`: AWS Batch; Spot fleet; ECR images; S3 staging
- `hpc`: Singularity; SLURM executor
- `bwa_mem2`: FALLBACK only — BWA-MEM2 + BQSR. Comment must say:
  "Use only if DRAGMAP hash table references are unavailable. BQSR
  IS run in this profile. This was the GATK Best Practice before
  DRAGEN-GATK was released in 2021."
- `pangenome`: enables pangenome arm; 32+ CPUs; 256 GB memory minimum

### Step 2.2 — Nextflow modules

Create `modules/<name>/main.nf` for every module. Every module MUST
have the full 8-line header comment block (§0.3) and inline comments
on every flag and parameter.

**`modules/dragmap/main.nf`** — PRIMARY ALIGNER

- Process: `DRAGMAP_ALIGN`
- Container: `nfcore/dragmap:1.3.0` (pin version)
- Comment block must explain:
  - Why DRAGMAP replaces BWA-MEM2 (GATK Best Practices since 2021)
  - That BQSR is NOT run after DRAGMAP in DRAGEN-GATK mode
  - The DRAGMAP reference format (hash table, not FASTA)
  - How to build the hash table (one-time setup)
  - Reference: Sci Reports 2022 PMID:36543800

**`modules/bwa_mem2/main.nf`** — FALLBACK ALIGNER

- Process: `BWA_MEM2_ALIGN`
- Must have a prominent comment: "FALLBACK ALIGNER — only activated
  with -profile bwa_mem2. The default pipeline uses DRAGMAP."
- Comment that BQSR IS required downstream when using this module

**`modules/gatk4/markduplicates/main.nf`**

- Comment: `--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500` for patterned
  flowcells (NovaSeq); default 100 is for unpatterned (HiSeq) only

**`modules/gatk4/haplotypecaller/main.nf`** — DRAGEN-GATK MODE

- Comment: `--dragen-mode true` activates BQD genotyping model
- Comment: DO NOT pass `--bqsr-recal-file` in DRAGEN-GATK mode
- Comment: `-ERC GVCF` required for joint genotyping
- Container: `broadinstitute/gatk:4.6.0.0`

**`modules/gatk4/applybqsr/main.nf`** — bwa_mem2 PROFILE ONLY

- Prominent comment: "ONLY used with -profile bwa_mem2. Must NOT
  be run in DRAGEN-GATK mode (incompatible with BQD model)."

**`modules/gatk4/vqsr/main.nf`**

- Two processes: `VQSR_SNP` and `VQSR_INDEL`
- Comment all training resources with prior values and why each
- Comment the per-ancestry mode: "when ancestry_vqsr=true, training
  resources are ancestry-stratified gnomAD v4.1 subsets"

**`modules/deepvariant/main.nf`** — v1.8.0

- Container: `google/deepvariant:1.8.0` (pin to 1.8.0)
- Comment v1.8.0 improvements over v1.6: SPRQ support, pangenome-aware
- Comment: "For TRIO samples, use modules/deeptrio/ instead"
- Comment the pileup image vs GATK graph assembly difference

**`modules/deeptrio/main.nf`** — NEW

- Container: `google/deepvariant:1.8.0`
- Comment: "Uses parent BAMs as additional image channels; improves
  de novo SNV sensitivity by ~15% vs per-sample DeepVariant"
- Reference: PMID:36050879

**`modules/vep/main.nf`** — v111

- Container: `ensemblorg/ensembl-vep:release_111`
- Comment the pick order: `mane_select,mane_plus_clinical,canonical`
- Comment each plugin with which ACMG evidence it supports:
  AlphaMissense (PP3/BP4 PRIMARY), SpliceAI (splicing PP3/BP4/BP7),
  Pangolin (tissue-specific splicing, supplements SpliceAI),
  dbNSFP_4.7 (bulk scores), gnomAD_v4.1 (PM2/BA1/BS1),
  ClinVar (PP5/BP6 — use with caution)

**`modules/alphamissense/main.nf`** — NEW

- Process: `ALPHAMISSENSE_LOOKUP`
- Container: biocontainers/htslib (for tabix)
- Comment: ClinGen SVI 2024 approved thresholds (≥0.564→PP3; ≤0.340→BP4)
- Reference: Cheng et al. 2023 Science PMID:37703350

**`modules/expansionhunter/main.nf`** — v5.0 catalog

- Comment the 60-locus catalog vs v4.0 (20 loci)
- Comment new loci: FGF14 (SCA27B), RFC1 (CANVAS)
- Comment: "TRGT is NOT used — it requires PacBio HiFi long reads"
- Create `modules/expansionhunter/LOCI.md` documenting all 60 loci
  with: disease name, gene, repeat motif, normal range, pathogenic
  threshold (from STRipy 2024 database)

**`modules/somalier/main.nf`**

- Comment the ancestry-to-VQSR-training-set mapping:
  AFR/AMR/EAS/EUR/SAS → gnomAD v4.1 ancestry-stratified subsets
- Comment the fallback for unassigned (admixed) samples

### Step 2.3 — Main GRCh38 pipeline (`pipelines/wgs_grch38.nf`)

Full file header docstring ≥ 300 words covering:

- DRAGEN-GATK mode and why BQSR is absent
- Ensemble calling strategy (INTERSECTION mode as default)
- Per-ancestry VQSR flow via somalier
- DeepTrio trigger condition (PED file + run_deeptrio=true)
- Mitochondrial and repeat expansion sub-workflows
- ASCII data flow diagram

Comment every channel transformation (`.map{}`, `.groupTuple()`,
`.branch{}`, `.join()`, `.combine()`).

Pipeline order:

```
READ_SAMPLE_SHEET → FASTP → DRAGMAP_ALIGN (or BWA_MEM2 if bwa_mem2 profile)
→ GATK4_MARKDUPLICATES → SAMTOOLS_FLAGSTAT → MOSDEPTH (coverage gate)
→ [if bwa_mem2 profile] GATK4_BQSR → GATK4_APPLYBQSR
→ SOMALIER (ancestry inference)
→ GATK4_HAPLOTYPECALLER (DRAGEN-GATK mode; scattered by interval)
→ GATK4_GENOTYPEGVCFS
→ VQSR_MODEL_SELECT (per-ancestry if ancestry_vqsr=true)
→ GATK4_VQSR (SNP + Indel)
→ [if run_deepvariant AND NOT trio] DEEPVARIANT (parallel)
  [if run_deeptrio AND trio]        DEEPTRIO (parallel)
→ ENSEMBLE_MERGE (INTERSECTION mode by default)
→ ALPHAMISSENSE_LOOKUP
→ VEP_ANNOTATE (v111; MANE Select; dbNSFP v4.7)
→ BAYESACMG_CLASSIFY (28 rules + ClinGen SVI 2024 + VCEP lookup)
→ PHENOPACKET_VALIDATE (if phenopacket_json provided)
→ EXOMISER_PRIORITISE (via Phenopackets v2)
→ INHERITANCE_FILTER
→ VARIANT_RANK
→ [if run_mito]       MITO_PIPELINE
→ [if run_expansions] EXPANSION_HUNTER
→ [if run_pgx]        CYRIUS_PGX
→ GENERATE_REPORT (HTML + PDF + JSON-LD audit trail)
→ CLINVAR_BEACON_INGEST
→ MULTIQC
```

### Step 2.4 — Ensemble caller (`modules/ensemble/`)

Write `modules/ensemble/main.nf` and `modules/ensemble/ensemble.py`:

`ensemble.py` module docstring must explain:

- Error profile complementarity: GATK4 (graph assembly) vs DeepVariant
  (pileup image CNN) have partially non-overlapping error profiles
- INTERSECTION mode (default): PASS if in BOTH callers → higher precision
- UNION mode (optional): PASS if in EITHER caller → higher sensitivity
- INFO fields added: `ENSEMBLE_CALLER` (GATK4/DV/BOTH) and `ENSEMBLE_MODE`
- Clinical rationale for INTERSECTION as default: reduces false positives
  in clinical reports

### Step 2.5 — GitHub Actions CI

**`.github/workflows/ci.yml`:**

- Trigger: push to main, pull_request to main
- Jobs with `name:` fields explaining each:
  1. `lint`: pre-commit all files
  2. `unit-tests`: pytest + coverage ≥ 90%
  3. `integration-test`: chr22 GIAB HG001; cache FASTQs 7 days; run
     Nextflow `-profile test`; assert hap.py SNP sensitivity ≥ 99.0%
  4. `docker-build`: build all three images; verify `--help` runs

Document acceptance criteria in the workflow file header comment:

```yaml
# Acceptance criteria (CI fails if these are not met):
# - SNP sensitivity ≥ 99.0% vs GIAB HG001 chr22 truth (hap.py)
# - pytest coverage ≥ 90%
# - All pre-commit hooks pass
# - All three Docker images build
```

**`.github/workflows/ci_pangenome.yml`:**

- Nightly at 02:00 UTC; cache HPRC chr22 graph 7 days
- Assert pangenome SNP sensitivity ≥ GRCh38 SNP sensitivity on HG001

### Step 2.6 — Commit

```
git add -A
git commit -m "feat(pipeline): DRAGEN-GATK WGS pipeline core

Aligner:
- DRAGMAP v1.3.0 primary (modules/dragmap/): BQSR SKIPPED in DRAGEN-GATK
  mode; BQD genotyping model in HC replaces it (documented extensively)
- BWA-MEM2 v2.2.1 fallback (modules/bwa_mem2/): bwa_mem2 profile only;
  BQSR IS run in this profile
- conf/dragen_gatk.config: run_bqsr=false with explanation

Variant calling:
- GATK 4.6.0.0 HaplotypeCaller: DRAGEN-GATK mode; --dragen-mode true
- GATK MarkDuplicates: OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 (NovaSeq)
- DeepVariant 1.8.0 parallel arm; DeepTrio 1.8.0 for trio samples
- Ensemble caller: INTERSECTION mode default; UNION optional

Annotation:
- VEP v111: MANE Select pick order; every plugin commented with
  the ACMG evidence it supports
- AlphaMissense module: ClinGen SVI thresholds documented (≥0.564/≤0.340)
- ExpansionHunter v5.0: 60-locus catalog; LOCI.md with all 60 loci

Infrastructure:
- somalier: ancestry inference → per-ancestry VQSR
- gnomAD v4.1 referenced throughout (not v4.0)
- GitHub Actions: chr22 CI (hap.py ≥99% SNP sensitivity); nightly pangenome
- 6 Architecture Decision Records in docs/adr/

Tool versions pinned: GATK 4.6.0.0 · DV 1.8.0 · VEP 111 · gnomAD v4.1"

git push origin main
```

---

## PART 3 — BAYESIAN ACMG CLASSIFIER (BayesACMG)

### Step 3.1 — Package setup (`bayesacmg/pyproject.toml`)

```toml
[project]
name = "bayesacmg"
version = "0.1.0"
description = """Bayesian ACMG/AMP variant classifier implementing ACGS
Best Practice Guidelines 2024 v1.2 (Durkie et al.) and ClinGen SVI 2024
recommendations. Implements all 28 ACMG/AMP criteria with: PM2 at Supporting
weight (ClinGen SVI 2024), AlphaMissense as primary PP3/BP4 predictor,
ClinGen SVI splicing subgroup framework (Walker et al. 2023), ACGS 2024 §6
mito-specific rules, VCEP gene-specific override lookup, novel ClinGen SVI
2024 evidence combinations (PVS1+PM2_Supporting=LP), and Bayesian
Dirichlet-Multinomial model with calibrated 95% credible intervals."""
requires-python = ">=3.12"
dependencies = [
    "pymc>=5.10",
    "arviz>=0.18",
    "numpy>=1.26",
    "scipy>=1.12",
    "pandas>=2.2",
    "pydantic>=2.6",
    "click>=8.1",
    "rich>=13",
    "httpx>=0.27",
    "tenacity>=8",
]
```

### Step 3.2 — Data models (`bayesacmg/src/bayesacmg/models.py`)

Full module docstring citing all four guidelines (Richards 2015,
Tavtigian 2020, ACGS 2024 v1.2, ClinGen SVI 2024).

```python
class EvidenceStrength(str, Enum):
    """
    ACMG/AMP evidence strength categories and point values.

    Point values: Tavtigian et al. 2020 PMID:32645316 (ClinGen SVI system).

    CRITICAL CHANGE FROM ACGS 2020: PM2 is now at SUPPORTING (1 pt).
    ClinGen SVI 2024 recommends applying PM2 at Supporting strength in
    most contexts. gnomAD v4.1 (807,162 individuals, April 2024) reveals
    high rates of ultra-rare variants, making rarity alone weaker evidence
    than assumed in 2015. Apply PM2 at SUPPORTING (1 pt) by default.
    VCEP specifications may override per-gene via vcep_client.py.
    Reference: ClinGen SVI PM2 guidance (2024)
    https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
    """
    VERY_STRONG       = "very_strong"       # PVS: 8 pts
    STRONG            = "strong"            # PS: 4 pts
    MODERATE          = "moderate"          # PM (not PM2): 2 pts
    SUPPORTING        = "supporting"        # PP + PM2: 1 pt
    STAND_ALONE       = "stand_alone"       # BA1: → Benign directly
    STRONG_BENIGN     = "strong_benign"     # BS: -4 pts
    SUPPORTING_BENIGN = "supporting_benign" # BP: -1 pt
```

Document every field of every model with inline comments.

### Step 3.3 — ACMG rules: pathogenic (`bayesacmg/src/bayesacmg/rules/pathogenic.py`)

Every rule function must cite:

1. Richards et al. 2015 PMID:25741868 (original rule)
2. ClinGen SVI recommendation (if updated)
3. ACGS 2024 v1.2 section (if UK practice differs)

**PVS1** — implement Abou Tayoun 2018 decision tree:

```python
def rule_pvs1(variant: VariantInput, transcript: TranscriptData,
              gene: GeneData) -> ACMGRule:
    """
    PVS1 — Null variant where LoF is a known disease mechanism.

    Guidelines:
    - Richards et al. 2015 PMID:25741868: original criterion definition
    - Abou Tayoun et al. 2018 PMID:30192042: PVS1 decision tree v1.0
    - Walker et al. 2023 PMID:36898414 (ClinGen SVI Splicing Subgroup):
      PVS1 for canonical splice site variants — refer to rules/splicing.py
    - ACGS 2024 v1.2 §5, Notes 1-3: UK-specific PVS1 application

    MANE Select requirement (ACGS 2024):
    PVS1 should be applied with reference to the MANE Select transcript.
    If the LoF only affects non-MANE transcripts, reduce strength by
    one level (e.g., Very Strong → Strong).
    MANE Select check via annotation.mane_select.adjust_pvs1_for_mane()
    """
```

**PM2** — MUST be at Supporting weight:

```python
def rule_pm2(variant: VariantInput) -> ACMGRule:
    """
    PM2 — Absent from/extremely rare in population databases.

    DEFAULT WEIGHT: SUPPORTING (1 point) — not Moderate (2 points).

    ClinGen SVI 2024 update: PM2 should be applied at Supporting strength
    in most contexts. gnomAD v4.1 (807,162 individuals, April 2024)
    reveals that many ultra-rare variants exist in the general population,
    making absence from population databases less distinctive evidence
    than assumed when the 2015 guidelines were written.

    When to apply:
    - Variant absent from gnomAD v4.1 (all populations combined)
    - OR variant present at extremely low frequency
      (AF < 0.0001 in most population groups; adjust per VCEP spec)

    Per-gene override:
    Some VCEP specifications allow PM2 at Moderate for very specific
    gene-disease pairs. Check vcep_client.py before applying this rule.

    References:
    - Richards et al. 2015 PMID:25741868: original PM2 definition
    - ClinGen SVI PM2 guidance (2024): https://clinicalgenome.org/tools/
    - ACGS 2024 v1.2 §5 Table 2, Appendix C: PM2 mini impact assessment
    """
```

**PP3** — AlphaMissense as PRIMARY predictor:

```python
def rule_pp3(variant: VariantInput,
             alphamissense_score: float | None,
             revel_score: float | None) -> ACMGRule:
    """
    PP3 — Multiple lines of computational evidence support pathogenicity.

    PRIMARY PREDICTOR: AlphaMissense (ClinGen SVI 2024 approved)
    AlphaMissense threshold: score ≥ 0.564 → PP3 (Supporting)
    Reference: Cheng et al. 2023 Science doi:10.1126/science.adg7492
               PMID:37703350

    ClinGen SVI 2024 approved tools for PP3:
    1. AlphaMissense: ≥ 0.564 (PRIMARY — use this first)
    2. REVEL: ≥ 0.7 (retained for comparison/validation)
    3. BayesDel: ≥ 0.13 (from dbNSFP v4.7)
    4. CADD: ≥ 25 (PHRED-scaled)

    IMPORTANT: For SPLICE SITE VARIANTS, do not use this function.
    Use rules/splicing.py rule_splicing_pp3_bp4_bp7() instead.
    The ClinGen SVI Splicing Subgroup (Walker et al. 2023) defines
    separate PP3/BP4/BP7 criteria for splice-impacting variants.

    Reference:
    - ClinGen SVI in silico tool recommendation (2024)
      https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
    """
```

Implement ALL 15 pathogenic rules (PVS1, PS1-4, PM1-6, PP1-5).

### Step 3.4 — Splicing framework (`bayesacmg/src/bayesacmg/rules/splicing.py`)

NEW FILE. Full module docstring citing Walker et al. 2023 PMID:36898414.

```python
"""
bayesacmg.rules.splicing
========================

ClinGen SVI Splicing Subgroup evidence framework for PP3, BP4, and BP7.

Reference:
    Walker et al. 2023 Am J Hum Genet 110:1046-1067
    doi:10.1016/j.ajhg.2023.06.002 PMID:36898414

Key thresholds (SpliceAI Δ score):
    ≥ 0.5 → strong splice impact → PP3 Strong
    ≥ 0.2 → moderate splice impact → PP3 Moderate
    < 0.1 → no splice impact predicted
    Synonymous variant + SpliceAI < 0.1 → BP7

When Pangolin and SpliceAI disagree:
    Use the more conservative (lower Δ score) estimate.
    Document the disagreement in the evidence_items field.

This function returns a list[ACMGRule] because PP3, BP4, and BP7
are mutually exclusive for the same variant+evidence — only one
can apply. The list will always contain exactly one ACMGRule.
"""
```

### Step 3.5 — Mito rules (`bayesacmg/src/bayesacmg/rules/mito.py`)

NEW FILE. Full module docstring citing ACGS 2024 v1.2 §6.

Key rules to implement:

- `rule_mito_haplogroup_defining()`: returns Benign if variant is
  haplogroup-defining per Haplogrep3 output (MUST run before any
  pathogenicity assessment)
- `rule_mito_pm2()`: mito-specific PM2 (different thresholds; check
  haplogroup first)
- `rule_mito_ba1()`: mito-specific BA1 (standard threshold does NOT
  apply to mito variants)
- `rule_mito_heteroplasmy_level()`: maps heteroplasmy % to clinical
  significance category per ACGS 2024 §6

### Step 3.6 — VCEP client (`bayesacmg/src/bayesacmg/vcep_client.py`)

NEW FILE. Query ClinGen CSpec registry before classifying each variant.

Module docstring must explain:

- Why VCEP specs are essential (BRCA1/2, RASopathy genes have
  different thresholds that OVERRIDE the general 28-rule framework)
- The CSpec registry API endpoint
- Caching strategy (24-hour TTL per gene)

### Step 3.7 — Novel combinations (`bayesacmg/src/bayesacmg/combinations.py`)

NEW FILE. ClinGen SVI 2024 combinations not in Richards 2015.

```python
"""
bayesacmg.combinations
======================

Novel ClinGen SVI 2024 evidence combinations extending Richards 2015.

Critical combination (addresses PM2 downgrade impact):
    PVS1 (Very Strong, 8 pts) + PM2_Supporting (1 pt) = 9 pts
    → Likely Pathogenic (threshold = 6 pts for LP; Bayesian PP ≥ 0.90)

Without this combination, novel LoF variants in LoF-mechanism genes
where the only secondary evidence is rarity would no longer reach LP
after the PM2 downgrade. This combination explicitly preserves LP
classification for such variants.

Reference: ClinGen SVI Working Group 2024 recommendations
https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
See also: ACGS 2024 v1.2 §5 Table 2 (UK implementation)
"""
```

### Step 3.8 — Bayesian model (`bayesacmg/src/bayesacmg/model.py`)

Module docstring ≥ 500 words. Must include:

- Scientific motivation for Bayesian over fixed point-scoring
- The Dirichlet-Multinomial model structure with diagram
- Prior choices and justification
- PM2 prior weight = Supporting (1 pt) per ClinGen SVI 2024
- Calibration strategy against ClinGen 500-variant set
- ECE acceptance criterion (< 0.05)

Critical inline comment in the model:

```python
# PM2 prior weight: ClinGen SVI 2024 updated PM2 from Moderate (2 pts)
# to Supporting (1 pt). The Dirichlet prior for PM2 evidence strength
# is therefore centred on the Supporting category, not Moderate.
# The concentration parameter is set to 5.0 for PM2 (vs 10.0 for more
# settled criteria like PVS1) to reflect ongoing community uncertainty
# about the "true" weight after the 2024 revision.
```

### Step 3.9 — Tests (`bayesacmg/tests/`)

pytest suite achieving ≥ 90% coverage:

- `test_rules_pathogenic.py`: one test per rule. Each test uses a real
  ClinVar RCV accession number documented in the test docstring.
  MUST include tests for:
  - PM2 returning SUPPORTING (1 pt) not MODERATE (2 pts)
  - PP3 using AlphaMissense score ≥ 0.564
  - BP4 using AlphaMissense score ≤ 0.340
  - PVS1 decision tree: frameshift, start-loss, canonical splice
  - Novel combination: PVS1 + PM2_Supporting = LP classification

- `test_rules_benign.py`: one test per benign rule
- `test_splicing.py`: all Walker et al. 2023 SpliceAI threshold branches
- `test_mito.py`: ACGS 2024 §6 mito rules
- `test_vcep_client.py`: mock ClinGen CSpec API
- `test_model.py`: Bayesian model unit tests
- `test_calibration.py`: ECE < 0.05 on fixture set

### Step 3.10 — Commit

```
git add -A
git commit -m "feat(bayesacmg): Bayesian ACMG classifier — ACGS 2024 + ClinGen SVI 2024

Guidelines implemented (all rules cite specific PMIDs and sections):
- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., Feb 2024)
  ratified 20 February 2024; replaces ACGS 2020
- ClinGen SVI 2024: PM2→Supporting (1 pt not 2 pt; PM2 docstring
  explains gnomAD v4.1 rationale); AlphaMissense primary PP3/BP4
  (thresholds ≥0.564/≤0.340; Cheng 2023 PMID:37703350 cited)
- Walker et al. 2023 PMID:36898414: splicing PP3/BP4/BP7 framework
  (rules/splicing.py): SpliceAI+Pangolin; BP7 for silent+no splice
- ACGS 2024 §6: mito-specific rules (rules/mito.py): haplogroup-first,
  mito PM2, mito BA1, heteroplasmy level mapping
- Novel ClinGen SVI 2024 combinations (combinations.py):
  PVS1+PM2_Supporting=LP; addresses PM2 downgrade impact
- VCEP lookup (vcep_client.py): ClinGen CSpec registry; 24hr cache

Package: bayesacmg 0.1.0 (PyPI-ready; JOSS paper draft included)
Tests: ≥90% coverage; real ClinVar RCV accessions in fixtures;
  PM2=Supporting verified; AlphaMissense thresholds verified"

git push origin main
```

---

## PART 4 — ANNOTATION STACK

### Step 4.1 — AlphaMissense client (`annotation/alphamissense_client.py`)

Full module docstring explaining:

- What AlphaMissense is and how it differs from REVEL/CADD
- ClinGen SVI 2024 approval and thresholds
- Download location and tabix indexing
- How scores map to PP3/BP4

```python
"""
annotation.alphamissense_client
================================

AlphaMissense (Cheng et al. 2023, Science doi:10.1126/science.adg7492
PMID:37703350) pathogenicity score lookup for missense variants.

ClinGen SVI 2024 approved AlphaMissense for PP3/BP4 evidence:
    AM ≥ 0.564 → PP3 (pathogenic supporting)
    AM ≤ 0.340 → BP4 (benign supporting)
    0.340 < AM < 0.564 → ambiguous

AlphaMissense covers all 71 million possible human missense variants.
It was fine-tuned from AlphaFold using population variant frequency
patterns to predict pathogenicity.

Why AlphaMissense over REVEL as primary:
    AlphaMissense outperforms EVE for 77% of ACMG genes (AUROC 0.91
    on ClinVar benchmark; Cheng et al. 2023). It is the most recently
    ClinGen SVI-approved tool and reflects the latest protein structure
    prediction technology. REVEL is retained as a secondary comparator.

Data: gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz (~2.7 GB)
Index: tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz
"""
```

### Step 4.2 — dbNSFP v4.7 client (`annotation/dbnsfp_client.py`)

```python
"""
annotation.dbnsfp_client
========================

Client for dbNSFP v4.7 (released May 2024).

dbNSFP v4.7 additions over v4.6:
- AlphaMissense scores (Cheng et al. 2023)
- ESM1b protein language model scores (Lin et al. 2023)
- Interpro domain annotations (updated)
- Updated MetaRNN, SIFT4G, MutPred2

In ClaritySeq, dbNSFP v4.7 is accessed via the VEP dbNSFP plugin
(configured in modules/vep/main.nf). This Python client provides
direct lookup for use in BayesACMG rule functions where VEP is not
available (e.g., for post-VEP score updates).
"""
```

### Step 4.3 — VEP runner (`annotation/vep_runner.py`)

Must document MANE Select prominently:

```python
"""
annotation.vep_runner
=====================

Runs VEP v111 and parses output into VariantInput objects.

MANE Select transcript priority:
    VEP is run with --pick_order mane_select,mane_plus_clinical,canonical
    This ensures MANE Select transcript is used for HGVSc/HGVSp notation.
    MANE (Matched Annotation from NCBI and EBI) Select is a jointly
    maintained transcript set for clinical genomics.
    Reference: Morales et al. 2022 Nature 604:310 PMID:35356062

    All downstream reports, audit trails, and FHIR resources use
    MANE Select HGVSc/HGVSp notation.
"""
```

### Step 4.4 — MANE Select client (`annotation/mane_select.py`)

NEW FILE:

```python
"""
annotation.mane_select
======================

MANE (Matched Annotation from NCBI and EBI) Select transcript utilities.

MANE Select provides one biologically-relevant transcript per protein-
coding gene, jointly maintained by NCBI RefSeq and EMBL-EBI Ensembl.
It is the recommended primary transcript for clinical variant reporting.

Reference: Morales et al. 2022 Nature 604:310 PMID:35356062

Functions:
- is_mane_select(transcript_id): check if MANE Select
- get_mane_select_for_gene(gene_symbol): return MANE Select NM/ENST
- adjust_pvs1_for_mane(variant, transcript): reduce PVS1 strength if
  LoF only affects non-MANE transcripts (ACGS 2024 v1.2 §5)
"""
```

### Step 4.5 — gnomAD v4.1 client (`annotation/gnomad_client.py`)

Must prominently note the version:

```python
"""
annotation.gnomad_client
========================

gnomAD v4.1 (released April 19, 2024) allele frequency client.

IMPORTANT: Always use v4.1, NOT v4.0.
gnomAD v4.0 had an allele number calculation bug that produced incorrect
allele frequencies for some variants. v4.1 fixes this bug.
v4.1 contains 807,162 individuals: 730,947 exomes + 76,215 genomes,
including 416,555 UK Biobank exomes not present in v3.

PM2 application (ClinGen SVI 2024):
    Apply PM2 at SUPPORTING (1 pt) not Moderate (2 pts).
    gnomAD v4.1's large size reveals that ultra-rare variants are common
    in the general population, making rarity less distinctive evidence.
"""
```

### Step 4.6 — Commit

```
git add -A
git commit -m "feat(annotation): full annotation stack with AlphaMissense primary

- VEP v111 runner: --pick_order mane_select,mane_plus_clinical,canonical;
  every plugin commented with ACMG evidence mapping
- AlphaMissense client: tabix lookup; thresholds ≥0.564/≤0.340 with
  Cheng 2023 PMID:37703350 cited; explains why primary over REVEL
- dbNSFP v4.7 client: documents v4.7 additions (AlphaMissense, ESM1b)
- MANE Select client: is_mane_select(); adjust_pvs1_for_mane() per
  ACGS 2024 v1.2 §5; Morales 2022 PMID:35356062 cited throughout
- gnomAD v4.1 client: v4.0 bug warning prominent in docstring;
  PM2 Supporting note with ClinGen SVI 2024 reference
- PanelApp client: green/amber; rate limiting; 100 req/min documented
- ClinGen Gene Curation + Dosage Sensitivity API client"

git push origin main
```

---

## PART 5 — CLINICAL REPORTING

### Step 5.1 — Phenopackets v2 (`phenopackets_input/`)

`schema_validator.py` — use phenopacket-tools + Python SDK:

```python
"""
phenopackets_input.schema_validator
=====================================

Validates Phenopackets v2.0 JSON using phenopacket-tools (strict
validation) and the phenopackets Python SDK.

Why phenopacket-tools in addition to the Python SDK:
    The Python SDK parses and creates Phenopackets but has less strict
    schema validation than the Java phenopacket-tools library.
    phenopacket-tools enforces additional constraints (valid OMIM IDs,
    HPO term existence, ontology version metadata) that the SDK misses.

`fhir_mapper.py` provides bidirectional Phenopackets v2 ↔ FHIR R4
mapping following the 2025 ontology-based rare disease CDM paper.
"""
```

### Step 5.2 — Report templates (ACGS 2024 compliant)

**`variants.html.j2`** — comment block must explain ALL ACGS 2024 additions:

```html
{# Columns (ACGS 2024 v1.2 compliant): Gene | MANE Select Transcript | HGVSc |
HGVSp | ACMG Class (ACGS 2024) | P(Pathogenic) [95% HDI] | gnomAD v4.1 AF |
AlphaMissense Score | ClinVar | Panel Status | Inheritance | HPO Relevance | Key
ACMG Evidence | Review Date (VUS only) MANE Select (ACGS 2024 v1.2 / Morales
2022 PMID:35356062): All HGVSc and HGVSp uses MANE Select transcript notation.
P(Pathogenic) [95% HDI] (NOVEL — BayesACMG): Shows Bayesian posterior
probability and 95% highest-density interval. Example: "LP — 78% [61%–91%]" This
is the first clinical WGS report template to show ACMG classification
uncertainty in a structured, quantified format. VUS Review Date (ACGS 2024 §9
REQUIREMENT): ACGS 2024 v1.2 §9: VUS variants should be reviewed within 2 years.
Each VUS displays "Review by: [date+2years]" stored in the reclassification
daemon schedule for automated follow-up. #}
```

**`mito.html.j2`** — ACGS 2024 §6 format:

```html
{# Follows ACGS 2024 v1.2 §6: "Classifying variants in the mitochondrial genome"
Required sections: 1. Haplogroup assignment (Haplogrep3) — MUST appear first
Haplogroup-defining variants excluded from pathogenic candidates 2.
Pathogenic/VUS mito variants with heteroplasmy % 3. MITOMAP disease associations
4. Tissue specificity note (blood vs affected tissue heteroplasmy) Citation in
report footer: "Mitochondrial variants classified per ACGS 2024 v1.2 §6 and ACGS
Best Practice Guidelines for Molecular Diagnosis of Mitochondrial Disease
(November 2020)" #}
```

### Step 5.3 — Report generator (`reporting/report_generator.py`)

ACGS 2024 additions in module docstring:

```python
"""
reporting.report_generator
===========================

Generates NHS GMS-style clinical reports (HTML + PDF + JSON-LD audit).

ACGS 2024 v1.2 compliance:
1. MANE Select transcript notation (§ Morales 2022 PMID:35356062)
2. VUS review date scheduling (§9): VUS → review by date+2yr;
   stored in reclassification daemon ClinVarSubmissionQueue table
3. ClinVar submission flag (§ Introduction): novel P/LP variants →
   written to pending_clinvar_submissions queue (NHS mandatory)
4. Classification scheme citation in footer:
   "Classified per ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al.)
    and ACMG/AMP guidelines (Richards et al. 2015 PMID:25741868)"

JSON-LD audit trail includes:
- gnomAD version: "v4.1" (NOT v4.0)
- VEP cache version: "111"
- AlphaMissense score file date
- BayesACMG version and prior weights (including PM2=Supporting)
- ACGS guidelines version: "2024 v1.2"
- ClinGen SVI recommendations date: "2024"
"""
```

### Step 5.4 — Commit

```
git add -A
git commit -m "feat(reporting): ACGS 2024 v1.2 compliant clinical reports

- Phenopackets v2 validator: phenopacket-tools + Python SDK (strict);
  FHIR R4 bidirectional mapper
- HTML templates — ACGS 2024 additions:
  - MANE Select transcript column (Morales 2022 PMID:35356062)
  - P(Pathogenic) 95% Bayesian HDI column (novel contribution)
  - VUS reclassification review date per ACGS 2024 §9
  - Mito section per ACGS 2024 §6 (haplogroup first; heteroplasmy %)
  - Classification scheme citation in footer (ACGS 2024 v1.2)
- Report generator: stores VUS review dates; writes novel P/LP to
  ClinVar pending queue (NHS mandatory per ACGS 2024 intro)
- Audit trail JSON-LD: gnomAD v4.1; VEP 111; AlphaMissense date;
  PM2=Supporting documented in audit"

git push origin main
```

---

## PART 6 — CLINVAR RECLASSIFICATION DAEMON

### Step 6.1 — Database models (`reclassification/models.py`)

Include `ClinVarSubmissionQueue` with docstring citing ACGS 2024:

```python
class ClinVarSubmissionQueue(Base):
    """
    Queue of novel P/LP variants awaiting ClinVar submission.

    NHS England mandate (ACGS 2024 v1.2 Introduction):
    'Submission of variants to ClinVar by NHS laboratories in England
    is now a requirement following completion of the information
    governance review process.'

    This table stores novel P/LP variants (not already in ClinVar)
    identified during pipeline runs. The clinvar_submitter.py daemon
    processes this queue and generates ClinVar API submission XML.
    """
```

### Step 6.2 — ClinVar diff (`reclassification/clinvar_diff.py`)

Module docstring citing ACGS 2024 §9:

```python
"""
reclassification.clinvar_diff
==============================

Weekly ClinVar release diff for variant reclassification detection.

Clinical context (ACGS 2024 v1.2, §9 — Reclassification):
    'Laboratories must have processes for periodic reclassification.
    For VUS variants, review within 2 years is recommended.'
    This daemon automates the detection step.

FHIR output: fhir_task.py generates FHIR R4 Task resources using
the Genomics Reporting IG v3.0.0 'recontact' task type.

gnomAD note: Use ClinVar builds from 2024+ which reference gnomAD v4.1.
Earlier ClinVar builds used gnomAD v2/v3 for BA1/BS1 evidence.
"""
```

### Step 6.3 — FHIR Task generator (`reclassification/fhir_task.py`)

Update to Genomics Reporting IG v3.0.0 + VRS:

```python
"""
reclassification.fhir_task
===========================

FHIR R4 Task resources for reclassification re-contact.

IG version: HL7 FHIR Genomics Reporting IG v3.0.0
URL: https://hl7.org/fhir/uv/genomics-reporting/STU3/

IG v3.0.0 features used:
- 'recontact' task type in Genomics Reporting IG code system
- DiagnosticImplication profile reference in Task.focus
- GA4GH VRS v2.0 variant identifiers in Task.input
  (from beacon_api.vrs_utils.compute_vrs_id())

GA4GH VRS identifiers ensure cross-system interoperability with other
GA4GH-compliant tools and databases (Beacon, ClinVar v2 API, etc.).
"""
```

### Step 6.4 — ClinVar submitter (`reclassification/clinvar_submitter.py`)

NEW FILE with NHS mandate context:

```python
"""
reclassification.clinvar_submitter
====================================

NHS-mandated ClinVar API submission for novel P/LP variants.

ACGS 2024 v1.2 Introduction:
'Submission of variants to ClinVar by NHS laboratories in England
is now a requirement following completion of the information
governance review process.'

This module generates and submits ClinVar API submissions:
- Triggered for variants classified as P or LP (BayesACMG PP ≥ 0.90)
  that are not already in ClinVar (no RCV/SCV accession)
- Uses ClinVar Submission API (NCBI): POST to submit.ncbi.nlm.nih.gov
- Requires NCBI_API_KEY and CLINVAR_SUBMISSION_ORG_ID from .env

Submission content per variant:
- Gene + MANE Select HGVSc (from MANE Select annotation)
- ACMG 5-tier classification (P/LP/VUS/LB/B)
- Evidence codes (ACMG/AMP criteria; ACGS 2024 applied)
- BayesACMG posterior probability (as submission note)
- gnomAD v4.1 allele frequency
- Phenotype from Phenopackets v2 if available
"""
```

### Step 6.5 — Commit

```
git add -A
git commit -m "feat(reclassification): ClinVar daemon + FHIR IG v3.0.0 + VRS + NHS submission

- ClinVar diff (clinvar_diff.py): ACGS 2024 §9 context; weekly FTP diff
- FHIR Task (fhir_task.py): Genomics Reporting IG v3.0.0; 'recontact'
  task type; GA4GH VRS v2.0 identifiers in Task.input
- ClinVar submitter (NEW, clinvar_submitter.py): NHS mandatory per ACGS
  2024 Introduction; NCBI API submission; MANE Select HGVSc notation
- ClinVarSubmissionQueue SQLAlchemy model: NHS mandate cited in docstring
- Celery beat: Monday 08:00 UTC; retry on FHIR errors; dead letter queue
- Tests: synthetic ClinVar XML; FHIR IG v3.0.0 schema validation;
  ClinVar API mock; VRS identifier correctness"

git push origin main
```

---

## PART 7 — GA4GH BEACON v2.1.1 API

### Step 7.1 — FastAPI app (`beacon_api/main.py`)

```python
"""
beacon_api.main
===============

GA4GH Beacon v2.1.1 API.

Specification: Beacon Framework v2.1.1 (released December 13, 2024)
v2.1.1 changes from v2.0: bug fixes; no breaking changes; better VRS
alignment. Reference: https://genomebeacons.org/news/

GA4GH VRS integration (v2.0):
All variant identifiers use VRS computed digests (24-character hashes
computed from the canonical allele tuple). Computed via vrs_utils.py.
Reference: Wagner et al. 2021 Cell Genomics PMID:35072137
VRS spec: https://vrs.ga4gh.org/
"""
```

### Step 7.2 — VRS utilities (`beacon_api/vrs_utils.py`)

NEW FILE:

```python
"""
beacon_api.vrs_utils
=====================

GA4GH Variant Representation Standard (VRS) v2.0 utilities.

VRS provides canonical, globally unique identifiers for genomic variants.
The same variant expressed in different coordinate systems (VCF, HGVS,
SPDI) produces the same VRS identifier, enabling cross-system matching.

VRS ID computation:
    Input: RefSeq contig + start (interbase) + end + ref + alt
    Output: "ga4gh:VA." + 24-char digest (SHA-512 truncated)

Used in:
- Beacon v2.1.1 /g_variants response (variant_id field)
- FHIR R4 Task.input (reclassification re-contact resources)
- ClinVar submission (future v2 API support)

Reference:
    Wagner et al. 2021 Cell Genomics doi:10.1016/j.xgen.2021.100027
    PMID:35072137
    VRS Python library: pip install ga4gh.vrs>=2.0.0
"""
```

### Step 7.3 — Commit

```
git add -A
git commit -m "feat(beacon): GA4GH Beacon v2.1.1 with VRS v2.0 identifiers

- Beacon v2.1.1 spec (December 2024): no breaking changes from v2.0;
  improved VRS alignment documented in main.py
- vrs_utils.py (NEW): GA4GH VRS v2.0 canonical identifiers
  (Wagner 2021 PMID:35072137); SHA-512 digest computation; used in
  /g_variants response and FHIR Task.input
- /g_variants: query by position, rsID, gene, HGVSc, or VRS ID
- /info: v2.1.1 BeaconInfoResponse with dataset metadata
- GA4GH Passports authentication; open/controlled tier granularity
- PostgreSQL 16 variant store; async SQLAlchemy 2.0
- Tests: v2.1.1 response schema; VRS ID correctness; Passports JWT"

git push origin main
```

---

## PART 8 — PGx, MITO, REPEAT EXPANSION, PANGENOME

### Step 8.1 — Cyrius PGx (`pgx/cyrius_runner.py`)

Module docstring must explain why GATK4 cannot detect CYP2D6 SVs,
what star alleles are, and what each metaboliser phenotype means
clinically (NM, IM, PM, UM).

### Step 8.2 — Mito sub-workflow (`pipelines/mito.nf`)

Header comment block citing ACGS 2024 §6 with the four key requirements:

1. Haplogrep3 BEFORE ACMG classification
2. Heteroplasmy as % (not genotype)
3. Haplogroup-defining variant exclusion
4. MITOMAP as primary disease database

### Step 8.3 — ExpansionHunter sub-workflow (`pipelines/expansions.nf`)

Header comment block must state:

- Using ExpansionHunter v5.0 catalog (60 loci)
- TRGT is NOT used (requires PacBio HiFi; out of scope)
- Reference: Current Protocols 2024 doi:10.1002/cpz1.70010
- List all 60 loci by disease group in `modules/expansionhunter/LOCI.md`

### Step 8.4 — Pangenome modules

In all pangenome modules, pin vg container:

```
// Container: quay.io/biocontainers/vg:1.56.0
// HPRC graph: v1.1 minigraph-cactus (47 genomes, GRCh38 coordinates)
// Reference: Liao et al. 2023 Nature 617:312 PMID:37165242
//            Sirén et al. 2021 Science 374:abg8871 PMID:34818024
```

### Step 8.5 — Commit

```
git add -A
git commit -m "feat(specialist): PGx + mito (ACGS 2024 §6) + expansions v5.0 + pangenome

- Cyrius CYP2D6: explains GATK4 SV limitation; PharmVar+CPIC integration
- Mito sub-workflow: ACGS 2024 §6 compliant; Haplogrep3 first; mito BA1
  does NOT use standard threshold; MITOMAP as primary disease database
- ExpansionHunter v5.0 (60 loci): TRGT NOT used (short reads only;
  documented with reference); FGF14/SCA27B and RFC1/CANVAS added;
  LOCI.md documents all 60 loci with STRipy 2024 thresholds
- vg giraffe: container pinned to vg:1.56.0; Liao 2023 + Sirén 2021 cited"

git push origin main
```

---

## PART 9 — INFRASTRUCTURE & DEVOPS

### Step 9.1 — Dockerfiles (Ubuntu 24.04)

**`docker/Dockerfile.pipeline`:**

```dockerfile
# Base: Ubuntu 24.04 LTS (Noble Numbat, released April 2024)
# Chosen over 22.04: LTS support to 2029; GLIBC 2.39 (GATK 4.6.0.0);
# updated OpenSSL 3.x; updated system Python (3.12 compatible)
FROM ubuntu:24.04 AS builder

# GATK 4.6.0.0 — MUST use this version
# 4.6.0.0 includes: DRAGEN-GATK mode improvements, BQD genotyping model,
# DragSTR model for STR-aware calling, improved small cohort VQSR
ARG GATK_VERSION=4.6.0.0

# DRAGMAP — primary aligner (DRAGEN-GATK Best Practices)
ARG DRAGMAP_VERSION=1.3.0

# DeepVariant is pulled as a separate Docker image in Nextflow modules
# NOT installed here (size: ~10 GB). See modules/deepvariant/main.nf.
```

Every RUN instruction commented explaining purpose and any pitfalls.
Every ARG version pin commented explaining when/why to update.

### Step 9.2 — Terraform (AWS provider v5; PostgreSQL 16; db.t4g.micro)

**`terraform/main.tf`:**

```hcl
# AWS Provider v5 (hashicorp/aws ~> 5.0)
# BREAKING CHANGE from v4: aws_s3_bucket_acl removed; use
# aws_s3_bucket_ownership_controls instead.
# Migration guide: registry.terraform.io/.../version-5-upgrade
terraform {
  required_providers {
    aws = { source = "hashicorp/aws", version = "~> 5.0" }
  }
}
```

**`terraform/rds.tf`:**

```hcl
# PostgreSQL 16 on db.t4g.micro (ARM Graviton2)
# PostgreSQL 16 improvements relevant here:
#   - Better JSONB query performance (ACMG evidence stored as JSONB)
#   - Improved logical replication (for Beacon database read replicas)
# db.t4g.micro vs db.t3.micro:
#   - 30% cheaper (~$12/month vs ~$16/month)
#   - AWS Graviton2; better price-performance for PostgreSQL workloads
resource "aws_db_instance" "clarityseq" {
  engine         = "postgres"
  engine_version = "16"       # Must be 16; not 14 or 15
  instance_class = "db.t4g.micro"
  ...
}
```

### Step 9.3 — Commit

```
git add -A
git commit -m "feat(infra): Docker Ubuntu 24.04 + GATK 4.6.0.0 + Terraform AWS v5

- Dockerfile.pipeline: Ubuntu 24.04 LTS (not 22.04); GATK 4.6.0.0
  (version rationale commented: BQD + DragSTR + DRAGEN-GATK mode);
  DRAGMAP v1.3.0 compiled from source; every RUN commented
- Dockerfile.beacon + Dockerfile.daemon: Ubuntu 24.04; Python 3.12
- Singularity .def: converted from pipeline Docker; %help section
- Terraform AWS provider v5 (~>5.0): v4→v5 breaking changes documented
  inline (S3 ACL removal, ownership controls required)
- RDS: PostgreSQL 16 on db.t4g.micro (30% cheaper than t3.micro;
  Graviton2; JSONB improvements for ACMG evidence storage)
- ElastiCache: Redis 7.x; VPC private subnets
- Budget alert: $100/month (aws_budgets_budget)"

git push origin main
```

---

## PART 10 — DOCUMENTATION

### Step 10.1 — All 6 Architecture Decision Records

Write all six ADRs in `docs/adr/`. Each follows:
Status → Context → Decision → Consequences → Alternatives considered

**ADR-001: Why DRAGMAP over BWA-MEM2**

- Context: GATK Best Practices changed; Broad production WARP uses DRAGMAP
- Decision: DRAGMAP primary; BWA-MEM2 fallback
- Key consequence: BQSR must be SKIPPED in DRAGMAP/DRAGEN-GATK mode

**ADR-005: Why AlphaMissense as primary PP3/BP4 predictor**

- Context: ClinGen SVI 2024 approved 4 tools; AM is newest and best-performing
- Decision: AlphaMissense primary (≥0.564/≤0.340); REVEL retained as comparator
- Reference: Cheng et al. 2023 Science PMID:37703350

**ADR-006: Why PM2 weight = Supporting (1 pt) not Moderate (2 pts)**

- Context: ClinGen SVI 2024 guidance; gnomAD v4.1 reveals ultra-rare variants
- Decision: PM2 = Supporting by default; VCEP override via vcep_client.py
- Reference: ClinGen SVI PM2 guidance 2024; ACGS 2024 v1.2 §5

### Step 10.2 — ACMG classification guide (`docs/guides/acmg_classification.md`)

NEW FILE documenting the full ACGS 2024 + ClinGen SVI 2024 implementation:

```markdown
# ACMG Variant Classification in ClaritySeq

## ACGS 2024 v1.2 + ClinGen SVI 2024 Implementation

### Guidelines implemented

- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., Feb 2024)
- Richards et al. 2015 PMID:25741868
- Tavtigian et al. 2020 PMID:32645316
- Walker et al. 2023 PMID:36898414 (splicing subgroup)

### Key change: PM2 is Supporting (not Moderate)

...

### AlphaMissense as primary PP3/BP4 predictor

Thresholds: ≥0.564 → PP3; ≤0.340 → BP4
...

### Splicing evidence framework (Walker et al. 2023)

...

### VCEP gene-specific overrides

...

### VUS review date (ACGS 2024 §9)

...

### ClinVar submission (NHS mandatory — ACGS 2024 Introduction)

...
```

### Step 10.3 — JOSS paper draft (`joss_paper/paper.md`)

```markdown
---
title: 'BayesACMG: a Bayesian ACMG/AMP variant classifier implementing
ACGS 2024 v1.2 with calibrated uncertainty quantification'
tags:
  - Python
  - genomics
  - variant classification
  - Bayesian statistics
  - ACGS 2024
  - ClinGen SVI
authors:
  - name: [Your Name]
    orcid: 0000-0000-0000-0000
bibliography: paper.bib
---

# Summary

[150 words: BayesACMG implements ACGS 2024 v1.2 guidelines with Bayesian
uncertainty; PM2 at Supporting per ClinGen SVI 2024; AlphaMissense primary]

# Statement of Need

[200 words: fixed ACMG scoring has no uncertainty; VUS counselling problem;
no existing tool implements ACGS 2024 with calibrated Bayesian uncertainty]

# Implementation

[200 words: 28 rules + splicing subgroup + mito §6 + VCEP lookup;
PyMC NUTS sampler; ClinGen 500-variant calibration; ECE < 0.05]

# References
```

### Step 10.4 — Commit

```
git add -A
git commit -m "docs: full documentation suite

- 6 ADRs: DRAGMAP (ADR-001), PyMC (ADR-002), Celery (ADR-003),
  vg giraffe (ADR-004), AlphaMissense/PP3-BP4 (ADR-005),
  PM2-Supporting (ADR-006) — each with context/decision/consequences
- ACMG classification guide: full ACGS 2024 v1.2 implementation docs
  (PM2 change, AlphaMissense, splicing, VCEP, VUS review, ClinVar)
- ClinVar submission guide: NHS mandatory workflow documentation
- Quickstart: DRAGMAP reference build; gnomAD v4.1 download;
  AlphaMissense TSV download and indexing
- Sphinx conf.py: autodoc + napoleon + mermaid + furo
- JOSS paper draft: BayesACMG with ACGS 2024 context"

git push origin main
```

---

## PART 11 — FINAL INTEGRATION & QUALITY PASS

### Step 11.1 — Full integration test (`tests/integration/test_full_pipeline.py`)

12 acceptance criteria — all must pass:

```python
"""
Integration test: full ClaritySeq pipeline on GIAB HG001 chr22.

All 12 criteria must pass for a green release. Each assertion documents
which clinical requirement it validates.

1. PASS VCF produced (DRAGEN-GATK + DeepVariant INTERSECTION ensemble)
2. BayesACMG posteriors sum to 1.0 (probability conservation)
3. Known BRCA2 pathogenic variant classified P or LP (clinical accuracy)
4. AlphaMissense scores present for all missense PASS variants
5. All missense HGVSc/HGVSp use MANE Select transcripts
6. PM2 evidence weight = Supporting (1 pt) in all classifications
7. HTML report is valid HTML5 (html5lib validation)
8. PDF report > 50 KB
9. Audit trail JSON-LD validates; contains gnomAD_version='v4.1'
10. Mito haplogroup assigned in report
11. hap.py SNP sensitivity ≥ 99.0% vs GIAB HG001 chr22 truth
12. hap.py SNP precision ≥ 99.5%
"""
```

### Step 11.2 — Pre-commit final run

```bash
pre-commit run --all-files
# Fix EVERY issue. No # noqa or # type: ignore without an explanation comment.
```

### Step 11.3 — Coverage assertion

```bash
pytest \
  --cov=bayesacmg --cov=beacon_api --cov=reclassification \
  --cov=annotation --cov=reporting --cov=pgx \
  --cov=phenopackets_input --cov=multi_ancestry \
  --cov-report=term-missing --cov-fail-under=90
```

### Step 11.4 — Final commit and tag

```
git add -A
git commit -m "feat: ClaritySeq v0.1.0 complete — ACGS 2024 + DRAGEN-GATK + 7 novel contributions

## Guidelines and tools (verified current as of June 2026)

Clinical guidelines:
- ACGS Best Practice Guidelines 2024 v1.2 (Durkie et al., Feb 2024)
- ClinGen SVI 2024: PM2→Supporting; AlphaMissense PP3/BP4 approved
- Walker et al. 2023 PMID:36898414: splicing PP3/BP4/BP7 framework
- ACGS 2024 §6: mito-specific classification
- ACGS 2024 §9: VUS review date scheduling
- ACGS 2024 Introduction: NHS ClinVar submission mandatory

Variant calling:
- DRAGMAP v1.3.0 (DRAGEN-GATK Best Practices; BQSR skipped)
- BWA-MEM2 fallback (bwa_mem2 profile; BQSR enabled in this mode)
- GATK 4.6.0.0 (DRAGEN-GATK mode; BQD genotyping model)
- DeepVariant 1.8.0 + DeepTrio 1.8.0

Annotation:
- gnomAD v4.1 (April 2024; allele number bug fixed from v4.0)
- AlphaMissense 2023 (ClinGen SVI primary PP3/BP4)
- dbNSFP v4.7 (AlphaMissense + ESM1b + REVEL + BayesDel)
- VEP v111 (MANE Select pick order; Morales 2022 PMID:35356062)
- ExpansionHunter v5.0 catalog (60 loci; adds FGF14, RFC1 AAGGG)

Standards:
- GA4GH Beacon v2.1.1 (December 2024)
- GA4GH VRS v2.0 identifiers throughout
- FHIR Genomics Reporting IG v3.0.0
- Phenopackets v2.0 + phenopacket-tools validation

Infrastructure:
- Python 3.12; Ubuntu 24.04 LTS; PostgreSQL 16; Redis 7.x
- Terraform AWS provider v5; db.t4g.micro (ARM, 30% cheaper)

## 7 novel contributions
1. Three-way benchmark: GRCh38 vs T2T-CHM13 vs HPRC pangenome
   (HG001/HG002/HG005; first Nextflow implementation)
2. BayesACMG: Bayesian ACMG classifier (ACGS 2024 + ClinGen SVI 2024)
   PM2=Supporting; AlphaMissense primary; splicing framework; VCEP lookup;
   novel PVS1+PM2_Supporting=LP combination; calibrated against ClinGen
   500-variant set (ECE < 0.05); published to PyPI
3. ClinVar reclassification daemon (Celery + FHIR IG v3.0.0 + VRS v2.0)
   with NHS-mandatory ClinVar submission (ACGS 2024)
4. GA4GH Beacon v2.1.1 + VRS v2.0 identifiers federation API
5. Multi-ancestry VQSR with somalier ancestry inference
6. Cyrius CYP2D6 SV integrated into clinical WGS report
7. ExpansionHunter v5.0 (60 loci) + Mutect2 mito (ACGS 2024 §6)

## Quality
- pytest coverage ≥ 90%; 12 integration test criteria
- All ACMG rules cite PMID/DOI; real ClinVar RCV in test fixtures
- PM2 = Supporting (1 pt) verified in integration test
- MANE Select transcripts verified in integration test
- gnomAD v4.1 (not v4.0) throughout"

git tag -a v0.1.0 -m "ClaritySeq v0.1.0

ACGS 2024 v1.2 · ClinGen SVI 2024 · DRAGEN-GATK · DeepVariant 1.8.0
gnomAD v4.1 · AlphaMissense (PP3/BP4) · VEP v111 (MANE Select)
ExpansionHunter v5.0 · Beacon v2.1.1 · VRS v2.0 · FHIR IG v3.0.0
PM2=Supporting · NHS ClinVar submission"

git push origin main --tags
```

---

## PART 12 — POST-BUILD VERIFICATION

```bash
# 1. All GitHub Actions workflows pass
gh run list --limit 10

# 2. Coverage ≥ 90%
pytest --cov --cov-report=term-missing 2>&1 | grep "TOTAL"

# 3. Docker images build correctly
docker build -f docker/Dockerfile.pipeline \
  --build-arg GATK_VERSION=4.6.0.0 \
  --build-arg DRAGMAP_VERSION=1.3.0 \
  -t clarityseq-pipeline:v0.1.0 .

# 4. BayesACMG installs and CLI works
pip install -e bayesacmg/[dev]
bayesacmg classify --help

# 5. Beacon API starts and returns v2.1.1
uvicorn beacon_api.main:app --reload &
curl -s http://localhost:8000/info | python3 -m json.tool | grep "apiVersion"
# Must show: "apiVersion": "2.1.1"
kill %1

# 6. Nextflow config validates
nextflow config -profile test 2>&1 | head -20

# 7. No genomic data committed
git ls-files | grep -E "\.(fastq|bam|vcf|cram|gvcf)(.gz)?$" | wc -l
# Must be 0

# 8. All Python files have module docstrings
python3 -c "
import ast, sys, pathlib
missing = []
for f in pathlib.Path('.').rglob('*.py'):
    skip = ['.venv','__pycache__','.mypy_cache','build','dist']
    if any(s in str(f) for s in skip): continue
    try:
        tree = ast.parse(f.read_text())
        if not (tree.body and isinstance(tree.body[0], ast.Expr) and
                isinstance(tree.body[0].value, ast.Constant)):
            missing.append(str(f))
    except: pass
if missing:
    print('Missing module docstrings:')
    for m in missing: print(' ', m)
    sys.exit(1)
print('All Python files have module docstrings ✓')
"

# 9. DRAGMAP is the default aligner (not BWA-MEM2)
grep -c "dragmap" nextflow.config
# Must be > 0

# 10. DeepVariant version is pinned to 1.8.0 everywhere
grep -rn "deepvariant" modules/ | grep -v "1.8.0" | grep -v "#"
# Must return nothing

# 11. gnomAD v4.1 referenced (not v4.0) in annotation modules
grep -rn "v4\.0\|4\.0" annotation/ bayesacmg/ | grep -i gnomad
# Must return nothing

# 12. PM2 weight is Supporting (1 pt), not Moderate (2 pts)
python3 -c "
from bayesacmg.rules.pathogenic import rule_pm2
from bayesacmg.models import VariantInput, EvidenceStrength
# Create a minimal variant absent from gnomAD
# ... (fixture setup)
result = rule_pm2(absent_variant)
assert result.strength == EvidenceStrength.SUPPORTING, \
    f'PM2 must be Supporting per ClinGen SVI 2024, got {result.strength}'
print('PM2 = Supporting (1 pt) ✓')
"

# 13. AlphaMissense threshold is correct
grep -rn "0\.564\|0\.340" annotation/alphamissense_client.py
# Must find both thresholds documented

# 14. Beacon returns v2.1.1
python3 -c "
import httpx
r = httpx.get('http://localhost:8000/info')
info = r.json()
assert '2.1.1' in info.get('apiVersion', ''), \
    'Beacon must return v2.1.1 per spec'
print('Beacon v2.1.1 ✓')
"
```

---

## NON-NEGOTIABLE RULES FOR CLAUDE CODE

1. **Never skip a commit.** Every numbered commit section executes exactly as
   written. Conventional Commits format required (auto-generates changelog).

2. **Never write placeholder functions.** Every function is fully implemented.
   If an external tool is needed, write complete logic and mark:
   `# TODO(clarityseq): requires <tool> vX.Y installed`

3. **PM2 is ALWAYS Supporting (1 pt).** If you write PM2 as Moderate (2 pts)
   anywhere, it violates ClinGen SVI 2024. Add this assertion to any PM2 test.

4. **BQSR is absent in DRAGEN-GATK mode.** If BQSR appears in the default
   DRAGMAP pipeline, it is WRONG. It only belongs in `-profile bwa_mem2`.

5. **gnomAD v4.1 only.** Any reference to gnomAD v4.0 is wrong (allele
   number bug). Search for v4.0 strings before each commit and fix them.

6. **MANE Select everywhere.** All HGVSc/HGVSp in code, tests, reports,
   examples must use MANE Select transcripts (not Ensembl canonical).

7. **AlphaMissense thresholds are 0.564 (PP3) and 0.340 (BP4).** These are
   the ClinGen SVI 2024 approved thresholds. Any other values are wrong.

8. **Beacon v2.1.1 (not v2.0).** The apiVersion field in /info must return
   "2.1.1". This was released December 2024.

9. **Version pins are sacred.** GATK 4.6.0.0, DeepVariant 1.8.0, VEP 111,
   Python 3.12, Ubuntu 24.04. Never use `latest` in Docker images.

10. **PMID/DOI in every threshold.** Every numeric cutoff cites the paper.

11. **Real ClinVar RCV accessions in tests.** No invented variant fixtures.

12. **ClinVar submission is mandatory.** `clinvar_submitter.py` must be
    implemented — not stubbed. NHS England mandate per ACGS 2024.

---

_ClaritySeq Claude Code Prompt — FINAL VERSION — June 2026_
_All tools verified against: ACGS 2024 v1.2 (Durkie et al., Feb 2024) ·_
_ClinGen SVI 2024 guidelines · GATK Best Practices (DRAGEN-GATK) ·_
_DeepVariant 1.8.0 release notes · gnomAD v4.1 release notes (April 2024) ·_
_Cheng et al. 2023 Science (AlphaMissense) · Morales 2022 Nature (MANE Select) ·_
_Walker et al. 2023 AJHG (splicing subgroup) · Beacon v2.1.1 spec (Dec 2024) ·_
_GA4GH VRS v2.0 · FHIR Genomics Reporting IG v3.0.0 · nf-core/sarek v3.8.1_
