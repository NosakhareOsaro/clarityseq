# ClaritySeq — Prerequisites

This document describes all tools required to run ClaritySeq and how to install them.
Run `make check-prereqs` to verify your environment.

---

## Required Tools

All tools listed below are **mandatory**. The pipeline will fail to start if any are missing.

### git ≥ 2.40

```bash
git --version
# Expected: git version 2.40.0 or later
```

**Install:**
- macOS: `brew install git`
- Ubuntu: `sudo apt-get install git`
- Official: <https://git-scm.com/downloads>

Why 2.40: ClaritySeq uses `git tag --sort=-creatordate` and `git log --follow` features
introduced in git 2.40.

---

### GitHub CLI (`gh`) — authenticated

```bash
gh auth status
# Expected: Logged in to github.com as <username>
```

**Install:**
- macOS: `brew install gh`
- Ubuntu: See <https://cli.github.com/>

**Authenticate:**
```bash
gh auth login
```

---

### Nextflow ≥ 24.x (DSL2)

```bash
nextflow -version
# Expected: nextflow version 24.04.x or later
```

**Install:**
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

DSL2 is enabled by `nextflow.enable.dsl = 2` in `nextflow.config`. Nextflow < 24.x
does not support all DSL2 features used by ClaritySeq.

---

### Docker ≥ 25.x

```bash
docker --version
# Expected: Docker version 25.x.x or later

docker run hello-world
# Expected: "Hello from Docker!" message
```

**Install:** <https://docs.docker.com/engine/install/>

**Post-install (Linux):** Add your user to the `docker` group:
```bash
sudo usermod -aG docker $USER
newgrp docker
```

---

### Python 3.12.x (exactly — not 3.11, not 3.13)

```bash
python3 --version
# Expected: Python 3.12.x

python3 -c "import sys; print(sys.version_info[:2] == (3, 12))"
# Expected: True
```

**Install:**
- macOS: `brew install python@3.12`
- Ubuntu: `sudo apt-get install python3.12 python3.12-venv python3.12-dev`
- pyenv: `pyenv install 3.12.3 && pyenv global 3.12.3`

**Why 3.12 exactly?**
- Python 3.12 introduced `sys.monitoring` (used by pytest-cov for accurate coverage)
- mypy strict mode has known false positives on 3.11 that are fixed in 3.12
- BayesACMG uses Python 3.12 pattern matching (`match/case`) for ACMG rule dispatch
- 3.13 is not yet stable for PyMC (the Bayesian model engine)

---

### conda or mamba

```bash
conda --version
# Expected: conda 24.x.x or later
# OR
mamba --version
```

**Install:**
- Miniforge (recommended): <https://github.com/conda-forge/miniforge>
- Miniconda: <https://docs.conda.io/en/latest/miniconda.html>

conda/mamba is used for tool environment management (VEP, GATK, etc.)
when not using Docker.

---

### Terraform ≥ 1.7

```bash
terraform --version
# Expected: Terraform v1.7.x or later
```

**Install:**
```bash
# macOS
brew tap hashicorp/tap
brew install hashicorp/tap/terraform

# Ubuntu
wget -O- https://apt.releases.hashicorp.com/gpg | sudo gpg --dearmor -o /usr/share/keyrings/hashicorp-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/hashicorp-archive-keyring.gpg] https://apt.releases.hashicorp.com $(lsb_release -cs) main" | sudo tee /etc/apt/sources.list.d/hashicorp.list
sudo apt update && sudo apt install terraform
```

Terraform 1.7+ is required for the AWS provider v5 features used in `terraform/`.

---

### Java ≥ 17 (required by GATK 4.6.0.0)

```bash
java -version
# Expected: openjdk version "17.x.x" or later

# GATK version check
gatk --version
# Expected: 4.6.0.0
```

**Install:**
- macOS: `brew install openjdk@17`
- Ubuntu: `sudo apt-get install openjdk-17-jdk`

GATK 4.6.0.0 requires Java 17. Java 11 will cause a runtime error.

---

## Optional Tools

These tools are not required for the standard pipeline but are needed for specific arms
or development workflows.

### Singularity or Apptainer ≥ 1.2

Required for `-profile hpc` (SLURM/HPC environments where Docker is not permitted).

```bash
singularity --version
# Expected: singularity-ce version 1.2.x or later
# OR
apptainer --version
```

**Install:**
- Apptainer (successor to Singularity): <https://apptainer.org/docs/user/main/quick_start.html>
- Most HPC systems have Singularity pre-installed; contact your HPC admin.

---

### Node.js ≥ 20 (LTS)

Required for Mermaid diagram generation (`make figures`).

```bash
node --version
# Expected: v20.x.x or later
```

**Install:**
- macOS: `brew install node@20`
- nvm: `nvm install 20 && nvm use 20`

After installing Node, install Mermaid CLI:
```bash
npm install -g @mermaid-js/mermaid-cli
```

---

## Bioinformatics Tool Versions

These tools are included in the Docker containers (`docker/Dockerfile.pipeline`).
You do not need to install them manually when using Docker profiles.

| Tool | Version | Purpose |
|------|---------|---------|
| DRAGMAP (dragen-os) | 1.3.0 | Primary aligner (DRAGEN-GATK) |
| BWA-MEM2 | 2.2.1 | Fallback aligner (-profile bwa_mem2 only) |
| GATK | 4.6.0.0 | HaplotypeCaller, MarkDuplicates, VQSR |
| DeepVariant | 1.8.0 | Parallel DL-based variant calling |
| VEP | 111 | Variant annotation (MANE Select) |
| samtools | 1.19 | BAM manipulation |
| fastp | 0.23.4 | FASTQ QC and adapter trimming |
| mosdepth | 0.3.6 | Coverage analysis |
| ExpansionHunter | 5.0 | STR expansion detection (60 loci) |
| Cyrius | 1.1.1 | CYP2D6 SV calling |
| Haplogrep3 | 3.2.1 | Mito haplogroup assignment |
| vg | Latest | Pangenome alignment (giraffe) |
| somalier | 0.2.19 | Ancestry inference and relatedness |

---

## Data Downloads

These large files must be downloaded manually before running the pipeline.
They are excluded from git (see `.gitignore`).

### GRCh38 FASTA

```bash
# GRCh38.p14 FASTA from Broad Google Cloud bucket
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta /data/
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai /data/
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict /data/
```

### DRAGMAP Hash Table

```bash
# Pre-built hg38 DRAGMAP hash table from Broad
gsutil -m cp -r gs://broad-references/hg38/v0/dragmap/ /data/dragmap_hg38/
# OR build from FASTA (takes ~2 hours on 32 CPUs):
# dragen-os --build-hash-table /data/dragmap_hg38/ --ht-reference /data/Homo_sapiens_assembly38.fasta
```

### gnomAD v4.1 (April 2024)

```bash
# gnomAD v4.1 genome VCF — use v4.1 specifically (see nextflow.config comment)
gsutil cp gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr*.vcf.bgz /data/gnomad_v4.1/
```

### AlphaMissense Scores

```bash
# AlphaMissense scores (Cheng et al. 2023, Science)
# ClinGen SVI approved for PP3/BP4 use in 2024
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz /data/
tabix -s 1 -b 2 -e 2 /data/AlphaMissense_hg38.tsv.gz
# File size: ~2.7 GB compressed
```

### VEP v111 Cache

```bash
# Install VEP v111 cache (GRCh38)
vep_install --NO_HTSLIB \
    -a c \
    -s homo_sapiens \
    -y GRCh38 \
    -c $HOME/.vep \
    --CACHE_VERSION 111
# Size: ~15 GB; download time: 30-60 minutes
```

---

## Environment Setup

```bash
# 1. Clone the repository
git clone https://github.com/clarityseq/clarityseq.git
cd clarityseq

# 2. Copy and fill environment variables
cp .env.example .env
# Edit .env with your values

# 3. Install Python dependencies
make install

# 4. Verify all prerequisites
make check-prereqs

# 5. Run the chr22 CI test
make ci-test
```

---

## NHS England Lab Requirements

For NHS England Genomic Medicine Service accreditation under ACGS 2024 v1.2:

- All data must remain in `eu-west-2` (AWS London) — do not use other regions
- ClinVar submission credentials (`NCBI_API_KEY`, `CLINVAR_SUBMISSION_ORG_ID`) are mandatory
- The reclassification daemon must be running continuously (see `make daemon-start`)
- Clinical reports must be generated in PDF format (`report_format = "pdf"`)
- Minimum sequencing coverage: 30× mean genome (enforced by `min_coverage = 30` in nextflow.config)
