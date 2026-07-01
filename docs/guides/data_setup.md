# Data Setup Guide

Download and prepare all required reference data before running ClaritySeq.

## 1. DRAGMAP hash table (PRIMARY reference)

DRAGMAP requires a pre-built hash table (not a FASTA file directly).

```bash
# Option A: Download pre-built from Broad (recommended)
aws s3 sync s3://broad-references/hg38/v0/dragmap/ /data/dragmap/

# Option B: Build from scratch (~4 hours; requires 64 GB RAM)
dragen-os --build-hash-table /data/GRCh38.fa --output-directory /data/dragmap/
```

Set in `nextflow.config`:
```groovy
params.dragmap_reference = "/data/dragmap"
```

## 2. GRCh38 FASTA (for BWA-MEM2 fallback and VEP)

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/GRCh38.fa.gz
```

## 3. gnomAD v4.1 (April 19, 2024)

**IMPORTANT**: Use v4.1, NOT v4.0. gnomAD v4.0 had an allele number calculation bug.

```bash
# Sites VCF only (do not download the full callset)
gsutil cp gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr*.vcf.bgz /data/gnomad/
```

## 4. AlphaMissense scores (~2.7 GB)

```bash
gsutil cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz /data/
# Tabix index (required for fast lookup)
tabix -s 1 -b 2 -e 2 /data/AlphaMissense_hg38.tsv.gz
```

Set in `.env`:
```
ALPHAMISSENSE_SCORES_PATH=/data/AlphaMissense_hg38.tsv.gz
```

## 5. VEP v111 cache

```bash
# Download cache (~15 GB)
vep_install -a cf -s homo_sapiens -y GRCh38 --CACHEDIR /data/vep_cache --CACHE_VERSION 111

# Download dbNSFP v4.7 plugin data
# (see VEP plugin documentation for dbNSFP setup)
```

## 6. ExpansionHunter catalog (v5.0, 60 loci)

The catalog is bundled in the ExpansionHunter Docker image. For custom additions, see `modules/expansionhunter/LOCI.md`.

## 7. HPRC pangenome graph (optional; for -profile pangenome)

```bash
# HPRC v1.1 minigraph-cactus graph (47 genomes)
aws s3 sync s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/ /data/hprc/
```

## Storage requirements

| Dataset | Size |
|---------|------|
| DRAGMAP hash table | ~60 GB |
| GRCh38 FASTA | ~3 GB |
| gnomAD v4.1 sites | ~400 GB |
| AlphaMissense | ~2.7 GB |
| VEP cache | ~15 GB |
| HPRC pangenome | ~50 GB |
| **Total** | **~530 GB** |

AWS S3 cost (eu-west-2): ~$12/month for 530 GB
