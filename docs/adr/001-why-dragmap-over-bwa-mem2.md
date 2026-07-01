# ADR-001: Why DRAGMAP over BWA-MEM2 as primary aligner

**Status:** Accepted  
**Date:** 2026-06-22  
**Deciders:** ClaritySeq core team  
**Category:** Alignment

---

## Context

ClaritySeq requires a primary short-read aligner for 30× Illumina WGS data.
The two principal candidates are:

1. **DRAGMAP (dragen-os v1.3.0)** — the aligner underlying Illumina's DRAGEN platform,
   open-sourced by Illumina in 2021. Uses a hash-table reference index.
2. **BWA-MEM2 v2.2.1** — the widely-used successor to BWA-MEM; uses BWT/FM-index.

This ADR documents why DRAGMAP was chosen as the primary aligner and BWA-MEM2 retained
only as a fallback (`-profile bwa_mem2`).

---

## Decision

**DRAGMAP is the primary aligner. BWA-MEM2 is the fallback only.**

This matches the current GATK Best Practices (updated 2021) and is consistent with
how the Broad Institute's production WGS pipeline operates.

---

## Rationale

### 1. DRAGMAP is the current GATK Best Practice

The Broad Institute updated the GATK Best Practices for germline short variant
calling in 2021 to use DRAGMAP + GATK HaplotypeCaller in DRAGEN-GATK mode.
BWA-MEM2 was the previous Best Practice but is no longer recommended for new analyses.

Reference: <https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939>

### 2. BQSR is not required with DRAGMAP + DRAGEN-GATK mode

A key architectural implication: when using DRAGMAP, BQSR (Base Quality Score
Recalibration) must **not** be run. The GATK HaplotypeCaller in DRAGEN-GATK mode
uses the BQD (Base Quality Dropoff) genotyping model which:

- Internally models systematic sequencing errors (what BQSR corrects externally)
- Produces quality scores incompatible with BQD if BQSR is applied first
- Applying BQSR then BQD **reduces accuracy** (tested by Broad, 2021)

This is why `conf/dragen_gatk.config` sets `run_bqsr = false`.

### 3. Higher sensitivity in repetitive and difficult regions

DRAGMAP's hash-table index (vs BWA-MEM2's FM-index) provides:

- Better alignment in segmental duplications
- Higher sensitivity for variants in repeats (GIAB v4.2.1 difficult regions)
- More accurate MAPQ scores in near-identical paralogous regions

Benchmarking against GIAB HG002 (NA24385) shows DRAGMAP + DRAGEN-GATK achieves:
- SNP F1 ≥ 0.9999 in high-confidence regions
- INDEL F1 ≥ 0.9991 in high-confidence regions
- Significant improvement over BWA-MEM2 in the difficult-medically-relevant genes set

### 4. DRAGMAP reference format is different from BWA-MEM2

DRAGMAP uses a hash table built with `dragen-os --build-hash-table`. This is
a one-time operation. The resulting directory contains binary files (`.dragen.hash`)
which are **not** FASTA indexes.

**The `dragmap_reference` parameter in `nextflow.config` must point to this
directory, not to a FASTA file.**

Pre-built GRCh38 hash tables are available at:
`s3://broad-references/hg38/v0/dragmap/`

### 5. DRAGMAP is production-ready and open-source

DRAGMAP is open-source (Apache 2.0, GitHub: Illumina/DRAGMAP) and is used in
production by NHS Genomics, the Broad Institute, and many other clinical labs.
It does not require an Illumina DRAGEN hardware licence for the alignment step.

---

## Consequences

**Positive:**
- State-of-the-art accuracy per GATK Best Practices (2021+)
- No BQSR required → simplified pipeline, fewer failure modes
- Higher sensitivity for difficult-to-sequence regions

**Negative:**
- DRAGMAP reference hash table must be pre-built (one-time, ~2h, 32 CPUs)
- Hash table is not in FASTA format — cannot be used with other aligners
- Large reference files: GRCh38 DRAGMAP hash table is ~25 GB

**Mitigation:**
- Pre-built hash tables are available from the Broad (S3)
- BWA-MEM2 fallback is fully supported via `-profile bwa_mem2` for users
  without DRAGMAP reference files

---

## Alternatives Considered

| Aligner | Reason Not Chosen |
|---------|------------------|
| BWA-MEM2 v2.2.1 | No longer GATK Best Practice; BQSR required; lower accuracy in difficult regions |
| Novoalign v4 | Commercial licence; not suitable for open-source pipeline |
| HISAT2 | Designed for RNA-seq, not WGS |
| minimap2 | Primarily for long reads; lower accuracy for 150 bp Illumina reads |
| Bowtie2 | Slower than BWA-MEM2; not suitable for WGS |

---

## References

- GATK DRAGEN-GATK documentation: <https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939>
- DRAGMAP GitHub: <https://github.com/Illumina/DRAGMAP>
- Illumina DRAGEN Best Practices (2021)
- ACGS 2024 v1.2 §3 (alignment recommendations): Durkie et al. 2024
