# ADR-004: Why vg giraffe over vg map for pangenome alignment

**Status:** Accepted  
**Date:** 2026-06-22  
**Deciders:** GenomeForge core team  
**Category:** Pangenome alignment

---

## Context

GenomeForge includes a pangenome alignment arm (`run_pangenome=true`) that aligns
Illumina WGS reads to the HPRC v1.1 pangenome graph. Two `vg` alignment algorithms
were evaluated:

1. **vg giraffe** — short-read aligner specifically designed for variation graphs;
   uses GBWT-based haplotype-aware alignment
2. **vg map** — the original vg aligner; general-purpose but slower

---

## Decision

**vg giraffe is the pangenome aligner for GenomeForge.**

---

## Rationale

### 1. vg giraffe is 4–12× faster than vg map

Sirén et al. 2021 (Nature Biotechnology, PMID:34385711) demonstrated:

- vg giraffe: ~600 million reads/hour on 48 CPU threads
- vg map: ~50 million reads/hour on 48 CPU threads

For 30× WGS (~900 million read pairs), this difference is critical:
- vg giraffe: ~3-4 hours
- vg map: ~18-36 hours (prohibitive for clinical use)

### 2. vg giraffe is more accurate for short reads

vg giraffe uses:
- **Minimizer index** (`.min` file) for fast seed finding
- **GBWT (Generalized BWT)** haplotype-aware graph traversal
- **Adaptive seed chaining** optimised for 150 bp reads

vg map uses the full variation graph GCSA2 index which is optimised for longer reads.
For 150 bp Illumina reads, giraffe consistently outperforms vg map in:
- Correctly mapped read rate
- Variant calling F1 (downstream GATK HC)
- Mapping quality calibration

Reference: Sirén et al. 2021 Nature Biotechnology doi:10.1038/s41587-021-00865-7

### 3. vg giraffe is the recommended algorithm for HPRC v1.1

The Human Pangenome Reference Consortium (HPRC) provides pre-built giraffe indexes
(`.gbz`, `.dist`, `.min`) optimised for the HPRC v1.1 graph. The HPRC documentation
explicitly recommends giraffe for short-read alignment.

These indexes are incompatible with vg map (which requires GCSA2 + XG indexes).
Pre-built giraffe indexes are available at:
`https://github.com/human-pangenomics/hpp_pangenome_resources`

### 4. vg giraffe outputs surjectable alignments

vg giraffe produces GAM (Graph Alignment/Map) files that can be surjected to a
linear GRCh38 BAM using `vg surject`. This is required for:
- GATK HaplotypeCaller (linear BAM required)
- Coverage analysis with mosdepth
- QC with samtools flagstat

### 5. Resource requirements

vg giraffe's minimum resource requirements for HPRC v1.1:
- **Memory**: ~100-200 GB RAM (the graph index is loaded into memory)
- **CPUs**: 32+ recommended (scales linearly to 128 CPUs)

These requirements are why `run_pangenome=false` by default and the `pangenome`
profile sets `max_memory = "256.GB"`.

---

## Consequences

**Positive:**
- 4–12× faster than vg map for 150 bp reads
- Better accuracy for short reads
- Native HPRC v1.1 index support
- Outputs can be surjected to GRCh38 for downstream tools

**Negative:**
- High memory requirements (~200 GB for HPRC v1.1 graph)
- Not suitable for smaller genomes or non-HPRC graphs
- giraffe indexes are specific to the pangenome build

**Mitigation:**
- `run_pangenome=false` by default (most clinical runs do not need pangenome)
- `-profile pangenome` enforces minimum resource requirements
- HPRC v1.1 indexes can be downloaded in ~2 hours on good bandwidth

---

## Alternatives Considered

| Aligner | Reason Not Chosen |
|---------|------------------|
| vg map | 4–12× slower; designed for long reads; GCSA2 indexes not provided by HPRC |
| vg mpmap | Multi-path aligner; more complex output; not needed for WGS |
| Minigraph | Does not use SNP variation; not suitable for variant-rich graphs |
| GraphAligner | No production support for HPRC graphs |
| PanAligner | Early-stage; not yet validated for clinical WGS |

---

## References

- Sirén et al. 2021: vg giraffe paper (PMID:34385711, Nature Biotechnology)
  doi:10.1038/s41587-021-00865-7
- HPRC v1.1 resources: <https://github.com/human-pangenomics/hpp_pangenome_resources>
- vg documentation: <https://github.com/vgteam/vg>
