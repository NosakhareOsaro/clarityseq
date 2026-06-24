#!/usr/bin/env nextflow
// ============================================================================
// GenomeForge — HPRC Pangenome Pipeline Arm
// ============================================================================
// Aligns short Illumina reads to the HPRC v1.1 minigraph-cactus pangenome
// graph using vg giraffe, outputs a GRCh38-coordinate BAM, and calls variants
// with GATK4 HaplotypeCaller + DeepVariant ensemble.
//
// IMPLEMENTATION STATUS: Partial — alignment step implemented; variant calling
// and benchmarking steps are stubs (see TODO comments below).
//
// RESOURCE REQUIREMENTS: ≥ 32 CPUs, ≥ 256 GB RAM
//   (HPRC v1.1 full GBZ graph ~22 GB; GBWT index ~20 GB; dist index ~4 GB)
//   Minimum disk: 500 GB for graph indexes + intermediate files
//
// ENABLE WITH: -profile pangenome  (sets params.run_pangenome = true)
// CI: Nightly only (see .github/workflows/ci_pangenome.yml)
//     Caches HPRC chr22 graph for 7 days to limit S3 egress costs.
//
// PANGENOME GRAPH: HPRC v1.1 minigraph-cactus (47 phased assemblies)
//   - 94 haplotypes from diverse global ancestries (AFR/AMR/EAS/EUR/SAS/OCE/MID)
//   - ~100 million variants vs GRCh38, including ~10 million SVs
//   - Download: https://github.com/human-pangenomics/hpp_pangenome_resources
//   Container: quay.io/biocontainers/vg:1.56.0
//
// vg GIRAFFE ALGORITHM:
//   Minimizer-seeded, GBWT-haplotype-guided graph alignment.
//   10–60× faster than vg map with equivalent accuracy.
//   Reference: Sirén et al. 2021 Science PMID:34818024
//
// PANGENOME BENEFIT:
//   Liao et al. 2023 Nature PMID:37165242 demonstrate:
//     - ~292 Mb of novel pangenome sequence not in GRCh38
//     - Reference bias reduced by 24% in non-EUR populations
//     - Improved SV calling sensitivity (novel insertions vs GRCh38)
//
// COMPARISON WITH GRCh38 ARM:
//   Nightly CI (ci_pangenome.yml) benchmarks pangenome SNP sensitivity vs
//   GRCh38 arm on GIAB HG001 chr22. Acceptance criterion: pangenome
//   sensitivity ≥ GRCh38 sensitivity (pangenome must not regress).
//
// REFERENCES:
//   Liao et al. 2023 Nature 617:312 PMID:37165242 (HPRC pangenome)
//   Sirén et al. 2021 Science 374:abg8871 PMID:34818024 (vg giraffe)
//   Eggertsson et al. 2017 Nat Genet PMID:28945250 (vg graph genome)
// ============================================================================

nextflow.enable.dsl = 2

include { VG_GIRAFFE_ALIGN   } from '../modules/vg_giraffe/main'
include { GATK4_MARKDUPLICATES } from '../modules/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX     } from '../modules/samtools/main'
include { MOSDEPTH_QC        } from '../modules/mosdepth/main'

workflow WGS_PANGENOME {
    take:
        ch_reads        // Channel: tuple(sample_id, fastq_1, fastq_2)
                        // Post-QC reads from FASTP_QC in the main pipeline
        graph_gbz       // Path: HPRC v1.1 GBZ graph + GBWT haplotype index (~22 GB)
        graph_dist      // Path: Distance index (.dist) for giraffe seeding (~4 GB)
        graph_min       // Path: Minimizer index (.min) for fast lookup (~20 GB)
        reference       // Path: GRCh38 FASTA + .fai (for BAM coordinate projection)
        reference_fai   // Path: GRCh38 .fai index

    main:
        // ── Step 1: Pangenome graph alignment (vg giraffe) ─────────────────
        // Aligns reads to HPRC v1.1 minigraph-cactus pangenome graph.
        // Output: coordinate-sorted BAM in GRCh38 linear coordinates.
        // Compatible with all downstream GATK4 and DeepVariant modules.
        // Container: quay.io/biocontainers/vg:1.56.0
        // Reference: Liao et al. 2023 PMID:37165242; Sirén et al. 2021 PMID:34818024
        VG_GIRAFFE_ALIGN(
            ch_reads,
            graph_gbz,
            graph_dist,
            graph_min,
            reference,
            reference_fai
        )

        // ── Step 2: Mark duplicates ────────────────────────────────────────
        // Same module as GRCh38 arm — BAM is now in GRCh38 coordinates.
        GATK4_MARKDUPLICATES(VG_GIRAFFE_ALIGN.out.bam)

        // ── Step 3: Index and QC ───────────────────────────────────────────
        SAMTOOLS_INDEX(GATK4_MARKDUPLICATES.out.bam)
        MOSDEPTH_QC(
            SAMTOOLS_INDEX.out.bam_bai,
            file("NO_FILE")  // No coverage BED — genome-wide
        )

        // ── TODO: Variant calling + ensemble merge (v2.0) ──────────────────
        // STEP 4: GATK4_HAPLOTYPECALLER on pangenome-aligned BAM
        // STEP 5: DEEPVARIANT_CALL (WGS model; pangenome BAM)
        // STEP 6: ENSEMBLE_MERGE (same INTERSECTION mode as GRCh38 arm)
        // STEP 7: VEP_ANNOTATE + ALPHAMISSENSE_LOOKUP
        // Tracked in: https://github.com/genomeforge/genomeforge/issues/PANGENOME-CALLING
        log.info "WGS_PANGENOME: alignment complete. Variant calling stub — see TODO."

    emit:
        // Pangenome-aligned, duplicate-marked BAM in GRCh38 coordinates
        bam            = GATK4_MARKDUPLICATES.out.bam
        bai            = SAMTOOLS_INDEX.out.bai
        // Stub: VCF channel empty until variant calling is implemented
        vcf            = Channel.empty()
        // Mosdepth coverage summary (pangenome alignment quality check)
        coverage       = MOSDEPTH_QC.out.summary
}
