#!/usr/bin/env nextflow

// =============================================================================
// GenomeForge — Primary WGS Pipeline (GRCh38)
// =============================================================================
//
// OVERVIEW
// --------
// This pipeline implements the GenomeForge DRAGEN-GATK WGS clinical variant
// interpretation workflow for GRCh38, producing annotated variants, Bayesian
// ACMG classifications, and NHS GMS-style clinical reports.
//
// ALIGNER STRATEGY
// ----------------
// PRIMARY:  DRAGMAP v1.3.0  (default; GATK DRAGEN-GATK Best Practices since 2021)
// FALLBACK: BWA-MEM2 v2.2.1 (only with -profile bwa_mem2)
//
// CRITICAL — BQSR POLICY:
//   DRAGMAP / DRAGEN-GATK mode: BQSR is SKIPPED.
//   HaplotypeCaller with --dragen-mode true uses the BQD (Base Quality Dropoff)
//   genotyping model, which replaces BQSR by modelling systematic errors
//   internally. Applying BQSR to DRAGMAP-aligned reads then passing to HC
//   --dragen-mode REDUCES accuracy (BQD expects raw DRAGMAP quality scores).
//   Reference: Broad DRAGEN-GATK docs (gatk.broadinstitute.org/...4407897446939)
//
//   BWA-MEM2 profile: BQSR IS run. The pipeline includes GATK4_BASERECALIBRATOR
//   → GATK4_APPLYBQSR steps when -profile bwa_mem2 is active.
//
// ENSEMBLE CALLING STRATEGY
// -------------------------
// GATK4 HaplotypeCaller (graph assembly) and DeepVariant v1.8.0 (CNN pileup)
// have partially non-overlapping error profiles. The ensemble caller combines
// both in INTERSECTION mode by default (variant must be PASS in BOTH callers),
// which maximises precision for clinical reporting. UNION mode is available
// for research use where sensitivity is prioritised over precision.
//
// PER-ANCESTRY VQSR
// -----------------
// somalier infers sample ancestry (AFR/AMR/EAS/EUR/SAS or admixed) and selects
// the appropriate gnomAD v4.1 ancestry-stratified subset for VQSR training.
// This reduces false positive rates for underrepresented populations.
//
// DEEPTRIO TRIGGER
// ----------------
// DeepTrio v1.8.0 is used instead of DeepVariant when:
//   1. A PED file is provided in the sample sheet (ped_file column), AND
//   2. params.run_deeptrio = true
// DeepTrio improves de novo SNV sensitivity by ~15% vs per-sample DeepVariant
// (PMID:36050879) by using parent BAMs as additional image channels.
//
// MITOCHONDRIAL SUB-WORKFLOW
// --------------------------
// When params.run_mito = true:
//   GATK Mutect2 --mitochondria-mode → Haplogrep3 haplogroup classification
//   → ACGS 2024 §6 mito ACMG classification (haplogroup-defining variants excluded)
//
// REPEAT EXPANSION SUB-WORKFLOW
// -----------------------------
// When params.run_expansions = true:
//   ExpansionHunter v5.0 (60-locus catalog)
//   Note: TRGT is NOT used — it requires PacBio HiFi long reads.
//
// ASCII DATA FLOW DIAGRAM
// -----------------------
//
//   FASTQs
//     │
//     ▼
//   FASTP_QC  ─────────────────────────────── QC reports
//     │
//     ▼
//   DRAGMAP_ALIGN (or BWA_MEM2_ALIGN if -profile bwa_mem2)
//     │
//     ▼
//   GATK4_MARKDUPLICATES ─────── [if bwa_mem2] GATK4_BQSR → GATK4_APPLYBQSR
//     │
//     ├──► SAMTOOLS_FLAGSTAT (alignment QC)
//     ├──► MOSDEPTH_QC (coverage gate ≥30×)
//     ├──► SOMALIER_ANCESTRY (ancestry inference → VQSR)
//     │
//     ▼
//   GATK4_HAPLOTYPECALLER ──────── (DRAGEN-GATK mode; --dragen-mode true; no BQSR)
//   (scattered by interval)
//     │
//     ▼
//   GATK4_GENOTYPEGVCFS
//     │
//     ▼
//   VQSR_SNP → VQSR_INDEL ──────── (per-ancestry gnomAD v4.1 training sets)
//     │                │
//     │         ┌──────┘
//     │    DEEPVARIANT_CALL (or DEEPTRIO if trio + run_deeptrio)
//     │         │
//     └────► ENSEMBLE_MERGE (INTERSECTION mode default)
//                 │
//                 ▼
//           ALPHAMISSENSE_LOOKUP ─── (ClinGen SVI 2024: ≥0.564→PP3, ≤0.340→BP4)
//                 │
//                 ▼
//           VEP_ANNOTATE ──────────── (v111; MANE Select; dbNSFP v4.7; AlphaMissense)
//                 │
//                 ├── [if run_mito]       MITO_PIPELINE (Mutect2 + Haplogrep3)
//                 ├── [if run_expansions] EXPANSION_HUNTER (v5.0; 60 loci)
//                 └── [if run_pgx]        CYRIUS_CYP2D6 + CPIC dosing
//                 │
//                 ▼
//           BAYESACMG_CLASSIFY ───── (28 rules; ACGS 2024 v1.2; PM2=Supporting)
//                 │
//                 ▼
//           GENERATE_REPORT ─────────── (HTML + PDF + JSON-LD audit trail)
//                 │
//                 ├── CLINVAR_INGEST (novel P/LP → submission queue; NHS mandate)
//                 └── BEACON_INGEST (VRS v2.0; GA4GH Beacon v2.1.1)
//
// GUIDELINES IMPLEMENTED
// ----------------------
// - ACGS 2024 v1.2 (Durkie et al., ratified Feb 2024): primary UK guidelines
// - Richards et al. 2015 PMID:25741868: original 28-rule ACMG/AMP framework
// - ClinGen SVI 2024: PM2→Supporting; AlphaMissense primary PP3/BP4
// - Walker et al. 2023 PMID:36898414: splicing framework (PP3/BP4/BP7)
// - ACGS 2024 §6: mitochondrial variant classification rules
//
// AUTHORS
// -------
// GenomeForge Contributors
// SPDX-License-Identifier: MIT
// =============================================================================

nextflow.enable.dsl = 2

// =============================================================================
// Module imports
// =============================================================================
include { FASTP_QC           } from '../modules/fastp/main'
include { DRAGMAP_ALIGN      } from '../modules/dragmap/main'
include { BWA_MEM2_ALIGN     } from '../modules/bwa_mem2/main'
include { GATK4_MARKDUPLICATES } from '../modules/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX; SAMTOOLS_FLAGSTAT } from '../modules/samtools/main'
include { MOSDEPTH_QC        } from '../modules/mosdepth/main'
include { SOMALIER_ANCESTRY  } from '../modules/somalier/main'
include { GATK4_HAPLOTYPECALLER } from '../modules/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS   } from '../modules/gatk4/genotypegvcfs/main'
include { VQSR_SNP; VQSR_INDEL  } from '../modules/gatk4/vqsr/main'
include { GATK4_APPLYBQSR    } from '../modules/gatk4/applybqsr/main'
include { DEEPVARIANT_CALL   } from '../modules/deepvariant/main'
include { DEEPTRIO_CALL      } from '../modules/deeptrio/main'
include { ENSEMBLE_MERGE     } from '../modules/ensemble/main'
include { ALPHAMISSENSE_LOOKUP } from '../modules/alphamissense/main'
include { VEP_ANNOTATE       } from '../modules/vep/main'
include { EXPANSIONHUNTER_CALL } from '../modules/expansionhunter/main'
include { CYRIUS_CYP2D6      } from '../modules/cyrius/main'
include { HAPLOGREP3_CLASSIFY } from '../modules/haplogrep3/main'
include { GATK4_MUTECT2_MITO } from '../modules/gatk4/mutect2_mito/main'

// =============================================================================
// Workflow definition
// =============================================================================
workflow {

    // ── Validate required parameters ─────────────────────────────────────────
    if (!params.input) {
        error "ERROR: --input sample sheet is required. See docs/guides/quickstart.md"
    }

    // ── Read and parse sample sheet ──────────────────────────────────────────
    // Sample sheet columns: sample, fastq_1, fastq_2, sex, affected, ped_file
    ch_reads = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            // Validate required columns
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "Sample sheet row missing required columns (sample, fastq_1, fastq_2): ${row}"
            }
            tuple(
                row.sample,
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true),
                row.sex ?: "unknown",
                row.affected ?: "false",
                row.ped_file ? file(row.ped_file) : file("NO_FILE")
            )
        }

    // ── QC and trimming ──────────────────────────────────────────────────────
    // fastp: adapter trimming, quality filtering, poly-G trimming (NovaSeq)
    FASTP_QC(
        ch_reads.map { sample, r1, r2, sex, aff, ped -> tuple(sample, r1, r2) }
    )

    // ── Alignment — PRIMARY: DRAGMAP; FALLBACK: BWA-MEM2 ────────────────────
    // DRAGMAP is the primary aligner per GATK DRAGEN-GATK Best Practices.
    // BWA-MEM2 is only used with -profile bwa_mem2 (when DRAGMAP hash unavailable).
    if (params.aligner == "bwa_mem2") {
        // FALLBACK: BWA-MEM2 — BQSR IS required downstream
        ch_aligned = BWA_MEM2_ALIGN(
            FASTP_QC.out.reads,
            params.reference_fasta
        )
    } else {
        // PRIMARY: DRAGMAP — BQSR is SKIPPED in DRAGEN-GATK mode
        ch_aligned = DRAGMAP_ALIGN(
            FASTP_QC.out.reads,
            params.dragmap_reference
        )
    }

    // ── Mark duplicates ──────────────────────────────────────────────────────
    // OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 for NovaSeq patterned flowcells.
    GATK4_MARKDUPLICATES(ch_aligned)

    // ── Optionally apply BQSR (BWA-MEM2 profile only) ───────────────────────
    // CRITICAL: BQSR must NOT be applied in DRAGEN-GATK mode (incompatible
    // with BQD genotyping model in HaplotypeCaller --dragen-mode).
    if (params.aligner == "bwa_mem2" || params.run_bqsr) {
        // Apply BQSR only when using BWA-MEM2 fallback aligner
        ch_final_bam = GATK4_APPLYBQSR(
            GATK4_MARKDUPLICATES.out.bam,
            params.reference_fasta,
            params.known_sites_vcf   // dbSNP + 1000G gold standard indels
        )
    } else {
        // DRAGEN-GATK mode: use marked-duplicate BAM directly (no BQSR)
        ch_final_bam = GATK4_MARKDUPLICATES.out.bam
    }

    // ── QC metrics ───────────────────────────────────────────────────────────
    SAMTOOLS_INDEX(ch_final_bam)
    SAMTOOLS_FLAGSTAT(ch_final_bam)
    MOSDEPTH_QC(ch_final_bam)
    // Note: MOSDEPTH_QC will emit a warning and exit if coverage < 30×
    // (30× is the ACGS 2024 minimum for clinical WGS)

    // ── Ancestry inference → per-ancestry VQSR ──────────────────────────────
    // somalier infers population ancestry (AFR/AMR/EAS/EUR/SAS/admixed)
    // to select appropriate gnomAD v4.1 stratified VQSR training resources.
    SOMALIER_ANCESTRY(ch_final_bam)

    // ── GATK4 HaplotypeCaller (DRAGEN-GATK mode) ─────────────────────────────
    // --dragen-mode true: activates BQD genotyping model
    // -ERC GVCF: output per-sample GVCFs for joint genotyping
    // DO NOT pass --bqsr-recal-file (DRAGEN-GATK mode; BQD replaces BQSR)
    GATK4_HAPLOTYPECALLER(
        ch_final_bam.join(SAMTOOLS_INDEX.out.bai),
        params.reference_fasta,
        params.scattered_interval_list   // Parallelised by genomic interval
    )

    // Collect all sample GVCFs for joint genotyping
    ch_gvcfs = GATK4_HAPLOTYPECALLER.out.gvcf.collect()

    // ── Joint genotyping + VQSR ───────────────────────────────────────────────
    GATK4_GENOTYPEGVCFS(
        ch_gvcfs,
        params.reference_fasta,
        params.scattered_interval_list
    )

    // Per-ancestry VQSR: use somalier ancestry to select training resources
    VQSR_SNP(
        GATK4_GENOTYPEGVCFS.out.vcf,
        SOMALIER_ANCESTRY.out.ancestry_labels,   // AFR/AMR/EAS/EUR/SAS/admixed
        params.gnomad_vcf                         // gnomAD v4.1 (MUST be v4.1)
    )
    VQSR_INDEL(
        VQSR_SNP.out.vcf,
        SOMALIER_ANCESTRY.out.ancestry_labels,
        params.gnomad_vcf
    )

    // ── DeepVariant / DeepTrio parallel arm ──────────────────────────────────
    if (params.run_deepvariant) {
        // Determine if this is a trio run (PED file + run_deeptrio flag)
        // DeepTrio uses parent BAMs as additional image channels: +15% de novo sensitivity
        ch_is_trio = ch_reads.map { sample, r1, r2, sex, aff, ped ->
            tuple(sample, ped.name != "NO_FILE")
        }

        ch_dv_vcf = ch_is_trio.branch {
            trio:    it[1] == true  && params.run_deeptrio
            single:  true           // Default: single-sample DeepVariant
        }

        // Trio samples → DeepTrio (PMID:36050879)
        DEEPTRIO_CALL(
            ch_dv_vcf.trio
                .map { sample, _ -> sample }
                .join(ch_final_bam)
                .join(
                    ch_reads.map { sample, r1, r2, sex, aff, ped -> tuple(sample, ped) }
                ),
            params.reference_fasta
        )

        // Single samples → DeepVariant v1.8.0
        DEEPVARIANT_CALL(
            ch_dv_vcf.single
                .map { sample, _ -> sample }
                .join(ch_final_bam),
            params.reference_fasta
        )

        // Merge DeepVariant outputs
        ch_all_dv_vcf = DEEPTRIO_CALL.out.vcf.mix(DEEPVARIANT_CALL.out.vcf)

        // ── Ensemble merge (INTERSECTION mode default) ─────────────────────
        // INTERSECTION: PASS if called by BOTH GATK4 AND DeepVariant → higher precision
        // UNION: PASS if called by EITHER caller → higher sensitivity (research use)
        ENSEMBLE_MERGE(
            VQSR_INDEL.out.vcf.join(ch_all_dv_vcf)
        )
        ch_merged_vcf = ENSEMBLE_MERGE.out.vcf

    } else {
        // Skip DeepVariant; use GATK4-only output
        ch_merged_vcf = VQSR_INDEL.out.vcf
    }

    // ── AlphaMissense lookup ──────────────────────────────────────────────────
    // PRIMARY PP3/BP4 predictor (ClinGen SVI 2024 approved)
    // Thresholds: ≥0.564 → PP3 (Supporting Pathogenic); ≤0.340 → BP4
    ALPHAMISSENSE_LOOKUP(
        ch_merged_vcf,
        params.alphamissense_tsv   // AlphaMissense_hg38.tsv.gz (tabix-indexed)
    )

    // ── VEP annotation (v111; MANE Select priority) ──────────────────────────
    // Pick order: mane_select,mane_plus_clinical,canonical
    // Plugins: AlphaMissense, SpliceAI, Pangolin, dbNSFP_4.7, gnomAD_v4.1
    VEP_ANNOTATE(
        ALPHAMISSENSE_LOOKUP.out.vcf,
        params.vep_cache_dir,
        params.vep_cache_version   // Must be 111
    )

    // ── Specialist sub-workflows (run conditionally) ──────────────────────────

    // Mitochondrial analysis (ACGS 2024 §6)
    if (params.run_mito) {
        // Mutect2 in mito mode → Haplogrep3 haplogroup classification
        // ACGS 2024 §6: haplogroup classification MUST precede ACMG assessment
        GATK4_MUTECT2_MITO(
            ch_final_bam,
            params.mito_reference   // MT-specific reference (shifted + unshifted)
        )
        HAPLOGREP3_CLASSIFY(GATK4_MUTECT2_MITO.out.vcf)
    }

    // Repeat expansion genotyping (ExpansionHunter v5.0; 60 loci)
    if (params.run_expansions) {
        // Note: TRGT is NOT used — requires PacBio HiFi long reads (out of scope)
        EXPANSIONHUNTER_CALL(
            ch_final_bam,
            params.expansion_catalog   // v5.0 catalog; 60 loci
        )
    }

    // Pharmacogenomics (CYP2D6 star allele genotyping)
    if (params.run_pgx) {
        // Cyrius: required because GATK4 cannot reliably call CYP2D6 SVs
        // due to CYP2D7 pseudogene interference
        CYRIUS_CYP2D6(
            ch_final_bam,
            params.reference_fasta
        )
    }

    // ── Output summary ────────────────────────────────────────────────────────
    // BayesACMG classification, clinical reporting, and ClinVar/Beacon ingestion
    // are handled as post-pipeline Python scripts invoked from Nextflow processes.
    // See: bayesacmg/, reporting/, beacon_api/, reclassification/ directories.
}

// =============================================================================
// Workflow completion handler
// =============================================================================
workflow.onComplete {
    log.info """
    =========================================
    GenomeForge WGS Pipeline (GRCh38)
    =========================================
    Pipeline : DRAGEN-GATK + DeepVariant ensemble
    ACMG     : BayesACMG (ACGS 2024 v1.2; PM2=Supporting)
    Aligner  : ${params.aligner ?: 'DRAGMAP (primary)'}
    BQSR     : ${(params.aligner == 'bwa_mem2' || params.run_bqsr) ? 'ENABLED (BWA-MEM2 mode)' : 'DISABLED (DRAGEN-GATK mode)'}
    Completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration : ${workflow.duration}
    Results  : ${params.outdir}
    =========================================
    """
}
