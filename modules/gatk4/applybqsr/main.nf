// ============================================================================
// Module: GATK4_APPLYBQSR
// Description: Applies pre-computed Base Quality Score Recalibration (BQSR)
//              tables to a BAM file, producing recalibrated base quality scores
//              that correct for systematic sequencing errors specific to the
//              instrument, cycle, and base context. This module is ONLY used
//              in the BWA-MEM2 alignment profile and must NOT be run when
//              using DRAGMAP + DRAGEN-GATK HaplotypeCaller mode.
// Guidelines: GATK Best Practices (BWA-MEM2 arm, pre-DRAGEN era);
//             ACGS Best Practice Guidelines v1.2 2024 §3.3 (fallback note)
// Inputs:  bam_recal_ch — tuple(sample_id, sorted.markdup.bam, recal.table)
//                         from GATK4_BASERECALIBRATOR (bwa_mem2 profile only)
//          reference    — GRCh38 FASTA + .fai + .dict
// Outputs: bqsr_bam_ch — tuple(sample_id, sample.bqsr.bam)
//          bqsr_bai_ch — tuple(sample_id, sample.bqsr.bam.bai)
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312
// Parameter rationale:
//   --bqsr-recal-file: the recalibration table from GATK4_BASERECALIBRATOR
//   --static-quantized-quals: quantise recalibrated scores to 10/20/30/40
//     (matches DRAGEN output range — enables direct comparison between arms)
//   --add-output-sam-program-record: record PG header for provenance tracking
//   --use-original-qualities: preserve OQ tag if present (for diagnostics)
// Version note: GATK 4.6.0.0. Pin to exact version.
// ============================================================================
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ONLY used with -profile bwa_mem2
// Must NOT be run in DRAGEN-GATK mode (incompatible with BQD model)
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// THIS MODULE IS ONLY ACTIVATED WHEN: -profile bwa_mem2
//
// WHY THIS MUST NOT RUN IN DRAGEN-GATK MODE:
// ===========================================
// In DRAGEN-GATK mode (--dragen-mode true in HaplotypeCaller), the BQD
// (Base Quality Dropoff) genotyping model corrects for systematic base errors
// INTERNALLY during variant calling. It expects RAW base quality scores from
// DRAGMAP alignment.
//
// If you run BQSR before DRAGEN-GATK HaplotypeCaller:
//   1. The BQD model receives modified (recalibrated) quality scores
//   2. BQD's internal recalibration is mis-calibrated against already-modified scores
//   3. Variant quality scores are DEGRADED — this reduces F1 on GIAB benchmarks
//
// BQSR is valid and required ONLY when:
//   - Aligner: BWA-MEM2 (-profile bwa_mem2)
//   - HaplotypeCaller: standard mode (NOT --dragen-mode)
//
// The nextflow.config dragen_gatk profile enforces: params.run_bqsr = false
// The bwa_mem2 profile sets: params.run_bqsr = true
// The pipeline conditionally includes this module via:
//   if (params.run_bqsr) { GATK4_APPLYBQSR(...) }
//

nextflow.enable.dsl = 2

process GATK4_APPLYBQSR {

    tag "${sample_id}"

    label 'process_medium'

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ONLY used with -profile bwa_mem2. Must NOT be run in DRAGEN-GATK mode.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/alignments/${sample_id}", mode: 'copy',
        pattern: "*.bqsr.bam*"

    input:
    // BAM from GATK4_MARKDUPLICATES + recalibration table from GATK4_BASERECALIBRATOR
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
    // Reference genome (required for ApplyBQSR interval-based processing)
    path reference
    path reference_fai
    path reference_dict

    output:
    // Recalibrated BAM — ready for GATK4_HAPLOTYPECALLER (standard mode, NOT --dragen-mode)
    tuple val(sample_id), path("${sample_id}.bqsr.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.bqsr.bam.bai"), emit: bai
    path "versions.yml",                                      emit: versions

    when:
    // Guard: this process only runs when run_bqsr=true (bwa_mem2 profile)
    params.run_bqsr == true

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 2)}g" : "14g"
    """
    # ── GATK4 ApplyBQSR ───────────────────────────────────────────────────
    # ONLY for -profile bwa_mem2. Must NOT run in DRAGEN-GATK mode.
    # Applies BQSR recalibration table to produce corrected base quality scores.
    gatk --java-options "-Xmx${avail_mem}" \\
        ApplyBQSR \\
        --input ${bam} \\
        --output ${sample_id}.bqsr.bam \\
        --reference ${reference} \\
        \\
        --bqsr-recal-file ${recal_table} \\
        # Recalibration table from GATK4_BASERECALIBRATOR.
        # Contains per-cycle, per-context, per-RG correction factors.
        \\
        --static-quantized-quals 10 \\
        --static-quantized-quals 20 \\
        --static-quantized-quals 30 \\
        --static-quantized-quals 40 \\
        # Quantise output quality scores to {10,20,30,40} Phred bins.
        # Reduces BAM file size and aligns with DRAGEN output quality levels
        # (enabling direct cross-arm comparison in ensemble caller).
        \\
        --add-output-sam-program-record \\
        # Adds @PG header record documenting BQSR parameters for audit trail
        \\
        --use-original-qualities \\
        # Preserve original (pre-BQSR) quality scores in OQ BAM tag.
        # Enables post-hoc inspection without re-running BQSR.
        \\
        --create-output-bam-index \\
        # Create BAI index in same pass (no separate samtools index call needed)
        \\
        --tmp-dir /tmp

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -o 'GATK v[0-9.]*' | sed 's/GATK v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.bqsr.bam ${sample_id}.bqsr.bam.bai versions.yml
    """
}
