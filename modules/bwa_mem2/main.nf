// ============================================================================
// Module: BWA_MEM2_ALIGN
// Description: FALLBACK ALIGNER — aligns paired-end WGS reads using BWA-MEM2
//              v2.2.1 against a GRCh38 FASTA index. This module is ONLY
//              activated when running with -profile bwa_mem2 (e.g., when a
//              DRAGMAP hash-table reference is unavailable). The default
//              ClaritySeq pipeline uses DRAGMAP. Unlike DRAGMAP, BWA-MEM2
//              alignment REQUIRES BQSR downstream — see GATK4_APPLYBQSR module.
// Guidelines: GATK Best Practices (pre-DRAGEN era, still valid for BWA);
//             Li & Durbin 2009 (BWA-MEM) extended to AVX512 in BWA-MEM2;
//             ACGS Best Practice Guidelines v1.2 2024 §3.1 (fallback note)
// Inputs:  reads_ch    — tuple(sample_id, [read1.fastq.gz, read2.fastq.gz])
//          reference   — path to GRCh38 FASTA (with BWA-MEM2 index files)
// Outputs: bam_ch      — tuple(sample_id, sample.sorted.bam)
//          bai_ch      — tuple(sample_id, sample.sorted.bam.bai)
// Container: nfcore/bwa-mem2:2.2.1
// Docs: https://github.com/bwa-mem2/bwa-mem2
//       https://gatk.broadinstitute.org/hc/en-us/articles/360035535912
// Parameter rationale:
//   -t: threads — use all available CPUs for maximum alignment throughput
//   -R: read group string — GATK requires @RG with SM, LB, PL, PU tags
//   -K 100000000: process 100 M bases per batch for deterministic output
//     across different thread counts (important for reproducible CI)
//   -Y: soft-clip supplementary alignments (vs hard-clip) — preserves
//     base quality scores in supplementary records for BQSR calculation
// Version note: v2.2.1 (released 2021) is the latest stable BWA-MEM2 release.
//   It uses AVX-512 SIMD for Smith-Waterman and is ~2x faster than BWA-MEM.
//   Pin to 2.2.1 explicitly — do NOT use 'latest' tag.
// ============================================================================
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// FALLBACK ALIGNER — only activated with -profile bwa_mem2
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// The DEFAULT ClaritySeq pipeline uses DRAGMAP (modules/dragmap/main.nf).
// Use this module ONLY when:
//   a) A DRAGMAP hash table has not been built for your reference, OR
//   b) You are running on infrastructure that cannot build the hash table
//
// HOW TO ACTIVATE:
//   nextflow run pipelines/wgs_grch38.nf -profile bwa_mem2 [other params]
//
// IMPORTANT — BQSR IS REQUIRED WHEN USING BWA-MEM2:
// ==================================================
// Unlike DRAGMAP+DRAGEN-GATK mode (which uses HaplotypeCaller's BQD model
// to correct systematic errors), BWA-MEM2-aligned reads require classical
// GATK Base Quality Score Recalibration (BQSR) before variant calling.
// When -profile bwa_mem2 is active, the pipeline automatically runs:
//   GATK4_BASERECALIBRATOR → GATK4_APPLYBQSR
// before GATK4_HAPLOTYPECALLER. The BQD model flag (--dragen-mode true)
// is NOT passed to HaplotypeCaller in this profile.
//
// This was the GATK Best Practice BEFORE DRAGEN-GATK was released in 2021.
// It remains valid and produces high-quality results, but DRAGMAP is preferred.
//

nextflow.enable.dsl = 2

process BWA_MEM2_ALIGN {

    tag "${sample_id}"

    label 'process_high'

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // FALLBACK ALIGNER — only activated with -profile bwa_mem2
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Pin to 2.2.1 — never use 'latest'. Image bundles bwa-mem2 2.2.1 +
    // samtools 1.19. Source: nfcore/bwa-mem2 Docker Hub.
    container 'nfcore/bwa-mem2:2.2.1'

    publishDir "${params.outdir}/alignments/${sample_id}", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? filename : null }

    input:
    // sample_id: string identifier propagated through all channels
    // reads: list [R1.fastq.gz, R2.fastq.gz] — paired-end required
    tuple val(sample_id), path(reads)
    // reference: GRCh38 FASTA — must have BWA-MEM2 index (.0123, .amb, .ann,
    //   .bwt.2bit.64, .pac) in same directory. Build once with:
    //   bwa-mem2 index GRCh38.fa
    path reference

    output:
    // Coordinate-sorted BAM — must go through BQSR before HaplotypeCaller
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: bam
    // BAI index co-located with BAM
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    // Construct @RG header string for GATK compatibility
    // SM: sample name — appears in VCF column header
    // LB: library — used by MarkDuplicates to avoid cross-library duplicate marking
    // PL: platform — GATK requires "ILLUMINA" (exact capitalisation)
    // PU: platform unit — flowcell.lane for optical duplicate distance calculation
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}_lib1\\tPL:ILLUMINA\\tPU:${sample_id}.flowcell.1"
    """
    # ── BWA-MEM2 alignment ────────────────────────────────────────────────
    # FALLBACK ALIGNER: BQSR MUST be run on the output BAM before HC.
    # bwa-mem2 mem: short-read aligner using FM-index (vs DRAGMAP hash-table)
    # -t: use all available CPUs
    # -R: embed read group header — required by all GATK tools downstream
    # -K 100000000: chunk size for deterministic multi-threaded output
    # -Y: soft-clip supplementary alignments to preserve BQ scores for BQSR
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        -K 100000000 \\
        -Y \\
        ${reference} \\
        ${read1} \\
        ${read2} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -m 4G \\
        -o ${sample_id}.sorted.bam \\
        -

    # ── Index ─────────────────────────────────────────────────────────────
    samtools index \\
        -@ ${task.cpus} \\
        ${sample_id}.sorted.bam

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai versions.yml
    """
}
