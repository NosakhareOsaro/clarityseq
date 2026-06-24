// ============================================================================
// Module: GATK4_MARKDUPLICATES
// Description: Marks PCR and optical duplicates in coordinate-sorted BAM files
//              using GATK4 MarkDuplicates (Picard). Duplicate marking is
//              required before HaplotypeCaller to prevent inflated read
//              support for false variant calls. Optical duplicates from
//              patterned flowcells (NovaSeq) are identified using a larger
//              pixel distance than the historical default.
// Guidelines: GATK Best Practices (MarkDuplicates);
//             ACGS Best Practice Guidelines v1.2 2024 §3.2;
//             Picard MarkDuplicates documentation
// Inputs:  bam_ch — tuple(sample_id, sorted.bam) from DRAGMAP_ALIGN
//                   or BWA_MEM2_ALIGN
// Outputs: md_bam_ch  — tuple(sample_id, markdup.bam) — duplicates flagged
//          md_bai_ch  — tuple(sample_id, markdup.bam.bai) — index
//          metrics_ch — per-sample duplication metrics (used by MultiQC)
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812
//       https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
// Parameter rationale:
//   --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500: CRITICAL for NovaSeq patterned
//     flowcells. Explanation below.
//   --CREATE_INDEX true: generate BAI alongside BAM in one pass (efficient)
//   --VALIDATION_STRINGENCY SILENT: suppress verbose SAM validation warnings
//     that do not affect duplicate marking correctness
//   --ASSUME_SORT_ORDER coordinate: explicitly assert sort order to skip
//     sort validation (BAM is guaranteed sorted by upstream DRAGMAP step)
//   --MAX_RECORDS_IN_RAM 4000000: hold 4M records in memory before spilling
//     to disk; tune based on available RAM vs I/O speed tradeoff
// Version note: GATK 4.6.0.0 (April 2024) is the latest stable release.
//   It includes updated Picard MarkDuplicates with improved flow-cell
//   duplicate detection for AVITI and Element instruments. Pin to 4.6.0.0.
// ============================================================================
//
// OPTICAL_DUPLICATE_PIXEL_DISTANCE — NovaSeq vs HiSeq:
// ======================================================
// --OPTICAL_DUPLICATE_PIXEL_DISTANCE controls how close two clusters must be
// on the flowcell surface to be considered optical duplicates (artefacts from
// the imaging system, not true PCR duplicates).
//
// PATTERNED FLOWCELLS (NovaSeq 6000, NovaSeq X, NextSeq 1000/2000):
//   Distance = 2500 pixels
//   Why: Patterned nanowell arrays have a fixed pitch of ~1000-2500 nm,
//   and optical duplicates from exclusion amplification can appear up to
//   2500 pixels apart. Using the default (100) would MISS the majority
//   of optical duplicates on NovaSeq, leading to over-counted read depth
//   and inflated false positive variant calls.
//   Reference: Illumina Technical Note "Optical Duplicate Rates" (2018)
//
// UNPATTERNED FLOWCELLS (HiSeq 2500, HiSeq 4000):
//   Distance = 100 pixels (Picard default)
//   Why: Random cluster distribution means optical duplicates are rare
//   beyond 100 pixels. Using 2500 on HiSeq would falsely flag legitimate
//   independent clusters as optical duplicates (under-calling true variants).
//
// GenomeForge defaults to 2500 (NovaSeq assumed). Change via:
//   params.optical_duplicate_pixel_distance = 100  # for HiSeq
//

nextflow.enable.dsl = 2

process GATK4_MARKDUPLICATES {

    tag "${sample_id}"

    label 'process_high_memory'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir [
        [
            path: "${params.outdir}/alignments/${sample_id}",
            mode: 'copy',
            pattern: "*.markdup.bam*"
        ],
        [
            path: "${params.outdir}/qc/markduplicates/${sample_id}",
            mode: 'copy',
            pattern: "*.metrics.txt"
        ]
    ]

    input:
    // sample_id: string identifier
    // bam: coordinate-sorted BAM from DRAGMAP_ALIGN or BWA_MEM2_ALIGN
    tuple val(sample_id), path(bam)

    output:
    // BAM with duplicate reads flagged (FLAG 0x400 set); not physically removed
    tuple val(sample_id), path("${sample_id}.markdup.bam"),     emit: bam
    // BAI index for random access by downstream GATK interval-scattered jobs
    tuple val(sample_id), path("${sample_id}.markdup.bam.bai"), emit: bai
    // Duplication metrics: % duplicates, estimated library size — for MultiQC
    tuple val(sample_id), path("${sample_id}.markdup.metrics.txt"), emit: metrics
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // pixel distance: 2500 for patterned flowcells (NovaSeq); 100 for HiSeq
    def pixel_distance = params.optical_duplicate_pixel_distance ?: 2500
    // JVM heap: reserve 2 GB for OS, give the rest to Picard
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 2)}g" : "14g"
    """
    # ── GATK4 MarkDuplicates ──────────────────────────────────────────────
    # Marks (does NOT remove) duplicate reads so downstream GATK tools
    # can ignore them without physically deleting data.
    #
    # OPTICAL_DUPLICATE_PIXEL_DISTANCE=${pixel_distance}:
    #   Set to 2500 for NovaSeq patterned flowcells (default 100 for HiSeq).
    #   See module header for full explanation.
    gatk --java-options "-Xmx${avail_mem} -XX:ParallelGCThreads=4" \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${sample_id}.markdup.bam \\
        \\
        --METRICS_FILE ${sample_id}.markdup.metrics.txt \\
        # Duplication rate metrics: consumed by MultiQC for QC reporting
        \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${pixel_distance} \\
        # CRITICAL: 2500 for NovaSeq patterned flowcells; 100 for HiSeq only.
        # Using 100 on NovaSeq MISSES most optical duplicates.
        \\
        --CREATE_INDEX true \\
        # Build BAI index in the same MarkDuplicates pass (no extra samtools step)
        \\
        --VALIDATION_STRINGENCY SILENT \\
        # Suppress non-critical SAM format warnings (does not affect correctness)
        \\
        --ASSUME_SORT_ORDER coordinate \\
        # Skip sort-order validation; BAM is guaranteed sorted by aligner step
        \\
        --MAX_RECORDS_IN_RAM 4000000 \\
        # Hold 4M records in RAM before spilling to temp files.
        # Increase on high-memory nodes to reduce I/O overhead.
        \\
        --TMP_DIR /tmp

    # ── Rename index to standard .bam.bai convention ──────────────────────
    # GATK writes <name>.bai but downstream tools expect <name>.bam.bai
    mv ${sample_id}.markdup.bai ${sample_id}.markdup.bam.bai 2>/dev/null || true

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -o 'GATK v[0-9.]*' | sed 's/GATK v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.markdup.bam \\
          ${sample_id}.markdup.bam.bai \\
          ${sample_id}.markdup.metrics.txt \\
          versions.yml
    """
}
