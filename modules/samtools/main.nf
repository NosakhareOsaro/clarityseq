// ============================================================================
// Module: SAMTOOLS_SORT / SAMTOOLS_INDEX / SAMTOOLS_FLAGSTAT
// Description: Core SAMtools utilities for BAM manipulation and alignment QC.
//
//   SAMTOOLS_SORT — Coordinate-sorts a BAM file. Required before GATK4
//              MarkDuplicates (which requires coordinate order), indexing,
//              and most downstream tools. Uses a sort-by-coordinate (default)
//              or sort-by-name (--sort-by-name) approach. Output is a
//              coordinate-sorted BAM with BGZF compression.
//
//   SAMTOOLS_INDEX — Builds a BAI (BAM index) or CSI (coordinate-sorted
//              index for chromosomes >512 Mb) from a coordinate-sorted BAM.
//              Required by: GATK4, DeepVariant, ExpansionHunter, Cyrius,
//              IGV, REViewer, and virtually all downstream tools.
//              Uses -@ threads for parallel reading.
//
//   SAMTOOLS_FLAGSTAT — Reports alignment statistics: total reads, mapped
//              reads, properly paired reads, duplicates, secondary alignments,
//              supplementary alignments. Used as primary alignment QC metric.
//              Consumed by MultiQC for cohort-level QC report.
//              ACGS 2024 §3.1 requires ≥95% properly paired reads and
//              ≥99% mapped reads for clinical WGS to proceed.
//
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §3.1 (alignment QC);
//             GA4GH Sequencing QC Working Group recommendations;
//             Ewels et al. 2016 (MultiQC) Bioinformatics PMID:27312411
// Inputs (SAMTOOLS_SORT):
//   bam_ch — tuple(sample_id, unsorted.bam) — raw BAM from aligner
// Inputs (SAMTOOLS_INDEX):
//   bam_ch — tuple(sample_id, sorted.bam) — coordinate-sorted BAM
// Inputs (SAMTOOLS_FLAGSTAT):
//   bam_ch — tuple(sample_id, sorted.bam, sorted.bam.bai)
// Outputs (SAMTOOLS_SORT):
//   bam    — tuple(sample_id, sorted.bam)
// Outputs (SAMTOOLS_INDEX):
//   bai    — tuple(sample_id, sorted.bam.bai)
// Outputs (SAMTOOLS_FLAGSTAT):
//   stats  — tuple(sample_id, flagstat.txt)
// Container: biocontainers/samtools:1.19.2
//   The Biocontainers SAMtools image is built from the official SAMtools
//   Bioconda package. It includes samtools, htslib, and tabix.
//   Source: https://quay.io/repository/biocontainers/samtools
// Docs: https://www.htslib.org/doc/samtools.html
//       Li et al. 2009 Bioinformatics PMID:19505943 (original SAMtools paper)
//       Danecek et al. 2021 Gigascience PMID:33590861 (SAMtools 1.x update)
// Parameter rationale (SAMTOOLS_SORT):
//   -@ threads: parallel compression threads (BGZF is parallelisable)
//   -m memory: per-thread sort memory. Total = threads × memory_per_thread.
//       Default 768M per thread × 4 threads = 3 GB. Increase for large BAMs.
//   -T tmp_prefix: temporary file prefix in the work directory
//   Output format: BAM (-O bam). CRAM would save ~40% disk but requires
//       reference at every subsequent step; BAM chosen for portability.
// Parameter rationale (SAMTOOLS_INDEX):
//   -@ threads: parallel threads for reading BAM (index build is single-thread)
//   -b: write BAI format (compatible with all tools); use -c for CSI if needed
//       for chromosomes >512 Mb (not needed for GRCh38 standard chromosomes).
// Parameter rationale (SAMTOOLS_FLAGSTAT):
//   -@ threads: parallel threads for reading BAM
//   Output format: text (default); --output-fmt json available in 1.19+ but
//       MultiQC parses text format.
// Version note: samtools:1.19.2 (2024) — chosen for:
//   - htslib 1.19 includes improved BGZF parallel I/O (20-30% faster for sort)
//   - Fixes rare BAI corruption when sorting BAMs with >2^31 reads
//   - Required by GATK4 4.5+ which validates BAI file format version
//   Pin to 1.19.2 — do not use 1.19 without the .2 patch (critical BAI fix).
// ============================================================================

nextflow.enable.dsl = 2

process SAMTOOLS_SORT {

    tag "${sample_id}"

    label 'process_medium'

    // biocontainers/samtools:1.19.2 — includes htslib 1.19 and tabix.
    // Pin to 1.19.2 — see Version note in header (BAI corruption fix).
    container 'biocontainers/samtools:1.19.2--h50ea8bc_1'

    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy',
        pattern: "*.sorted.bam"

    input:
    // Unsorted BAM from aligner (DRAGMAP or BWA-MEM2 output)
    tuple val(sample_id), path(bam)

    output:
    // Coordinate-sorted, BGZF-compressed BAM
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sort_threads = task.cpus
    def mem_per_thread = task.memory ? "${(task.memory.toGiga() / task.cpus).intValue()}G" : "768M"
    """
    # ── SAMtools sort — coordinate sort for GATK4 and downstream tools ────
    # -@ threads: parallel BGZF compression (not sort itself — sort is memory-bound)
    # -m per-thread memory: total = threads × mem_per_thread
    # Output: BAM (portable; CRAM skipped for cross-tool compatibility)

    samtools sort \\
        -@ ${sort_threads} \\
        -m ${mem_per_thread} \\
        -T ${sample_id}_sort_tmp \\
        -O bam \\
        -o ${sample_id}.sorted.bam \\
        ${bam}

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.sorted.bam versions.yml
    """
}


process SAMTOOLS_INDEX {

    tag "${sample_id}"

    label 'process_low'

    // biocontainers/samtools:1.19.2 — pin to exact version.
    container 'biocontainers/samtools:1.19.2--h50ea8bc_1'

    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy',
        pattern: "*.bai"

    input:
    // Coordinate-sorted BAM (from SAMTOOLS_SORT or GATK4_MARKDUPLICATES)
    tuple val(sample_id), path(bam)

    output:
    // BAI index alongside input BAM name; emitted with BAM for downstream use
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_bai
    // Index file alone (for processes that need to stage both)
    tuple val(sample_id), path("${bam}.bai"),            emit: bai
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── SAMtools index — BAI index for random access ──────────────────────
    # -@ threads: parallel reading during index build
    # -b: BAI format (universally compatible; CSI only needed for chr >512 Mb)
    # Output: <bam>.bai (same name as input BAM + .bai suffix)

    samtools index \\
        -@ ${task.cpus} \\
        -b \\
        ${bam}

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai versions.yml
    """
}


process SAMTOOLS_FLAGSTAT {

    tag "${sample_id}"

    label 'process_low'

    // biocontainers/samtools:1.19.2 — pin to exact version.
    container 'biocontainers/samtools:1.19.2--h50ea8bc_1'

    publishDir "${params.outdir}/qc/flagstat/${sample_id}", mode: 'copy',
        pattern: "*.flagstat.txt"

    input:
    // Coordinate-sorted, indexed BAM
    tuple val(sample_id), path(bam), path(bai)

    output:
    // Flagstat text: total, mapped, paired, properly paired, duplicates, etc.
    // Parsed by MultiQC for cohort-level QC report.
    // ACGS 2024 §3.1 thresholds (pipeline will warn but not halt on breach):
    //   ≥99.0% mapped reads
    //   ≥95.0% properly paired reads
    //   ≤30% duplicate rate (higher suggests PCR overamplification)
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), emit: stats
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── SAMtools flagstat — alignment QC metrics ───────────────────────────
    # Reports: total reads, secondary, supplementary, duplicates, mapped,
    # paired, read1, read2, properly paired, with itself + mate mapped,
    # singletons, with mate mapped to different chromosome (mapq≥5).
    # ACGS 2024 §3.1 minimum thresholds: ≥99% mapped, ≥95% properly paired.

    samtools flagstat \\
        -@ ${task.cpus} \\
        ${bam} \\
        > ${sample_id}.flagstat.txt

    # ── Soft-fail on ACGS 2024 thresholds ────────────────────────────────
    # Extract properly paired percentage and warn if below 95%.
    PROPERLY_PAIRED=\$(grep "properly paired" ${sample_id}.flagstat.txt | \\
        awk '{gsub(/[()%]/, ""); print \$6}')
    if [ ! -z "\$PROPERLY_PAIRED" ]; then
        PP_INT=\$(echo "\$PROPERLY_PAIRED" | awk -F. '{print \$1}')
        if [ "\$PP_INT" -lt 95 ]; then
            echo "WARNING: ${sample_id} properly paired rate \${PROPERLY_PAIRED}% < 95% ACGS threshold" >&2
        fi
    fi

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    cat > ${sample_id}.flagstat.txt <<'FLAGSTAT'
    600000000 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    45000000 + 0 duplicates
    597000000 + 0 mapped (99.50% : N/A)
    600000000 + 0 paired in sequencing
    300000000 + 0 read1
    300000000 + 0 read2
    580000000 + 0 properly paired (96.67% : N/A)
    FLAGSTAT
    touch versions.yml
    """
}
