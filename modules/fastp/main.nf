// ============================================================================
// Module: FASTP_QC
// Description: Performs adapter trimming, quality filtering, and poly-G tail
//              trimming on paired-end Illumina WGS reads using fastp v0.23.4.
//              Poly-G trimming is critical for NovaSeq data because the
//              two-channel chemistry encodes dark (no signal) cycles as G,
//              producing artefactual poly-G tails that impair alignment and
//              inflate false variant calls near read ends. Generates JSON and
//              HTML QC reports consumed by MultiQC.
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §2 (pre-alignment QC);
//             Illumina NovaSeq poly-G technical note (2019);
//             Chen et al. 2018 iScience PMID:30423086 (fastp publication)
// Inputs:  reads_ch — tuple(sample_id, [read1.fastq.gz, read2.fastq.gz])
// Outputs: trimmed_reads — tuple(sample_id, [trimmed_R1.fq.gz, trimmed_R2.fq.gz])
//          json_report   — per-sample fastp JSON for MultiQC aggregation
//          html_report   — per-sample fastp HTML report
// Container: biocontainers/fastp:0.23.4
// Docs: https://github.com/OpenGene/fastp
//       https://doi.org/10.1016/j.isci.2018.11.021
// Parameter rationale:
//   --detect_adapter_for_pe: auto-detect adapters from PE overlap (no adapter
//     sequence needed); robust to mixed library preparations
//   --trim_poly_g: enable poly-G tail trimming (critical for NovaSeq);
//     default minimum length 10 nt for poly-G detection
//   --poly_g_min_len 10: require ≥10 consecutive G bases to trigger trimming
//     (avoids over-trimming legitimate G-rich regions)
//   --trim_poly_x: also trim poly-X (poly-A, poly-C, poly-T) as a precaution
//   --qualified_quality_phred 20: Q<20 bases fail quality filter (1% error)
//   --unqualified_percent_limit 40: discard reads with >40% unqualified bases
//   --length_required 36: discard reads shorter than 36 bp after trimming
//     (36 bp is the practical minimum for unique alignment to GRCh38)
//   --thread: use available CPUs; fastp scales well to 16 threads
//   --json / --html: structured output for MultiQC consumption
// Version note: fastp 0.23.4 (Nov 2023) is the latest stable release.
//   v0.23.x improves adapter overlap detection vs v0.22.x and adds
//   support for newer Illumina adapter sequences. Pin to 0.23.4.
// ============================================================================

nextflow.enable.dsl = 2

process FASTP_QC {

    tag "${sample_id}"

    label 'process_medium'

    // Pin to exact version. biocontainers/fastp:0.23.4 is built from the
    // official BioContainers registry and contains only fastp + dependencies.
    container 'biocontainers/fastp:0.23.4--h5f740d0_3'

    publishDir [
        [
            path: "${params.outdir}/qc/fastp/${sample_id}",
            mode: 'copy',
            pattern: "*.{json,html}"
        ],
        [
            path: "${params.outdir}/reads/trimmed/${sample_id}",
            mode: 'copy',
            pattern: "*.trimmed.fastq.gz",
            // Only save trimmed reads if params.save_trimmed is true
            // (default: false — saves disk; trimmed reads consumed immediately)
            enabled: params.save_trimmed ?: false
        ]
    ]

    input:
    // sample_id: string identifier for naming output files and log messages
    // reads: list [R1.fastq.gz, R2.fastq.gz] — must be exactly two files
    tuple val(sample_id), path(reads)

    output:
    // Adapter-trimmed, quality-filtered paired reads for downstream alignment
    tuple val(sample_id), path("${sample_id}.trimmed_{R1,R2}.fastq.gz"), emit: reads
    // JSON report: machine-readable stats consumed by MultiQC
    path "${sample_id}.fastp.json",                                       emit: json
    // HTML report: human-readable QC summary per sample
    path "${sample_id}.fastp.html",                                       emit: html
    path "versions.yml",                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    # ── fastp: adapter trimming + QC filtering ────────────────────────────
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${sample_id}.trimmed_R1.fastq.gz \\
        --out2 ${sample_id}.trimmed_R2.fastq.gz \\
        \\
        --detect_adapter_for_pe \\
        # Auto-detect adapters by PE overlap — works for TruSeq, Nextera, etc.
        # No need to specify adapter sequences manually.
        \\
        --trim_poly_g \\
        # CRITICAL for NovaSeq: two-channel chemistry encodes dark cycles as G,
        # causing artefactual poly-G tails that mis-map near contig boundaries
        # and inflate false positive variant calls at read ends.
        \\
        --poly_g_min_len 10 \\
        # Require ≥10 consecutive Gs to trigger poly-G trimming.
        # This avoids over-trimming genuinely G-rich exonic regions.
        \\
        --trim_poly_x \\
        # Also trim poly-A/C/T tails (less common but seen in some libraries)
        \\
        --qualified_quality_phred 20 \\
        # Base quality threshold: Q<20 (>1% error probability) flagged as low quality
        \\
        --unqualified_percent_limit 40 \\
        # Discard read if >40% of bases are below Q20 threshold
        \\
        --length_required 36 \\
        # Minimum read length after trimming: 36 bp is the floor for
        # reliable unique alignment to GRCh38 (mappability track based
        # on 36-mers covers ~95% of the genome)
        \\
        --thread ${task.cpus} \\
        # Use all available CPUs; fastp is multi-threaded for I/O and processing
        \\
        --json ${sample_id}.fastp.json \\
        # Structured JSON for MultiQC aggregation across samples
        \\
        --html ${sample_id}.fastp.html
        # Human-readable HTML report with duplication and GC bias plots

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.trimmed_R1.fastq.gz \\
          ${sample_id}.trimmed_R2.fastq.gz \\
          ${sample_id}.fastp.json \\
          ${sample_id}.fastp.html \\
          versions.yml
    """
}
