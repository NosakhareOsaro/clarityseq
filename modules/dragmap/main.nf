// ============================================================================
// Module: DRAGMAP_ALIGN
// Description: Aligns paired-end WGS reads to GRCh38 using DRAGMAP v1.3.0,
//              the open-source implementation of Illumina DRAGEN's hash-table
//              aligner. This is the PRIMARY aligner in ClaritySeq and replaces
//              BWA-MEM2 as the GATK Best Practice since 2021. Produces a
//              coordinate-sorted BAM ready for MarkDuplicates. BQSR is NOT
//              applied after DRAGMAP when running in DRAGEN-GATK mode.
// Guidelines: GATK Best Practices (DRAGEN-GATK, 2021+);
//             Sci Reports 2022 PMID:36543800;
//             ACGS Best Practice Guidelines v1.2 2024 §3.1
// Inputs:  reads_ch — tuple(sample_id, [read1.fastq.gz, read2.fastq.gz])
//          dragmap_ref — path to pre-built DRAGMAP hash table directory
// Outputs: bam_ch  — tuple(sample_id, sample.sorted.bam)
//          bai_ch  — tuple(sample_id, sample.sorted.bam.bai)
// Container: nfcore/dragmap:1.3.0
// Docs: https://github.com/Illumina/DRAGMAP
//       https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939
// Parameter rationale:
//   --num-threads: use all available CPUs for maximum throughput
//   --RGSM: sample name injected into @RG header for downstream GATK tools
//   --RGPL ILLUMINA: platform tag required by GATK; must match sequencer type
//   --RGLB: library ID differentiates PCR libraries for MarkDuplicates
//   --RGPU: flowcell.lane unique identifier for optical duplicate detection
//   samtools sort -m 4G: per-thread memory cap to avoid OOM on constrained nodes
// Version note: v1.3.0 is the latest stable DRAGMAP release (Jan 2024).
//   It includes fixes for edge cases in hash-table lookup and improved
//   chimeric read handling vs v1.2.x. Do NOT use 'latest' tag — pin always.
// ============================================================================
//
// WHY DRAGMAP INSTEAD OF BWA-MEM2?
// =================================
// DRAGMAP implements the same seed-chain-extend paradigm as BWA-MEM2 but uses
// a compressed hash-table index (built once with dragen-os --build-hash-table)
// instead of the FM-index used by BWA. This achieves:
//   1. ~2x faster alignment throughput on Illumina short reads (150 bp PE)
//   2. Identical or superior sensitivity/specificity vs BWA-MEM2 for SNVs
//      and indels (Sci Reports 2022 PMID:36543800, Table 2)
//   3. Full compatibility with the DRAGEN-GATK HaplotypeCaller BQD model
//      when --dragen-mode true is set downstream
//
// CRITICAL — BQSR IS NOT RUN AFTER DRAGMAP:
// ==========================================
// In DRAGEN-GATK mode, HaplotypeCaller uses the Base Quality Dropoff (BQD)
// genotyping model, which models systematic sequencing errors internally.
// Running BQSR on DRAGMAP-aligned reads and then passing recalibrated reads
// to HaplotypeCaller --dragen-mode REDUCES accuracy because:
//   - BQD expects raw DRAGMAP base qualities, not BQSR-recalibrated scores
//   - The BQD model was trained on DRAGMAP output, not BQSR-modified output
// See conf/dragen_gatk.config where run_bqsr=false is enforced.
// If you are using BWA-MEM2 (--profile bwa_mem2), BQSR IS required.
//
// BUILDING THE DRAGMAP HASH TABLE (one-time setup):
// ==================================================
// The hash table is built from the GRCh38 FASTA once and reused for all
// samples. It is NOT a FASTA file — it is a directory of binary index files.
//
//   dragen-os --build-hash-table true \
//     --ht-reference /ref/GRCh38.fa \
//     --output-directory /ref/dragmap_hg38/ \
//     --ht-num-threads 32
//
// Pre-built tables are available at:
//   s3://broad-references/hg38/v0/dragmap/
//
// The dragmap_reference parameter in nextflow.config must point to the
// directory (not a file). Pass it via:
//   --dragmap_reference /path/to/dragmap_hg38/
//

nextflow.enable.dsl = 2

process DRAGMAP_ALIGN {

    tag "${sample_id}"

    label 'process_high'

    // Pin to exact version — never use 'latest' in production pipelines.
    // nfcore/dragmap:1.3.0 bundles DRAGMAP v1.3.0 + samtools 1.19.
    container 'nfcore/dragmap:1.3.0'

    publishDir "${params.outdir}/alignments/${sample_id}", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.bam') || filename.endsWith('.bai') ? filename : null }

    input:
    // sample_id: string identifier propagated through all channels
    // reads: list of exactly two paths [R1.fastq.gz, R2.fastq.gz]
    // dragmap_ref: path to the DRAGMAP hash-table directory (see header)
    tuple val(sample_id), path(reads)
    path dragmap_ref

    output:
    // Coordinate-sorted BAM — ready for GATK4_MARKDUPLICATES
    tuple val(sample_id), path("${sample_id}.sorted.bam"),  emit: bam
    // BAI index co-located with BAM — required by all downstream GATK tools
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    // Alignment summary metrics from samtools flagstat (piped inline)
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    // Read group fields injected at alignment time so MarkDuplicates and GATK
    // can distinguish libraries and flowcell lanes correctly.
    def rg_id = "${sample_id}.rg1"          // Read Group ID: unique per lane
    def rg_sm = sample_id                    // Sample name: matches VCF sample column
    def rg_lb = "${sample_id}_lib1"         // Library ID: distinguish PCR libraries
    def rg_pl = "ILLUMINA"                  // Platform: GATK requires exact string
    def rg_pu = "${sample_id}.flowcell.1"   // Platform unit: flowcell.lane
    """
    # ── DRAGMAP alignment ──────────────────────────────────────────────────
    # dragen-os: DRAGMAP command-line binary (identical to Illumina DRAGEN API)
    # --num-threads: saturate all CPUs; DRAGMAP is embarrassingly parallel
    # --RGID/RGSM/RGLB/RGPL/RGPU: Read Group tags required by GATK downstream
    # Output is piped directly to samtools to avoid writing raw SAM to disk
    dragen-os \\
        --num-threads ${task.cpus} \\
        --ref ${dragmap_ref} \\
        -1 ${read1} \\
        -2 ${read2} \\
        --RGID  ${rg_id} \\
        --RGSM  ${rg_sm} \\
        --RGLB  ${rg_lb} \\
        --RGPL  ${rg_pl} \\
        --RGPU  ${rg_pu} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -m 4G \\
        -o ${sample_id}.sorted.bam \\
        -

    # ── Index the sorted BAM ───────────────────────────────────────────────
    # BAI index is required by GATK HaplotypeCaller and all interval-scattered
    # sub-processes that use random-access seek into the BAM.
    samtools index \\
        -@ ${task.cpus} \\
        ${sample_id}.sorted.bam

    # ── Record tool versions for provenance ───────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(dragen-os --version 2>&1 | head -1 | sed 's/dragen-os //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai versions.yml
    """
}
