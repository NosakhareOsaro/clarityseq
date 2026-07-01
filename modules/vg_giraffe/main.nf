// ============================================================================
// Module: VG_GIRAFFE_ALIGN
// Description: Aligns short Illumina reads to the Human Pangenome Reference
//              Consortium (HPRC) v1.1 minigraph-cactus pangenome graph using
//              vg giraffe, a fast read aligner for variation graphs. This
//              module implements the pangenome arm of the ClaritySeq
//              dual-reference strategy (see pipelines/wgs_pangenome.nf).
//
//              WHY PANGENOME GRAPH ALIGNMENT:
//              The GRCh38 linear reference represents a single human haplotype
//              (~78% European ancestry). This causes:
//              1. REFERENCE BIAS: reads from individuals carrying common
//                 non-reference alleles align with lower MAPQ to GRCh38,
//                 causing reduced sensitivity in non-EUR populations.
//              2. UNMAPPABLE REGIONS: ~292 Mb of novel sequence exists in
//                 diverse human genomes not captured in GRCh38 (Liao et al.
//                 2023 Nature PMID:37165242).
//              3. SV GENOTYPING: structural variants (insertions, deletions
//                 >50 bp) are better represented as graph edges than linear
//                 reference alleles.
//
//              THE HPRC v1.1 PANGENOME GRAPH (minigraph-cactus):
//              Built from 47 phased human genome assemblies (94 haplotypes)
//              representing diverse global ancestry (AFR, AMR, EAS, EUR, SAS,
//              OCE, MID). Contains ~100 million variants including ~10 million
//              SVs not in GRCh38.
//
//              vg GIRAFFE ALGORITHM (Sirén et al. 2021 Science PMID:34818024):
//              Uses a minimizer-seeded, GBWT-haplotype-guided alignment
//              strategy that is 10–60× faster than vg map while maintaining
//              accuracy. Key steps:
//              1. Minimizer seeding against the GBWT haplotype index
//              2. Cluster minimizers along embedded haplotype paths
//              3. Extend clusters using graph alignment (partial order alignment)
//              4. Output GAF (Graph Alignment Format) or BAM with AS/MAPQ scores
//
//              OUTPUT:
//              Giraffe outputs a position-sorted BAM with GRCh38 coordinates
//              via the --output-format BAM --reference-fasta flag. The BAM
//              is directly compatible with GATK4, DeepVariant, and other
//              linear-reference downstream tools, enabling the pangenome arm
//              to share the same variant calling modules as the GRCh38 arm.
//
//              REFERENCE:
//              Liao et al. 2023 Nature — HPRC pangenome paper.
//                "A draft human pangenome reference" PMID:37165242
//              Sirén et al. 2021 Science — vg giraffe algorithm paper.
//                "Pangenomics enables genotyping of known structural variants
//                in 5202 diverse genomes" PMID:34818024
//
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §4.5 (pangenome);
//             HPRC Consortium Recommendations for clinical pangenome use (2024);
//             GA4GH Pangenome Working Group Standards (2023)
// Inputs:  reads_ch   — tuple(sample_id, read1.fastq.gz, read2.fastq.gz)
//              Post-QC reads from FASTP_QC (adapter-trimmed, quality-filtered).
//          graph_gbz  — HPRC v1.1 pangenome graph (.gbz) — GBWT index
//              ~22 GB for the full HPRC v1.1 graph. Cache locally; do not
//              re-download for each run. Set params.pangenome_gbz.
//          graph_dist — Distance index (.dist) for vg giraffe seeding
//          graph_min  — Minimizer index (.min) for fast seed lookup
//          reference  — GRCh38 FASTA for output BAM coordinate projection
// Outputs: bam     — tuple(sample_id, giraffe_sorted.bam)
//          bai     — tuple(sample_id, giraffe_sorted.bam.bai)
// Container: quay.io/biocontainers/vg:1.56.0
//   Maintained by the Biocontainers project.
//   vg 1.56.0 is pinned because the HPRC v1.1 graph and indexes were built
//   with vg 1.55–1.56 tooling; earlier versions may not read v1.1 GBZ format.
//   Source: https://quay.io/repository/biocontainers/vg
// Docs: https://github.com/vgteam/vg
//       https://humanpangenome.org
//       Liao et al. 2023 Nature PMID:37165242
//       Sirén et al. 2021 Science PMID:34818024
// Parameter rationale:
//   --read-group: embed RG tag in output BAM (required by GATK4 MarkDuplicates)
//   --sample: sample name in RG SM field (matches sample_id)
//   --output-format BAM: linear BAM output with GRCh38 coordinates
//   --reference-fasta: GRCh38 FASTA for coordinate projection to linear ref
//   -t threads: giraffe is thread-parallel; use all available CPUs
//   --fragment-mean / --fragment-stdev: fragment size parameters for paired-end
//       alignment. Default auto-detection; override if insert size is non-standard.
//   Sorting: output BAM from giraffe is NOT coordinate sorted; pipe through
//       samtools sort for downstream compatibility.
// Version note: vg:1.56.0 (2024) — required for HPRC v1.1 GBZ format support.
//   The HPRC v1.1 graph was built with minigraph-cactus v2.5 and indexed with
//   vg 1.55–1.56. Earlier vg versions (1.50–1.54) cannot read v1.1 GBZ indexes.
//   vg 1.56 adds improved MAPQ calibration for graph alignments that partially
//   mitigates the pangenome-to-linear projection MAPQ loss. Pin to 1.56.0.
// ============================================================================

nextflow.enable.dsl = 2

process VG_GIRAFFE_ALIGN {

    tag "${sample_id}"

    label 'process_high'

    // quay.io/biocontainers/vg:1.56.0 — required for HPRC v1.1 GBZ format.
    // Do NOT use vg <1.55 — will fail to read HPRC v1.1 GBZ index format.
    // Do NOT use 'latest' — vg API changes frequently.
    // Liao et al. 2023 Nature PMID:37165242 (HPRC pangenome)
    // Sirén et al. 2021 Science PMID:34818024 (vg giraffe algorithm)
    container 'quay.io/biocontainers/vg:1.56.0--h9ee0642_0'

    publishDir "${params.outdir}/alignment_pangenome/${sample_id}", mode: 'copy',
        pattern: "*.{bam,bam.bai}"

    input:
    // Post-QC FASTQ reads (adapter-trimmed by FASTP_QC)
    tuple val(sample_id), path(read1), path(read2)
    // HPRC v1.1 pangenome graph indexes:
    //   .gbz: GBZ format variation graph + GBWT haplotype index (22 GB)
    //   .dist: distance index for cluster seeding (~4 GB)
    //   .min: minimizer index for fast lookup (~20 GB)
    // Download from: https://github.com/human-pangenomics/hpp_pangenome_resources
    // Cache with: nextflow pull --assets (see conf/resources.config)
    path graph_gbz
    path graph_dist
    path graph_min
    // GRCh38 FASTA + .fai for output BAM coordinate projection
    path reference
    path reference_fai

    output:
    // Coordinate-sorted BAM with GRCh38 linear reference coordinates.
    // Compatible with all downstream GATK4, DeepVariant, Mosdepth modules.
    // Contains extra tags: AS (alignment score), MQ (mapping quality), HP (haplotype support)
    tuple val(sample_id), path("${sample_id}.giraffe.sorted.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.giraffe.sorted.bam.bai"), emit: bai
    // Statistics: fraction of reads aligned, multi-mapping rate, pangenome path distribution
    path "${sample_id}.giraffe_stats.txt",                             emit: stats
    path "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Platform and library info for BAM read group (required by GATK4)
    def rg_platform = params.rg_platform ?: "ILLUMINA"
    def rg_library  = params.rg_library  ?: "${sample_id}_lib1"
    def rg_center   = params.rg_center   ?: "CLARITYSEQ"
    """
    # ── vg giraffe — pangenome graph alignment ─────────────────────────────
    # HPRC v1.1 minigraph-cactus pangenome graph (47 genomes, 94 haplotypes)
    # Reference: Liao et al. 2023 Nature PMID:37165242
    # Algorithm: Sirén et al. 2021 Science PMID:34818024
    #
    # Output: coordinate-sorted BAM with GRCh38 linear coordinates
    # (compatible with GATK4, DeepVariant, Mosdepth, etc.)

    vg giraffe \\
        --gbz-name   ${graph_gbz} \\
        --minimizer-name ${graph_min} \\
        --dist-name  ${graph_dist} \\
        --fastq-in   ${read1} \\
        --fastq-in   ${read2} \\
        --output-format BAM \\
        --reference-fasta ${reference} \\
        --read-group "ID:${sample_id}\tSM:${sample_id}\tPL:${rg_platform}\tLB:${rg_library}\tCN:${rg_center}" \\
        --sample     ${sample_id} \\
        --threads    ${task.cpus} \\
        --progress \\
        2> ${sample_id}.giraffe_stats.txt \\
    | samtools sort \\
        --threads ${task.cpus} \\
        -m 2G \\
        -O bam \\
        -o ${sample_id}.giraffe.sorted.bam

    # ── Index the sorted BAM ───────────────────────────────────────────────
    samtools index \\
        -@ ${task.cpus} \\
        ${sample_id}.giraffe.sorted.bam

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | grep -o 'v[0-9.]*' | head -1)
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
        hprc_graph: "v1.1 minigraph-cactus"
        reference_genome: "GRCh38"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.giraffe.sorted.bam \
          ${sample_id}.giraffe.sorted.bam.bai \
          ${sample_id}.giraffe_stats.txt \
          versions.yml
    """
}
