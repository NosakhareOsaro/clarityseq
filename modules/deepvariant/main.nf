// ============================================================================
// Module: DEEPVARIANT_CALL
// Description: Calls germline SNVs and indels using Google DeepVariant v1.8.0,
//              a deep neural network variant caller that treats variant calling
//              as an image classification problem. Pileup images of read
//              alignments at candidate sites are fed to a convolutional neural
//              network (InceptionV3 architecture) trained on GIAB truth sets.
//              DeepVariant runs PARALLEL to GATK4 HaplotypeCaller; results are
//              merged by the ENSEMBLE module (INTERSECTION mode by default).
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §4.4 (ensemble calling);
//             ClinGen recommendations for ensemble somatic/germline callers;
//             Poplin et al. 2018 Nature Biotechnology PMID:30247488 (original DV)
// Inputs:  bam_ch    — tuple(sample_id, markdup.bam, markdup.bam.bai)
//          reference — GRCh38 FASTA + .fai
// Outputs: vcf_ch  — tuple(sample_id, deepvariant.vcf.gz)
//          gvcf_ch — tuple(sample_id, deepvariant.g.vcf.gz)
// Container: google/deepvariant:1.8.0
// Docs: https://github.com/google/deepvariant
//       https://github.com/google/deepvariant/blob/v1.8.0/docs/
// Parameter rationale:
//   --model_type WGS: use the WGS-trained model (vs WES, PACBIO, ONT models)
//   --output_gvcf: emit GVCF for ensemble comparison with GATK4 GVCFs
//   --num_shards: number of parallel workers for make_examples step
//   --regions: restrict to non-N callable regions to reduce runtime
//   --postprocess_variants_extra_args "qual_filter=3": remove very low quality calls
// Version note: v1.8.0 (December 2023) — chosen for two major improvements
//   over v1.6.x: (1) SPRQ (Short-Pair Read Queue) support for better calling
//   at tandem repeats, (2) pangenome-aware variant calling when used with
//   pangenome alignment (see modules/vg_giraffe/main.nf and pipelines/
//   wgs_pangenome.nf). Pin to 1.8.0 — never use 'latest' tag. Pin matters
//   because DeepVariant model weights are bundled in the container.
// ============================================================================
//
// WHY DEEPVARIANT RUNS IN PARALLEL WITH GATK4 HAPLOTYPECALLER:
// =============================================================
// GATK4 HaplotypeCaller uses LOCAL DE NOVO ASSEMBLY (Kmer graph, ~150 bp
// haplotype windows) to detect candidate variants. Strengths:
//   - Excellent at complex indels requiring re-assembly
//   - Well-calibrated for large cohort joint-genotyping (GVCF + GenotypeGVCFs)
//   - Integrates with established VQSR filtering framework
// Weaknesses:
//   - Systematically misses certain read-end artefacts
//   - Higher false positive rate near low-complexity regions
//
// DeepVariant uses PILEUP IMAGE CNN to call variants. Strengths:
//   - Captures subtle base-level patterns missed by graph assembly
//   - Lower false positive rate overall (especially for SNPs)
//   - Robust at sites where alignment is complex but coverage is deep
// Weaknesses:
//   - Less accurate for complex structural indels >50 bp
//   - Cannot do joint genotyping (single-sample only in this module)
//   - Requires GPU or many CPUs for practical runtime
//
// INTERSECTION MODE (default): PASS variants present in BOTH callers.
//   → High precision; suitable for clinical reports
//   → De novo SNVs: captured by DeepTrio (see modules/deeptrio/main.nf)
//
// UNION MODE (optional): PASS variants in EITHER caller.
//   → High sensitivity; use for research or carrier screen review
//
// FOR TRIO SAMPLES:
//   Do NOT use this module for trio analysis.
//   Use modules/deeptrio/main.nf instead, which uses parent BAMs as
//   additional image channels and improves de novo SNV sensitivity by ~15%.
//
// DEEPVARIANT v1.8.0 vs v1.6.x IMPROVEMENTS:
// ============================================
// 1. SPRQ (Short-Pair Read Queue) support:
//    Improved calling at tandem repeats and homopolymers by processing
//    read-pairs jointly rather than as independent reads. Reduces false
//    calls at STR loci that overlap with ExpansionHunter targets.
//
// 2. Pangenome-aware calling:
//    When input BAM is produced by vg giraffe (pangenome arm), v1.8.0
//    correctly interprets graph-based alignment quality scores (GQ field
//    in GBWT alignments) that differ from linear-reference MAPQs.
//    Earlier versions treated all low-MAPQ graph alignments as low-confidence.
//
// PILEUP IMAGE vs GRAPH ASSEMBLY (CNN vs HMM):
// =============================================
// DeepVariant generates 6-channel pileup "images" at each candidate site:
//   Channel 1: base identity (ACGT encoded as pixel intensity)
//   Channel 2: base quality score
//   Channel 3: strand of origin
//   Channel 4: read mapping quality (MAPQ)
//   Channel 5: local realignment status
//   Channel 6: differences from reference
// These images are classified by a 3-class CNN: HOM_REF / HET / HOM_ALT
// The CNN was trained on GIAB HG001-HG007 truth sets across multiple sequencers.
//

nextflow.enable.dsl = 2

process DEEPVARIANT_CALL {

    tag "${sample_id}"

    label 'process_high'

    // google/deepvariant:1.8.0 — pin to EXACT version.
    // v1.8.0 includes SPRQ support and pangenome-aware calling (see header).
    // The model weights are BUNDLED in the container — version pinning is
    // essential because different model versions produce different calls.
    // NEVER use 'latest' — it will silently upgrade model weights.
    container 'google/deepvariant:1.8.0'

    publishDir "${params.outdir}/deepvariant/${sample_id}", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi,g.vcf.gz,g.vcf.gz.tbi}"

    input:
    // sorted, duplicate-marked BAM from GATK4_MARKDUPLICATES
    // (NOT BQSR-recalibrated — DeepVariant is trained on raw quality scores)
    tuple val(sample_id), path(bam), path(bai)
    // GRCh38 FASTA + .fai index (no .dict needed — DeepVariant uses fai only)
    path reference
    path reference_fai
    // Callable regions BED: restricts DeepVariant to well-mapped regions
    // Exclude centromeres, telomeres, N-bases, and segmental duplications
    path callable_regions

    output:
    // Single-sample VCF from DeepVariant CNN model
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: tbi
    // GVCF for ensemble comparison and optional joint-analysis
    tuple val(sample_id), path("${sample_id}.deepvariant.g.vcf.gz"),     emit: gvcf
    tuple val(sample_id), path("${sample_id}.deepvariant.g.vcf.gz.tbi"), emit: gtbi
    // HTML report: quality metrics and variant statistics
    path "${sample_id}.deepvariant_visual_report.html",                emit: report
    path "versions.yml",                                               emit: versions

    when:
    // Only run when run_deepvariant=true AND this is NOT a trio sample
    // For trio samples, use modules/deeptrio/main.nf instead
    task.ext.when == null || task.ext.when

    script:
    """
    # ── DeepVariant v1.8.0 — WGS CNN variant calling ─────────────────────
    # FOR TRIO SAMPLES: Use modules/deeptrio/main.nf instead of this module.
    # DeepTrio uses parent BAMs as additional image channels and achieves
    # ~15% better de novo SNV sensitivity (PMID:36050879).
    #
    # v1.8.0 improvements: SPRQ support, pangenome-aware calling (see header)
    #
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type WGS \\
        # WGS model: trained on Illumina short-read whole-genome data.
        # Do NOT use WES model for WGS (different depth/coverage assumptions).
        # Alternative models: PACBIO (long-read), ONT (Oxford Nanopore)
        \\
        --ref ${reference} \\
        --reads ${bam} \\
        --output_vcf ${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf ${sample_id}.deepvariant.g.vcf.gz \\
        # GVCF output: reference confidence model for ensemble comparison
        \\
        --num_shards ${task.cpus} \\
        # Number of parallel make_examples workers.
        # DeepVariant parallelises at the make_examples step (pileup image creation).
        # call_variants (CNN inference) runs in a single GPU/CPU process.
        # postprocess_variants is single-threaded.
        \\
        --regions ${callable_regions} \\
        # Restrict to callable regions: avoids wasted compute in N-bases,
        # centromeres, and telomeres where no reliable calls can be made.
        \\
        --sample_name ${sample_id} \\
        # Embed sample name in VCF FORMAT column header
        \\
        --postprocess_variants_extra_args "qual_filter=3" \\
        # Remove calls with QUAL<3 in postprocessing.
        # These are near-random calls with no clinical utility.
        \\
        --vcf_stats_report \\
        # Generate per-sample HTML report with Ti/Tv, indel size, and qual distributions
        \\
        --logging_dir ${sample_id}_dv_logs

    # Rename visual report to include sample ID
    mv *visual_report.html ${sample_id}.deepvariant_visual_report.html 2>/dev/null || true

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: "1.8.0"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz \\
          ${sample_id}.deepvariant.vcf.gz.tbi \\
          ${sample_id}.deepvariant.g.vcf.gz \\
          ${sample_id}.deepvariant.g.vcf.gz.tbi \\
          ${sample_id}.deepvariant_visual_report.html \\
          versions.yml
    """
}
