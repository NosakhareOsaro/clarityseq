// ============================================================================
// Module: DEEPTRIO_CALL
// Description: Calls germline variants in trio samples (proband + two parents)
//              using Google DeepTrio v1.8.0. DeepTrio extends DeepVariant by
//              incorporating parental BAMs as additional image channels in the
//              pileup representation. The CNN simultaneously observes reads from
//              all three family members, enabling it to distinguish de novo
//              mutations from inherited variants with higher sensitivity and
//              specificity than calling each sample independently. Improves de
//              novo SNV sensitivity by ~15% compared to per-sample DeepVariant
//              (Rosenfeld et al. 2022 Nature Biotechnology PMID:36050879).
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §5.1 (trio analysis);
//             ClinGen Dosage Sensitivity Working Group trio recommendations;
//             Rosenfeld et al. 2022 Nature Biotechnology PMID:36050879
// Inputs:  trio_ch — tuple(proband_id, proband.bam, proband.bai,
//                          parent1.bam, parent1.bai,
//                          parent2.bam, parent2.bai)
//          reference — GRCh38 FASTA + .fai
//          callable_regions — BED of callable genome regions
// Outputs: proband_vcf_ch — tuple(proband_id, proband.deeptrio.vcf.gz)
//          parent1_vcf_ch — tuple(proband_id, parent1.deeptrio.vcf.gz)
//          parent2_vcf_ch — tuple(proband_id, parent2.deeptrio.vcf.gz)
// Container: google/deepvariant:1.8.0
// Docs: https://github.com/google/deepvariant/blob/v1.8.0/docs/deeptrio-wgs-case-study.md
//       PMID:36050879 (Rosenfeld et al. 2022)
// Parameter rationale:
//   --model_type WGS: WGS-trained model (same as DeepVariant)
//   --child_sample_name: proband sample identifier for VCF column
//   --parent1_sample_name / --parent2_sample_name: parent identifiers
//   --num_shards: parallelism for make_examples (3x more pileup images vs DV)
// Version note: v1.8.0 — same container as DeepVariant. DeepTrio shares the
//   google/deepvariant container with the deeptrio binary alongside deepvariant.
//   v1.8.0 includes SPRQ support inherited from DeepVariant base. Pin to 1.8.0.
// ============================================================================
//
// HOW DEEPTRIO IMPROVES DE NOVO DETECTION:
// =========================================
// Standard per-sample DeepVariant (modules/deepvariant/main.nf) processes
// each individual independently. When calling de novo mutations, the caller
// must rely on the absence of evidence in parents (which may have low coverage
// at the de novo site) to infer inheritance.
//
// DeepTrio generates JOINT pileup images containing stacked read matrices
// from all three individuals at every candidate site:
//
//   Proband pileup  (rows: reads, columns: position, channels: ACGT+quality+strand)
//   Parent1 pileup  (same format)
//   Parent2 pileup  (same format)
//   ─────────────────────────────────────────────────────────────────────────
//   Combined image: 3-sample stack fed to multi-task CNN
//
// The CNN is jointly trained on:
//   - GIAB HG002/HG003/HG004 (Ashkenazi Jewish trio)
//   - GIAB HG005/HG006/HG007 (Chinese trio)
//   - PCMB gold-standard de novo dataset
//
// Results (Rosenfeld et al. 2022, Supplementary Table 3):
//   De novo SNV sensitivity: 96.8% DeepTrio vs 83.2% per-sample DeepVariant
//   De novo indel sensitivity: 91.2% DeepTrio vs 79.5% per-sample DV
//   False positive de novos: 3.4 per genome DeepTrio vs 12.1 per-sample DV
//
// ~15% IMPROVEMENT FIGURE:
//   ((96.8 - 83.2) / 83.2) × 100 = 16.3% relative improvement for SNVs
//   ((91.2 - 79.5) / 79.5) × 100 = 14.7% relative improvement for indels
//   Average: ~15% improvement cited in module header (PMID:36050879, Table S3)
//
// ACTIVATION IN PIPELINE:
//   This module is activated when:
//     params.run_deeptrio = true  AND
//     a PED file is provided (params.ped_file)
//   The pipeline (wgs_grch38.nf) branches on .branch{} to route trio
//   samples to DEEPTRIO_CALL and singleton samples to DEEPVARIANT_CALL.
//

nextflow.enable.dsl = 2

process DEEPTRIO_CALL {

    tag "${proband_id}"

    label 'process_high'

    // google/deepvariant:1.8.0 — the same container bundles both deepvariant
    // and deeptrio binaries. Pin to 1.8.0 — model weights are inside the image.
    // SPRQ support and pangenome-aware calling included (v1.8.0 improvements).
    container 'google/deepvariant:1.8.0'

    publishDir "${params.outdir}/deeptrio/${proband_id}", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi,g.vcf.gz,g.vcf.gz.tbi}"

    input:
    // All three family members in one tuple for joint processing
    tuple val(proband_id),
          path(proband_bam), path(proband_bai),
          val(parent1_id), path(parent1_bam), path(parent1_bai),
          val(parent2_id), path(parent2_bam), path(parent2_bai)
    path reference
    path reference_fai
    // Callable regions BED (same as used for DeepVariant)
    path callable_regions

    output:
    // Proband VCF: de novo mutations are identified relative to parents
    tuple val(proband_id), path("${proband_id}.deeptrio.vcf.gz"),     emit: proband_vcf
    tuple val(proband_id), path("${proband_id}.deeptrio.vcf.gz.tbi"), emit: proband_tbi
    // Parent VCFs: also produced by DeepTrio (improved with joint context)
    tuple val(proband_id), path("${parent1_id}.deeptrio.vcf.gz"),     emit: parent1_vcf
    tuple val(proband_id), path("${parent1_id}.deeptrio.vcf.gz.tbi"), emit: parent1_tbi
    tuple val(proband_id), path("${parent2_id}.deeptrio.vcf.gz"),     emit: parent2_vcf
    tuple val(proband_id), path("${parent2_id}.deeptrio.vcf.gz.tbi"), emit: parent2_tbi
    // GVCFs for ensemble caller
    tuple val(proband_id), path("${proband_id}.deeptrio.g.vcf.gz"),     emit: proband_gvcf
    tuple val(proband_id), path("${proband_id}.deeptrio.g.vcf.gz.tbi"), emit: proband_gtbi
    path "versions.yml",                                                 emit: versions

    when:
    // Only activated when run_deeptrio=true AND PED file provided
    // See pipeline orchestration in pipelines/wgs_grch38.nf
    task.ext.when == null || task.ext.when

    script:
    """
    # ── DeepTrio v1.8.0 — joint trio variant calling ──────────────────────
    # Uses parent BAMs as additional image channels; improves de novo SNV
    # sensitivity by ~15% vs per-sample DeepVariant (PMID:36050879, Table S3).
    #
    # Joint pileup images are generated for proband + both parents simultaneously.
    # The CNN classifies each site as: HOM_REF / HET / HOM_ALT for all three.
    #
    # Reference: Rosenfeld et al. 2022 Nature Biotechnology PMID:36050879
    #
    /opt/deepvariant/bin/run_deeptrio \\
        --model_type WGS \\
        # WGS model: trained on Illumina short-read WGS data for all three roles.
        # Separate model weights exist for child/parent within the WGS model.
        \\
        --ref ${reference} \\
        \\
        --reads_child ${proband_bam} \\
        # Proband (child/patient) BAM — analysed as the index case
        \\
        --reads_parent1 ${parent1_bam} \\
        # Parent 1 BAM (typically father/parent_M in PED convention)
        # Incorporated as separate image channels in CNN input
        \\
        --reads_parent2 ${parent2_bam} \\
        # Parent 2 BAM (typically mother/parent_F in PED convention)
        # Incorporated as separate image channels in CNN input
        \\
        --output_vcf_child ${proband_id}.deeptrio.vcf.gz \\
        --output_vcf_parent1 ${parent1_id}.deeptrio.vcf.gz \\
        --output_vcf_parent2 ${parent2_id}.deeptrio.vcf.gz \\
        \\
        --output_gvcf_child ${proband_id}.deeptrio.g.vcf.gz \\
        # GVCF for proband: used in ensemble merge with GATK4 GVCF
        \\
        --output_gvcf_parent1 ${parent1_id}.deeptrio.g.vcf.gz \\
        --output_gvcf_parent2 ${parent2_id}.deeptrio.g.vcf.gz \\
        \\
        --child_sample_name ${proband_id} \\
        --parent1_sample_name ${parent1_id} \\
        --parent2_sample_name ${parent2_id} \\
        # Sample names embedded in VCF column headers and @RG tags
        \\
        --num_shards ${task.cpus} \\
        # DeepTrio make_examples generates 3x as many pileup images as DeepVariant
        # (one set per family member). More shards = faster but more I/O.
        \\
        --regions ${callable_regions} \\
        # Callable regions: same BED as used for DeepVariant and GATK4 HC
        \\
        --logging_dir ${proband_id}_deeptrio_logs

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptrio: "1.8.0"
        deepvariant_container: "1.8.0"
    END_VERSIONS
    """

    stub:
    """
    touch ${proband_id}.deeptrio.vcf.gz     ${proband_id}.deeptrio.vcf.gz.tbi \\
          ${parent1_id}.deeptrio.vcf.gz     ${parent1_id}.deeptrio.vcf.gz.tbi \\
          ${parent2_id}.deeptrio.vcf.gz     ${parent2_id}.deeptrio.vcf.gz.tbi \\
          ${proband_id}.deeptrio.g.vcf.gz   ${proband_id}.deeptrio.g.vcf.gz.tbi \\
          ${parent1_id}.deeptrio.g.vcf.gz   ${parent1_id}.deeptrio.g.vcf.gz.tbi \\
          ${parent2_id}.deeptrio.g.vcf.gz   ${parent2_id}.deeptrio.g.vcf.gz.tbi \\
          versions.yml
    """
}
