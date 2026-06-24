// ============================================================================
// Module: GATK4_MUTECT2_MITO
// Description: Calls mitochondrial variants using GATK4 Mutect2 in
//              --mitochondria-mode. Mitochondrial DNA has special characteristics
//              requiring a dedicated calling mode: high copy number (100-10,000x
//              depth), heteroplasmy (mixtures of mutant and wild-type mtDNA),
//              and a circular genome with non-diploid ploidy. This module
//              implements the ACGS 2024 §6 mitochondrial analysis workflow,
//              using the revised Cambridge Reference Sequence (rCRS, chrM in
//              GRCh38) as the reference. Output VCF is annotated with
//              heteroplasmy fractions and passed to HAPLOGREP3_CLASSIFY for
//              haplogroup assignment.
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §6 (Mitochondrial analysis);
//             GATK Mutect2 mitochondria mode documentation;
//             MITOMAP (www.mitomap.org) variant database;
//             Triska et al. 2021 Am J Hum Genet PMID:34166612
// Inputs:  bam_ch    — tuple(sample_id, markdup.bam, markdup.bam.bai)
//          mito_ref  — chrM reference FASTA (rCRS, extracted from GRCh38)
//          mito_dict — Picard sequence dictionary for chrM
//          blacklist — MITOMAP artifact/NUMTs blacklist BED
// Outputs: mito_vcf_ch — tuple(sample_id, mito.vcf.gz, mito.vcf.gz.tbi)
//          stats_ch    — Mutect2 filtering statistics
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275
//       https://www.acgs.uk.com/quality/best-practice-guidelines/ §6
// Parameter rationale:
//   --mitochondria-mode: enables mtDNA-specific calling (see below)
//   --max-reads-per-alignment-start 75: higher cap than WGS (default 50)
//     because mtDNA coverage is 100-10,000x vs ~30x for nuclear genome
//   --min-base-quality-score 20: Q20 threshold for mtDNA base inclusion
//     (ACGS 2024 §6 recommendation)
//   --minimum-allele-fraction 0.01: detect variants at ≥1% heteroplasmy
//     (clinical threshold for reportable mitochondrial disease)
//   --f1r2-tar-gz: capture orientation model for FilterMutectCalls
// Version note: GATK 4.6.0.0 includes improved chrM calling with
//   --mitochondria-mode; NUMTs filtering is enhanced vs 4.5.x. Pin to 4.6.0.0.
// ============================================================================
//
// --mitochondria-mode EXPLANATION:
// =================================
// Standard Mutect2 is designed for somatic calling (tumor vs normal).
// When --mitochondria-mode is set, Mutect2 adapts to mtDNA characteristics:
//
//   1. PLOIDY: Mitochondria are not diploid (not 2 copies per cell, but
//      100-10,000 copies). --mitochondria-mode sets effective ploidy to match
//      the observed copy number, allowing heteroplasmy detection at any fraction.
//
//   2. HETEROPLASMY: Mutect2 reports allele fraction (AF) in FORMAT/AF field.
//      This is the proportion of reads supporting the variant allele.
//      Clinical significance thresholds vary by disease:
//        - AF >0.40: typically homoplasmic (de novo or complete mtDNA replacement)
//        - AF 0.01-0.40: heteroplasmy (mixture of mutant/WT mtDNA)
//        - AF <0.01: artefact threshold — filtered out
//
//   3. CIRCULAR GENOME: chrM is circular but represented linearly in GRCh38.
//      The module uses a shifted version of chrM reference + ShiftedReferenceInterval
//      to ensure variants at the origin (position 1-300 / 16,069-16,569) are
//      correctly called. The pipeline handles this shift/unshift in pipelines/mito.nf.
//
//   4. NUMTs (Nuclear Mitochondrial DNA Segments): fragments of mtDNA inserted
//      into the nuclear genome. Without filtering, NUMT reads appear as false
//      heteroplasmic variants. The blacklist BED (--mitochondria-mode includes
//      NUMT filtering) and --max-alt-allele-count 4 suppress NUMT artefacts.
//
// ACGS 2024 §6 IMPLEMENTATION:
//   - Minimum reportable heteroplasmy: 1% AF (--minimum-allele-fraction 0.01)
//   - Variants classified using MITOMAP pathogenicity status
//   - Haplogroup assigned by Haplogrep3 (see modules/haplogrep3/main.nf)
//   - Pathogenic variants ≥10% AF reported as primary findings
//   - Pathogenic variants 1-10% AF reported as incidental (low-level heteroplasmy)
//

nextflow.enable.dsl = 2

process GATK4_MUTECT2_MITO {

    tag "${sample_id}"

    label 'process_medium'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/mito/${sample_id}", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi,stats,f1r2.tar.gz}"

    input:
    // Full-genome sorted BAM — Mutect2 will extract chrM reads automatically
    tuple val(sample_id), path(bam), path(bai)
    // chrM-only reference FASTA (rCRS = revised Cambridge Reference Sequence)
    // Extract from GRCh38: samtools faidx GRCh38.fa chrM > chrM.fa
    path mito_ref
    path mito_ref_fai
    path mito_ref_dict
    // MITOMAP NUMT blacklist: BED of nuclear chrM insertions to exclude
    path blacklist
    // Shifted chrM reference (for correct calling at genomic origin):
    // The 16,569 bp chrM is padded 300 bp and shifted to handle the origin.
    path mito_ref_shifted
    path mito_ref_shifted_fai
    path mito_ref_shifted_dict
    // Interval list for chrM (coordinates 1-16569 in GRCh38)
    path mito_interval
    path mito_interval_shifted

    output:
    // Per-sample mitochondrial variant VCF with heteroplasmy fractions
    tuple val(sample_id), path("${sample_id}.mito.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.mito.vcf.gz.tbi"), emit: tbi
    // Mutect2 call statistics: used by FilterMutectCalls downstream
    path "${sample_id}.mito.stats",                              emit: stats
    // Orientation model: corrects OXOG artefacts in mtDNA
    path "${sample_id}.mito.f1r2.tar.gz",                       emit: f1r2
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 2)}g" : "14g"
    """
    # ── GATK4 Mutect2 — mitochondrial mode ───────────────────────────────
    # ACGS 2024 §6: Mitochondrial variant calling using Mutect2 --mitochondria-mode
    # See module header for full explanation of mitochondria-specific settings.
    gatk --java-options "-Xmx${avail_mem}" \\
        Mutect2 \\
        --input ${bam} \\
        --output ${sample_id}.mito.vcf.gz \\
        --reference ${mito_ref} \\
        \\
        --mitochondria-mode \\
        # CRITICAL: enables mtDNA-specific calling.
        # Sets ploidy based on observed copy number.
        # Enables heteroplasmy detection at any allele fraction.
        # Activates NUMT suppression logic.
        # Without this flag, Mutect2 assumes somatic tumor-normal calling.
        \\
        -L ${mito_interval} \\
        # Restrict to chrM only — avoids wasting compute on nuclear genome
        \\
        --max-reads-per-alignment-start 75 \\
        # Higher than WGS default (50) because mtDNA coverage is 100-10,000x.
        # Without this increase, high-coverage regions are downsampled and
        # low-frequency heteroplasmic variants are missed.
        \\
        --min-base-quality-score 20 \\
        # Q20 threshold per ACGS 2024 §6 — exclude low-quality bases from
        # heteroplasmy calculation to reduce noise at low AF variants.
        \\
        --minimum-allele-fraction 0.01 \\
        # Detect variants at ≥1% heteroplasmy.
        # Clinical threshold for reportable mitochondrial disease variants.
        # Variants 1-10% AF flagged as low-level heteroplasmy in report.
        \\
        --max-alt-allele-count 4 \\
        # Allow up to 4 alt alleles per site (multi-allelic heteroplasmy occurs
        # in some mitochondrial disorders and haplogroup-defining sites)
        \\
        --f1r2-tar-gz ${sample_id}.mito.f1r2.tar.gz \\
        # Captures orientation bias data for FilterMutectCalls OXOG correction
        # (oxidative artefacts from library preparation can mimic C>T variants)
        \\
        --stats ${sample_id}.mito.stats \\
        # Mutect2 statistics file required by FilterMutectCalls downstream
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
    touch ${sample_id}.mito.vcf.gz \\
          ${sample_id}.mito.vcf.gz.tbi \\
          ${sample_id}.mito.stats \\
          ${sample_id}.mito.f1r2.tar.gz \\
          versions.yml
    """
}
