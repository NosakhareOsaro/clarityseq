// ============================================================================
// Module: VQSR_SNP / VQSR_INDEL
// Description: Variant Quality Score Recalibration (VQSR) for germline SNPs
//              and indels. VQSR trains a Gaussian Mixture Model on known-true
//              variant sites (training resources) to learn the multivariate
//              annotation profile of real variants, then applies a sensitivity
//              tranche system to filter. Two separate processes are used because
//              SNPs and indels have different error profiles and require
//              different training resources with different prior weights.
//              When params.ancestry_vqsr=true, training resources are
//              ancestry-stratified gnomAD v4.1 subsets (AFR/AMR/EAS/EUR/SAS)
//              inferred by SOMALIER_ANCESTRY upstream. This improves calibration
//              for non-European individuals where pan-ethnic gnomAD can have
//              inflated AF at population-private variants.
// Guidelines: GATK Best Practices (VQSR);
//             ACGS Best Practice Guidelines v1.2 2024 §4.3;
//             gnomAD v4.1 (released April 19, 2024)
// Inputs:  vcf_ch          — merged cohort VCF from GATK4_GENOTYPEGVCFS
//          training_resources — paths to gnomAD/HapMap/Mills VCF resources
// Outputs: recal_ch    — recalibration table (.recal)
//          tranches_ch — tranche sensitivity table (.tranches)
//          rscript_ch  — R script for diagnostic plots
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360036510892
//       https://gnomad.broadinstitute.org/news/2024-05-gnomad-v4-1-updates/
// Parameter rationale: see inline comments on each --resource flag
// Version note: GATK 4.6.0.0; gnomAD v4.1 (April 2024; use v4.1 only).
// ============================================================================
//
// TRAINING RESOURCES AND PRIOR VALUES:
// =====================================
// Each training resource has a "prior" value representing the expected log10
// probability that a variant in that resource is a true positive. Higher prior
// = greater weight in training the GMM model.
//
// SNP resources:
//   HapMap 3.3.b37 lifted to hg38    prior=15.0  (highest: near-perfect truth)
//   1000G Omni2.5 lifted to hg38     prior=12.0  (high-confidence chip sites)
//   1000G Phase1 SNPs lifted to hg38 prior=10.0  (broad population variation)
//   dbSNP b156 hg38                  prior=2.0   (lowest: many false positives)
//   gnomAD v4.1 sites (non-neuro)    prior=6.0   (large, well-curated)
//
// Indel resources:
//   Mills & 1000G gold standard      prior=12.0  (curated indel truth set)
//   Axiom Exome Plus v1              prior=10.0  (exome chip indels)
//   dbSNP b156 indels                prior=2.0   (noisy; low prior)
//   gnomAD v4.1 indels (non-neuro)   prior=6.0
//
// GNOMAD v4.1 NOTE:
//   IMPORTANT: Use gnomAD v4.1 (released April 2024); prior releases had an
//   allele number (AN) bug causing underestimated frequencies in certain
//   callsets. gnomAD v4.1 corrects AN values for all variants.
//   v4.1 contains 807,162 individuals (730,947 exomes + 76,215 genomes
//   including 416,555 UK Biobank samples).
//   Reference: https://gnomad.broadinstitute.org/news/2024-05-gnomad-v4-1-updates/
//
// PER-ANCESTRY VQSR (when ancestry_vqsr=true):
// ==============================================
// Standard VQSR uses pan-ethnic gnomAD as training resource. This can
// mis-calibrate for non-European samples where:
//   - Population-private variants appear at low AF in pan-ethnic gnomAD
//   - The GMM model may classify these as low-quality sites
//
// When ancestry_vqsr=true, the pipeline uses SOMALIER_ANCESTRY output
// to select ancestry-specific gnomAD v4.1 subsets:
//   AFR → gnomAD v4.1 African/African American subset
//   AMR → gnomAD v4.1 Latino/Admixed American subset
//   EAS → gnomAD v4.1 East Asian subset
//   EUR → gnomAD v4.1 European (non-Finnish) subset
//   SAS → gnomAD v4.1 South Asian subset
//
// Fallback for admixed/unassigned samples:
//   If SOMALIER assigns probability <0.90 to any single ancestry,
//   the sample falls back to pan-ethnic gnomAD v4.1 (conservative choice).
//   This is logged with a WARNING in the MultiQC report.
//

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Process 1: VQSR for SNPs
// ─────────────────────────────────────────────────────────────────────────────
process VQSR_SNP {

    tag "SNP_VQSR"

    label 'process_high_memory'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/vqsr/snp", mode: 'copy'

    input:
    // Merged multi-sample cohort VCF from GenotypeGVCFs (all intervals merged)
    tuple path(vcf), path(tbi)
    // Reference genome
    path reference
    path reference_fai
    path reference_dict
    // Training resources for SNP VQSR (ancestry-stratified if ancestry_vqsr=true)
    path hapmap          // HapMap 3.3 hg38 liftover  — prior 15
    path hapmap_tbi
    path omni            // 1000G Omni 2.5 hg38        — prior 12
    path omni_tbi
    path onekg_snps      // 1000G Phase1 SNPs hg38     — prior 10
    path onekg_snps_tbi
    path dbsnp           // dbSNP b156 hg38            — prior 2
    path dbsnp_tbi
    path gnomad          // gnomAD v4.1 non-neuro hg38 — prior 6
    path gnomad_tbi
    // Ancestry label (string) from SOMALIER_ANCESTRY: AFR/AMR/EAS/EUR/SAS/mixed
    val ancestry

    output:
    // Recalibration table: SNP-mode scores for ApplyVQSR
    path "snp.recal",     emit: recal
    // Tranche sensitivity table: maps VQSLOD to Ti/Tv and sensitivity
    path "snp.tranches",  emit: tranches
    // R script for VQSR diagnostic plots (Ti/Tv vs tranche)
    path "snp.plots.R",   emit: rscript
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 4)}g" : "28g"
    // Log ancestry selection for audit trail
    def ancestry_note = ancestry ?: "mixed (pan-ethnic fallback)"
    """
    echo "INFO: SNP VQSR using ancestry=${ancestry_note} gnomAD v4.1 training set" >&2
    # If ancestry_vqsr=true and ancestry is known, gnomAD path is the
    # ancestry-specific subset selected by the pipeline orchestrator.
    # If ancestry=mixed, gnomAD path is the pan-ethnic gnomAD v4.1 non-neuro.

    gatk --java-options "-Xmx${avail_mem}" \\
        VariantRecalibrator \\
        --variant ${vcf} \\
        --output snp.recal \\
        --tranches-file snp.tranches \\
        --rscript-file snp.plots.R \\
        --reference ${reference} \\
        \\
        --mode SNP \\
        # SNP mode: trains on Ti/Tv ratio; different model from INDEL mode
        \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \\
        # HapMap 3.3: highest-confidence SNP truth set for training.
        # prior=15.0 means P(true variant) ≈ 10^15 / (10^15 + 1) ≈ 1.000
        # known=false: does not affect dbSNP ID lookup
        # training=true: include in GMM model training
        # truth=true: count as concordant truth sites in tranche calculation
        \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \\
        # 1000G Omni 2.5 chip genotypes: high-confidence but not infallible.
        # truth=false: not used for sensitivity calculation (chip artefacts present)
        # prior=12.0: very high confidence but slightly below HapMap
        \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekg_snps} \\
        # 1000G Phase1 high-confidence SNPs: broad population variation.
        # prior=10.0: confident but population-level (not single-sample truth)
        \\
        --resource:dbsnp,known=true,training=false,training=false,truth=false,prior=2.0 ${dbsnp} \\
        # dbSNP b156: used for rsID annotation only (known=true).
        # training=false: NOT used to train the model (contains many false positives)
        # prior=2.0: minimal weight — presence in dbSNP is weak evidence
        \\
        --resource:gnomad,known=false,training=true,truth=false,prior=6.0 ${gnomad} \\
        # gnomAD v4.1 non-neuro sites (807,162 individuals).
        # Ancestry-stratified subset selected by SOMALIER_ANCESTRY if ancestry_vqsr=true.
        # IMPORTANT: v4.1 ONLY — prior releases had AN calculation bug (see header).
        # prior=6.0: moderate confidence — large but not perfectly curated
        \\
        --use-allele-specific-annotations \\
        # Use AS_* annotations (AS_QD, AS_FS, AS_MQ) for multi-allelic handling.
        # Requires -G AS_StandardAnnotation in HaplotypeCaller (already set).
        \\
        -an QD \\
        # QualByDepth: variant confidence / depth — key discriminant for real variants
        \\
        -an MQ \\
        # RMS Mapping Quality: high MQ = reads map well to variant region
        \\
        -an MQRankSum \\
        # Rank sum test: compares MQ of ref-supporting vs alt-supporting reads
        \\
        -an ReadPosRankSum \\
        # Rank sum: compares position of variant within reads (end-of-read artefacts)
        \\
        -an FS \\
        # FisherStrand: strand bias as Phred score; artefacts are strand-biased
        \\
        -an SOR \\
        # StrandOddsRatio: improved strand bias metric over FS for high-depth data
        \\
        --truth-sensitivity-tranche 100.0 \\
        --truth-sensitivity-tranche 99.95 \\
        --truth-sensitivity-tranche 99.9 \\
        --truth-sensitivity-tranche 99.5 \\
        --truth-sensitivity-tranche 99.0 \\
        --truth-sensitivity-tranche 95.0 \\
        # Six tranches from 95% to 100% sensitivity.
        # Clinical pipeline applies PASS filter at the 99.5% tranche:
        # balances sensitivity (captures rare pathogenic SNVs) vs
        # precision (avoids false positive clinical reports).
        \\
        --max-gaussians 6 \\
        # Number of Gaussian components in mixture model.
        # 6 is GATK recommendation for cohorts ≥30 samples.
        # Use 4 for smaller cohorts if model fails to converge.
        \\
        --tmp-dir /tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -o 'GATK v[0-9.]*' | sed 's/GATK v//')
    END_VERSIONS
    """

    stub:
    """
    touch snp.recal snp.tranches snp.plots.R versions.yml
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Process 2: VQSR for Indels
// ─────────────────────────────────────────────────────────────────────────────
process VQSR_INDEL {

    tag "INDEL_VQSR"

    label 'process_high_memory'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/vqsr/indel", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path reference
    path reference_fai
    path reference_dict
    // Training resources for Indel VQSR
    path mills          // Mills & 1000G gold standard indels — prior 12
    path mills_tbi
    path axiom          // Axiom Exome Plus v1 — prior 10
    path axiom_tbi
    path dbsnp          // dbSNP b156 indels — prior 2
    path dbsnp_tbi
    path gnomad         // gnomAD v4.1 indels — prior 6
    path gnomad_tbi
    val ancestry

    output:
    path "indel.recal",    emit: recal
    path "indel.tranches", emit: tranches
    path "indel.plots.R",  emit: rscript
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 4)}g" : "28g"
    def ancestry_note = ancestry ?: "mixed (pan-ethnic fallback)"
    """
    echo "INFO: Indel VQSR using ancestry=${ancestry_note} gnomAD v4.1 training set" >&2

    gatk --java-options "-Xmx${avail_mem}" \\
        VariantRecalibrator \\
        --variant ${vcf} \\
        --output indel.recal \\
        --tranches-file indel.tranches \\
        --rscript-file indel.plots.R \\
        --reference ${reference} \\
        \\
        --mode INDEL \\
        # INDEL mode: uses different annotation profile from SNP mode.
        # Indels have lower Ti/Tv signal; model relies more on QD, FS, SOR.
        \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \\
        # Mills & 1000G gold standard indels (hg38 liftover).
        # Curated set of high-confidence insertion/deletion variants.
        # truth=true: highest sensitivity target for indel calling.
        # prior=12.0: very high confidence — manually validated indels.
        \\
        --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${axiom} \\
        # Axiom Exome Plus genotyping array: captures common coding indels.
        # prior=10.0: high confidence for common indels in exonic regions.
        # Not used as truth (array-based, not sequencing-validated).
        \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \\
        # dbSNP indels: rsID annotation only.
        # training=false, truth=false: too noisy for model training.
        # prior=2.0: minimal weight.
        \\
        --resource:gnomad,known=false,training=true,truth=false,prior=6.0 ${gnomad} \\
        # gnomAD v4.1 indels (ancestry-stratified if ancestry_vqsr=true).
        # Includes 76,215 genome samples with deep coverage for indel calling.
        # prior=6.0: moderate confidence.
        \\
        --use-allele-specific-annotations \\
        \\
        -an QD \\
        -an DP \\
        # DP (depth): included for indels — low depth = low confidence.
        # Not used for SNPs because aggregate DP is less informative there.
        \\
        -an FS \\
        -an SOR \\
        -an ReadPosRankSum \\
        -an MQRankSum \\
        \\
        --truth-sensitivity-tranche 100.0 \\
        --truth-sensitivity-tranche 99.9 \\
        --truth-sensitivity-tranche 99.0 \\
        --truth-sensitivity-tranche 95.0 \\
        --truth-sensitivity-tranche 90.0 \\
        # Indels use 90%-100% tranches.
        # Clinical filter applied at 99.0% indel sensitivity tranche.
        # Indels are harder to model than SNPs; 99% is more conservative
        # than the 99.5% used for SNPs.
        \\
        --max-gaussians 4 \\
        # Fewer Gaussians than SNPs: indel space is smaller and less well-sampled.
        # 4 Gaussians avoids overfitting on small cohorts.
        \\
        --tmp-dir /tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -o 'GATK v[0-9.]*' | sed 's/GATK v//')
    END_VERSIONS
    """

    stub:
    """
    touch indel.recal indel.tranches indel.plots.R versions.yml
    """
}
