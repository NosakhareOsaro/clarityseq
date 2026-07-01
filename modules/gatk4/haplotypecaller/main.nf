// ============================================================================
// Module: GATK4_HAPLOTYPECALLER
// Description: Calls germline SNVs and indels using GATK4 HaplotypeCaller in
//              DRAGEN-GATK mode. Operates in GVCF mode (-ERC GVCF) for joint
//              genotyping across cohorts via GenomicsDBImport/GenotypGVCFs.
//              Scattered by interval list for parallelism. This is the primary
//              variant caller in the DRAGEN-GATK arm of ClaritySeq.
// Guidelines: GATK DRAGEN-GATK Best Practices (2021+);
//             ACGS Best Practice Guidelines v1.2 2024 §4;
//             Van der Auwera & O'Connor "Genomics in the Cloud" (2020) Ch.11
// Inputs:  bam_interval_ch — tuple(sample_id, markdup.bam, markdup.bam.bai,
//                            interval_list) — one channel item per interval
//          reference       — path to GRCh38 FASTA + .fai + .dict
// Outputs: gvcf_ch   — tuple(sample_id, interval_id, sample.g.vcf.gz)
//          tbi_ch    — tuple(sample_id, interval_id, sample.g.vcf.gz.tbi)
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939
//       https://gatk.broadinstitute.org/hc/en-us/articles/360037225632
// Parameter rationale:
//   --dragen-mode true: activates the BQD genotyping model (see below)
//   -ERC GVCF: emit reference confidence model for joint genotyping
//   -G StandardAnnotation: standard INFO/FORMAT annotations for VCF
//   -G AS_StandardAnnotation: allele-specific annotations (AS_QD, AS_FS, etc.)
//     required for VQSR allele-specific mode
//   --sample-ploidy 2: diploid human; change to 1 for mitochondria (separate module)
//   -L: limit to current interval (scattered execution)
//   --native-pair-hmm-threads: use all available CPUs for PairHMM
//   --max-reads-per-alignment-start 50: cap reads per position to prevent
//     excessive runtime at highly duplicated loci
// Version note: GATK 4.6.0.0 includes improved DRAGEN-GATK BQD model
//   calibration for gnomAD v4.1 allele frequencies. Pin to 4.6.0.0.
// ============================================================================
//
// DRAGEN-GATK MODE EXPLAINED:
// ============================
// --dragen-mode true activates the Base Quality Dropout (BQD) genotyping
// model in HaplotypeCaller. BQD models the pattern of base quality dropoff
// at read ends that characterises systematic Illumina sequencing errors.
// This replaces the function of BQSR (Base Quality Score Recalibration):
//
//   Traditional GATK best practice:
//     ALIGN (BWA) → BQSR (learn error model) → APPLYBQSR → HC
//
//   DRAGEN-GATK best practice:
//     ALIGN (DRAGMAP) → HC --dragen-mode true  ← NO BQSR STEP
//
// WHY NO BQSR IN DRAGEN-GATK MODE:
//   1. HaplotypeCaller --dragen-mode infers base-quality errors using BQD
//      internally during genotyping — it does not rely on pre-calibrated
//      quality scores from BQSR.
//   2. Applying BQSR to DRAGMAP output then passing to HC --dragen-mode
//      REDUCES accuracy because BQD expects raw DRAGMAP quality scores.
//   3. Broad benchmarks on GIAB truth sets show DRAGEN-GATK mode without
//      BQSR achieves F1 ≥ 0.999 for SNVs on NA12878 chr22.
//
// DO NOT PASS --bqsr-recal-file IN THIS MODULE.
// If you need BQSR (BWA-MEM2 profile), see modules/gatk4/applybqsr/main.nf
// and ensure --dragen-mode is NOT set in that profile.
//
// -ERC GVCF — JOINT GENOTYPING MODE:
//   HaplotypeCaller emits a GVCF (Genomic VCF) that includes reference
//   confidence intervals in addition to variant sites. This enables:
//   - Joint genotyping across multiple samples via GenomicsDBImport + GenotypeGVCFs
//   - More accurate genotyping at polymorphic sites seen in other cohort members
//   - Batch-mode processing: add new samples without re-running existing ones
//
// INTERVAL SCATTERING:
//   This module receives one interval_list per process invocation. The pipeline
//   creates ~50 scattered intervals from the GRCh38 callable region BED, then
//   groupTuple() merges per-sample GVCFs before GenomicsDBImport.
//

nextflow.enable.dsl = 2

process GATK4_HAPLOTYPECALLER {

    tag "${sample_id}:${interval_id}"

    label 'process_high'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/gvcfs/${sample_id}", mode: 'copy',
        pattern: "*.g.vcf.gz*"

    input:
    // One tuple per scattered interval: sample BAM + index + interval file
    tuple val(sample_id), path(bam), path(bai), val(interval_id), path(interval_list)
    // GRCh38 reference: FASTA + .fai (samtools index) + .dict (Picard dict)
    path reference
    path reference_fai
    path reference_dict

    output:
    // Per-interval GVCF for joint genotyping (grouped by sample downstream)
    tuple val(sample_id), val(interval_id), path("${sample_id}.${interval_id}.g.vcf.gz"),     emit: gvcf
    // TBI index for random access
    tuple val(sample_id), val(interval_id), path("${sample_id}.${interval_id}.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 2)}g" : "14g"
    """
    # ── GATK4 HaplotypeCaller — DRAGEN-GATK mode ──────────────────────────
    #
    # KEY FLAGS:
    # --dragen-mode true: activates BQD genotyping model (replaces BQSR)
    # -ERC GVCF: emit reference confidence blocks for joint genotyping
    # DO NOT pass --bqsr-recal-file in this module (DRAGEN-GATK mode)
    #
    gatk --java-options "-Xmx${avail_mem} -XX:ParallelGCThreads=4" \\
        HaplotypeCaller \\
        --input ${bam} \\
        --output ${sample_id}.${interval_id}.g.vcf.gz \\
        --reference ${reference} \\
        \\
        -L ${interval_list} \\
        # Restrict to current scatter interval for parallelism.
        # The pipeline scatters across ~50 GRCh38 callable-region intervals.
        \\
        --emit-ref-confidence GVCF \\
        # -ERC GVCF: required for joint genotyping. Emits reference confidence
        # blocks at non-variant sites so cohort genotyping can distinguish
        # "no call" (missing data) from "homozygous reference".
        \\
        --dragen-mode true \\
        # DRAGEN-GATK mode: activates BQD (Base Quality Dropoff) genotyping
        # model. This REPLACES BQSR — do NOT also run GATK4_APPLYBQSR on
        # DRAGMAP-aligned reads before passing to this process.
        \\
        -G StandardAnnotation \\
        # Standard annotations: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum
        # Used by VQSR training and ClinVar/ACMG interpretation
        \\
        -G AS_StandardAnnotation \\
        # Allele-specific variants of standard annotations (AS_QD, AS_FS, etc.)
        # Required for VQSR allele-specific mode (--AS in VQSR module)
        \\
        -G StandardHCAnnotation \\
        # HaplotypeCaller-specific: PGT, PID (phasing), RPA (repeat alleles)
        \\
        --sample-ploidy 2 \\
        # Diploid human genome. Use 1 for mitochondria (handled by separate
        # modules/gatk4/mutect2_mito/main.nf — NOT this module).
        \\
        --native-pair-hmm-threads ${task.cpus} \\
        # PairHMM likelihood calculation is parallelised with this thread count
        \\
        --max-reads-per-alignment-start 50 \\
        # Cap reads per position to prevent excessive runtime at repeat loci
        # or very high-coverage regions (e.g., mitochondria captured in WGS)
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
    touch ${sample_id}.${interval_id}.g.vcf.gz \\
          ${sample_id}.${interval_id}.g.vcf.gz.tbi \\
          versions.yml
    """
}
