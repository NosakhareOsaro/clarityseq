// ============================================================================
// Module: GATK4_GENOTYPEGVCFS
// Description: Performs joint genotyping across all samples in a cohort using
//              GATK4 GenotypeGVCFs. Takes a GenomicsDB workspace (created by
//              GenomicsDBImport from per-sample GVCFs) and produces a
//              multi-sample VCF with genotype calls at all variant sites.
//              Joint genotyping improves sensitivity at low-frequency variants
//              and provides more accurate allele frequency estimates than
//              single-sample calling.
// Guidelines: GATK Best Practices (joint genotyping pipeline);
//             ACGS Best Practice Guidelines v1.2 2024 §4.2;
//             Poplin et al. 2018 Nature Genetics PMID:29431740
// Inputs:  genomicsdb_ch — tuple(interval_id, path to GenomicsDB workspace)
//          reference     — GRCh38 FASTA + .fai + .dict
//          dbsnp         — dbSNP VCF for rs-ID annotation (b156 or later)
// Outputs: vcf_ch — tuple(interval_id, genotyped.vcf.gz)
//          tbi_ch — tuple(interval_id, genotyped.vcf.gz.tbi)
// Container: broadinstitute/gatk:4.6.0.0
// Docs: https://gatk.broadinstitute.org/hc/en-us/articles/360036711071
// Parameter rationale:
//   --dbsnp: annotate with rsIDs; required for VQSR training resource matching
//   -G StandardAnnotation: standard site-level annotations in output VCF
//   -G AS_StandardAnnotation: allele-specific annotations for VQSR AS mode
//   --only-output-calls-starting-in-intervals: restrict output to current
//     interval boundaries (prevents duplicate records at interval edges)
//   --heterozygosity 0.001: prior for human germline heterozygosity (1/1000)
//     per Li 2011 genome-wide estimate; used in genotype likelihood calculation
//   --indel-heterozygosity 0.000125: GATK default prior for indel rates
// Version note: GATK 4.6.0.0. Includes improved GenomicsDB consolidation
//   and allele-specific annotation accuracy. Pin to 4.6.0.0.
// ============================================================================

nextflow.enable.dsl = 2

process GATK4_GENOTYPEGVCFS {

    tag "${interval_id}"

    label 'process_high'

    // GATK 4.6.0.0 — pin exact version, never 'latest'
    container 'broadinstitute/gatk:4.6.0.0'

    publishDir "${params.outdir}/genotyped", mode: 'copy',
        pattern: "*.genotyped.vcf.gz*"

    input:
    // GenomicsDB workspace containing GVCFs for all samples at this interval
    // Created by GATK4_GENOMICSDBIMPORT (upstream in wgs_grch38.nf)
    tuple val(interval_id), path(genomicsdb)
    // Reference genome files (FASTA + fai + dict required by GATK)
    path reference
    path reference_fai
    path reference_dict
    // dbSNP VCF: used to annotate variant sites with rsIDs
    // Use dbSNP b156 (hg38) or later; required for VQSR resource matching
    path dbsnp
    path dbsnp_tbi

    output:
    // Per-interval genotyped VCF with all samples and all variant sites
    tuple val(interval_id), path("${interval_id}.genotyped.vcf.gz"),     emit: vcf
    tuple val(interval_id), path("${interval_id}.genotyped.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = task.memory ? "${(task.memory.toGiga() - 2)}g" : "14g"
    """
    # ── GATK4 GenotypeGVCFs — joint genotyping ────────────────────────────
    # Reads all GVCFs from the GenomicsDB workspace and emits genotype
    # calls for every sample at every polymorphic site in the cohort.
    gatk --java-options "-Xmx${avail_mem} -XX:ParallelGCThreads=4" \\
        GenotypeGVCFs \\
        --variant gendb://${genomicsdb} \\
        # Input: GenomicsDB workspace URI (gendb:// prefix required)
        \\
        --output ${interval_id}.genotyped.vcf.gz \\
        --reference ${reference} \\
        \\
        --dbsnp ${dbsnp} \\
        # Annotate variants with rsIDs from dbSNP.
        # Required for VQSR to match training resource variants by rsID.
        # Use b156 (hg38) — matches GATK resource bundle.
        \\
        -G StandardAnnotation \\
        # QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum — standard VQSR features
        \\
        -G AS_StandardAnnotation \\
        # Allele-specific: AS_QD, AS_FS, AS_MQ — required for VQSR AS mode
        \\
        --only-output-calls-starting-in-intervals \\
        # Prevents duplicate records at interval boundaries when scattering.
        # Each interval emits only variants whose POS falls within the interval.
        \\
        --heterozygosity 0.001 \\
        # Prior probability of heterozygosity at any site (Li 2011 estimate).
        # Used in the Bayesian genotype likelihood calculation.
        \\
        --indel-heterozygosity 0.000125 \\
        # Prior for indel rate (GATK default); lower than SNV rate
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
    touch ${interval_id}.genotyped.vcf.gz \\
          ${interval_id}.genotyped.vcf.gz.tbi \\
          versions.yml
    """
}
