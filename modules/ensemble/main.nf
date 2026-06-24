// ============================================================================
// Module: ENSEMBLE_MERGE
// Description: Merges GATK4 HaplotypeCaller (post-VQSR) and DeepVariant v1.8.0
//              VCFs into a single ensemble call set using modules/ensemble/
//              ensemble.py. In INTERSECTION mode (default), only variants called
//              PASS by BOTH callers are retained, maximising precision for
//              clinical reporting. In UNION mode, variants called PASS by EITHER
//              caller are retained, maximising sensitivity for research use.
//              ENSEMBLE_CALLER and ENSEMBLE_MODE INFO fields are added to each
//              record so downstream tools can filter or stratify by origin.
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §4.4 (ensemble calling);
//             ClinGen Variant Curation Expert Panel recommendations;
//             GA4GH Benchmarking Team (Zook et al. 2019 Nat Biotechnol
//             PMID:30936564)
// Inputs:  ensemble_ch — tuple(sample_id, gatk4_vcf, deepvariant_vcf)
//              sample_id     : cohort-unique sample identifier (string)
//              gatk4_vcf     : post-VQSR GATK4 SNP+indel VCF (plain or .gz)
//              deepvariant_vcf: DeepVariant v1.8.0 VCF (plain or .gz)
// Outputs: vcf — tuple(sample_id, ensemble.vcf.gz)
//          tbi — tuple(sample_id, ensemble.vcf.gz.tbi)
// Container: python:3.12-slim
//   Rationale: ensemble.py uses only the Python standard library (argparse,
//   gzip, pathlib, typing). A minimal python:3.12-slim image is sufficient —
//   no heavy bioinformatics container required. bgzip is installed via apt
//   at runtime so the output VCF is tabix-compatible.
//   python:3.12-slim is chosen over python:3.12-alpine because slim images
//   use glibc (required for compiled htslib binaries) while alpine uses musl.
// Docs: modules/ensemble/ensemble.py (inline docstring and comments)
//       https://github.com/google/deepvariant/blob/v1.8.0/docs/
//       https://gatk.broadinstitute.org/hc/en-us/articles/360035531112
// Parameter rationale:
//   --mode INTERSECTION: clinical default — see ensemble.py module docstring
//                        and ACGS 2024 §4.4.
//   bgzip -c: compress to stdout piped to file; avoids double-write.
//   tabix -p vcf: index the bgzip-compressed output for downstream processes
//                 that require random access (VEP, bcftools).
// Version note: python:3.12-slim — Python 3.12 LTS; slim variant avoids
//   unnecessary packages (test suite, tk, etc.) reducing image size ~40%.
//   Pin to 3.12-slim not 3-slim so the Python minor version is reproducible.
// ============================================================================

nextflow.enable.dsl = 2

process ENSEMBLE_MERGE {

    tag "${sample_id}"

    label 'process_low'

    // python:3.12-slim — standard library only; bgzip/tabix installed inline.
    // Pin to 3.12-slim — see Version note in header above.
    container 'python:3.12-slim'

    publishDir "${params.outdir}/ensemble/${sample_id}", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi}"

    input:
    // Joined channel from VQSR_INDEL.out.vcf and DEEPVARIANT_CALL.out.vcf
    // Both must be for the same sample_id (joined with .join() operator).
    tuple val(sample_id), path(gatk_vcf), path(dv_vcf)

    output:
    // Ensemble VCF: bgzip-compressed, tabix-indexed
    tuple val(sample_id), path("${sample_id}.ensemble.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.ensemble.vcf.gz.tbi"), emit: tbi
    // Plain-text log from ensemble.py (written to stderr, captured here)
    path "${sample_id}.ensemble.log",                                emit: log
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // params.ensemble_mode defaults to INTERSECTION if not set.
    def mode = params.ensemble_mode ?: "INTERSECTION"
    """
    # ── Install bgzip/tabix from htslib (only dependency beyond stdlib) ─────
    apt-get update -qq && apt-get install -y --no-install-recommends \
        tabix \
    > /dev/null 2>&1

    # ── Copy ensemble.py from the module directory into the work dir ─────────
    # Nextflow stages files referenced in 'path' inputs automatically.
    # ensemble.py is referenced via params.ensemble_script (set in nextflow.config).
    cp ${moduleDir}/ensemble.py ./ensemble.py

    # ── Run ensemble merge ───────────────────────────────────────────────────
    python3 ensemble.py \\
        --gatk-vcf  ${gatk_vcf} \\
        --dv-vcf    ${dv_vcf} \\
        --output    ${sample_id}.ensemble.vcf \\
        --mode      ${mode} \\
        2> ${sample_id}.ensemble.log

    # ── Compress and index output VCF ────────────────────────────────────────
    bgzip -c ${sample_id}.ensemble.vcf > ${sample_id}.ensemble.vcf.gz
    tabix -p vcf ${sample_id}.ensemble.vcf.gz

    # ── Versions ─────────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        bgzip: \$(bgzip --version | head -1 | awk '{print \$NF}')
        tabix: \$(tabix --version | head -1 | awk '{print \$NF}')
        ensemble_mode: "${mode}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.ensemble.vcf.gz \
          ${sample_id}.ensemble.vcf.gz.tbi \
          ${sample_id}.ensemble.log \
          versions.yml
    """
}
