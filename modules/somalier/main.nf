// ============================================================================
// Module: SOMALIER_RELATE / SOMALIER_ANCESTRY
// Description: Performs sample quality-control checks using Somalier v0.2.19:
//
//   SOMALIER_RELATE  — pairwise relatedness estimation and sample-swap
//              detection. Extracts genotypes at ~17,000 high-quality SNP sites
//              and computes kinship coefficients between all sample pairs.
//              Flags unexpected duplicates (kinship≥0.354), unexpected
//              first-degree relatives (0.177–0.354), and sample swaps
//              (expected family structure not matched). ACGS 2024 mandates
//              relatedness QC for all clinical WGS to detect laboratory
//              mix-ups before variant reporting.
//
//   SOMALIER_ANCESTRY — projects sample onto a principal-component space
//              derived from 1000 Genomes superpopulations (AFR/AMR/EAS/EUR/SAS)
//              and assigns the most likely ancestry label (or "admixed" when no
//              superpopulation is clearly dominant). Ancestry label is consumed
//              by VQSR to select the per-ancestry training resource weights
//              (see conf/dragen_gatk.config for ancestry→VQSR mapping).
//
//              Per-ancestry VQSR mapping (GATK4 training resource weights):
//                AFR: higher weight on AFR-specific 1000G + gnomAD-AFR sites
//                AMR: balanced 1000G AMR + gnomAD-AMR; reduced EUR weight
//                EAS: EAS-specific 1000G + gnomAD-EAS; reduced hapmap3 weight
//                EUR: standard GATK4 VQSR resources (1000G + hapmap3)
//                SAS: SAS-specific 1000G + gnomAD-SAS sites
//                admixed: use all-population resources at equal weight (safest
//                         default when ancestry is uncertain)
//              Rationale: VQSR trained exclusively on EUR data has lower
//              sensitivity and specificity in non-EUR populations (Chen et al.
//              2019 Am J Hum Genet PMID:31564432). Per-ancestry resources
//              improve Ti/Tv and indel size distributions in non-EUR samples.
//
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §3.2 (sample QC),
//             §4.3 (per-ancestry VQSR);
//             GA4GH Genomic Data Infrastructure recommendations (2023)
// Inputs (SOMALIER_RELATE):
//   bam_ch     — tuple(sample_id, bam, bai) — all samples in the cohort
//   sites_vcf  — Somalier sites VCF for GRCh38 (~17,000 SNP positions)
//              Download: https://github.com/brentp/somalier#sites-files
//   reference  — GRCh38 FASTA + .fai
// Inputs (SOMALIER_ANCESTRY):
//   somalier_ch — directory of *.somalier files from SOMALIER_RELATE
//   labels_tsv  — 1000G labels file (sample→superpopulation)
//              Bundled in brentp/somalier:v0.2.19 at /opt/1kg-samples.tsv
//   pc_tsv      — 1000G PCA file (pre-computed principal components)
//              Bundled in brentp/somalier:v0.2.19 at /opt/1kg-pcs.tsv
// Outputs (SOMALIER_RELATE):
//   pairs       — tuple(cohort_id, pairs.tsv) — pairwise kinship table
//   samples     — tuple(cohort_id, samples.tsv) — per-sample QC metrics
//   html        — tuple(cohort_id, relate.html)
// Outputs (SOMALIER_ANCESTRY):
//   ancestry    — tuple(cohort_id, ancestry.tsv) [sample, predicted_ancestry]
//   pc_plot     — tuple(cohort_id, ancestry.html)
// Container: brentp/somalier:v0.2.19
//   Official Somalier container maintained by Brent Pedersen (author).
//   Includes bundled 1000G PCA reference files.
//   Source: https://hub.docker.com/r/brentp/somalier
// Docs: https://github.com/brentp/somalier
//       Pedersen et al. 2020 Genome Med PMID:32664994
// Parameter rationale:
//   somalier extract: default sites VCF covers LD-pruned bi-allelic SNPs with
//       MAF>0.05 in all 5 superpopulations. --sample-prefix not used (sample
//       names read from BAM SM tag).
//   somalier relate --ped: when PED file provided, expected vs observed
//       kinship comparisons are made. Flag --infer when PED absent.
//   somalier ancestry --n-pcs 5: 5 PCs sufficient to assign superpopulation;
//       additional PCs capture within-population structure not needed for VQSR.
//   --labels: 1000G superpopulation labels required for PC-projection ancestry.
// Version note: v0.2.19 (2023) — latest stable release as of 2024-Q4.
//   v0.2.18 → v0.2.19 fixes ancestry misassignment in admixed samples when
//   PC-projection is borderline (within 0.5σ of multiple centroids).
//   Pin to v0.2.19 — the bundled 1000G PCA coordinates must match the
//   software version.
// ============================================================================

nextflow.enable.dsl = 2

process SOMALIER_RELATE {

    tag "${cohort_id}"

    label 'process_medium'

    // brentp/somalier:v0.2.19 — includes bundled 1000G reference files.
    // Pin to exact version — PCA coordinates bundled in image are version-specific.
    container 'brentp/somalier:v0.2.19'

    publishDir "${params.outdir}/somalier/relate", mode: 'copy',
        pattern: "*.{tsv,html}"

    input:
    // All per-sample BAMs collected for cohort-level relatedness check
    tuple val(cohort_id), path(bams), path(bais)
    // Somalier sites VCF for GRCh38 (17k SNP positions)
    path sites_vcf
    // GRCh38 FASTA + index
    path reference
    path reference_fai
    // PED file (optional — pass 'NO_FILE' if unavailable)
    path ped_file

    output:
    // Pairwise kinship table: columns sample_a, sample_b, kinship, n, hom_concordance, etc.
    tuple val(cohort_id), path("relate/${cohort_id}.pairs.tsv"),   emit: pairs
    // Per-sample QC: depth, het/hom ratio, predicted sex, etc.
    tuple val(cohort_id), path("relate/${cohort_id}.samples.tsv"), emit: samples
    // Interactive HTML relatedness plot
    tuple val(cohort_id), path("relate/${cohort_id}.html"),        emit: html
    // Individual *.somalier binary files (for ancestry step)
    path "somalier_files/*.somalier",                              emit: somalier_files
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ped_arg = (ped_file.name != "NO_FILE") ? "--ped ${ped_file}" : "--infer"
    """
    # ── Step 1: Extract somalier sites from each BAM ─────────────────────
    # Extracts read counts at ~17,000 high-quality SNP positions.
    # Output: one *.somalier binary file per sample.
    mkdir -p somalier_files

    for bam in ${bams.join(' ')}; do
        somalier extract \\
            --sites        ${sites_vcf} \\
            --fasta        ${reference} \\
            --out-dir      somalier_files/ \\
            "\${bam}"
    done

    # ── Step 2: Pairwise relatedness calculation ──────────────────────────
    # Computes kinship coefficient (Φ) for all sample pairs.
    # Kinship thresholds (ACGS 2024 §3.2):
    #   Φ ≥ 0.354: duplicate / monozygotic twins (flag as possible swap)
    #   Φ 0.177–0.354: first-degree relatives (parent-child, full sibling)
    #   Φ 0.088–0.177: second-degree relatives
    #   Φ < 0.088: unrelated
    mkdir -p relate

    somalier relate \\
        somalier_files/*.somalier \\
        ${ped_arg} \\
        --output-prefix relate/${cohort_id}

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier version 2>&1 | grep -o 'v[0-9.]*')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somalier_files relate
    touch somalier_files/${cohort_id}.somalier \
          relate/${cohort_id}.pairs.tsv \
          relate/${cohort_id}.samples.tsv \
          relate/${cohort_id}.html \
          versions.yml
    """
}


process SOMALIER_ANCESTRY {

    tag "${cohort_id}"

    label 'process_low'

    // brentp/somalier:v0.2.19 — 1000G PCA reference bundled in image.
    container 'brentp/somalier:v0.2.19'

    publishDir "${params.outdir}/somalier/ancestry", mode: 'copy',
        pattern: "*.{tsv,html}"

    input:
    // cohort_id + directory of *.somalier files from SOMALIER_RELATE
    tuple val(cohort_id), path(somalier_files)
    // 1000G superpopulation labels: columns IID, SuperPop (AFR/AMR/EAS/EUR/SAS)
    // Bundled in container at /opt/1kg-samples.tsv; override with params.somalier_labels
    path labels_tsv
    // 1000G pre-computed PCA file (PC1..PC10 for each 1000G sample)
    // Bundled in container at /opt/1kg-pcs.tsv; override with params.somalier_pcs
    path pc_tsv

    output:
    // ancestry.tsv: columns sample, predicted_ancestry (AFR/AMR/EAS/EUR/SAS/admixed),
    //               PC1..PC10, and confidence scores per superpopulation
    // Consumed by VQSR_SNP and VQSR_INDEL to select per-ancestry training resources.
    // Per-ancestry VQSR mapping (see module header for full rationale):
    //   AFR → higher weight on AFR-specific 1000G + gnomAD-AFR
    //   AMR → balanced 1000G AMR + gnomAD-AMR; reduced EUR weight
    //   EAS → EAS-specific 1000G + gnomAD-EAS; reduced hapmap3 weight
    //   EUR → standard GATK4 VQSR resources (1000G + hapmap3)
    //   SAS → SAS-specific 1000G + gnomAD-SAS
    //   admixed → all-population resources at equal weight
    tuple val(cohort_id), path("${cohort_id}.ancestry.tsv"), emit: ancestry
    // HTML interactive PCA plot: query samples projected onto 1000G PCs
    tuple val(cohort_id), path("${cohort_id}.ancestry.html"), emit: pc_plot
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── Ancestry inference via 1000G PC projection ────────────────────────
    # Projects query samples onto PCs computed from 1000G reference panel.
    # Assigns superpopulation (AFR/AMR/EAS/EUR/SAS) or "admixed" when the
    # nearest centroid is within 0.5σ of two or more superpopulations.
    #
    # --n-pcs 5: five PCs are sufficient to distinguish five 1000G superpops.
    # Additional PCs capture within-population substructure not needed for VQSR.

    somalier ancestry \\
        ${somalier_files.join(' ')} \\
        --labels  ${labels_tsv} \\
        --pc-file ${pc_tsv} \\
        --n-pcs   5 \\
        --output-prefix ${cohort_id}

    # Rename output to expected file names
    mv ${cohort_id}.somalier-ancestry.tsv  ${cohort_id}.ancestry.tsv  2>/dev/null || true
    mv ${cohort_id}.somalier-ancestry.html ${cohort_id}.ancestry.html 2>/dev/null || true

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier version 2>&1 | grep -o 'v[0-9.]*')
    END_VERSIONS
    """

    stub:
    """
    touch ${cohort_id}.ancestry.tsv \
          ${cohort_id}.ancestry.html \
          versions.yml
    """
}
