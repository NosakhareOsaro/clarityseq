// ============================================================================
// Module: ALPHAMISSENSE_LOOKUP
// Description: Performs tabix-indexed lookup of pre-computed AlphaMissense
//              pathogenicity scores for missense variants in a VCF. AlphaMissense
//              is a deep learning model from Google DeepMind that predicts the
//              clinical pathogenicity of amino acid substitutions using
//              protein structure information from AlphaFold2. Scores are added
//              to the VCF INFO field and evaluated against ClinGen SVI 2024
//              approved thresholds for ACMG PP3/BP4 criteria assignment.
// Guidelines: ClinGen SVI Pathogenicity Calibration (2024 recommendations);
//             ACMG/AMP 2015 variant classification PP3/BP4 criteria;
//             Cheng et al. 2023 Science PMID:37703350 (AlphaMissense paper)
// Inputs:  vcf_ch            — tuple(sample_id, vep.vcf.gz) from VEP_ANNOTATE
//          alphamissense_tsv — pre-computed AlphaMissense TSV (tabix-indexed)
//                             Download: gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz
//                             Size: ~2.7 GB; must be tabix-indexed
// Outputs: am_vcf_ch — tuple(sample_id, am_annotated.vcf.gz) — AM scores in INFO
//          am_tbi_ch — tuple(sample_id, am_annotated.vcf.gz.tbi)
// Container: biocontainers/htslib:1.19
// Docs: https://github.com/google-deepmind/alphamissense
//       https://doi.org/10.1126/science.adg7492 (Cheng et al. 2023 Science)
//       https://clinicalgenome.org/working-groups/sequence-variant-interpretation/
// Parameter rationale:
//   tabix: perform random-access lookup into pre-computed score TSV
//   The tabix approach avoids re-running AlphaMissense (which requires
//   AlphaFold2 protein structure inputs) for every pipeline run.
//   The pre-computed TSV covers all possible missense variants in canonical
//   transcripts across the human proteome.
// Version note: htslib 1.19 (Jan 2024) includes tabix improvements for
//   large TSV lookups. Pin to 1.19 — use biocontainers registry for
//   reproducibility. AlphaMissense scores themselves are database-versioned
//   separately from this container.
// ============================================================================
//
// ALPHAMISSENSE SCORES AND CLINGEN SVI 2024 THRESHOLDS:
// =======================================================
// AlphaMissense predicts pathogenicity on a continuous scale [0, 1]:
//   Score → Predicted class → ACMG evidence code (ClinGen SVI 2024):
//
//   ≥ 0.564 → "likely_pathogenic" → PP3 (pathogenic supporting evidence)
//   0.341 to 0.563 → "ambiguous"  → NO ACMG code assigned
//   ≤ 0.340 → "likely_benign"    → BP4 (benign supporting evidence)
//
// These thresholds were calibrated by ClinGen SVI in 2024 to achieve:
//   - Specificity ≥ 0.99 for PP3 (≥0.564 threshold)
//   - Specificity ≥ 0.99 for BP4 (≤0.340 threshold)
// Against ClinVar pathogenic/benign variants with ≥2-star review status.
//
// INFO FIELDS ADDED TO VCF:
//   AM_SCORE: AlphaMissense score (0-1, two decimal precision)
//   AM_CLASS: likely_pathogenic / ambiguous / likely_benign
//   AM_PP3: 1 if score ≥ 0.564; 0 otherwise
//   AM_BP4: 1 if score ≤ 0.340; 0 otherwise
//
// REFERENCE:
//   Cheng et al. 2023 Science "Accurate proteome-wide missense variant
//   effect prediction with AlphaMissense" PMID:37703350
//   Model coverage: 71 million missense variants across 19,233 human proteins
//   Concordance with ClinVar: 99% (P/LP sites) and 97% (B/LB sites)
//
// NOTE: AlphaMissense scores are ALSO provided by the VEP AlphaMissense plugin
// (modules/vep/main.nf). This dedicated module provides standalone lookup for:
//   1. Samples where VEP is not run (research-only pipelines)
//   2. Cross-validation of VEP plugin output
//   3. Bulk batch lookup without full VEP annotation overhead
//

nextflow.enable.dsl = 2

process ALPHAMISSENSE_LOOKUP {

    tag "${sample_id}"

    label 'process_low'

    // biocontainers/htslib:1.19 — provides tabix and bgzip.
    // htslib 1.19 (Jan 2024): tabix improvements for large indexed files.
    // Pin to 1.19 — never 'latest'. Container from BioContainers registry.
    container 'biocontainers/htslib:1.19--h81da01d_1'

    publishDir "${params.outdir}/annotation/${sample_id}", mode: 'copy',
        pattern: "*.am_annotated.vcf.gz*"

    input:
    // VEP-annotated VCF (or ensemble VCF if VEP not run)
    tuple val(sample_id), path(vcf), path(tbi)
    // AlphaMissense pre-computed scores TSV (tabix-indexed)
    // Format: chr  pos  ref  alt  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
    // Download from: gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz (~2.7 GB)
    // Build tabix index: tabix -s 1 -b 2 -e 2 -c '#' AlphaMissense_hg38.tsv.gz
    path alphamissense_tsv
    path alphamissense_tbi

    output:
    // VCF with AM_SCORE, AM_CLASS, AM_PP3, AM_BP4 added to INFO column
    tuple val(sample_id), path("${sample_id}.am_annotated.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.am_annotated.vcf.gz.tbi"), emit: tbi
    // Statistics: count of PP3, BP4, ambiguous assignments
    path "${sample_id}.am_stats.txt",                                    emit: stats
    path "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ClinGen SVI 2024 approved thresholds
    def pp3_threshold = 0.564   // ≥ this score → PP3 (likely pathogenic)
    def bp4_threshold = 0.340   // ≤ this score → BP4 (likely benign)
    """
    # ── AlphaMissense tabix lookup + INFO field annotation ────────────────
    # ClinGen SVI 2024 approved thresholds:
    #   score ≥ ${pp3_threshold} → PP3 (likely_pathogenic)
    #   score ≤ ${bp4_threshold} → BP4 (likely_benign)
    #   ${bp4_threshold} < score < ${pp3_threshold} → ambiguous (no ACMG code)
    # Reference: Cheng et al. 2023 PMID:37703350

    # ── Add VCF header definitions for new INFO fields ─────────────────────
    bcftools view -h ${vcf} > header.txt
    cat >> header.txt << 'HEADER'
##INFO=<ID=AM_SCORE,Number=1,Type=Float,Description="AlphaMissense pathogenicity score (0-1); Cheng et al. 2023 PMID:37703350">
##INFO=<ID=AM_CLASS,Number=1,Type=String,Description="AlphaMissense class: likely_pathogenic / ambiguous / likely_benign">
##INFO=<ID=AM_PP3,Number=0,Type=Flag,Description="AlphaMissense score >= 0.564; supports PP3 per ClinGen SVI 2024">
##INFO=<ID=AM_BP4,Number=0,Type=Flag,Description="AlphaMissense score <= 0.340; supports BP4 per ClinGen SVI 2024">
HEADER

    # ── Python lookup script: annotate each missense variant ──────────────
    python3 << 'PYEOF'
import gzip
import sys
import subprocess
import re

def am_lookup(chrom, pos, ref, alt, am_tsv):
    """Look up AlphaMissense score using tabix for a single variant."""
    cmd = ["tabix", am_tsv, f"{chrom}:{pos}-{pos}"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split("\\n"):
            if not line:
                continue
            fields = line.split("\\t")
            # TSV columns: chr, pos, ref, alt, uniprot_id, transcript_id,
            #              protein_variant, am_pathogenicity, am_class
            if len(fields) >= 9 and fields[2] == ref and fields[3] == alt:
                return float(fields[7]), fields[8]
    except subprocess.CalledProcessError:
        pass
    return None, None

# Process VCF
pp3_count = 0
bp4_count = 0
ambiguous_count = 0
no_score_count = 0

with gzip.open("${vcf}", "rt") as fin, open("${sample_id}.am_annotated.unsorted.vcf", "w") as fout:
    for line in fin:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                # Insert new INFO headers before the column header line
                with open("header.txt") as hf:
                    for hline in hf:
                        if not hline.startswith("#CHROM"):
                            fout.write(hline)
            fout.write(line)
            continue
        fields = line.rstrip("\\n").split("\\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info = fields[7]
        # Only look up missense variants (VEP annotated missense_variant in INFO)
        if "missense_variant" in info:
            score, am_class = am_lookup(chrom, pos, ref, alt, "${alphamissense_tsv}")
            if score is not None:
                am_class = am_class.strip()
                info += f";AM_SCORE={score:.3f};AM_CLASS={am_class}"
                if score >= ${pp3_threshold}:
                    info += ";AM_PP3"
                    pp3_count += 1
                elif score <= ${bp4_threshold}:
                    info += ";AM_BP4"
                    bp4_count += 1
                else:
                    ambiguous_count += 1
            else:
                no_score_count += 1
        fields[7] = info
        fout.write("\\t".join(fields) + "\\n")

# Write statistics
with open("${sample_id}.am_stats.txt", "w") as sf:
    sf.write(f"# AlphaMissense annotation statistics\\n")
    sf.write(f"# ClinGen SVI 2024 thresholds: PP3>={${pp3_threshold}}, BP4<={${bp4_threshold}}\\n")
    sf.write(f"PP3_count\\t{pp3_count}\\n")
    sf.write(f"BP4_count\\t{bp4_count}\\n")
    sf.write(f"Ambiguous_count\\t{ambiguous_count}\\n")
    sf.write(f"No_AM_score\\t{no_score_count}\\n")
PYEOF

    # ── Compress and index output VCF ─────────────────────────────────────
    bgzip -c ${sample_id}.am_annotated.unsorted.vcf > ${sample_id}.am_annotated.vcf.gz
    tabix -p vcf ${sample_id}.am_annotated.vcf.gz

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
        alphamissense_db: "AlphaMissense_hg38.tsv.gz (2023)"
        pp3_threshold: "${pp3_threshold}"
        bp4_threshold: "${bp4_threshold}"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.am_annotated.vcf.gz \\
          ${sample_id}.am_annotated.vcf.gz.tbi \\
          ${sample_id}.am_stats.txt \\
          versions.yml
    """
}
