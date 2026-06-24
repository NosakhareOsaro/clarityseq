// ============================================================================
// Module: HAPLOGREP3_CLASSIFY
// Description: Classifies the mitochondrial haplogroup of a sample from a
//              mitochondrial VCF using HaploGrep3. Haplogroup classification
//              assigns each sample to a node in the Phylotree mitochondrial
//              phylogenetic tree (Build 17, van Oven & Kayser 2009), producing
//              a haplogroup label (e.g. H1a1, L3e2b, A2) and a quality score.
//
//              CLINICAL RELEVANCE — MUST RUN BEFORE ACMG CLASSIFICATION:
//              Mitochondrial variant pathogenicity is strongly haplogroup-
//              dependent. The ACGS Mito Guidelines 2024 §6 mandate haplogroup
//              classification before applying ACMG criteria to mito variants:
//
//              1. PRIVATE vs PHYLOGENETIC variants:
//                 Many mtDNA positions are haplogroup-defining (phylogenetic).
//                 A variant that appears "rare" in gnomAD may be common within
//                 its haplogroup. BayesACMG uses the haplogroup label to query
//                 gnomAD mtDNA haplogroup-stratified frequencies (BA1/BS1 rules).
//
//              2. HETEROPLASMY vs HOMOPLASMASY context:
//                 Haplogroup-defining variants are typically homoplasmic.
//                 A candidate pathogenic variant at a haplogroup-defining
//                 position that is heteroplasmic requires additional scrutiny
//                 (may be a somatic rather than germline de novo mutation).
//
//              3. CO-PHYLOGENETIC VARIANTS (benign modifier):
//                 Variants on the same phylogenetic branch as a known
//                 pathogenic variant are often co-inherited benign modifiers.
//                 HaploGrep3 output flags variants as "phylogenetic" vs
//                 "private" to support this distinction.
//
//              4. MITOMAP vs PHYLOTREE reconciliation:
//                 MITOMAP lists variants as pathogenic/benign. ACGS §6
//                 requires reconciling MITOMAP status with phylogenetic
//                 position — HaploGrep3 haplogroup output is the input for
//                 this reconciliation step in BayesACMG rules/mito.py.
//
// Guidelines: ACGS Mitochondrial DNA Variant Interpretation Guidelines 2024 §6;
//             MITOMAP (https://www.mitomap.org);
//             Phylotree Build 17 (van Oven & Kayser 2009 Hum Mutat PMID:18853457);
//             Weissensteiner et al. 2021 (HaploGrep3) Nucleic Acids Res PMID:33963836
// Inputs:  mito_vcf_ch — tuple(sample_id, mito.vcf.gz) from GATK4_MUTECT2_MITO
//              Mitochondrial VCF produced by Mutect2 running in
//              --mitochondria-mode (chrM only, with heteroplasmy FILTER removed).
//              Must be filtered to PASS variants only before HaploGrep3 input.
// Outputs: haplogroup  — tuple(sample_id, haplogrep.txt) [haplogroup + score]
//          phylo_vcf   — tuple(sample_id, haplogrep.phylo.vcf) [annotated VCF]
// Container: quay.io/biocontainers/haplogrep:3.2.1
//   Maintained by the Biocontainers project.
//   HaploGrep3 is a Java application; the biocontainers image bundles OpenJDK 17
//   and the HaploGrep3 JAR. Phylotree Build 17 is embedded in the JAR.
//   Source: https://quay.io/repository/biocontainers/haplogrep
// Docs: https://github.com/seppinho/haplogrep-cmd
//       https://haplogrep.i-med.ac.at
// Parameter rationale:
//   --format vcf: input is VCF format (vs FASTA/hsd legacy formats)
//   --out txt: plaintext haplogroup report (tab-delimited: SampleID, Haplogroup,
//       Quality, Range, Not_Found_Polys, Found_Polys, Remaining_Polys, AAC_In_Seqs)
//   --extend-report: emit extended per-variant annotation columns including
//       "phylogenetic" / "private" classification. Required by BayesACMG mito.py.
//   --metric kosambi: Kosambi distance metric for haplogroup quality score.
//       Alternative: kulczynski2 (similar results for high-coverage mito).
//   --chip false: not an Affymetrix chip input; do NOT set --chip for WGS VCF.
// Version note: haplogrep:3.2.1 — HaploGrep3 (2021+) replaces HaploGrep2.
//   Key improvement: native VCF input with heteroplasmy support. HaploGrep2
//   required converting VCF to HSD format (lossy — dropped heteroplasmy info).
//   HaploGrep3 retains heteroplasmy level (AF field from Mutect2) and uses it
//   when classifying borderline haplogroup assignments. Pin to 3.2.1.
// ============================================================================

nextflow.enable.dsl = 2

process HAPLOGREP3_CLASSIFY {

    tag "${sample_id}"

    label 'process_low'

    // quay.io/biocontainers/haplogrep:3.2.1 — HaploGrep3 + OpenJDK 17.
    // Phylotree Build 17 is embedded in the JAR; no external reference needed.
    // Pin to 3.2.1 — haplogroup assignments depend on bundled Phylotree version.
    container 'quay.io/biocontainers/haplogrep:3.2.1'

    publishDir "${params.outdir}/haplogrep3/${sample_id}", mode: 'copy',
        pattern: "*.{txt,vcf}"

    input:
    // Mitochondrial VCF from GATK4_MUTECT2_MITO (chrM only, PASS-filtered).
    // ACGS 2024 §6: haplogroup classification must precede ACMG mito classification.
    tuple val(sample_id), path(mito_vcf)

    output:
    // Tab-delimited haplogroup report: SampleID, Haplogroup, Quality, ...
    // Consumed by BayesACMG rules/mito.py for haplogroup-stratified frequency lookup.
    tuple val(sample_id), path("${sample_id}.haplogrep.txt"),      emit: haplogroup
    // VCF annotated with phylogenetic position (phylogenetic vs private per variant)
    // Required for BayesACMG BS1/BA1 haplogroup-stratified variant classification.
    tuple val(sample_id), path("${sample_id}.haplogrep.vcf"),      emit: phylo_vcf
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── HAPLOGREP 3.2.1 — Mitochondrial haplogroup classification ─────────
    # MUST run before ACMG mito variant classification (ACGS 2024 §6).
    # See module header for clinical rationale.
    #
    # --format vcf: direct VCF input (HaploGrep3 natively parses VCF).
    #    No HSD conversion needed (unlike HaploGrep2).
    # --extend-report: adds "phylogenetic" / "private" column per variant.
    #    BayesACMG mito.py reads this to distinguish haplogroup-defining
    #    variants from candidate pathogenic private variants.

    haplogrep classify \\
        --in     ${mito_vcf} \\
        --format vcf \\
        --out    ${sample_id}.haplogrep.txt \\
        --extend-report \\
        --metric kosambi \\
        --chip   false

    # ── Generate VCF-format annotated output ───────────────────────────────
    # Produces a VCF with HAPLOGROUP and PHYLO_STATUS INFO fields appended.
    haplogrep classify \\
        --in     ${mito_vcf} \\
        --format vcf \\
        --out    ${sample_id}.haplogrep.vcf \\
        --output-format vcf \\
        --extend-report \\
        --metric kosambi \\
        --chip   false

    # ── Validate output (fail early if haplogroup assignment failed) ───────
    if [ ! -s "${sample_id}.haplogrep.txt" ]; then
        echo "ERROR: HaploGrep3 produced empty output for ${sample_id}" >&2
        exit 1
    fi

    # ── Log the assigned haplogroup for pipeline traceability ─────────────
    echo "=== HaploGrep3 haplogroup assignment for ${sample_id} ===" >&2
    cat "${sample_id}.haplogrep.txt" >&2

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep: \$(haplogrep --version 2>&1 | grep -o '[0-9.]*' | head -1)
        phylotree_build: "17"
    END_VERSIONS
    """

    stub:
    """
    echo -e "SampleID\tHaplogroup\tQuality\n${sample_id}\tH1a1\t0.99" > ${sample_id}.haplogrep.txt
    touch ${sample_id}.haplogrep.vcf versions.yml
    """
}
