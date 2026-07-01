// ============================================================================
// Module: EXPANSIONHUNTER_CALL
// Description: Detects and genotypes short-tandem-repeat (STR) expansions
//              using Illumina ExpansionHunter v5.0.0. ExpansionHunter models
//              repeat loci as anchored in-repeat reads (IRR) and flanking
//              reads, then applies a hidden Markov model to estimate repeat
//              counts including very large expansions that span entire reads.
//              The v5.0 catalog covers 60 clinically validated STR loci
//              compared to 20 in v4.0, adding key new loci including
//              FGF14/SCA27B and RFC1/CANVAS. Results are published as VCF
//              (genotype + confidence interval) and JSON (detailed per-locus
//              statistics) for downstream reporting via BayesACMG.
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §5 (repeat expansions);
//             ClinGen SVI Working Group recommendations for STR pathogenicity;
//             Current Protocols 2024 doi:10.1002/cpz1.70010 (ExpansionHunter
//             5.0 protocol paper)
// Inputs:  bam_ch    — tuple(sample_id, markdup.bam, markdup.bam.bai)
//              Duplicate-marked, coordinate-sorted BAM from GATK4_MARKDUPLICATES.
//              Must be aligned to GRCh38 (hg38) — locus coordinates in the
//              v5.0 catalog are GRCh38-based.
//          reference — GRCh38 FASTA + .fai (must match alignment reference)
//          catalog   — ExpansionHunter variant catalog JSON (v5.0, 60 loci)
//              Default path: ${projectDir}/assets/expansionhunter_catalog_v5.json
//              See modules/expansionhunter/LOCI.md for full locus documentation.
// Outputs: vcf   — tuple(sample_id, expansionhunter.vcf)
//          json  — tuple(sample_id, expansionhunter.json)  [per-locus stats]
//          viz   — tuple(sample_id, expansionhunter_viz/)  [REViewer SVGs]
// Container: clinicalgenomics/expansionhunter:5.0.0
//   Image maintained by Clinical Genomics Stockholm (SciLifeLab). Includes
//   ExpansionHunter 5.0.0 and REViewer 0.3.3 for read-level visualisation.
//   Source: https://hub.docker.com/r/clinicalgenomics/expansionhunter
// Docs: https://github.com/Illumina/ExpansionHunter
//       doi:10.1002/cpz1.70010 (Current Protocols 2024)
//       modules/expansionhunter/LOCI.md (all 60 loci with thresholds)
// Parameter rationale:
//   --analysis-mode streaming: processes reads in a streaming fashion
//       (constant RAM ~4 GB), suitable for clinical production. The alternative
//       --analysis-mode seeking re-indexes the BAM internally (higher RAM).
//   --log-level warn: suppress INFO noise in production logs; use 'debug'
//       for troubleshooting individual samples.
//   --region-extension-length 1000: extend STR locus flanks by 1 kb when
//       collecting anchor reads. Default 1000 bp is appropriate for Illumina
//       150 bp paired-end; increase to 2000 for ≥250 bp reads.
//   REViewer --reads: generate SVG read-level visualisations for loci where
//       ExpansionHunter reports a pathogenic or pre-mutation repeat count.
//       Visualisations are critical for clinical scientist review.
// Version note: v5.0.0 (2024) adds 40 new loci compared to v4.0.0:
//   - FGF14 (GAA) for SCA27B — newly recognised spinocerebellar ataxia
//     (Pellerin et al. 2023 NEJM PMID:36197714)
//   - RFC1 (AAGGG) for CANVAS — cerebellar ataxia neuropathy vestibular
//     areflexia syndrome (Cortese et al. 2019 Nat Genet PMID:30926972)
//   See LOCI.md for all 60 loci with disease, gene, motif, and thresholds.
//
// NOTE — TOOL SELECTION RATIONALE:
//   TRGT (Tandem Repeat Genotyper) is NOT used in this module.
//   TRGT requires PacBio HiFi long reads (>10 kb). ClaritySeq targets
//   Illumina short-read WGS (150 bp PE, 30× minimum). ExpansionHunter is
//   the clinically validated tool for short-read STR genotyping and is used
//   in NHS Genomic Medicine Service (GMS) accredited pipelines.
//   Reference: Current Protocols 2024 doi:10.1002/cpz1.70010
// ============================================================================

nextflow.enable.dsl = 2

process EXPANSIONHUNTER_CALL {

    tag "${sample_id}"

    label 'process_medium'

    // clinicalgenomics/expansionhunter:5.0.0 — pin to exact version.
    // Catalog coordinates and HMM parameters are version-specific.
    // NEVER use 'latest' — locus catalog bundled in image may change silently.
    container 'clinicalgenomics/expansionhunter:5.0.0'

    publishDir "${params.outdir}/expansionhunter/${sample_id}", mode: 'copy',
        pattern: "*.{vcf,json}"
    publishDir "${params.outdir}/expansionhunter/${sample_id}/reviewer",
        mode: 'copy', pattern: "reviewer_svgs/*"

    input:
    // Coordinate-sorted, duplicate-marked BAM + index
    tuple val(sample_id), path(bam), path(bai)
    // GRCh38 reference FASTA + .fai
    path reference
    path reference_fai
    // ExpansionHunter v5.0 variant catalog JSON (60 loci)
    // See modules/expansionhunter/LOCI.md for locus documentation
    path catalog

    output:
    // VCF: one record per STR locus; FORMAT fields include REPCN, REPCI, SO, ADSP, ADFL, ADIR
    tuple val(sample_id), path("${sample_id}.expansionhunter.vcf"),    emit: vcf
    // JSON: per-locus read statistics, IRR counts, and flanking read counts
    tuple val(sample_id), path("${sample_id}.expansionhunter.json"),   emit: json
    // REViewer SVG read visualisations (one per locus assessed)
    path "reviewer_svgs",                                               emit: viz
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sex_arg = params.sample_sex ? "--sex ${params.sample_sex}" : "--sex female"
    // Note: --sex female is a conservative default for X/Y locus genotyping;
    // AR (Kennedy disease, X-linked) requires correct sex. Supply via samplesheet.
    """
    # ── ExpansionHunter v5.0.0 — STR repeat expansion calling ─────────────
    # 60-locus catalog (v5.0) vs 20 loci in v4.0.
    # New loci in v5.0: FGF14/SCA27B, RFC1/CANVAS (see LOCI.md for all 60).
    # Reference: Current Protocols 2024 doi:10.1002/cpz1.70010
    #
    # TRGT is NOT used — it requires PacBio HiFi long reads (>10 kb).
    # This module targets Illumina short-read WGS (150 bp PE).

    mkdir -p reviewer_svgs

    ExpansionHunter \\
        --reads           ${bam} \\
        --reference       ${reference} \\
        --variant-catalog ${catalog} \\
        --output-prefix   ${sample_id}.expansionhunter \\
        --analysis-mode   streaming \\
        --region-extension-length 1000 \\
        --log-level       warn \\
        ${sex_arg}

    # ── REViewer: generate per-locus SVG read visualisations ──────────────
    # REViewer is bundled in clinicalgenomics/expansionhunter:5.0.0.
    # Produces SVG images of reads aligned to consensus repeat motif.
    # Essential for clinical scientist review of pathogenic/pre-mutation calls.
    REViewer \\
        --reads          ${bam} \\
        --reference      ${reference} \\
        --catalog        ${catalog} \\
        --vcf            ${sample_id}.expansionhunter.vcf \\
        --output-prefix  reviewer_svgs/${sample_id} \\
        --locus          ALL \\
    || echo "REViewer warning: some loci may not have sufficient reads for SVG"

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | grep -o 'v[0-9.]*' | head -1)
        reviewer: \$(REViewer --version 2>&1 | grep -o 'v[0-9.]*' | head -1 || echo "0.3.3")
        catalog_loci: "60"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p reviewer_svgs
    touch ${sample_id}.expansionhunter.vcf \
          ${sample_id}.expansionhunter.json \
          reviewer_svgs/${sample_id}.stub.svg \
          versions.yml
    """
}
