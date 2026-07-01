#!/usr/bin/env nextflow
// ============================================================================
// ClaritySeq — Short Tandem Repeat Expansion Sub-Workflow
// ============================================================================
// Genotypes repeat expansions using ExpansionHunter v5.0 (60-locus catalog).
//
// CATALOG: ExpansionHunter v5.0 — 60 loci (see modules/expansionhunter/LOCI.md)
//   v5.0 vs v4.0: 60 loci (v5.0) vs 20 loci (v4.0).
//   Key NEW loci in v5.0:
//     - FGF14 (GAA): SCA27B — Pellerin et al. 2023 NEJM PMID:36197714
//     - RFC1 (AAGGG): CANVAS — Cortese et al. 2019 Nat Genet PMID:30926972
//   Reference: Dolzhenko et al. 2024 Current Protocols doi:10.1002/cpz1.70010
//
// TOOL SELECTION — WHY EXPANSIONHUNTER AND NOT TRGT:
//   TRGT (Tandem Repeat Genotyper by PacBio) is NOT used because:
//     - TRGT requires PacBio HiFi long reads (>10 kb continuous reads)
//     - ClaritySeq targets Illumina short-read WGS (150 bp PE)
//     - ExpansionHunter is validated for short-read STR calling and is used
//       in NHS GMS accredited pipelines for clinical STR detection
//   Reference: Dolzhenko et al. 2024 Current Protocols doi:10.1002/cpz1.70010
//
// THRESHOLDS: STRipy database 2024 release (https://stripy.org/)
//   Normal/pathogenic thresholds per locus documented in:
//   modules/expansionhunter/LOCI.md
//
// GUIDELINES: ACGS Best Practice Guidelines v1.2 2024 §5 (repeat expansions)
//             ClinGen SVI Working Group STR pathogenicity recommendations
//
// OUTPUTS:
//   VCF: one record per locus; FORMAT REPCN (repeat count), REPCI (confidence
//        interval), SO (spanning origin), ADSP/ADFL/ADIR (allele depth fields)
//   JSON: per-locus read-level statistics for clinical scientist review
//   REViewer SVGs: read-level visualisations for pathogenic/pre-mutation calls
//
// ============================================================================

nextflow.enable.dsl = 2

include { EXPANSIONHUNTER_CALL } from '../modules/expansionhunter/main'

workflow EXPANSION_PIPELINE {
    take:
        ch_bam             // Channel: tuple(sample_id, bam, bai)
        reference          // Path: GRCh38 FASTA + .fai (must match alignment reference)
        reference_fai      // Path: GRCh38 FASTA .fai index
        expansion_catalog  // Path: ExpansionHunter v5.0 variant catalog JSON (60 loci)
                           // Default: ${projectDir}/assets/expansionhunter_catalog_v5.json
                           // See modules/expansionhunter/LOCI.md for all 60 loci

    main:
        // ExpansionHunter v5.0: streaming mode, 60-locus catalog
        // Container: clinicalgenomics/expansionhunter:5.0.0
        // See modules/expansionhunter/main.nf for full parameter documentation
        // Note: TRGT is NOT used — requires PacBio HiFi long reads
        EXPANSIONHUNTER_CALL(
            ch_bam,
            reference,
            reference_fai,
            expansion_catalog
        )

    emit:
        // Per-locus VCF with REPCN, REPCI, SO, ADSP, ADFL, ADIR fields
        vcf           = EXPANSIONHUNTER_CALL.out.vcf
        // Per-locus JSON with detailed read-level statistics
        json          = EXPANSIONHUNTER_CALL.out.json
        // REViewer SVG visualisations (one per locus — for clinical review)
        reviewer_svgs = EXPANSIONHUNTER_CALL.out.viz
}
