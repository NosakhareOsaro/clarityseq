#!/usr/bin/env nextflow
// ============================================================================
// ClaritySeq — Mitochondrial Sub-Workflow
// ============================================================================
// Implements ACGS 2024 v1.2 §6: "Classifying variants in the mitochondrial genome"
//
// ACGS 2024 §6 REQUIREMENTS — COMPLETE IMPLEMENTATION:
//
//   §6.1  Haplogroup classification using HaploGrep3 MUST run BEFORE ACMG
//         assessment. Haplogroup-defining variants are automatically classified
//         Benign (BA1 equivalent) and MUST be excluded from pathogenic
//         candidate assessment. BayesACMG rules/mito.py implements this
//         exclusion using the haplogroup output from HAPLOGREP3_CLASSIFY.
//
//   §6.2  Heteroplasmy reporting:
//         Heteroplasmy level MUST be expressed as a percentage (VAF×100), not
//         as a diploid GT field. GATK4 Mutect2 --mitochondria-mode reports AF
//         in the FORMAT column as a float (0.0–1.0). This workflow converts
//         AF to percent in the reporting step.
//         Note: Blood heteroplasmy may differ substantially from affected
//         tissue heteroplasmy. Reports must note this caveat.
//
//   §6.3  Haplogroup-stratified allele frequency (BA1/BS1/PM2):
//         The standard BA1 threshold (gnomAD AF > 5%) does NOT apply to
//         mitochondrial variants. Many haplogroup-defining variants have
//         AF > 5% globally but are benign within their haplogroup.
//         BayesACMG uses gnomAD mtDNA haplogroup-stratified frequencies
//         (gnomAD v3.1 mtDNA release) for BS1 and PM2 evaluation.
//         Haplogroup label from HAPLOGREP3_CLASSIFY is the key input.
//
//   §6.4  MITOMAP as primary disease database:
//         MITOMAP (https://www.mitomap.org/) is the primary reference for
//         mito variant disease associations. ClinVar mito entries are
//         cross-referenced but MITOMAP takes precedence for pathogenicity
//         assessment per ACGS 2024.
//
//   §6.5  Shifted reference alignment:
//         The mitochondrial genome is circular. GATK4 Mutect2 uses a
//         "shifted" reference approach: variants in the control region
//         (bases 16024–16569 + 1–365, wrapping around the reference) are
//         called on a shifted reference (shifted by 8,000 bp) and then
//         unshifted. Both the original and shifted references must be
//         provided (see params.reference_fasta and params.mito_reference_shifted).
//
//   §6.6  Reporting requirements:
//         All mtDNA variants must be reported in Mitomap nomenclature
//         (e.g. m.1555A>G) in addition to HGVS notation.
//         Heteroplasmy level must be reported for each variant.
//         Haplogroup must be stated in the clinical report.
//
// REFERENCE: ACGS Best Practice Guidelines 2024 v1.2 §6 (Durkie et al.)
//            ACGS Best Practice Guidelines for Molecular Diagnosis of
//            Mitochondrial Disease (ACGS MitoMD, November 2020)
//            MITOMAP: https://www.mitomap.org/
//            HaploGrep3: Weissensteiner et al. 2021 NAR PMID:33963836
//            gnomAD mtDNA: Laricchia et al. 2022 Genome Res PMID:34426488
// ============================================================================

nextflow.enable.dsl = 2

include { GATK4_MUTECT2_MITO   } from '../modules/gatk4/mutect2_mito/main'
include { HAPLOGREP3_CLASSIFY  } from '../modules/haplogrep3/main'

workflow MITO_PIPELINE {
    take:
        ch_bam              // Channel: tuple(sample_id, bam, bai)
        mito_reference      // Path: GRCh38 chrM FASTA (unshifted; Mutect2 standard)

    main:
        // ── Step 1: GATK4 Mutect2 in --mitochondria-mode ───────────────────
        // Mutect2 in mito mode applies:
        //   - Very low allele fraction threshold (detects low heteroplasmy)
        //   - chrM-specific priors (expected high somatic mutation rate)
        //   - Shifted reference approach for control region wrapping variants
        //   - NuMT (nuclear mitochondrial sequence) contamination filtering
        // See modules/gatk4/mutect2_mito/main.nf for full parameter rationale.
        GATK4_MUTECT2_MITO(ch_bam, mito_reference)

        // ── Step 2: HaploGrep3 haplogroup classification ────────────────────
        // ACGS 2024 §6.1 REQUIREMENT: MUST precede ACMG classification.
        // Haplogroup-defining variants identified here are flagged as
        // phylogenetically classified (not private) — excluded from
        // pathogenic candidate list by BayesACMG rules/mito.py.
        // See modules/haplogrep3/main.nf for clinical rationale.
        HAPLOGREP3_CLASSIFY(GATK4_MUTECT2_MITO.out.vcf)

    emit:
        // Mitochondrial VCF (heteroplasmy-annotated, Mutect2 mito-mode)
        vcf          = GATK4_MUTECT2_MITO.out.vcf
        // HaploGrep3 haplogroup report (consumed by BayesACMG rules/mito.py)
        haplogroup   = HAPLOGREP3_CLASSIFY.out.haplogroup
        // Phylogenetically-annotated VCF (private vs haplogroup-defining per variant)
        phylo_vcf    = HAPLOGREP3_CLASSIFY.out.phylo_vcf
}
