// ============================================================================
// Module: VEP_ANNOTATE
// Description: Annotates ensemble-merged VCF with functional consequences,
//              population frequencies, pathogenicity predictions, and clinical
//              significance using Ensembl Variant Effect Predictor (VEP) v111.
//              VEP is the primary annotation engine in GenomeForge, providing
//              transcript-level consequences using MANE Select as the canonical
//              clinical transcript. Multiple plugins are loaded to support ACMG
//              variant classification criteria (PP3, BP4, PM2, BA1, BS1, etc.).
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §5 (annotation);
//             ACMG/AMP 2015 variant interpretation guidelines PMID:25741868;
//             ClinGen SVI recommendations 2024;
//             Morales et al. 2022 Nature PMID:35356062 (MANE Select)
// Inputs:  vcf_ch      — tuple(sample_id, ensemble.vcf.gz) from ENSEMBLE_MERGE
//          vep_cache   — path to VEP offline cache v111 (GRCh38)
//          plugins_dir — path to VEP plugins directory
//          reference   — GRCh38 FASTA (for HGMD lookups)
//          alphamissense_tsv — AlphaMissense scores TSV (tabix-indexed)
// Outputs: annotated_vcf_ch — tuple(sample_id, annotated.vcf.gz)
//          vep_stats_ch     — VEP annotation statistics HTML
// Container: ensemblorg/ensembl-vep:release_111
// Docs: https://www.ensembl.org/info/docs/tools/vep/index.html
//       https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
// Parameter rationale:
//   --pick_order: transcript selection priority order (see below)
//   --format vcf / --vcf: VCF input and output
//   --everything: enables most annotation fields (see inline notes for overrides)
//   --offline: use local cache; no network required
//   --cache_version 111: must match the cache directory version
//   --assembly GRCh38: must match reference genome build
// Version note: VEP v111 (2024) — chosen because:
//   1. v111 includes MANE Select v1.3 annotations (January 2024 release)
//   2. v111 integrates dbNSFP v4.7 (latest as of Jan 2024)
//   3. v111 improves LOFTEE annotations for loss-of-function variant calls
//   4. v111 is compatible with gnomAD v4.1 AF annotations
//   Do NOT use v110 or earlier — MANE Select annotations differ significantly.
// ============================================================================
//
// TRANSCRIPT PICK ORDER (--pick_order):
// =======================================
// VEP can annotate variants against multiple transcripts. --pick selects ONE
// canonical annotation per variant using a priority order. GenomeForge uses:
//
//   mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,rank,length
//
//   1. mane_select: MANE Select transcript — the single transcript per gene
//      agreed upon by NCBI and Ensembl as the default clinical transcript.
//      Per ACGS 2024 v1.2 and Morales et al. 2022 (PMID:35356062), MANE
//      Select must be used as the primary transcript for clinical reporting.
//      In VEP v111, ~99% of genes with clinical significance have MANE Select.
//
//   2. mane_plus_clinical: MANE Plus Clinical — disease-specific transcripts
//      added when a clinically relevant isoform differs from MANE Select.
//      Used for BRCA1 exon skipping isoforms, TTN cardiac transcripts, etc.
//
//   3. canonical: Ensembl canonical transcript (fallback if no MANE Select)
//      Used for genes not yet covered by MANE (rare, decreasing over time)
//
//   4-8: appris, tsl, biotype, rank, length — further tiebreakers for
//      genes with multiple equally-canonical transcripts
//
// PLUGINS AND ACMG EVIDENCE CODES:
// ==================================
//
//   AlphaMissense (primary: PP3/BP4):
//     ClinGen SVI 2024 approved AlphaMissense for PP3 and BP4 evidence.
//     Thresholds: score ≥0.564 → PP3 (likely pathogenic missense);
//                 score ≤0.340 → BP4 (likely benign missense)
//     0.340 < score < 0.564 → ambiguous (no code assigned)
//     Cheng et al. 2023 Science PMID:37703350
//     Note: GenomeForge also runs modules/alphamissense/main.nf independently
//     for tabix lookup; this plugin provides the same scores via VEP output.
//
//   SpliceAI (splicing PP3/BP4/BP7):
//     Delta scores for acceptor/donor gain/loss per transcript.
//     Thresholds (Jaganathan et al. 2019 Cell PMID:30661751):
//       DS_AG/DS_AL/DS_DG/DS_DL ≥0.5 → strong splicing impact → PP3/BP4
//       DS_* 0.1-0.5 → moderate impact → relevant for BP7 adjustment
//       DS_* <0.1 → no predicted splice impact → supports BP7
//
//   Pangolin (tissue-specific splicing, supplements SpliceAI):
//     Predicts tissue-specific splice-site usage changes.
//     Where SpliceAI gives a combined score, Pangolin provides per-tissue
//     resolution (e.g., distinguishes brain-specific from ubiquitous splicing).
//     Used as supplementary evidence alongside SpliceAI for BP7 coding variants.
//     Zeng et al. 2022 Genome Biology PMID:35240905
//
//   dbNSFP_4.7 (bulk in-silico predictions):
//     Provides 50+ in-silico prediction scores in one plugin. Key scores used:
//     CADD_phred (PP3: >25; BP4: <15), REVEL (PP3: >0.75; BP4: <0.15),
//     SIFT (BP4 supporting: tolerated), PolyPhen-2 (PP3 supporting: probably_damaging)
//     Note: individual scores have lower evidential weight than AlphaMissense;
//     use as supporting evidence, not standalone PP3/BP4.
//
//   gnomAD_v4.1 (PM2/BA1/BS1):
//     Population allele frequency filtering — critical for variant classification.
//     BA1 (stand-alone benign): AF >0.05 in any gnomAD v4.1 population
//     BS1 (strong benign): 0.005 < AF ≤ 0.05 in relevant population
//     PM2 (moderate pathogenic): absent from gnomAD v4.1 (or extremely rare)
//     IMPORTANT: Use ancestry-matched population AF where possible:
//       (gnomAD AFR/AMR/EAS/EUR/SAS subset from SOMALIER_ANCESTRY output)
//     Reference: gnomAD v4.1 released April 2024; 807,162 individuals.
//     v4.0 had AN bug — do NOT use v4.0 annotations.
//
//   ClinVar (PP5/BP6 — use with caution):
//     PP5 (supporting pathogenic): ClinVar P/LP with ≥2-star review status
//     BP6 (supporting benign): ClinVar B/LB with ≥2-star review status
//     CAUTION: ClinVar entries below 2-star have high inter-lab disagreement.
//     Never use PP5/BP6 with conflict status or 1-star entries.
//     Always verify ClinVar entries independently before clinical reporting.
//     Landrum et al. 2016 Nucleic Acids Res PMID:26582918
//

nextflow.enable.dsl = 2

process VEP_ANNOTATE {

    tag "${sample_id}"

    label 'process_high'

    // ensemblorg/ensembl-vep:release_111 — pin to release_111 exactly.
    // VEP version determines plugin API compatibility and cache format.
    // Cache must be v111 (GRCh38) to match this container.
    container 'ensemblorg/ensembl-vep:release_111'

    publishDir "${params.outdir}/annotation/${sample_id}", mode: 'copy',
        pattern: "*.{vcf.gz,vcf.gz.tbi,html}"

    input:
    // Ensemble-merged VCF (GATK4 ∩ DeepVariant) from ENSEMBLE_MERGE
    tuple val(sample_id), path(vcf), path(tbi)
    // VEP offline cache: must be v111 GRCh38 version
    // Download: perl /opt/vep/src/ensembl-vep/INSTALL.pl -a c -s homo_sapiens -y GRCh38 --CACHE_VERSION 111
    path vep_cache_dir
    // VEP plugins directory: AlphaMissense, SpliceAI, Pangolin, dbNSFP, gnomAD, ClinVar
    path plugins_dir
    // GRCh38 FASTA (required for HGMD and some plugin lookups)
    path reference
    path reference_fai
    // AlphaMissense tabix-indexed TSV (from params.alphamissense_tsv)
    path alphamissense_tsv
    path alphamissense_tbi
    // SpliceAI pre-computed scores (download from Illumina Basespace)
    path spliceai_snv
    path spliceai_snv_tbi
    path spliceai_indel
    path spliceai_indel_tbi
    // dbNSFP v4.7 database (from https://sites.google.com/site/jpopgen/dbNSFP)
    path dbnsfp
    path dbnsfp_tbi
    // gnomAD v4.1 VCF for population AF annotation
    path gnomad_vcf
    path gnomad_tbi
    // ClinVar VCF (latest clinvar.vcf.gz from NCBI FTP)
    path clinvar_vcf
    path clinvar_tbi

    output:
    // Fully annotated VCF with all plugin fields in INFO column
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz.tbi"), emit: tbi
    // VEP annotation statistics HTML
    path "${sample_id}.vep_stats.html",                        emit: stats
    // TSV summary for MultiQC plugin
    path "${sample_id}.vep_summary.txt",                       emit: summary
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cache_version = params.vep_cache_version ?: 111
    """
    # ── VEP v111 annotation ───────────────────────────────────────────────
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.vcf.gz \\
        --stats_file ${sample_id}.vep_stats.html \\
        \\
        --vcf \\
        # Output VCF format (vs default VEP tab-delimited output).
        # Required for downstream ACMG classifier and ClinVar beacon ingest.
        \\
        --compress_output bgzip \\
        # bgzip-compress output VCF for space efficiency and tabix indexing
        \\
        --offline \\
        # Use local cache only — no network requests.
        # Cache must be at vep_cache_dir/homo_sapiens/${cache_version}_GRCh38/
        \\
        --cache \\
        --dir_cache ${vep_cache_dir} \\
        --cache_version ${cache_version} \\
        # VEP cache version must match container version (111)
        \\
        --assembly GRCh38 \\
        # Genome assembly: must be GRCh38 (hg38).
        # T2T arm uses separate annotation after liftover (pipelines/wgs_t2t.nf).
        \\
        --species homo_sapiens \\
        \\
        --fasta ${reference} \\
        # Reference FASTA for lookup in HGMD and some plugin operations
        \\
        --pick \\
        # Select ONE consequence per variant (the most severe / most canonical)
        # Prevents per-transcript information explosion in clinical VCF.
        \\
        --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,rank,length \\
        # TRANSCRIPT PRIORITY (most → least preferred):
        # 1. mane_select: MANE Select transcript (ACGS 2024 §5, Morales 2022)
        # 2. mane_plus_clinical: disease-specific clinical isoforms
        # 3. canonical: Ensembl canonical (fallback if no MANE)
        # 4-8: appris/tsl/biotype/rank/length tiebreakers
        \\
        --everything \\
        # Enable all annotation fields: HGVS, Existing_variation, AF, CLIN_SIG,
        # DOMAINS, EXON, INTRON, HGVSc, HGVSp, Codons, Amino_acids, etc.
        \\
        --hgvs \\
        # HGVS notation (c. and p.) — required for clinical report variant naming
        \\
        --hgvsg \\
        # HGVS genomic notation (g.) — for variant databases (ClinVar, LOVD)
        \\
        --symbol \\
        # Gene symbol (HGNC) in output — required for clinical report
        \\
        --biotype \\
        # Transcript biotype (protein_coding, retained_intron, etc.)
        \\
        --canonical \\
        # Flag canonical transcripts (those not selected by --pick)
        \\
        --mane \\
        # Annotate MANE Select and MANE Plus Clinical transcripts
        \\
        --numbers \\
        # Exon/intron numbers (e.g., exon 3/24) for clinical reporting
        \\
        --sift b \\
        # SIFT score + prediction (T=tolerated / D=deleterious)
        # Supporting evidence for BP4 (tolerated) or PP3 (deleterious)
        \\
        --polyphen b \\
        # PolyPhen-2 HVAR score + prediction for missense variants
        # Supporting evidence for BP4 (benign) or PP3 (possibly/probably damaging)
        \\
        --fork ${task.cpus} \\
        # Parallelise VEP across input VCF chunks
        \\
        --plugin AlphaMissense,file=${alphamissense_tsv} \\
        # AlphaMissense: PRIMARY evidence for PP3/BP4 (ClinGen SVI 2024 approved).
        # Thresholds: ≥0.564 → PP3; ≤0.340 → BP4; ambiguous otherwise.
        # Cheng et al. 2023 Science PMID:37703350
        \\
        --plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel},cutoff=0.1 \\
        # SpliceAI: splicing PP3/BP4/BP7 evidence.
        # DS_AG/AL/DG/DL ≥0.5 → strong splicing impact → PP3 or BP4
        # DS_* 0.1-0.5 → moderate impact → consider for BP7
        # Jaganathan et al. 2019 Cell PMID:30661751
        \\
        --plugin Pangolin,${plugins_dir}/Pangolin/pangolin.py \\
        # Pangolin: tissue-specific splice prediction (supplement to SpliceAI).
        # Useful for BP7: if SpliceAI is borderline and Pangolin shows no impact
        # in relevant tissues, BP7 (silent/intronic with no splice impact) is stronger.
        # Zeng et al. 2022 Genome Biology PMID:35240905
        \\
        --plugin dbNSFP,${dbnsfp},CADD_phred,REVEL_score,SIFT4G_score,MutationAssessor_score,MetaSVM_score,MetaLR_score,M-CAP_score,VEST4_score,MutPred_score,PrimateAI_score \\
        # dbNSFP v4.7: bulk in-silico prediction scores in one lookup.
        # Key scores and PP3/BP4 thresholds:
        #   CADD_phred: >25=PP3, <15=BP4 (Kircher et al. 2014 Nature Genetics)
        #   REVEL: >0.75=PP3, <0.15=BP4 (Ioannidis et al. 2016 Am J Hum Genet)
        # These are SUPPORTING evidence — lower weight than AlphaMissense.
        \\
        --plugin gnomADc,${gnomad_vcf} \\
        # gnomAD v4.1 population allele frequencies for PM2/BA1/BS1:
        #   BA1 (stand-alone benign): AF >0.05 in any population
        #   BS1 (strong benign): 0.005 < AF ≤ 0.05 in matched population
        #   PM2 (moderate pathogenic): absent or extremely rare in gnomAD
        # IMPORTANT: Use gnomAD v4.1 (April 2024). v4.0 had AN bug.
        \\
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
        # ClinVar: PP5/BP6 supporting evidence — USE WITH CAUTION.
        # PP5: ClinVar P/LP with ≥2-star review (criteria_provided_multiple_submitters)
        # BP6: ClinVar B/LB with ≥2-star review
        # NEVER use single-star or conflict-status ClinVar entries for PP5/BP6.
        # Always verify ClinVar assertions independently before clinical reporting.
        \\
        --no_stats \\
        # Generate separate stats file (--stats_file above) rather than inline

    # ── Index the output VCF ──────────────────────────────────────────────
    tabix -p vcf ${sample_id}.vep.vcf.gz

    # ── Generate TSV summary for MultiQC ──────────────────────────────────
    bcftools stats ${sample_id}.vep.vcf.gz | grep "^SN" > ${sample_id}.vep_summary.txt

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: \$(vep --version 2>&1 | sed 's/ensembl-vep : //')
        htslib: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.vep.vcf.gz \\
          ${sample_id}.vep.vcf.gz.tbi \\
          ${sample_id}.vep_stats.html \\
          ${sample_id}.vep_summary.txt \\
          versions.yml
    """
}
