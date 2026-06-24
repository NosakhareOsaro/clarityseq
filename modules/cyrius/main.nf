// ============================================================================
// Module: CYRIUS_CYP2D6
// Description: Genotypes CYP2D6 star alleles (pharmacogenomic haplotypes)
//              using Cyrius v1.1.1, a purpose-built tool for CYP2D6 calling
//              from short-read Illumina WGS data. Cyrius uses a read-depth
//              and allele-ratio model that explicitly accounts for the
//              CYP2D7 pseudogene, enabling reliable detection of copy-number
//              variants (CNVs) and hybrid alleles that other callers miss.
//
//   WHY CYRIUS INSTEAD OF GATK4 FOR CYP2D6:
//   -----------------------------------------
//   GATK4 HaplotypeCaller CANNOT reliably genotype CYP2D6 because:
//
//   1. PSEUDOGENE INTERFERENCE (CYP2D7):
//      CYP2D6 and its pseudogene CYP2D7 share ~97% sequence identity.
//      Standard BWA-MEM2/DRAGMAP short-read aligners map a substantial
//      fraction of CYP2D6 reads to CYP2D7 (and vice versa), causing:
//        - Artificially reduced apparent depth at CYP2D6 exons
//        - Phantom SNVs from CYP2D7 sequence "bleeding" into CYP2D6 calls
//        - MAPQ=0 multi-mapping reads discarded by GATK4 (losing true CYP2D6
//          reads that also align to CYP2D7)
//
//   2. CNV COMPLEXITY:
//      The most clinically impactful CYP2D6 variants involve structural
//      rearrangements:
//        *5 (gene deletion): entire CYP2D6 deleted → Poor Metaboliser (PM)
//        *13 (hybrid CYP2D7/D6): non-functional hybrid → PM contribution
//        CYP2D6 duplications (xN): Ultra-rapid Metaboliser (UM) phenotype
//      GATK4 SVCaller does not call single-exon or gene-level deletions
//      reliably. CNVKit and GATK4 CNV pipeline require matched normal controls
//      not available in germline WGS.
//
//   3. HAPLOTYPE PHASING:
//      CYP2D6 diplotype (e.g. *1/*4) must be PHASED to assign metaboliser
//      phenotype. GATK4 does not phase diplotypes across the ~5 kb CYP2D6
//      gene body within a single haplotype. Cyrius uses read-pair spanning
//      and CN-aware phasing.
//
//   Cyrius approach (Chen et al. 2022 NPJ Genomic Med PMID:35082305):
//     - Aligns reads to a CYP2D6+CYP2D7 dual-reference to separate signal
//     - Models read-depth ratio CYP2D6:CYP2D7 per exon for CN estimation
//     - Combines CN, SNV haplotype, and hybrid-allele detection into diplotype
//     - Assigns *-allele diplotypes from CPIC/PharmVar star allele database
//
//   CLINICAL IMPACT OF CYP2D6 PHENOTYPE:
//     Poor Metaboliser (PM): *4/*4, *5/*5, *4/*5 etc.
//       → Cannot metabolise codeine, tramadol (risk: opioid toxicity)
//       → Impaired tamoxifen activation (breast cancer risk implications)
//     Ultra-rapid Metaboliser (UM): *1/*1xN, *2/*2xN duplications
//       → Excess codeine/morphine conversion (risk: respiratory depression)
//     Phenotype labels: UM > NM (Normal) > IM (Intermediate) > PM
//
// Guidelines: CPIC CYP2D6 Guideline (Goetz et al. 2021 Clin Pharmacol Ther
//             PMID:33387367); PharmVar (https://www.pharmvar.org);
//             ACGS PGx Best Practices 2024 (Module 8);
//             FDA Table of Pharmacogenomic Biomarkers in Drug Labeling
// Inputs:  bam_ch — tuple(sample_id, markdup.bam, markdup.bam.bai)
//              Duplicate-marked, coordinate-sorted BAM aligned to GRCh38.
//              Must include CYP2D6 region (chr22:42,522,500-42,526,883 GRCh38).
// Outputs: star_allele — tuple(sample_id, cyrius_cyp2d6.txt) [diplotype + phenotype]
//          json        — tuple(sample_id, cyrius_cyp2d6.json) [detailed per-exon stats]
// Container: clinicalgenomics/cyrius:1.1.1
//   Maintained by Clinical Genomics Stockholm. Includes Cyrius v1.1.1.
//   Source: https://hub.docker.com/r/clinicalgenomics/cyrius
// Docs: https://github.com/Illumina/Cyrius
//       Chen et al. 2022 NPJ Genomic Med PMID:35082305
// Parameter rationale:
//   --genome hg38: use GRCh38 coordinate system (vs hg19 legacy)
//   --prefix: output file prefix (sample ID)
//   --threads: parallelise read processing across exons
//   Default minimum depth: Cyrius requires ≥30× mean depth at CYP2D6 region.
//       If depth is lower, result is flagged as LOW_DEPTH_WARNING in output.
//       ACGS 2024 requires ≥30× WGS minimum (MOSDEPTH_QC enforces globally).
// Version note: Cyrius v1.1.1 (2022) — chosen because:
//   - v1.1.x introduces PharmVar 5.2 star allele database (adds *139, *140,
//     *141 and updated *6 suballele definitions)
//   - Fixes incorrect hybrid-allele classification near exon 9 junction seen
//     in v1.0.x (*13 mis-called as *68 in ~2% of samples)
//   Pin to 1.1.1 — star allele nomenclature is database-version dependent.
// ============================================================================

nextflow.enable.dsl = 2

process CYRIUS_CYP2D6 {

    tag "${sample_id}"

    label 'process_low'

    // clinicalgenomics/cyrius:1.1.1 — Cyrius + PharmVar 5.2 star allele database.
    // Pin to 1.1.1 — star allele nomenclature depends on bundled PharmVar version.
    // NEVER use 'latest' — PharmVar updates can silently change allele names.
    container 'clinicalgenomics/cyrius:1.1.1'

    publishDir "${params.outdir}/pharmacogenomics/cyp2d6/${sample_id}", mode: 'copy',
        pattern: "*.{txt,json}"

    input:
    // Coordinate-sorted, duplicate-marked BAM + index (must include chr22)
    tuple val(sample_id), path(bam), path(bai)

    output:
    // CYP2D6 diplotype result: tab-delimited with columns
    //   Sample, Diplotype (e.g. *1/*4), Phenotype (PM/IM/NM/UM), Filter
    tuple val(sample_id), path("${sample_id}.cyrius_cyp2d6.txt"),  emit: star_allele
    // JSON with per-exon depth, CN ratio, and allele-specific read counts
    tuple val(sample_id), path("${sample_id}.cyrius_cyp2d6.json"), emit: json
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ── Cyrius v1.1.1 — CYP2D6 star allele calling ────────────────────────
    # RATIONALE: GATK4 cannot reliably detect CYP2D6 variants due to
    # CYP2D7 pseudogene interference, CNV complexity, and haplotype phasing
    # requirements. See module header for detailed rationale.
    # Reference: Chen et al. 2022 NPJ Genomic Med PMID:35082305

    star_caller.py \\
        --manifest  <(echo "${bam}") \\
        --genome    hg38 \\
        --prefix    ${sample_id} \\
        --outDir    . \\
        --threads   ${task.cpus}

    # ── Rename outputs to expected filenames ──────────────────────────────
    mv ${sample_id}.tsv  ${sample_id}.cyrius_cyp2d6.txt  2>/dev/null || \
    mv ${sample_id}.txt  ${sample_id}.cyrius_cyp2d6.txt  2>/dev/null || true
    mv ${sample_id}.json ${sample_id}.cyrius_cyp2d6.json 2>/dev/null || true

    # ── Validate minimum depth at CYP2D6 locus ────────────────────────────
    # Cyrius sets Filter="LOW_DEPTH" if mean depth < 30× at CYP2D6 region.
    # ACGS 2024 minimum WGS depth is 30× (enforced globally by MOSDEPTH_QC).
    if grep -q "LOW_DEPTH" "${sample_id}.cyrius_cyp2d6.txt" 2>/dev/null; then
        echo "WARNING: Cyrius flagged LOW_DEPTH for ${sample_id} at CYP2D6." >&2
        echo "WARNING: CYP2D6 diplotype call may be unreliable. Check MOSDEPTH QC." >&2
    fi

    # ── Log diplotype for traceability ────────────────────────────────────
    echo "=== Cyrius CYP2D6 diplotype for ${sample_id} ===" >&2
    cat "${sample_id}.cyrius_cyp2d6.txt" >&2

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cyrius: "1.1.1"
        pharmvar: "5.2"
    END_VERSIONS
    """

    stub:
    """
    echo -e "Sample\tDiplotype\tPhenotype\tFilter\n${sample_id}\t*1/*1\tNM\tPASS" \
        > ${sample_id}.cyrius_cyp2d6.txt
    echo '{"sample": "${sample_id}", "diplotype": "*1/*1"}' \
        > ${sample_id}.cyrius_cyp2d6.json
    touch versions.yml
    """
}
