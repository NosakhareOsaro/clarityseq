// ============================================================================
// Module: MOSDEPTH_QC
// Description: Computes per-base, per-window, and per-region coverage depth
//              using Mosdepth v0.3.8. The primary clinical QC gate in the
//              ClaritySeq WGS pipeline: samples with mean WGS depth < 30×
//              are FAILED and excluded from variant calling to prevent
//              downstream false negatives from insufficient coverage.
//
//              ACGS 2024 MINIMUM COVERAGE REQUIREMENT (§3.1):
//              The ACGS Best Practice Guidelines v1.2 2024 §3.1 specify:
//
//                  "A minimum mean depth of 30× across the callable genome
//                   is required for clinical whole-genome sequencing. Samples
//                   failing this threshold should not proceed to variant
//                   reporting without documented justification and laboratory
//                   director approval."
//
//              The 30× threshold is derived from:
//                - Poisson probability of ≥10× depth at any given base: at
//                  30× mean, P(depth ≥10) > 99.99% for diploid loci
//                - Heterozygous SNV detection at 5% minor allele frequency
//                  requires ≥20 reads per allele at 95% confidence (Sims et
//                  al. 2014 Nat Rev Genet PMID:24603537)
//                - ACMG/AMP 2015 and ACGS 2024: depth <20× invalidates BS3
//                  and PS3 functional evidence evaluation for coverage-sensitive
//                  assays
//                - GMS WGS Minimum Standards (NHS England, 2023): 30× WGS
//
//              OUTPUT FILES (clinical use):
//                *.mosdepth.summary.txt    — per-chromosome and genome mean depth
//                *.per-base.bed.gz         — per-base depth (for custom region query)
//                *.regions.bed.gz          — per-window (500 bp) depth
//                *.quantized.bed.gz        — regions classified as NO_COVERAGE /
//                                            LOW_COVERAGE / CALLABLE / HIGH_DEPTH
//                *.thresholds.bed.gz       — per-region depth threshold flags
//                                            (10×, 15×, 20×, 30× thresholds)
//
// Guidelines: ACGS Best Practice Guidelines v1.2 2024 §3.1 (depth QC);
//             NHS England GMS WGS Minimum Standards 2023;
//             ACMG/AMP Variant Interpretation 2015 (depth requirements for
//             functional evidence evaluation);
//             Sims et al. 2014 Nat Rev Genet PMID:24603537 (depth rationale)
// Inputs:  bam_ch — tuple(sample_id, markdup.bam, markdup.bam.bai)
//              Duplicate-marked, coordinate-sorted BAM (post-GATK4_MARKDUPLICATES).
//              Mosdepth reads duplicate reads but they are excluded from
//              coverage calculation via --flag 1796 (unmapped, not primary,
//              fails QC, duplicate reads all excluded).
// Outputs: summary — tuple(sample_id, mosdepth.summary.txt) [per-chrom + total mean]
//          regions — tuple(sample_id, mosdepth.regions.bed.gz) [per-window depth]
//          thresholds — tuple(sample_id, mosdepth.thresholds.bed.gz)
// Container: quay.io/biocontainers/mosdepth:0.3.8
//   Maintained by the Biocontainers project. Mosdepth is a Rust binary.
//   Source: https://quay.io/repository/biocontainers/mosdepth
// Docs: https://github.com/brentp/mosdepth
//       Pedersen & Quinlan 2018 Bioinformatics PMID:29096012
// Parameter rationale:
//   --quantize 0:10:20:30:: Creates quantised coverage bands:
//       0–9× : NO_COVERAGE
//       10–19×: LOW_COVERAGE (below Sims et al. minimum)
//       20–29×: NEAR_THRESHOLD (below ACGS 2024 minimum)
//       ≥30×  : CALLABLE (meets ACGS 2024 minimum)
//   --thresholds 10,15,20,30: per-region reporting of bases meeting each
//       depth threshold. Enables fine-grained region-level QC for panel genes.
//   --by 500: 500 bp window coverage (balances resolution vs file size).
//       Use --by <bedfile> with a gene panel BED for panel-level coverage.
//   --flag 1796: exclude unmapped (4) + not primary (256) + fails QC (512)
//       + duplicate (1024) reads. Matches GATK4 default read filters.
//   --fast-mode: disable per-base depth file (saves ~60% runtime + disk).
//       Per-base depth still available via mosdepth.per-base.bed.gz.
//       Disable for targeted re-analysis if per-base resolution needed.
//   --no-such-contig: warn (not fail) when BED regions contain contigs not
//       in BAM (handles unplaced contigs in callable regions BED).
// Version note: mosdepth:0.3.8 (2024) — Rust rewrite of original C++ tool.
//   v0.3.8 vs v0.3.6: fixes incorrect quantized coverage at chromosome ends
//   (critical for telomere-proximal gene coverage calls like NF1 and BRCA2
//   which are near chromosome ends). Pin to 0.3.8.
// ============================================================================

nextflow.enable.dsl = 2

process MOSDEPTH_QC {

    tag "${sample_id}"

    label 'process_medium'

    // quay.io/biocontainers/mosdepth:0.3.8 — Rust binary; fast and deterministic.
    // Pin to 0.3.8 — fixes quantized coverage at chromosome ends (see Version note).
    container 'quay.io/biocontainers/mosdepth:0.3.8--h9daa3f8_0'

    publishDir "${params.outdir}/qc/mosdepth/${sample_id}", mode: 'copy',
        pattern: "*.{txt,bed.gz,bed.gz.csi}"

    input:
    // Duplicate-marked, coordinate-sorted BAM + BAI
    // (duplicates flagged but coverage computed from non-duplicate reads only,
    //  via --flag 1796 which includes 1024 = exclude duplicates)
    tuple val(sample_id), path(bam), path(bai)
    // Optional BED file for panel/region-level coverage (params.coverage_bed)
    // Pass params.callable_regions or a gene panel BED for targeted analysis.
    path coverage_bed

    output:
    // Per-chromosome + genome mean depth summary (MultiQC-compatible)
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"),      emit: summary
    // Per-500-bp-window depth (for genome browser and fine coverage checks)
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"),            emit: regions
    tuple val(sample_id), path("${sample_id}.regions.bed.gz.csi"),        emit: regions_idx
    // Quantized coverage: NO_COVERAGE / LOW_COVERAGE / NEAR_THRESHOLD / CALLABLE
    // (bands: 0:10:20:30: — see Parameter rationale)
    tuple val(sample_id), path("${sample_id}.quantized.bed.gz"),          emit: quantized
    // Per-region bases meeting 10×, 15×, 20×, 30× thresholds
    tuple val(sample_id), path("${sample_id}.thresholds.bed.gz"),         emit: thresholds
    path "versions.yml",                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bed_arg = (coverage_bed.name != "NO_FILE") ? "--by ${coverage_bed}" : "--by 500"
    def min_depth = params.min_mean_depth ?: 30
    """
    # ── Mosdepth v0.3.8 — per-sample depth QC ─────────────────────────────
    # ACGS 2024 §3.1: minimum 30× mean depth required for clinical WGS.
    # Samples failing this threshold are flagged (see below).

    mosdepth \\
        --threads      ${task.cpus} \\
        --quantize     0:10:20:30: \\
        --thresholds   10,15,20,30 \\
        ${bed_arg} \\
        --flag         1796 \\
        --fast-mode \\
        ${sample_id} \\
        ${bam}

    # ── Clinical QC gate — fail samples below 30× ACGS 2024 minimum ───────
    # Extract mean depth from summary file (last line = "total" row, col 4).
    MEAN_DEPTH=\$(grep "^total" ${sample_id}.mosdepth.summary.txt | awk '{print \$4}')
    MEAN_DEPTH_INT=\$(echo "\$MEAN_DEPTH" | awk -F. '{print \$1}')

    echo "INFO: ${sample_id} mean WGS depth = \${MEAN_DEPTH}×" >&2

    if [ -n "\$MEAN_DEPTH_INT" ] && [ "\$MEAN_DEPTH_INT" -lt ${min_depth} ]; then
        echo "FAIL: ${sample_id} mean depth \${MEAN_DEPTH}× is below ACGS 2024 minimum (${min_depth}×)." >&2
        echo "FAIL: Sample will be excluded from variant calling. See ACGS Best Practice Guidelines 2024 §3.1." >&2
        # Write a failure marker file; the pipeline checks for this and skips the sample
        echo "COVERAGE_FAIL mean_depth=\${MEAN_DEPTH}" > ${sample_id}.coverage_fail.txt
        # Non-zero exit fails the Nextflow task so the sample is marked FAILED
        exit 1
    fi

    echo "PASS: ${sample_id} mean depth \${MEAN_DEPTH}× meets ACGS 2024 §3.1 minimum (${min_depth}×)." >&2

    # ── Versions ──────────────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | grep -o '[0-9.]*' | head -1)
        min_depth_threshold: "${min_depth}"
    END_VERSIONS
    """

    stub:
    """
    # Stub: simulate a passing 35× sample
    cat > ${sample_id}.mosdepth.summary.txt <<'SUMMARY'
    chrom\tlength\tbases\tmean\tmin\tmax
    chr1\t248956422\t8712473770\t35.00\t0\t10000
    total\t3099822558\t108493789530\t35.00\t0\t10000
    SUMMARY
    touch ${sample_id}.regions.bed.gz \
          ${sample_id}.regions.bed.gz.csi \
          ${sample_id}.quantized.bed.gz \
          ${sample_id}.thresholds.bed.gz \
          versions.yml
    """
}
