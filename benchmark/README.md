# benchmark/

Variant calling benchmarking against GIAB truth sets.

## Acceptance criteria (CI)

- **SNP sensitivity**: ≥ 99.0% vs GIAB HG001 chr22 truth (hap.py)
- **Indel sensitivity**: ≥ 98.5% vs GIAB HG001 chr22 truth (hap.py)
- **Pangenome arm**: SNP sensitivity ≥ GRCh38 arm (nightly CI only)

## Truth sets used

| Sample | GIAB ID | Ethnicity | Usage |
|--------|---------|-----------|-------|
| HG001 (NA12878) | GIAB v4.2.1 | CEU (European) | Primary CI benchmark |
| HG002 (NA24385) | GIAB v4.2.1 | AJ (Ashkenazi Jewish) | Three-way benchmark |
| HG005 (NA24631) | GIAB v4.2.1 | CHN (Chinese) | Three-way benchmark |

## Three-way benchmark

`results/` contains pre-run benchmark outputs comparing:
1. GRCh38 (DRAGMAP + GATK4 + DeepVariant ensemble)
2. T2T-CHM13 v2.0 (with CrossMap liftover)
3. HPRC pangenome (vg giraffe v1.1)

## Tools

- **hap.py** (Illumina): primary benchmarking tool
- **RTG vcfeval**: secondary benchmarking
- **ga4gh-benchmarking-tools**: GA4GH-compatible metrics

## Running

```bash
# Single sample
python benchmark/run_hap_py.py --vcf results/HG001.vcf.gz --truth GIAB/HG001_v4.2.1.vcf.gz

# Parse results
python benchmark/parse_results.py --input results/HG001_hap_py/ --assert-snp-sensitivity 99.0
```
