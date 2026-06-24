# pgx/

Pharmacogenomics (PGx) pipeline components.

## CYP2D6 star allele genotyping (Cyrius)

**Why Cyrius instead of GATK4?**
CYP2D6 is a gene with extensive structural variation (copy number variants, gene conversions with the *CYP2D7* pseudogene). Standard short-read aligners and GATK4 HaplotypeCaller cannot reliably distinguish CYP2D6 from CYP2D7 — alignments are frequently misassigned, producing false variant calls. Cyrius uses a tailored algorithm that models CYP2D6 SV structure explicitly.

## Star allele nomenclature

Star alleles (e.g., `*1`, `*2`, `*4`) represent haplotypes with known functional effects. Diplotype (pair of star alleles) predicts metaboliser phenotype:

| Phenotype | Abbreviation | Clinical implication |
|-----------|-------------|----------------------|
| Normal Metaboliser | NM | Standard dosing |
| Intermediate Metaboliser | IM | Reduce dose for some drugs |
| Poor Metaboliser | PM | Avoid certain drugs; alternative required |
| Ultrarapid Metaboliser | UM | Increase dose or use alternative |

## Data sources

- **PharmVar**: Star allele definitions — https://www.pharmvar.org/
- **CPIC**: Drug dosing recommendations — https://cpicpgx.org/
- **Cyrius**: Container `clinicalgenomics/cyrius:1.1.1`
