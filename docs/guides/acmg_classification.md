# ACMG Variant Classification in ClaritySeq

## Guidelines implemented

ClaritySeq implements ACMG/AMP classification under four overlapping guideline documents:

| Guideline | Reference | Notes |
|-----------|-----------|-------|
| ACMG/AMP 2015 | Richards et al. 2015 PMID:25741868 | Original 28-rule framework |
| Bayesian scoring | Tavtigian et al. 2020 PMID:32645316 | Point-based Bayesian system |
| Splicing subgroup | Walker et al. 2023 PMID:36898414 | PP3/BP4/BP7 for splice variants |
| ACGS 2024 v1.2 | Durkie et al. Feb 2024 | UK Best Practice; replaces ACGS 2020 |
| ClinGen SVI 2024 | ClinGen working group | PM2→Supporting; AlphaMissense approved |

---

## Key change: PM2 is Supporting (not Moderate)

**Before (ACGS 2020)**: PM2 applied at **Moderate** weight (2 points).

**After (ClinGen SVI 2024 + ACGS 2024)**: PM2 applied at **Supporting** weight (1 point).

### Why?

gnomAD v4.1 (released April 19, 2024) contains 807,162 individuals — the largest population database ever. Analysis of this dataset revealed that:

- Ultra-rare variants (AF < 0.0001) are far more common than assumed in 2015
- Absence from (or extreme rarity in) population databases is less distinctive evidence than it was when gnomAD contained ~100,000 individuals
- Many ultra-rare variants in dominant disease genes turn out to be population-specific common variants at low global frequency

### Impact

PM2 downgrade from Moderate (2 pts) to Supporting (1 pt) affects variant classifications near the LP/VUS threshold. ClaritySeq implements the **novel ClinGen SVI 2024 combination** to address this:

```
PVS1 (8 pts) + PM2_Supporting (1 pt) = 9 pts → Likely Pathogenic
```

Without this combination rule, novel LoF variants in LoF-mechanism genes (where rarity is the only secondary evidence) would fall below the LP threshold (≥6 pts) after the PM2 downgrade.

**Implementation**: `bayesacmg/src/bayesacmg/rules/pathogenic.py:rule_pm2()`
**Combination**: `bayesacmg/src/bayesacmg/combinations.py`

### VCEP override

Some VCEP gene-specific specifications allow PM2 at Moderate for specific gene-disease pairs. Check `bayesacmg/src/bayesacmg/vcep_client.py` before applying the default Supporting weight.

---

## AlphaMissense as primary PP3/BP4 predictor

**ClinGen SVI 2024 approved four in silico tools** for PP3/BP4 evidence. ClaritySeq uses AlphaMissense as the **primary** tool:

| Tool | PP3 threshold | BP4 threshold | Status in ClaritySeq |
|------|--------------|--------------|----------------------|
| **AlphaMissense** | **≥ 0.564** | **≤ 0.340** | **PRIMARY** |
| REVEL | ≥ 0.7 | ≤ 0.15 | Secondary comparator |
| BayesDel | ≥ 0.13 | ≤ -0.18 | From dbNSFP v4.7 |
| CADD PHRED | ≥ 25 | ≤ 15 | From dbNSFP v4.7 |

### Why AlphaMissense is primary

- AlphaMissense (Cheng et al. 2023 Science PMID:37703350) was fine-tuned from AlphaFold using population variant frequency patterns
- AUROC 0.91 on ClinVar benchmark; outperforms EVE for 77% of ACMG genes
- Covers all 71 million possible human missense variants
- Most recently ClinGen SVI-approved tool

### Ambiguous range

AlphaMissense scores between 0.340 and 0.564 are **ambiguous** — neither PP3 nor BP4 is applied. This range represents genuine uncertainty in the model.

**Implementation**: `bayesacmg/src/bayesacmg/rules/pathogenic.py:rule_pp3()`
**AlphaMissense client**: `annotation/alphamissense_client.py`

---

## Splicing evidence framework (Walker et al. 2023)

For canonical splice site variants and deep intronic variants, the ClinGen SVI Splicing Subgroup framework applies (Walker et al. 2023 PMID:36898414):

| SpliceAI Δ score | Evidence |
|-----------------|---------|
| ≥ 0.5 | Strong splice impact → PP3 Strong |
| ≥ 0.2 | Moderate splice impact → PP3 Moderate |
| 0.1–0.2 | Weak / inconclusive |
| < 0.1 | No splice impact predicted |
| Synonymous + < 0.1 | → BP7 (no splice impact for silent variant) |

**When Pangolin and SpliceAI disagree**: Use the more conservative (lower Δ score) estimate. Document the disagreement in `evidence_items`.

**Implementation**: `bayesacmg/src/bayesacmg/rules/splicing.py`

---

## MANE Select transcript notation

All HGVSc and HGVSp notation in ClaritySeq uses **MANE Select** transcripts.

MANE (Matched Annotation from NCBI and EBI) Select provides one biologically-relevant transcript per protein-coding gene, jointly maintained by NCBI RefSeq and EMBL-EBI Ensembl.

Reference: Morales et al. 2022 Nature 604:310 PMID:35356062

**VEP pick order**: `mane_select,mane_plus_clinical,canonical`

If a LoF variant only affects non-MANE transcripts, PVS1 strength is reduced by one level (e.g., Very Strong → Strong) per ACGS 2024 v1.2 §5.

**Implementation**: `annotation/mane_select.py:adjust_pvs1_for_mane()`

---

## Mitochondrial variant classification (ACGS 2024 §6)

Mitochondrial variants follow ACGS 2024 v1.2 §6 — a separate set of rules from nuclear variants:

1. **Haplogroup classification first** (Haplogrep3) — haplogroup-defining variants are automatically Benign; must be excluded before any pathogenicity assessment
2. **Heteroplasmy as %** (not genotype) — heteroplasmy level maps to clinical significance per ACGS 2024 §6
3. **Mito-specific BA1**: The standard BA1 threshold (AF > 5%) does NOT apply to mtDNA variants. Use mito-specific thresholds.
4. **MITOMAP** as primary disease database (not ClinVar alone)

**Implementation**: `bayesacmg/src/bayesacmg/rules/mito.py`

---

## VCEP gene-specific overrides

ClinGen Variant Curation Expert Panels (VCEPs) publish gene-specific specifications that override the general 28-rule framework for specific genes (e.g., BRCA1/2, RASopathy genes, CDH1).

ClaritySeq queries the **ClinGen CSpec registry API** before classifying each variant. If a VCEP specification exists for the gene, those thresholds and rule modifications take precedence.

Common VCEP overrides:
- BRCA1/2: Modified PM2 threshold; additional evidence codes
- RASopathy genes: Modified functional evidence requirements

**Implementation**: `bayesacmg/src/bayesacmg/vcep_client.py`
**API**: https://cspec.genome.network/cspec/api/svi/

---

## VUS review dates (ACGS 2024 §9)

ACGS 2024 v1.2 §9 requires: *"For VUS variants, review within 2 years is recommended."*

ClaritySeq implements this by:
1. Every VUS classification triggers creation of a `VUSReviewSchedule` record (variant_id, review_due = classification_date + 2 years)
2. The ClinVar reclassification daemon checks these dates weekly
3. The clinical report displays: **"Review by: [date]"** for each VUS
4. When a VUS approaches its review date, the daemon sends a notification (Slack/email)

**Implementation**: `reclassification/models.py:VUSReviewSchedule`, `reclassification/daemon.py`

---

## ClinVar submission (NHS mandatory — ACGS 2024)

ACGS 2024 v1.2 Introduction: *"Submission of variants to ClinVar by NHS laboratories in England is now a requirement following completion of the information governance review process."*

ClaritySeq automates ClinVar submission for novel P/LP variants:
1. Variant classified as P or LP (BayesACMG posterior probability ≥ 0.90)
2. Not already in ClinVar (no RCV/SCV accession)
3. Added to `ClinVarSubmissionQueue` table
4. `clinvar_submitter.py` generates ClinVar API submission XML
5. Submission includes: MANE Select HGVSc/HGVSp, ACMG classification, evidence codes (ACGS 2024 applied), BayesACMG posterior probability (as submission note), gnomAD v4.1 AF

**Implementation**: `reclassification/clinvar_submitter.py`
**Prerequisites**: `NCBI_API_KEY` and `CLINVAR_SUBMISSION_ORG_ID` in `.env`

---

## Bayesian posterior probabilities

ClaritySeq reports a **Bayesian posterior probability P(Pathogenic)** alongside the ACMG 5-tier classification, with a **95% highest-density interval (HDI)**.

Example report entry:
```
BRCA1 NM_007294.4:c.5266dupC (p.Gln1756ProfsTer74)
Class: Pathogenic | P(Path) = 99% [96%–100%] | gnomAD v4.1: absent
```

This Bayesian uncertainty quantification is a novel contribution of ClaritySeq — no other clinical WGS platform reports calibrated ACMG classification uncertainty in this format.

**Implementation**: `bayesacmg/src/bayesacmg/model.py`
**Calibration**: ECE < 0.05 on ClinGen 500-variant set (`calibration/`)
