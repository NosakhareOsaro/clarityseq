# ADR-005: Why AlphaMissense as Primary PP3/BP4 Predictor

**Status**: Accepted
**Date**: 2026-06-23
**Supersedes**: N/A (new decision)

---

## Context

The ClinGen Sequence Variant Interpretation (SVI) working group approved four computational tools for PP3/BP4 evidence in 2024:
1. AlphaMissense (Cheng et al. 2023, Science PMID:37703350)
2. REVEL (Ioannidis et al. 2016 PMID:27666373)
3. BayesDel (Feng et al. 2017 PMID:28324717)
4. CADD PHRED (Kircher et al. 2014 PMID:24487276)

ClaritySeq must select a primary tool and specify clear thresholds.

## Decision

**AlphaMissense is the primary PP3/BP4 predictor in ClaritySeq.**

Thresholds (ClinGen SVI 2024 calibration):
- AM score ≥ 0.564 → PP3 (Supporting Pathogenic)
- AM score ≤ 0.340 → BP4 (Supporting Benign)
- 0.340 < score < 0.564 → ambiguous (no evidence applied)

REVEL, BayesDel, and CADD are retained in dbNSFP v4.7 output as secondary comparators but do not determine PP3/BP4 classification unless AlphaMissense is unavailable (indels, variants outside model coverage).

## Rationale

### Why AlphaMissense over REVEL?

1. **Performance**: AlphaMissense achieves AUROC 0.91 on the ClinVar P/LP vs B/LB benchmark; outperforms EVE for 77% of ACMG disease genes (Cheng et al. 2023).

2. **Coverage**: AlphaMissense scores all 71 million possible human missense variants. REVEL requires pre-computed dbNSFP scores and has gaps in rare amino acid substitutions.

3. **Technology**: AlphaMissense uses the AlphaFold2 protein structure backbone fine-tuned with population frequency patterns. It represents the current state of the art in missense pathogenicity prediction.

4. **Recency**: AlphaMissense was approved by ClinGen SVI in 2024 — the most recent tool approval. REVEL was published in 2016 and reflects older training data and models.

5. **ClinGen SVI recommendation**: The 2024 guidance explicitly notes AlphaMissense as the recommended primary tool when available.

### Why not use REVEL as primary?

REVEL was the de facto primary tool in ClaritySeq's predecessor pipelines (ACGS 2020 era). It remains valid under ClinGen SVI 2024 but is outperformed by AlphaMissense on every benchmark measure. Retaining REVEL as primary would be a missed opportunity to use the most accurate approved tool.

### Ambiguous range

The 0.340–0.564 range covers approximately 18% of all missense variants — these are genuinely uncertain by the model. Applying any evidence code in this range would introduce noise. The conservative choice (no evidence) is correct here per Bayesian principles: uncertain predictor evidence contributes minimal information.

### For splice-impacting variants

This rule does NOT apply to splice variants. For variants with significant splice impact (canonical splice sites, deep intronic), use the Walker et al. 2023 (PMID:36898414) splicing framework in `bayesacmg/src/bayesacmg/rules/splicing.py`. SpliceAI ≥ 0.5 → PP3 Strong; SpliceAI ≥ 0.2 → PP3 Moderate.

## Consequences

### Positive
- Uses the most accurate ClinGen SVI 2024-approved tool
- Clear, well-calibrated thresholds from the SVI working group
- Complete coverage of missense variant space

### Negative
- AlphaMissense does not score indels or splice variants (handled by separate modules)
- ~2.7 GB data file must be downloaded and tabix-indexed (see `docs/guides/data_setup.md`)
- Cannot be used without the pre-downloaded TSV file (no API available)

## Alternatives considered

| Tool | AUROC (ClinVar) | ClinGen SVI 2024 status | Decision |
|------|-----------------|------------------------|----------|
| AlphaMissense | 0.91 | Approved (2024) | **PRIMARY** |
| REVEL | 0.87 | Approved (retained) | Secondary comparator |
| BayesDel | 0.85 | Approved | From dbNSFP v4.7 |
| CADD PHRED | 0.82 | Approved | From dbNSFP v4.7 |

## References

- Cheng et al. 2023 Science doi:10.1126/science.adg7492 PMID:37703350
- ClinGen SVI in silico calibration memo v2024-05
- ACGS 2024 v1.2 §5 (computational evidence)
