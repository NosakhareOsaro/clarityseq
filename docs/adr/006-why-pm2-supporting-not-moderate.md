# ADR-006: Why PM2 Weight = Supporting (1 pt) Not Moderate (2 pts)

**Status**: Accepted
**Date**: 2026-06-23
**Supersedes**: Previous implementation using PM2 at Moderate (2 pts)

---

## Context

PM2 (absent from/extremely rare in population databases) was defined in Richards et al. 2015 (PMID:25741868) as a Moderate pathogenicity criterion (2 points in the Tavtigian 2020 scoring system). ClaritySeq must decide which weight to apply.

The ClinGen SVI Working Group published 2024 guidance recommending that PM2 be applied at **Supporting strength (1 pt)** in most contexts, citing evidence from gnomAD v4.1.

## Decision

**PM2 is applied at Supporting weight (1 point) by default in ClaritySeq.**

VCEP gene-specific specifications may override this to Moderate (2 pts) for specific gene-disease pairs. Check `bayesacmg/src/bayesacmg/vcep_client.py` before every PM2 application.

## Rationale

### gnomAD v4.1 evidence

gnomAD v4.1 (released April 19, 2024) contains 807,162 individuals:
- 730,947 exomes (including 416,555 UK Biobank exomes not in v3)
- 76,215 genomes

Analysis of gnomAD v4.1 revealed that ultra-rare variants (AF < 0.0001) are far more common in the general population than assumed in 2015 when gnomAD contained ~100,000 individuals. Many variants classified as "absent from population databases" in gnomAD v2/v3 turned out to be population-specific ultra-rare variants in gnomAD v4.1's larger, more diverse dataset.

This means **rarity alone is weaker evidence of pathogenicity** than previously thought. Applying 2 points for PM2 overstates the evidence.

### ACGS 2024 v1.2 §5 alignment

ACGS 2024 v1.2 Appendix C includes a "PM2 mini impact assessment" noting the gnomAD v4.1 findings and recommending Supporting strength as the default implementation for UK laboratories.

### Impact on novel LoF variants

The PM2 downgrade from Moderate (2 pts) to Supporting (1 pt) reduces the score for a typical novel LoF variant:

**Before (ACGS 2020)**:
- PVS1 (8) + PM2_Moderate (2) = 10 pts → Pathogenic

**After (ACGS 2024)**:
- PVS1 (8) + PM2_Supporting (1) = 9 pts → Still Likely Pathogenic (threshold ≥ 6 pts)

ClaritySeq implements the ClinGen SVI 2024 explicit combination: `PVS1 + PM2_Supporting = LP` at 9 pts, preserving LP classification for novel LoF variants in LoF-mechanism genes where rarity is the only secondary evidence. See `bayesacmg/src/bayesacmg/combinations.py`.

## Consequences

### Positive
- Aligned with ClinGen SVI 2024 and ACGS 2024 v1.2
- Reduces over-classification of ultra-rare variants as P/LP based on rarity alone
- More accurately reflects the evidential weight of population frequency data

### Negative
- Breaks compatibility with pipelines implementing ACGS 2020 (PM2 at Moderate)
- May require re-classification of existing VUS variants that were LP due to PM2 at Moderate
- Requires the `PVS1+PM2_Supporting=LP` combination rule to preserve LP for typical novel LoF variants

### Migration note

Existing classifications using PM2 at Moderate should be flagged for reclassification review. The ClinVar reclassification daemon (`reclassification/daemon.py`) can be configured to flag such variants.

## Alternatives considered

### Keep PM2 at Moderate (ACGS 2020)
Rejected: Contradicts ClinGen SVI 2024 guidance and ACGS 2024 v1.2. Would misrepresent current scientific consensus. gnomAD v4.1 data specifically warrants the downgrade.

### Apply PM2 at Moderate when gnomAD AF = 0 (completely absent)
Partially considered: The gnomAD v4.1 rationale applies even to completely absent variants — they may exist in unsequenced populations. Supporting is appropriate even for complete absence.

### VCEP-dependent weight only (no default)
Rejected: Most genes do not have VCEP specifications. A clear default is needed. Supporting is the appropriate default.

## Implementation

```python
# bayesacmg/src/bayesacmg/rules/pathogenic.py
def rule_pm2(variant: VariantInput) -> ACMGRule:
    # CRITICAL: return SUPPORTING (1 pt), NOT MODERATE (2 pts)
    return ACMGRule(
        rule_id="PM2",
        strength=EvidenceStrength.SUPPORTING,  # ClinGen SVI 2024
        ...
    )
```

## References

- ClinGen SVI PM2 guidance (2024): https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
- ACGS 2024 v1.2 §5 Table 2 and Appendix C (Durkie et al.)
- gnomAD v4.1: Chen et al. 2024 Nature doi:10.1038/s41586-024-07701-9
- Richards et al. 2015 PMID:25741868 (original PM2 definition)
- Tavtigian et al. 2020 PMID:32645316 (point-scoring system)
