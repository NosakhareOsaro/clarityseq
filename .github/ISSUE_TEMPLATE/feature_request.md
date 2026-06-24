---
name: Feature request
about: Suggest a new feature or enhancement for GenomeForge
title: "[FEAT] "
labels: enhancement
assignees: ""
---

## Feature description

<!-- What would you like added or changed? Be specific — what is the new capability? -->

## Use case / motivation

<!-- Who would use this feature, and in what clinical or research context?
     Is this a requirement from an accreditation body, NHS policy, or ACGS guideline?
     Is this blocking a clinical use case? -->

## ACMG/Clinical guideline context

<!-- If this feature implements or updates a specific guideline, complete this section.
     All clinical changes MUST cite their evidence source. -->

- Guideline / evidence source:
  `ACGS 2024 v1.2 §__ / ClinGen SVI 2024 (date) / Richards 2015 PMID:25741868 / Other:`
- Current behaviour in GenomeForge:
- Desired behaviour:
- ACMG rule(s) affected: `PVS1 / PM2 / PP3 / BA1 / BS1 / ...`
- Does this change any ACMG rule weight? `Y / N`
  If yes: what weight changes, and what is the evidence level (PS / PM / PP / BA / BS / BP)?
- Does this affect variant pathogenicity calls? `Y / N`
  If yes: a re-validation against GIAB HG001 and/or ClinVar gold-star variants will be required.

## Proposed implementation

<!-- Optional but encouraged: suggest how this could be implemented.
     For Nextflow modules: which container? which tool version? where in the data flow?
     For Python: which module/class? approximate interface? -->

**Component affected:**
- [ ] New Nextflow module (`modules/<tool>/main.nf`)
- [ ] Existing Nextflow module update
- [ ] BayesACMG rule (`bayesacmg/src/bayesacmg/rules/`)
- [ ] Annotation client (`annotation/`)
- [ ] Reporting (`reporting/`)
- [ ] ClinVar submission (`reclassification/`)
- [ ] CI / benchmarking
- [ ] Documentation only

**Estimated implementation complexity:**
- [ ] Small (< 1 day, < 100 lines)
- [ ] Medium (1–3 days, 1–2 files)
- [ ] Large (> 3 days, multiple files or new container)

## Acceptance criteria

<!-- How will we know this feature is correctly implemented?
     Define testable, measurable criteria. -->

- [ ]
- [ ]
- [ ]

## Alternatives considered

<!-- What alternatives were evaluated and why were they rejected? -->

## Priority

- [ ] P1 — Critical: ACGS 2024 compliance gap or clinical safety issue
- [ ] P2 — High: affects classification accuracy or pipeline correctness
- [ ] P3 — Medium: quality of life, performance, or new guideline implementation
- [ ] P4 — Low: nice to have; no clinical urgency

## References

<!-- Papers, guidelines, or GitHub issues that are relevant -->
-
