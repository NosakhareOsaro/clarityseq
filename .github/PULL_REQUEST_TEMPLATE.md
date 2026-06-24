## Summary

<!-- Describe what this PR does and why. For clinical changes, explain the
     patient-safety rationale and cite the relevant guideline. -->

## Type of change

- [ ] Bug fix (non-breaking fix for incorrect behaviour)
- [ ] New feature (new module, tool, or capability)
- [ ] ACMG/clinical guideline update (fill out the Clinical Compliance section below)
- [ ] Nextflow module (new or updated DSL2 module in modules/)
- [ ] Pipeline (change to pipelines/*.nf)
- [ ] Refactoring / non-breaking cleanup (no behaviour change)
- [ ] Documentation only
- [ ] Infrastructure / CI / Docker

## Linked issues

<!-- Reference any GitHub issues this PR closes or relates to -->
Closes #

## Changes made

<!-- Bullet-point list of the specific changes in this PR -->
-
-

## Clinical compliance (complete for any change affecting variant classification or reporting)

<!-- CRITICAL: changes that affect classification output require all checkboxes below -->

- [ ] PM2 weight remains at **Supporting (1 pt)** per ClinGen SVI 2024
      (NEVER change to Moderate without a new ClinGen SVI publication and team consensus)
- [ ] AlphaMissense PP3 threshold ≥ **0.564**, BP4 threshold ≤ **0.340** per ClinGen SVI 2024
- [ ] Any new/changed guideline threshold cites the specific PMID or guideline section
- [ ] MANE Select transcripts used for all HGVSc/HGVSp notation (ACGS 2024 §4.1)
- [ ] VUS review date scheduling logic is unaffected, or if changed: reviewed by clinical lead
      (ACGS 2024 §9: VUS must be re-evaluated every 18 months minimum)
- [ ] Haplogroup classification still precedes mito ACMG assessment (ACGS 2024 §6)
- [ ] Per-ancestry VQSR mapping unchanged, or if changed: re-validated against GIAB truth sets

## Pipeline/module changes (complete for Nextflow changes)

- [ ] All new Nextflow processes have the required 12-line header comment block
- [ ] Container image pinned to exact version tag (no `latest`)
- [ ] `versions.yml` emitted from every new process
- [ ] `stub:` block present in every new process
- [ ] Module imports verified against actual process names in the .nf file
- [ ] `publishDir` paths follow the pattern `${params.outdir}/<tool>/<sample_id>`

## Tests

- [ ] Unit tests pass: `pytest bayesacmg/tests/ annotation/tests/ --cov-fail-under=90`
- [ ] ensemble.py tests pass: `pytest modules/ensemble/`
- [ ] Pre-commit hooks pass: `pre-commit run --all-files`
- [ ] If modifying pipeline: `nextflow run pipelines/wgs_grch38.nf -profile test,stub`
- [ ] If modifying Docker: image builds, passes `--help`, and Hadolint finds no errors
- [ ] If modifying VQSR: re-run hap.py benchmark and confirm SNP sensitivity ≥ 99.0%

## Code quality checklist

- [ ] Type annotations on all new Python functions (`-> ReturnType`)
- [ ] Google-style docstrings on all new public functions
- [ ] Every ACMG rule function cites: PMID:25741868 + ClinGen SVI section + ACGS 2024 §
- [ ] Every threshold/cutoff has an inline comment citing the paper/guideline
- [ ] No Python file exceeds 500 lines (split into sub-modules if needed)
- [ ] Tool versions pinned (no `latest` in Nextflow modules or Dockerfiles)
- [ ] No secrets or credentials committed (detect-secrets pre-commit hook passes)
- [ ] `CHANGELOG.md` updated (if user-facing change)

## Reviewer guidance

<!-- Tell the reviewer where to focus their attention -->
Key files to review:
-
-

Potential risks:
-
