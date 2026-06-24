---
title: 'GenomeForge: A DRAGEN-GATK Pangenome-Aware Clinical WGS Platform with
  Bayesian ACMG Classification (ACGS 2024 v1.2) and Automated ClinVar
  Reclassification Monitoring'
tags:
  - Python
  - genomics
  - whole-genome sequencing
  - variant classification
  - Bayesian statistics
  - ACGS 2024
  - ClinGen SVI
  - DRAGEN-GATK
  - pangenome
  - FHIR
  - GA4GH Beacon
authors:
  - name: "[Your Name]"
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: "[Your Institution]"
    index: 1
date: 2026-06-23
bibliography: paper.bib
---

# Summary

GenomeForge is a production-grade, research-novel whole-genome sequencing (WGS)
clinical variant interpretation platform implementing seven novel contributions
not present in any existing public pipeline. The platform integrates DRAGMAP +
GATK4 DRAGEN-GATK variant calling with DeepVariant v1.8.0 ensemble calling,
HPRC pangenome alignment (vg giraffe), T2T-CHM13 v2.0 alignment, VEP v111
annotation with AlphaMissense as the primary pathogenicity predictor, and a
Bayesian ACMG/AMP classifier (BayesACMG) implementing all 28 ACMG/AMP criteria
under ACGS Best Practice Guidelines 2024 v1.2 [@durkie2024] and ClinGen SVI 2024
recommendations. BayesACMG provides calibrated 95% Bayesian credible intervals
on variant classifications — enabling quantified uncertainty communication for
Variants of Uncertain Significance (VUS) counselling — a capability absent from
all existing clinical WGS platforms. The platform produces NHS GMS-style clinical
reports with automated VUS reclassification monitoring (Celery daemon), a GA4GH
Beacon v2.1.1 federation API with VRS v2.0 identifiers, and NHS-mandated ClinVar
submission workflow per ACGS 2024.

# Statement of Need

Clinical WGS interpretation pipelines face three unresolved challenges. First,
fixed point-based ACMG scoring [@richards2015] provides no uncertainty
quantification — a variant scoring 5 of 6 points for Likely Pathogenic receives
the same classification as one scoring 10 points, despite vastly different
evidence strength. This creates difficulties in VUS counselling where patients
and clinicians need to understand classification certainty. Second, existing
pipelines implement ACGS 2020 guidelines and have not been updated to ACGS 2024
v1.2 [@durkie2024], which introduces critical changes: PM2 at Supporting weight
(not Moderate) reflecting gnomAD v4.1's [@chen2024] revelation that ultra-rare
variants are far more common than assumed; and AlphaMissense [@cheng2023] as the
ClinGen SVI-approved primary PP3/BP4 predictor replacing REVEL. Third, no
existing public platform implements the NHS-mandated ClinVar submission workflow
introduced in ACGS 2024 Introduction.

GenomeForge addresses all three challenges while adding four further novel
contributions: three-way benchmarking (GRCh38 vs T2T-CHM13 vs HPRC pangenome)
on GIAB truth sets; GA4GH Beacon v2.1.1 federation with VRS v2.0; Phenopackets
v2 + Exomiser 14 phenotype-driven prioritisation; and FHIR R4 Genomics Reporting
IG v3.0.0 recontact task generation.

# Seven Novel Contributions

1. **Calibrated Bayesian ACMG uncertainty** (BayesACMG): Dirichlet-Multinomial
   model calibrated on ClinGen curations (ECE < 0.05); 95% HDI on all
   classifications; PM2 prior centred on Supporting per ClinGen SVI 2024.

2. **ACGS 2024 v1.2 full implementation**: PM2→Supporting; AlphaMissense primary
   PP3/BP4 (≥0.564/≤0.340); Walker 2023 splicing framework; ACGS 2024 §6 mito
   rules; NHS ClinVar submission mandate; VUS 2-year review scheduling.

3. **Three-way reference benchmark**: Systematic SNV/indel sensitivity comparison
   across GRCh38, T2T-CHM13 v2.0, and HPRC pangenome on HG001/HG002/HG005.

4. **GA4GH Beacon v2.1.1 federation**: VRS v2.0 identifiers; GA4GH Passports;
   compatible with European Beacon Network and EGA registration.

5. **Phenopackets v2 + Exomiser 14 integration**: Phenotype-driven variant
   prioritisation directly from clinical Phenopackets v2 input; FHIR R4
   Genomics Reporting IG v3.0.0 recontact Task output.

6. **DRAGEN-GATK + DeepVariant ensemble**: First open-source pipeline combining
   DRAGMAP (DRAGEN Best Practices since 2021; BQSR skipped per BQD model) with
   DeepVariant v1.8.0 (SPRQ support; pangenome-aware) in both INTERSECTION and
   UNION ensemble modes.

7. **Complete NHS GMS reporting**: MANE Select transcripts; Bayesian HDI column;
   VUS review dates; mito section (ACGS 2024 §6); JSON-LD audit trail with
   gnomAD/VEP/AlphaMissense provenance; ClinVar submission queue.

# Implementation

BayesACMG implements all 28 ACMG/AMP criteria [@richards2015] using the
Tavtigian et al. 2020 [@tavtigian2020] point-scoring framework within a
Dirichlet-Multinomial Bayesian model built with PyMC [@abril2023]. The PM2 prior
is set at Supporting strength (1 point; concentration parameter 5.0, reflecting
community uncertainty after the 2024 revision) rather than the Moderate weight
(2 points) used in pre-2024 implementations. For PP3/BP4, AlphaMissense
[@cheng2023] is queried first (thresholds ≥0.564/≤0.340 per ClinGen SVI 2024);
REVEL, BayesDel, and CADD from dbNSFP v4.7 are retained as secondary comparators.
The Walker et al. 2023 [@walker2023] splicing framework is implemented as a
separate `rules/splicing.py` module: SpliceAI Δ ≥ 0.5 → PP3 Strong; ≥ 0.2 →
PP3 Moderate; synonymous + Δ < 0.1 → BP7. Mitochondrial variants follow ACGS
2024 §6: Haplogrep3 classification runs before ACMG assessment, haplogroup-
defining variants are automatically classified Benign, and separate mito-specific
BA1/PM2 thresholds apply.

The pipeline uses DRAGMAP v1.3.0 as the primary aligner (GATK Best Practices
since 2021), with BQSR intentionally omitted in DRAGEN-GATK mode: HaplotypeCaller
4.6.0.0's BQD (Base Quality Dropoff) model replaces BQSR by modelling systematic
errors internally, and applying BQSR to DRAGMAP-aligned reads reduces accuracy
[@broad2022]. BWA-MEM2 with BQSR is retained as an explicit fallback
(`-profile bwa_mem2`) for settings where DRAGMAP hash tables are unavailable.

Ancestry-stratified VQSR uses somalier ancestry inference to select the
appropriate gnomAD v4.1 population subset (AFR/AMR/EAS/EUR/SAS) for training.
The ClinVar reclassification daemon runs weekly (Celery + Redis 7), downloading
the ClinVar FTP release, diffing against stored variant classifications, and
generating FHIR R4 Task resources (Genomics Reporting IG v3.0.0) for recontact
where reclassification is detected. VUS variants receive automated review date
scheduling (date + 2 years) per ACGS 2024 §9.

# References
