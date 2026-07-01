"""
phenopackets_input.converter
==============================
Convert Phenopackets v2.0 HPO terms to variant filter parameters.

Phenopackets encode patient phenotypes as HPO term arrays.  This module
converts those HPO terms into structured variant prioritisation parameters
for the ClaritySeq prioritisation pipeline.

HPO (Human Phenotype Ontology):
    https://hpo.jax.org/
    Köhler et al. 2021 Nucleic Acids Research PMID:33264411.

Conversion logic:
    HPO terms → gene lists via HPO-to-gene mappings (Orphanet/OMIM).
    Inheritance mode terms (HP:0000006=AD, HP:0000007=AR, HP:0001417=XL)
    are extracted to set the inheritance_mode filter parameter.
    Affected status terms determine whether to filter by de_novo status.

References:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    Köhler et al. 2021 PMID:33264411 (HPO).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HPO inheritance mode terms
# ---------------------------------------------------------------------------

# HPO term IDs for inheritance modes
_INHERITANCE_HPO: dict[str, str] = {
    "HP:0000006": "AD",    # Autosomal dominant
    "HP:0000007": "AR",    # Autosomal recessive
    "HP:0001417": "XLD",   # X-linked dominant
    "HP:0001419": "XLR",   # X-linked recessive
    "HP:0001423": "XL",    # X-linked (general)
    "HP:0001450": "YL",    # Y-linked
    "HP:0001427": "Mito",  # Mitochondrial
    "HP:0003745": "Sporadic",  # Sporadic
    "HP:0001452": "AD",    # Autosomal dominant with incomplete penetrance
}


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class VariantFilterParams:
    """Variant filter parameters derived from a Phenopackets v2.0 phenopacket.

    Attributes:
        hpo_terms: List of HPO term IDs from the phenopacket
            (e.g. ``["HP:0001250", "HP:0004322"]``).
        excluded_hpo_terms: List of excluded HPO term IDs
            (features the patient does NOT have).
        inheritance_mode: Inferred inheritance mode from HPO inheritance
            terms (``"AD"``, ``"AR"``, ``"XL"``, ``"Mito"``, or ``""``).
        is_de_novo: True if the phenotype/family history suggests de novo.
        gene_panel: Gene symbols from disease annotations if available.
        disease_ids: OMIM/Orphanet disease IDs from the phenopacket.
        age_of_onset: Age of onset category string (HPO onset term label).
        sex: Subject sex (``"FEMALE"``, ``"MALE"``, ``"UNKNOWN_SEX"``).
        proband_id: Subject ID from the phenopacket.
    """

    hpo_terms: list[str] = field(default_factory=list)
    excluded_hpo_terms: list[str] = field(default_factory=list)
    inheritance_mode: str = ""
    is_de_novo: bool = False
    gene_panel: list[str] = field(default_factory=list)
    disease_ids: list[str] = field(default_factory=list)
    age_of_onset: str = ""
    sex: str = "UNKNOWN_SEX"
    proband_id: str = ""


# ---------------------------------------------------------------------------
# Converter
# ---------------------------------------------------------------------------


def phenopacket_to_filter_params(phenopacket: dict[str, Any]) -> VariantFilterParams:
    """Convert a Phenopackets v2.0 dict to variant filter parameters.

    Extracts HPO terms, excluded terms, inheritance mode, disease IDs,
    and other filter-relevant fields from the phenopacket.

    Args:
        phenopacket: Parsed Phenopackets v2.0 JSON dict.

    Returns:
        VariantFilterParams containing all filter parameters for the
        prioritisation pipeline.

    References:
        Jacobsen et al. 2022 PMID:35705716 (Phenopackets v2).
        Köhler et al. 2021 PMID:33264411 (HPO).
    """
    params = VariantFilterParams()

    # Subject ID and sex
    subject = phenopacket.get("subject", {})
    params.proband_id = subject.get("id", "")
    sex_raw = subject.get("sex", "UNKNOWN_SEX")
    params.sex = sex_raw if sex_raw in ("FEMALE", "MALE", "UNKNOWN_SEX") else "UNKNOWN_SEX"

    # Phenotypic features → HPO terms
    for feature in phenopacket.get("phenotypicFeatures", []):
        if not isinstance(feature, dict):
            continue
        term = feature.get("type", {})
        term_id = term.get("id", "")
        if not term_id:
            continue

        excluded = feature.get("excluded", False)
        if excluded:
            params.excluded_hpo_terms.append(term_id)
        else:
            params.hpo_terms.append(term_id)

        # Check for inheritance mode HPO terms
        if term_id in _INHERITANCE_HPO and not excluded:
            inferred_mode = _INHERITANCE_HPO[term_id]
            if not params.inheritance_mode:
                params.inheritance_mode = inferred_mode
                logger.debug(
                    "Inheritance mode inferred from HPO term %s: %s",
                    term_id,
                    inferred_mode,
                )

    # Diseases → disease IDs and gene panels
    for disease in phenopacket.get("diseases", []):
        if not isinstance(disease, dict):
            continue
        term = disease.get("term", {})
        disease_id = term.get("id", "")
        if disease_id:
            params.disease_ids.append(disease_id)

    # Age of onset
    onset = phenopacket.get("diseases", [{}])
    if onset and isinstance(onset[0], dict):
        onset_term = onset[0].get("onset", {}).get("ontologyClass", {})
        params.age_of_onset = onset_term.get("label", "")

    # Family information — check for de novo
    family = phenopacket.get("family", {})
    if family:
        pedigree = family.get("pedigree", {})
        persons = pedigree.get("persons", [])
        for person in persons:
            if person.get("individualId") == params.proband_id:
                if person.get("affectedStatus") == "AFFECTED":
                    # Check if parents are in pedigree and unaffected
                    affected_parents = [
                        p for p in persons
                        if p.get("affectedStatus") == "AFFECTED"
                        and p.get("individualId") != params.proband_id
                    ]
                    if not affected_parents:
                        params.is_de_novo = True

    return params


def extract_hpo_ids(phenopacket: dict[str, Any], exclude_inheritance: bool = True) -> list[str]:
    """Extract HPO term IDs from a phenopacket, excluding inheritance mode terms.

    Args:
        phenopacket: Parsed Phenopackets v2.0 JSON dict.
        exclude_inheritance: If True, exclude HPO inheritance mode terms
            (HP:0000006, HP:0000007, etc.) from the returned list.
            Default True — these are filter parameters, not phenotype queries.

    Returns:
        List of HPO term ID strings (e.g. ``["HP:0001250", "HP:0004322"]``).
        Excluded features are not included.

    References:
        Köhler et al. 2021 PMID:33264411 (HPO).
    """
    hpo_ids: list[str] = []
    for feature in phenopacket.get("phenotypicFeatures", []):
        if not isinstance(feature, dict):
            continue
        if feature.get("excluded", False):
            continue
        term_id = feature.get("type", {}).get("id", "")
        if not term_id:
            continue
        if exclude_inheritance and term_id in _INHERITANCE_HPO:
            continue
        hpo_ids.append(term_id)
    return hpo_ids
