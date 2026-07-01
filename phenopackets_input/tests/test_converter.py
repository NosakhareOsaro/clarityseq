"""
phenopackets_input.tests.test_converter
=========================================
pytest tests for Phenopackets v2.0 → variant filter parameter conversion.

Tests cover:
    - phenopacket_to_filter_params: HPO extraction, inheritance mode,
      disease IDs, sex, proband ID, de_novo detection.
    - extract_hpo_ids: HPO list extraction with/without inheritance terms.
    - VariantFilterParams dataclass defaults.

References:
    Jacobsen et al. 2022 Nature Biotechnology PMID:35705716 (Phenopackets v2).
    Köhler et al. 2021 PMID:33264411 (HPO).
"""

from __future__ import annotations

import pytest

from phenopackets_input.converter import (
    VariantFilterParams,
    extract_hpo_ids,
    phenopacket_to_filter_params,
)


# ---------------------------------------------------------------------------
# Minimal Phenopackets fixture
# ---------------------------------------------------------------------------

MINIMAL_PHENOPACKET = {
    "id": "PP-001",
    "subject": {"id": "PATIENT-001", "sex": "FEMALE"},
    "phenotypicFeatures": [
        {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
        {"type": {"id": "HP:0004322", "label": "Short stature"}, "excluded": False},
    ],
    "diseases": [
        {
            "term": {"id": "OMIM:615846", "label": "CDKL5 Deficiency Disorder"},
        }
    ],
    "metaData": {"created": "2024-01-01T00:00:00Z", "createdBy": "test"},
}

AD_PHENOPACKET = {
    "id": "PP-AD-001",
    "subject": {"id": "PATIENT-002", "sex": "MALE"},
    "phenotypicFeatures": [
        {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
        {"type": {"id": "HP:0000006", "label": "Autosomal dominant inheritance"}, "excluded": False},
    ],
    "diseases": [],
    "metaData": {},
}

AR_PHENOPACKET = {
    "id": "PP-AR-001",
    "subject": {"id": "PATIENT-003", "sex": "FEMALE"},
    "phenotypicFeatures": [
        {"type": {"id": "HP:0000007", "label": "Autosomal recessive inheritance"}, "excluded": False},
    ],
    "diseases": [],
    "metaData": {},
}


# ---------------------------------------------------------------------------
# VariantFilterParams dataclass
# ---------------------------------------------------------------------------


class TestVariantFilterParams:
    """Tests for VariantFilterParams defaults."""

    def test_default_sex_is_unknown(self) -> None:
        """Default sex is UNKNOWN_SEX."""
        params = VariantFilterParams()
        assert params.sex == "UNKNOWN_SEX"

    def test_default_inheritance_mode_empty(self) -> None:
        """Default inheritance_mode is empty string."""
        params = VariantFilterParams()
        assert params.inheritance_mode == ""

    def test_default_is_de_novo_false(self) -> None:
        """Default is_de_novo is False."""
        params = VariantFilterParams()
        assert params.is_de_novo is False

    def test_default_gene_panel_empty(self) -> None:
        """Default gene_panel is empty list."""
        params = VariantFilterParams()
        assert params.gene_panel == []


# ---------------------------------------------------------------------------
# phenopacket_to_filter_params tests
# ---------------------------------------------------------------------------


class TestPhenopacketToFilterParams:
    """Tests for phenopacket_to_filter_params()."""

    def test_hpo_terms_extracted(self) -> None:
        """Non-excluded HPO terms are extracted into hpo_terms."""
        params = phenopacket_to_filter_params(MINIMAL_PHENOPACKET)
        assert "HP:0001250" in params.hpo_terms
        assert "HP:0004322" in params.hpo_terms

    def test_excluded_hpo_terms_separated(self) -> None:
        """Excluded HPO terms go into excluded_hpo_terms, not hpo_terms."""
        pp = {
            "id": "PP-002",
            "subject": {"id": "P-002", "sex": "MALE"},
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
                {"type": {"id": "HP:0000252", "label": "Microcephaly"}, "excluded": True},
            ],
            "diseases": [],
        }
        params = phenopacket_to_filter_params(pp)
        assert "HP:0001250" in params.hpo_terms
        assert "HP:0000252" not in params.hpo_terms
        assert "HP:0000252" in params.excluded_hpo_terms

    def test_sex_female_extracted(self) -> None:
        """FEMALE sex is extracted from subject."""
        params = phenopacket_to_filter_params(MINIMAL_PHENOPACKET)
        assert params.sex == "FEMALE"

    def test_sex_male_extracted(self) -> None:
        """MALE sex is extracted from subject."""
        params = phenopacket_to_filter_params(AD_PHENOPACKET)
        assert params.sex == "MALE"

    def test_unknown_sex_falls_back(self) -> None:
        """Unrecognised sex string falls back to UNKNOWN_SEX."""
        pp = dict(MINIMAL_PHENOPACKET)
        pp["subject"] = {"id": "P", "sex": "NOT_SPECIFIED"}
        params = phenopacket_to_filter_params(pp)
        assert params.sex == "UNKNOWN_SEX"

    def test_proband_id_extracted(self) -> None:
        """Proband ID is extracted from subject.id."""
        params = phenopacket_to_filter_params(MINIMAL_PHENOPACKET)
        assert params.proband_id == "PATIENT-001"

    def test_disease_ids_extracted(self) -> None:
        """Disease OMIM IDs are extracted into disease_ids."""
        params = phenopacket_to_filter_params(MINIMAL_PHENOPACKET)
        assert "OMIM:615846" in params.disease_ids

    def test_autosomal_dominant_inheritance_detected(self) -> None:
        """HP:0000006 → inheritance_mode='AD'."""
        params = phenopacket_to_filter_params(AD_PHENOPACKET)
        assert params.inheritance_mode == "AD"

    def test_autosomal_recessive_inheritance_detected(self) -> None:
        """HP:0000007 → inheritance_mode='AR'."""
        params = phenopacket_to_filter_params(AR_PHENOPACKET)
        assert params.inheritance_mode == "AR"

    def test_inheritance_term_included_in_hpo_terms(self) -> None:
        """Inheritance mode HPO term is also included in hpo_terms."""
        params = phenopacket_to_filter_params(AD_PHENOPACKET)
        assert "HP:0000006" in params.hpo_terms

    def test_empty_phenopacket_returns_defaults(self) -> None:
        """Empty phenopacket returns default VariantFilterParams."""
        params = phenopacket_to_filter_params({})
        assert params.hpo_terms == []
        assert params.inheritance_mode == ""
        assert params.sex == "UNKNOWN_SEX"

    def test_de_novo_detected_from_pedigree(self) -> None:
        """De novo detected when proband is affected and parents are unaffected."""
        pp = {
            "id": "PP-003",
            "subject": {"id": "PROBAND", "sex": "MALE"},
            "phenotypicFeatures": [],
            "diseases": [],
            "family": {
                "pedigree": {
                    "persons": [
                        {"individualId": "PROBAND", "affectedStatus": "AFFECTED"},
                        {"individualId": "FATHER", "affectedStatus": "UNAFFECTED"},
                        {"individualId": "MOTHER", "affectedStatus": "UNAFFECTED"},
                    ]
                }
            },
        }
        params = phenopacket_to_filter_params(pp)
        assert params.is_de_novo is True

    def test_not_de_novo_when_parent_affected(self) -> None:
        """Not de novo when a parent is also affected."""
        pp = {
            "id": "PP-004",
            "subject": {"id": "PROBAND", "sex": "FEMALE"},
            "phenotypicFeatures": [],
            "diseases": [],
            "family": {
                "pedigree": {
                    "persons": [
                        {"individualId": "PROBAND", "affectedStatus": "AFFECTED"},
                        {"individualId": "FATHER", "affectedStatus": "AFFECTED"},
                        {"individualId": "MOTHER", "affectedStatus": "UNAFFECTED"},
                    ]
                }
            },
        }
        params = phenopacket_to_filter_params(pp)
        assert params.is_de_novo is False

    def test_first_disease_onset_extracted(self) -> None:
        """Age of onset is extracted from the first disease entry."""
        pp = {
            "id": "PP-005",
            "subject": {"id": "P", "sex": "MALE"},
            "phenotypicFeatures": [],
            "diseases": [
                {
                    "term": {"id": "OMIM:123456", "label": "Test disease"},
                    "onset": {
                        "ontologyClass": {"id": "HP:0003623", "label": "Neonatal onset"}
                    },
                }
            ],
        }
        params = phenopacket_to_filter_params(pp)
        assert params.age_of_onset == "Neonatal onset"

    def test_feature_missing_type_id_skipped(self) -> None:
        """Features without type.id are silently skipped."""
        pp = {
            "id": "PP-006",
            "subject": {"id": "P", "sex": "MALE"},
            "phenotypicFeatures": [
                {"type": {}, "excluded": False},
                {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
            ],
            "diseases": [],
        }
        params = phenopacket_to_filter_params(pp)
        assert params.hpo_terms == ["HP:0001250"]


# ---------------------------------------------------------------------------
# extract_hpo_ids tests
# ---------------------------------------------------------------------------


class TestExtractHpoIds:
    """Tests for extract_hpo_ids()."""

    def test_extracts_non_excluded_hpo(self) -> None:
        """Returns non-excluded HPO terms."""
        ids = extract_hpo_ids(MINIMAL_PHENOPACKET)
        assert "HP:0001250" in ids
        assert "HP:0004322" in ids

    def test_excludes_excluded_features(self) -> None:
        """Excluded features are not returned."""
        pp = {
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}, "excluded": False},
                {"type": {"id": "HP:0000252"}, "excluded": True},
            ]
        }
        ids = extract_hpo_ids(pp)
        assert "HP:0001250" in ids
        assert "HP:0000252" not in ids

    def test_excludes_inheritance_terms_by_default(self) -> None:
        """Inheritance mode HPO terms are excluded when exclude_inheritance=True."""
        ids = extract_hpo_ids(AD_PHENOPACKET, exclude_inheritance=True)
        assert "HP:0000006" not in ids
        assert "HP:0001250" in ids

    def test_includes_inheritance_terms_when_not_excluded(self) -> None:
        """Inheritance mode HPO terms are included when exclude_inheritance=False."""
        ids = extract_hpo_ids(AD_PHENOPACKET, exclude_inheritance=False)
        assert "HP:0000006" in ids

    def test_empty_phenopacket_returns_empty(self) -> None:
        """Empty phenopacket returns empty list."""
        ids = extract_hpo_ids({})
        assert ids == []

    def test_feature_without_type_id_skipped(self) -> None:
        """Feature without type.id is skipped."""
        pp = {
            "phenotypicFeatures": [
                {"type": {}, "excluded": False},
            ]
        }
        ids = extract_hpo_ids(pp)
        assert ids == []

    def test_non_dict_feature_skipped(self) -> None:
        """A non-dict entry in phenotypicFeatures is skipped without error."""
        pp = {
            "phenotypicFeatures": [
                "not-a-dict",
                {"type": {"id": "HP:0001250"}, "excluded": False},
            ]
        }
        ids = extract_hpo_ids(pp)
        assert ids == ["HP:0001250"]


# ---------------------------------------------------------------------------
# Malformed / non-dict entry handling
# ---------------------------------------------------------------------------


class TestMalformedEntries:
    """Tests for graceful handling of non-dict list entries."""

    def test_non_dict_phenotypic_feature_skipped(self) -> None:
        """A non-dict entry in phenotypicFeatures is skipped in filter params."""
        pp = {
            "id": "PP-007",
            "subject": {"id": "P", "sex": "MALE"},
            "phenotypicFeatures": [
                "not-a-dict",
                {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
            ],
            "diseases": [],
        }
        params = phenopacket_to_filter_params(pp)
        assert params.hpo_terms == ["HP:0001250"]

    def test_non_dict_disease_skipped(self) -> None:
        """A non-dict entry in diseases is skipped without raising."""
        pp = {
            "id": "PP-008",
            "subject": {"id": "P", "sex": "MALE"},
            "phenotypicFeatures": [],
            "diseases": [
                "not-a-dict",
                {"term": {"id": "OMIM:100100", "label": "Test disease"}},
            ],
        }
        params = phenopacket_to_filter_params(pp)
        assert params.disease_ids == ["OMIM:100100"]
