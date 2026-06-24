"""Tests for the ClinVar VCF diff module.

Uses synthetic ClinVar VCF data to verify:
- VCF parsing correctness.
- Reclassification event detection between two VCF releases.
- Cross-reference of local variants against ClinVar data.
- Edge cases: new variants, missing fields, low review stars.
"""

from __future__ import annotations

import gzip
import tempfile
from datetime import date, datetime
from pathlib import Path
from typing import Generator

import pytest

from reclassification.clinvar_diff import (
    ClinVarRecord,
    _iter_vcf_records,
    _parse_clinvar_date,
    _parse_info,
    _parse_review_stars,
    _requires_recontact,
    diff_variants,
    find_reclassified_variants,
)
from reclassification.models import ClinicalSignificance


# ---------------------------------------------------------------------------
# Synthetic ClinVar VCF content helpers
# ---------------------------------------------------------------------------

# Minimal VCF header
VCF_HEADER = """\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ClinVar allele ID">
##INFO=<ID=CLNDATE,Number=.,Type=String,Description="Date last evaluated">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _make_vcf_line(
    chrom: str = "chr17",
    pos: int = 43094692,
    accession: str = "RCV000112345",
    ref: str = "G",
    alt: str = "A",
    clnsig: str = "Pathogenic",
    clnrevstat: str = "criteria_provided,_single_submitter",
    alleleid: int = 12345,
    clndate: str = "2024-01-15",
) -> str:
    """Build a single synthetic ClinVar VCF data line."""
    info = (
        f"CLNSIG={clnsig};"
        f"CLNREVSTAT={clnrevstat};"
        f"ALLELEID={alleleid};"
        f"CLNDATE={clndate}"
    )
    return f"{chrom}\t{pos}\t{accession}\t{ref}\t{alt}\t.\t.\t{info}\n"


def _write_vcf(lines: list[str], compressed: bool = False) -> Path:
    """Write VCF content to a temporary file and return its Path."""
    suffix = ".vcf.gz" if compressed else ".vcf"
    tmp = tempfile.NamedTemporaryFile(
        suffix=suffix, delete=False, mode="wb" if compressed else "w"
    )
    content = VCF_HEADER + "".join(lines)
    if compressed:
        tmp.write(gzip.compress(content.encode()))
    else:
        tmp.write(content)
    tmp.flush()
    return Path(tmp.name)


# ---------------------------------------------------------------------------
# Tests: _parse_info
# ---------------------------------------------------------------------------


class TestParseInfo:
    """Tests for the INFO field parser."""

    def test_simple_key_value(self):
        info = _parse_info("CLNSIG=Pathogenic;ALLELEID=123")
        assert info["CLNSIG"] == "Pathogenic"
        assert info["ALLELEID"] == "123"

    def test_flag_field(self):
        """Flag fields without '=' should be stored as 'true'."""
        info = _parse_info("CLNSIG=VUS;FLAG")
        assert info["FLAG"] == "true"

    def test_empty_string(self):
        info = _parse_info("")
        assert info == {}

    def test_multiple_values_with_pipe(self):
        """Pipe-separated values should be preserved as-is by _parse_info."""
        info = _parse_info("CLNSIG=Pathogenic|Likely_pathogenic")
        assert info["CLNSIG"] == "Pathogenic|Likely_pathogenic"


# ---------------------------------------------------------------------------
# Tests: _parse_review_stars
# ---------------------------------------------------------------------------


class TestParseReviewStars:
    """Tests for ClinVar review star rating parsing."""

    def test_single_submitter(self):
        assert _parse_review_stars("criteria_provided,_single_submitter") == 1

    def test_expert_panel(self):
        assert _parse_review_stars("reviewed_by_expert_panel") == 3

    def test_practice_guideline(self):
        assert _parse_review_stars("practice_guideline") == 4

    def test_no_assertion(self):
        assert _parse_review_stars("no_assertion_provided") == 0

    def test_conflicting(self):
        assert _parse_review_stars("criteria_provided,_conflicting_interpretations") == 1

    def test_multiple_submitters(self):
        assert _parse_review_stars(
            "criteria_provided,_multiple_submitters,_no_conflicts"
        ) == 2

    def test_unknown_value(self):
        assert _parse_review_stars("something_new") == 0

    def test_pipe_separated(self):
        """Should parse only the first pipe-separated value."""
        result = _parse_review_stars(
            "reviewed_by_expert_panel|criteria_provided,_single_submitter"
        )
        assert result == 3


# ---------------------------------------------------------------------------
# Tests: _parse_clinvar_date
# ---------------------------------------------------------------------------


class TestParseClinvarDate:
    """Tests for ClinVar date string parsing."""

    def test_standard_format(self):
        result = _parse_clinvar_date("2024-01-15")
        assert result == date(2024, 1, 15)

    def test_slash_format(self):
        result = _parse_clinvar_date("2024/01/15")
        assert result == date(2024, 1, 15)

    def test_missing_dot(self):
        assert _parse_clinvar_date(".") is None

    def test_empty_string(self):
        assert _parse_clinvar_date("") is None

    def test_none_string(self):
        assert _parse_clinvar_date(None) is None  # type: ignore

    def test_invalid_format(self):
        # Should return None without raising
        assert _parse_clinvar_date("not-a-date") is None


# ---------------------------------------------------------------------------
# Tests: _requires_recontact
# ---------------------------------------------------------------------------


class TestRequiresRecontact:
    """Tests for clinical recontact determination logic."""

    def test_vus_to_pathogenic_requires_recontact(self):
        assert _requires_recontact(
            ClinicalSignificance.VUS, ClinicalSignificance.PATHOGENIC
        ) is True

    def test_pathogenic_to_benign_requires_recontact(self):
        assert _requires_recontact(
            ClinicalSignificance.PATHOGENIC, ClinicalSignificance.BENIGN
        ) is True

    def test_likely_pathogenic_to_vus_requires_recontact(self):
        assert _requires_recontact(
            ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.VUS
        ) is True

    def test_benign_to_likely_benign_no_recontact(self):
        """Benign to likely benign is clinically insignificant."""
        assert _requires_recontact(
            ClinicalSignificance.BENIGN, ClinicalSignificance.LIKELY_BENIGN
        ) is False

    def test_vus_to_conflicting(self):
        """VUS to conflicting warrants recontact."""
        assert _requires_recontact(
            ClinicalSignificance.VUS, ClinicalSignificance.CONFLICTING
        ) is True

    def test_pathogenic_to_likely_pathogenic_recontact(self):
        """Both actionable — recontact still applies."""
        assert _requires_recontact(
            ClinicalSignificance.PATHOGENIC, ClinicalSignificance.LIKELY_PATHOGENIC
        ) is True


# ---------------------------------------------------------------------------
# Tests: _iter_vcf_records
# ---------------------------------------------------------------------------


class TestIterVcfRecords:
    """Tests for VCF record iteration and parsing."""

    def test_basic_record_parsing(self):
        vcf_path = _write_vcf([
            _make_vcf_line(
                chrom="chr17", pos=43094692, ref="G", alt="A",
                clnsig="Pathogenic",
                clnrevstat="criteria_provided,_single_submitter",
                alleleid=12345,
                clndate="2024-06-01",
            )
        ])
        records = list(_iter_vcf_records(vcf_path))
        assert len(records) == 1
        r = records[0]
        assert r.chrom == "chr17"
        assert r.pos == 43094692
        assert r.ref == "G"
        assert r.alt == "A"
        assert r.clnsig == ClinicalSignificance.PATHOGENIC
        assert r.review_stars == 1
        assert r.clinvar_date == date(2024, 6, 1)
        assert r.variation_id == 12345

    def test_chromosome_prefix_added(self):
        """Chromosome without 'chr' prefix should have it added."""
        vcf_path = _write_vcf([
            _make_vcf_line(chrom="17", pos=100, clnsig="Benign")
        ])
        records = list(_iter_vcf_records(vcf_path))
        assert records[0].chrom == "chr17"

    def test_gzipped_vcf_parsed(self):
        """Gzipped VCF should be transparently decompressed."""
        vcf_path = _write_vcf(
            [_make_vcf_line(clnsig="Likely_pathogenic")],
            compressed=True,
        )
        records = list(_iter_vcf_records(vcf_path))
        assert len(records) == 1
        assert records[0].clnsig == ClinicalSignificance.LIKELY_PATHOGENIC

    def test_unknown_clnsig_skipped(self):
        """Records with unrecognised CLNSIG should be skipped."""
        vcf_path = _write_vcf([
            _make_vcf_line(clnsig="drug_response")  # Not in CLNSIG_MAP
        ])
        records = list(_iter_vcf_records(vcf_path))
        # drug_response is not in CLNSIG_MAP → clnsig is None → skipped
        # (the iterator only yields records where clnsig is not None)
        assert len(records) == 0

    def test_missing_clnsig_skipped(self):
        """Records without CLNSIG INFO field should be skipped."""
        line = "chr1\t100\tRCV000001\tA\tT\t.\t.\tALLELEID=1\n"
        vcf_path = _write_vcf([line])
        records = list(_iter_vcf_records(vcf_path))
        assert len(records) == 0

    def test_multiple_records(self):
        vcf_path = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic"),
            _make_vcf_line(pos=200, clnsig="Benign"),
            _make_vcf_line(pos=300, clnsig="Uncertain_significance"),
        ])
        records = list(_iter_vcf_records(vcf_path))
        assert len(records) == 3
        assert records[0].clnsig == ClinicalSignificance.PATHOGENIC
        assert records[1].clnsig == ClinicalSignificance.BENIGN
        assert records[2].clnsig == ClinicalSignificance.VUS

    def test_variant_key_format(self):
        vcf_path = _write_vcf([
            _make_vcf_line(chrom="chr7", pos=117548628, ref="CT", alt="C")
        ])
        records = list(_iter_vcf_records(vcf_path))
        assert records[0].key == "chr7:117548628:CT:C"


# ---------------------------------------------------------------------------
# Tests: diff_variants
# ---------------------------------------------------------------------------


class TestDiffVariants:
    """Tests for the ClinVar VCF diff function."""

    def test_detects_vus_to_pathogenic(self):
        """VUS → Pathogenic reclassification should be detected."""
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Uncertain_significance", alleleid=1)
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic", alleleid=1,
                           clndate="2024-12-09")
        ])
        events = diff_variants(old_vcf, new_vcf)
        assert len(events) == 1
        assert events[0].old_class == ClinicalSignificance.VUS.value
        assert events[0].new_class == ClinicalSignificance.PATHOGENIC.value
        assert events[0].recontact_required is True

    def test_no_change_not_detected(self):
        """Variants with the same classification should not produce events."""
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic")
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic")
        ])
        events = diff_variants(old_vcf, new_vcf)
        assert len(events) == 0

    def test_new_variant_not_detected(self):
        """Variants only in new VCF (not in old) should not be reclassifications."""
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic")
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic"),
            _make_vcf_line(pos=200, clnsig="Benign"),  # New variant
        ])
        events = diff_variants(old_vcf, new_vcf)
        assert len(events) == 0

    def test_low_star_rating_filtered(self):
        """Reclassifications with insufficient review stars should be filtered."""
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Uncertain_significance",
                           clnrevstat="no_assertion_provided")
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic",
                           clnrevstat="no_assertion_provided")
        ])
        # min_stars=1 by default; no_assertion_provided = 0 stars → filtered
        events = diff_variants(old_vcf, new_vcf, min_stars=1)
        assert len(events) == 0

    def test_benign_to_likely_benign_no_recontact(self):
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Benign")
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Likely_benign")
        ])
        events = diff_variants(old_vcf, new_vcf)
        assert len(events) == 1
        assert events[0].recontact_required is False

    def test_multiple_reclassifications(self):
        old_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Uncertain_significance", alleleid=1),
            _make_vcf_line(pos=200, clnsig="Pathogenic", alleleid=2),
            _make_vcf_line(pos=300, clnsig="Benign", alleleid=3),
        ])
        new_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic", alleleid=1),
            _make_vcf_line(pos=200, clnsig="Likely_pathogenic", alleleid=2),
            _make_vcf_line(pos=300, clnsig="Benign", alleleid=3),  # Unchanged
        ])
        events = diff_variants(old_vcf, new_vcf)
        assert len(events) == 2


# ---------------------------------------------------------------------------
# Tests: find_reclassified_variants
# ---------------------------------------------------------------------------


class TestFindReclassifiedVariants:
    """Tests for local catalogue vs ClinVar cross-reference."""

    def _make_local_variant(
        self,
        variant_id: str = "v001",
        chrom: str = "chr17",
        pos: int = 43094692,
        ref: str = "G",
        alt: str = "A",
        classification: str = "Uncertain significance",
        report_date: str = "2023-06-01",
    ) -> dict:
        return {
            "variant_id": variant_id,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "classification": classification,
            "report_date": report_date,
        }

    def test_detects_local_reclassification(self):
        """Local VUS reclassified to P in ClinVar should be detected."""
        local = [self._make_local_variant(classification="Uncertain significance")]
        clinvar_vcf = _write_vcf([
            _make_vcf_line(
                chrom="chr17", pos=43094692, ref="G", alt="A",
                clnsig="Pathogenic",
                clnrevstat="criteria_provided,_single_submitter",
            )
        ])
        events = find_reclassified_variants(local, clinvar_vcf)
        assert len(events) == 1
        assert events[0].old_class == "Uncertain significance"
        assert events[0].new_class == "Pathogenic"
        assert events[0].recontact_required is True
        assert events[0].variant_id == "v001"

    def test_matching_classification_not_detected(self):
        """Variant with same classification in local and ClinVar — no event."""
        local = [self._make_local_variant(classification="Pathogenic")]
        clinvar_vcf = _write_vcf([
            _make_vcf_line(clnsig="Pathogenic")
        ])
        events = find_reclassified_variants(local, clinvar_vcf)
        assert len(events) == 0

    def test_variant_not_in_clinvar_skipped(self):
        """Local variants absent from ClinVar should not produce events."""
        local = [
            self._make_local_variant(pos=99999999)  # Not in synthetic ClinVar VCF
        ]
        clinvar_vcf = _write_vcf([
            _make_vcf_line(pos=43094692, clnsig="Pathogenic")
        ])
        events = find_reclassified_variants(local, clinvar_vcf)
        assert len(events) == 0

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            find_reclassified_variants([], Path("/nonexistent/path/clinvar.vcf.gz"))

    def test_multiple_local_variants(self):
        """Multiple local variants reclassified simultaneously."""
        local = [
            self._make_local_variant(
                variant_id="v001", pos=100, classification="Uncertain significance"
            ),
            self._make_local_variant(
                variant_id="v002", pos=200, classification="Pathogenic"
            ),
        ]
        clinvar_vcf = _write_vcf([
            _make_vcf_line(pos=100, clnsig="Pathogenic"),    # v001: VUS→P
            _make_vcf_line(pos=200, clnsig="Likely_benign"), # v002: P→LB
        ])
        events = find_reclassified_variants(local, clinvar_vcf)
        assert len(events) == 2
        assert all(e.recontact_required for e in events)

    def test_unknown_local_classification_skipped(self):
        """Local variants with unrecognised classification are skipped gracefully."""
        local = [
            self._make_local_variant(classification="not-a-real-classification")
        ]
        clinvar_vcf = _write_vcf([
            _make_vcf_line(clnsig="Pathogenic")
        ])
        # Should not raise — should silently skip with a warning
        events = find_reclassified_variants(local, clinvar_vcf)
        assert len(events) == 0
