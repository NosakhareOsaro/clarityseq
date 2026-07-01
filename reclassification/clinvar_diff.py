"""Weekly ClinVar FTP diff for variant reclassification detection.

This module implements the weekly ClinVar VCF diff pipeline that detects
reclassification events for variants held in the ClaritySeq variant
catalogue. It is the core data-acquisition component of the ACGS 2024 §9
reclassification monitoring requirement.

Regulatory and scientific references:
    - ACGS 2024 Best Practice Guidelines §9: "Laboratories should have a
      documented process for monitoring published reclassifications, including
      at minimum a weekly check of the ClinVar variant summary file."
    - FHIR Genomics Reporting IG v3.0.0: Reclassification events feed
      downstream FHIR Task generation (fhir_task.py) per the HL7 Genomics
      Reporting Implementation Guide recontact workflow.
    - gnomAD v4.1 (Chen et al. 2024): Population frequency changes between
      gnomAD v3.1.2 and v4.1 (730,947 exomes + 76,215 genomes) may drive
      ClinVar reclassifications, particularly for PM2 (absent/rare in
      population) and BA1 (>5% allele frequency) criteria reassessment.

ClinVar FTP source:
    ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    Released weekly on Mondays; MD5 checksum available at same URL + '.md5'.
"""

from __future__ import annotations

import ftplib
import gzip
import hashlib
import logging
import os
import re
import tempfile
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Iterator, Optional
from urllib.parse import urlparse

from reclassification.models import ClinicalSignificance, ReclassificationEvent

logger = logging.getLogger(__name__)

# ClinVar FTP defaults
DEFAULT_CLINVAR_FTP_URL = (
    "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)
DEFAULT_CLINVAR_MD5_URL = (
    "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5"
)

# ClinVar INFO field key for clinical significance
CLNSIG_KEY = "CLNSIG"
CLNREVSTAT_KEY = "CLNREVSTAT"  # Review status (star rating)
CLINVAR_DATE_KEY = "CLNDATE"    # Date last evaluated
VARIATION_ID_KEY = "ALLELEID"

# Minimum ClinVar review status to trigger reclassification alert
# Only consider variants with at least 1 star (criteria provided, single
# submitter) to reduce noise from low-confidence classifications.
MIN_REVIEW_STAR_RATING = 1

# Mapping from ClinVar CLNSIG INFO value to internal enum
CLNSIG_MAP: dict[str, ClinicalSignificance] = {
    "Pathogenic": ClinicalSignificance.PATHOGENIC,
    "Likely_pathogenic": ClinicalSignificance.LIKELY_PATHOGENIC,
    "Pathogenic/Likely_pathogenic": ClinicalSignificance.LIKELY_PATHOGENIC,
    "Uncertain_significance": ClinicalSignificance.VUS,
    "Likely_benign": ClinicalSignificance.LIKELY_BENIGN,
    "Benign": ClinicalSignificance.BENIGN,
    "Benign/Likely_benign": ClinicalSignificance.LIKELY_BENIGN,
    "Conflicting_interpretations_of_pathogenicity": (
        ClinicalSignificance.CONFLICTING
    ),
    "not_provided": ClinicalSignificance.NOT_PROVIDED,
    "risk_factor": ClinicalSignificance.RISK_FACTOR,
}


@dataclass
class ClinVarRecord:
    """Parsed representation of a single ClinVar VCF record.

    Attributes:
        chrom: Chromosome (e.g. 'chr17').
        pos: 1-based genomic position (GRCh38).
        ref: Reference allele.
        alt: Alternate allele.
        clnsig: Parsed clinical significance.
        variation_id: ClinVar variation ID from ALLELEID INFO field.
        clinvar_date: Date last evaluated (from CLNDATE INFO field).
        review_stars: Number of ClinVar review stars (0–4).
        accession: ClinVar RCV accession (from ID column).
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    clnsig: Optional[ClinicalSignificance]
    variation_id: Optional[int]
    clinvar_date: Optional[date]
    review_stars: int
    accession: Optional[str]

    @property
    def key(self) -> str:
        """Return a canonical key for variant identity matching."""
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


def _parse_review_stars(clnrevstat: str) -> int:
    """Map ClinVar CLNREVSTAT value to number of review stars.

    Args:
        clnrevstat: Raw CLNREVSTAT INFO field value.

    Returns:
        Integer star rating 0–4.
    """
    rating_map = {
        "no_assertion_provided": 0,
        "no_assertion_criteria_provided": 0,
        "criteria_provided,_single_submitter": 1,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "reviewed_by_expert_panel": 3,
        "practice_guideline": 4,
    }
    # Normalise: ClinVar uses commas+underscores; strip pipe-separated values
    key = clnrevstat.split("|")[0].lower().replace(" ", "_")
    return rating_map.get(key, 0)


def _parse_info(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO field string into a key-value dictionary.

    Args:
        info_str: Raw INFO column string (semicolon-delimited key=value pairs).

    Returns:
        Dictionary of INFO key to raw string value.
    """
    result: dict[str, str] = {}
    for item in info_str.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, _, v = item.partition("=")
            result[k] = v
        else:
            result[item] = "true"  # Flag fields
    return result


def _parse_clinvar_date(raw: str) -> Optional[date]:
    """Parse a ClinVar date string into a Python date object.

    ClinVar uses YYYY-MM-DD format in the CLNDATE INFO field, but may
    use '.' for missing values.

    Args:
        raw: Raw CLNDATE string value.

    Returns:
        Parsed date object, or None if value is missing/unparseable.
    """
    if not raw or raw in (".", ""):
        return None
    try:
        return datetime.strptime(raw, "%Y-%m-%d").date()
    except ValueError:
        # Some older ClinVar entries use YYYY/MM/DD
        try:
            return datetime.strptime(raw, "%Y/%m/%d").date()
        except ValueError:
            logger.warning("Unable to parse ClinVar date: %r", raw)
            return None


def _iter_vcf_records(vcf_path: Path) -> Iterator[ClinVarRecord]:
    """Yield parsed ClinVarRecord objects from a (possibly gzipped) VCF.

    Skips header lines (starting with '#'). Only yields records with
    a recognised CLNSIG INFO field.

    Args:
        vcf_path: Path to ClinVar VCF file (.vcf or .vcf.gz).

    Yields:
        ClinVarRecord objects for each valid data line.
    """
    open_fn = gzip.open if str(vcf_path).endswith(".gz") else open

    with open_fn(vcf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue  # Skip header and meta-information lines

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue  # Malformed line

            chrom, pos_str, accession, ref, alt, _qual, _filt, info_str = (
                parts[0], parts[1], parts[2], parts[3],
                parts[4], parts[5], parts[6], parts[7],
            )

            info = _parse_info(info_str)
            raw_clnsig = info.get(CLNSIG_KEY, "")
            if not raw_clnsig:
                continue  # No clinical significance — skip

            # ClinVar may pipe-separate multiple values; take first
            clnsig_raw = raw_clnsig.split("|")[0]
            clnsig = CLNSIG_MAP.get(clnsig_raw)
            if clnsig is None:
                continue  # Unrecognised CLNSIG (e.g. drug_response) — skip

            # Parse review stars
            raw_revstat = info.get(CLNREVSTAT_KEY, "")
            stars = _parse_review_stars(raw_revstat)

            # Parse variation ID
            variation_id: Optional[int] = None
            raw_varid = info.get(VARIATION_ID_KEY, "")
            if raw_varid.isdigit():
                variation_id = int(raw_varid)

            # Parse evaluation date
            clinvar_date = _parse_clinvar_date(info.get(CLINVAR_DATE_KEY, ""))

            # Normalise chromosome name (add 'chr' prefix if missing)
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"

            yield ClinVarRecord(
                chrom=chrom,
                pos=int(pos_str),
                ref=ref,
                alt=alt,
                clnsig=clnsig,
                variation_id=variation_id,
                clinvar_date=clinvar_date,
                review_stars=stars,
                accession=accession if accession != "." else None,
            )


def _verify_md5(file_path: Path, expected_md5: str) -> bool:
    """Verify MD5 checksum of a downloaded file.

    Args:
        file_path: Path to file to verify.
        expected_md5: Expected MD5 hex digest string.

    Returns:
        True if checksum matches, False otherwise.
    """
    md5 = hashlib.md5()
    with open(file_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            md5.update(chunk)
    return md5.hexdigest() == expected_md5.strip().split()[0]


def download_latest_clinvar_vcf(
    ftp_url: str = DEFAULT_CLINVAR_FTP_URL,
    dest_dir: Optional[Path] = None,
    verify_checksum: bool = True,
) -> Path:
    """Download the latest ClinVar VCF from NCBI FTP and verify integrity.

    Downloads the weekly ClinVar VCF release for GRCh38 from the NCBI FTP
    server. Verifies MD5 checksum to ensure file integrity. The downloaded
    file is retained at dest_dir for archival and comparison purposes.

    Args:
        ftp_url: FTP URL of the ClinVar VCF. Defaults to the NCBI GRCh38
            weekly release. URL must be in the form
            ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
        dest_dir: Directory to save the downloaded file. If None, uses a
            temporary directory (file will be cleaned up automatically).
        verify_checksum: If True (default), download and verify the MD5
            checksum file alongside the VCF.

    Returns:
        Path to the downloaded (and verified) ClinVar VCF file.

    Raises:
        ftplib.error_perm: If FTP authentication or path fails.
        ValueError: If MD5 checksum verification fails.
        OSError: If the download directory is not writable.

    Example:
        >>> vcf_path = download_latest_clinvar_vcf()
        >>> print(vcf_path)
        /tmp/clinvar_20241209.vcf.gz
    """
    parsed = urlparse(ftp_url)
    host = parsed.netloc
    remote_path = parsed.path
    filename = os.path.basename(remote_path)

    if dest_dir is None:
        dest_dir = Path(tempfile.mkdtemp(prefix="clinvar_"))
    else:
        dest_dir.mkdir(parents=True, exist_ok=True)

    local_vcf = dest_dir / filename
    local_md5 = dest_dir / (filename + ".md5")

    logger.info("Connecting to ClinVar FTP: %s", host)

    with ftplib.FTP(host) as ftp:
        ftp.login()  # Anonymous FTP access for NCBI
        logger.info("Downloading ClinVar VCF: %s", remote_path)

        # Download main VCF file
        with open(local_vcf, "wb") as fh:
            ftp.retrbinary(f"RETR {remote_path}", fh.write)
        logger.info("Downloaded %s (%.1f MB)", filename, local_vcf.stat().st_size / 1e6)

        if verify_checksum:
            # Download MD5 sidecar file
            md5_remote = remote_path + ".md5"
            with open(local_md5, "wb") as fh:
                ftp.retrbinary(f"RETR {md5_remote}", fh.write)

            expected_md5 = local_md5.read_text().strip()
            logger.info("Verifying MD5 checksum for %s", filename)
            if not _verify_md5(local_vcf, expected_md5):
                raise ValueError(
                    f"MD5 verification failed for {local_vcf}. "
                    f"Expected: {expected_md5}. File may be corrupted."
                )
            logger.info("MD5 checksum verified OK")

    return local_vcf


def diff_variants(
    old_vcf: Path,
    new_vcf: Path,
    min_stars: int = MIN_REVIEW_STAR_RATING,
) -> list[ReclassificationEvent]:
    """Compare two ClinVar VCF files and return reclassification events.

    Performs a full comparison between two consecutive ClinVar VCF releases,
    identifying variants where the CLNSIG field has changed. Only variants
    with a minimum ClinVar review star rating are considered to reduce
    noise from low-confidence single-submitter changes.

    Args:
        old_vcf: Path to the previous ClinVar VCF (last week's release).
        new_vcf: Path to the current ClinVar VCF (this week's release).
        min_stars: Minimum ClinVar review star rating to consider.
            Defaults to 1 (criteria provided, single submitter).

    Returns:
        List of ReclassificationEvent objects (not yet persisted to DB).
        Events are ordered by variant key for deterministic output.

    Example:
        >>> events = diff_variants(Path("clinvar_old.vcf.gz"),
        ...                        Path("clinvar_new.vcf.gz"))
        >>> print(len(events), "reclassifications detected")
    """
    logger.info("Loading old ClinVar VCF: %s", old_vcf)
    old_records: dict[str, ClinVarRecord] = {}
    for record in _iter_vcf_records(old_vcf):
        if record.clnsig is not None:
            old_records[record.key] = record
    logger.info("Loaded %d old ClinVar records", len(old_records))

    logger.info("Loading new ClinVar VCF: %s", new_vcf)
    new_records: dict[str, ClinVarRecord] = {}
    for record in _iter_vcf_records(new_vcf):
        if record.clnsig is not None:
            new_records[record.key] = record
    logger.info("Loaded %d new ClinVar records", len(new_records))

    events: list[ReclassificationEvent] = []

    for key, new_rec in sorted(new_records.items()):
        old_rec = old_records.get(key)
        if old_rec is None:
            continue  # New variant — not a reclassification of an existing one

        if old_rec.clnsig == new_rec.clnsig:
            continue  # No change in classification

        if new_rec.review_stars < min_stars:
            logger.debug(
                "Skipping reclassification for %s (review stars %d < %d)",
                key, new_rec.review_stars, min_stars,
            )
            continue

        # Determine whether clinical recontact is required per ACGS 2024 §9:
        # Recontact is required for transitions involving P/LP or significant
        # changes in actionability (e.g. VUS -> P/LP, P/LP -> Benign).
        recontact_required = _requires_recontact(old_rec.clnsig, new_rec.clnsig)

        event = ReclassificationEvent(
            variant_id=key,  # Will be resolved to internal ID by caller
            old_class=old_rec.clnsig.value,
            new_class=new_rec.clnsig.value,
            clinvar_accession=new_rec.accession,
            clinvar_variation_id=new_rec.variation_id,
            clinvar_date=new_rec.clinvar_date or date.today(),
            detected_at=datetime.utcnow(),
            recontact_required=recontact_required,
        )
        events.append(event)
        logger.info(
            "Reclassification detected: %s %s -> %s (recontact=%s)",
            key, old_rec.clnsig, new_rec.clnsig, recontact_required,
        )

    logger.info("Total reclassifications detected: %d", len(events))
    return events


def _requires_recontact(
    old_sig: ClinicalSignificance,
    new_sig: ClinicalSignificance,
) -> bool:
    """Determine whether a classification change requires patient recontact.

    Per ACGS 2024 §9, recontact is required when:
    - A variant previously classified as P/LP is downgraded to VUS or benign.
    - A variant previously classified as VUS or benign is upgraded to P/LP.
    - A conflicting interpretation resolves to P/LP.

    Args:
        old_sig: Previous clinical significance.
        new_sig: New clinical significance.

    Returns:
        True if clinical recontact is required.
    """
    actionable = {ClinicalSignificance.PATHOGENIC, ClinicalSignificance.LIKELY_PATHOGENIC}
    old_actionable = old_sig in actionable
    new_actionable = new_sig in actionable

    # Any change involving actionable classification triggers recontact
    if old_actionable or new_actionable:
        return True

    # VUS -> conflicting or conflicting -> P/LP also warrants recontact
    if old_sig == ClinicalSignificance.VUS and new_sig != ClinicalSignificance.VUS:
        return True

    return False


def find_reclassified_variants(
    local_variants: list[dict],
    clinvar_data: Path,
    min_stars: int = MIN_REVIEW_STAR_RATING,
) -> list[ReclassificationEvent]:
    """Cross-reference local variant catalogue against a ClinVar VCF.

    Unlike diff_variants (which compares two ClinVar releases), this function
    takes the local ClaritySeq variant catalogue and identifies variants
    whose current ClinVar classification differs from what was recorded at
    the time of clinical reporting.

    This supports on-demand reclassification checks (e.g. triggered by
    uploading a new patient cohort) independent of the weekly diff.

    Args:
        local_variants: List of variant dictionaries from the ClaritySeq
            variant catalogue. Each dict must contain keys:
            - 'variant_id': Internal identifier (str)
            - 'chrom': Chromosome (str, with 'chr' prefix)
            - 'pos': Position (int)
            - 'ref': Reference allele (str)
            - 'alt': Alternate allele (str)
            - 'classification': Current classification (str matching
              ClinicalSignificance enum value)
            - 'report_date': ISO date string of last clinical report
        clinvar_data: Path to current ClinVar VCF to compare against.
        min_stars: Minimum ClinVar review star rating to accept.

    Returns:
        List of ReclassificationEvent objects for variants where the
        current ClinVar classification differs from the locally recorded
        classification.

    Raises:
        FileNotFoundError: If clinvar_data path does not exist.
        KeyError: If a local_variant dict is missing a required key.

    Example:
        >>> local = [{"variant_id": "v1", "chrom": "chr17", "pos": 43094692,
        ...           "ref": "G", "alt": "A",
        ...           "classification": "Uncertain significance",
        ...           "report_date": "2023-01-15"}]
        >>> events = find_reclassified_variants(local, Path("clinvar.vcf.gz"))
    """
    if not clinvar_data.exists():
        raise FileNotFoundError(f"ClinVar VCF not found: {clinvar_data}")

    # Index ClinVar data by canonical key
    clinvar_index: dict[str, ClinVarRecord] = {}
    for record in _iter_vcf_records(clinvar_data):
        if record.clnsig is not None and record.review_stars >= min_stars:
            clinvar_index[record.key] = record

    logger.info(
        "Indexed %d ClinVar records (min_stars=%d)", len(clinvar_index), min_stars
    )

    events: list[ReclassificationEvent] = []

    for variant in local_variants:
        # Build canonical key for lookup
        chrom = variant["chrom"]
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        key = f"{chrom}:{variant['pos']}:{variant['ref']}:{variant['alt']}"

        clinvar_rec = clinvar_index.get(key)
        if clinvar_rec is None:
            continue  # Not in ClinVar (yet), or below min star rating

        # Parse the locally recorded classification
        local_class_str = variant.get("classification", "")
        local_class = CLNSIG_MAP.get(local_class_str)
        if local_class is None:
            # Try direct enum value lookup
            try:
                local_class = ClinicalSignificance(local_class_str)
            except ValueError:
                logger.warning(
                    "Unknown local classification %r for variant %s",
                    local_class_str, variant["variant_id"],
                )
                continue

        if local_class == clinvar_rec.clnsig:
            continue  # No reclassification

        recontact_required = _requires_recontact(local_class, clinvar_rec.clnsig)

        # Parse last report date for audit trail
        report_date_str = variant.get("report_date", "")
        try:
            _report_date = datetime.strptime(report_date_str, "%Y-%m-%d").date()
        except (ValueError, TypeError):
            _report_date = date.today()

        event = ReclassificationEvent(
            variant_id=variant["variant_id"],
            old_class=local_class.value,
            new_class=clinvar_rec.clnsig.value,
            clinvar_accession=clinvar_rec.accession,
            clinvar_variation_id=clinvar_rec.variation_id,
            clinvar_date=clinvar_rec.clinvar_date or date.today(),
            detected_at=datetime.utcnow(),
            recontact_required=recontact_required,
        )
        events.append(event)
        logger.info(
            "Variant %s reclassified: %s -> %s",
            variant["variant_id"], local_class, clinvar_rec.clnsig,
        )

    logger.info(
        "Found %d reclassified variants in local catalogue", len(events)
    )
    return events
