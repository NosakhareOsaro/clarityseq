"""SQLAlchemy ORM models for the variant reclassification daemon.

This module defines the database schema for tracking ClinVar reclassification
events, NHS-mandated ClinVar submissions, patient-variant linkage, and
VUS review scheduling.

Regulatory references:
    - ACGS 2024 Best Practice Guidelines §9: "Laboratories should have a
      process for the monitoring of reclassification of variants in ClinVar
      and should recontact patients where clinically appropriate."
    - ACGS 2024 Introduction: NHS Genomic Medicine Service mandates ClinVar
      submission for all P/LP variants identified in diagnostic WGS testing.
    - ACGS 2024 §9: VUS review should occur every 2 years to incorporate
      new evidence, including gnomAD v4.1 population frequency updates.
"""

from __future__ import annotations

import enum
from datetime import date, datetime
from typing import Optional

from sqlalchemy import (
    Date,
    DateTime,
    Enum,
    ForeignKey,
    Index,
    String,
    Text,
    UniqueConstraint,
    func,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    """Declarative base for all reclassification ORM models."""

    pass


class ClinicalSignificance(str, enum.Enum):
    """ClinVar clinical significance classifications.

    Follows ClinVar controlled vocabulary; maps to ACMG/AMP 2015 five-tier
    classification scheme as adopted by ACGS 2024.
    """

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"
    VUS = "Uncertain significance"
    LIKELY_BENIGN = "Likely benign"
    BENIGN = "Benign"
    CONFLICTING = "Conflicting interpretations"
    NOT_PROVIDED = "Not provided"
    RISK_FACTOR = "risk factor"


class SubmissionStatus(str, enum.Enum):
    """ClinVar submission lifecycle states."""

    PENDING = "pending"
    SUBMITTED = "submitted"
    PROCESSING = "processing"
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    ERROR = "error"


class ReclassificationEvent(Base):
    """Records a single variant reclassification detected from ClinVar diff.

    Populated by the weekly ClinVar FTP diff process (clinvar_diff.py).
    Each row represents a transition of a variant from one clinical
    significance tier to another, as detected by comparing consecutive
    ClinVar VCF releases.

    Attributes:
        id: Auto-incremented primary key.
        variant_id: Internal GenomeForge variant identifier (FK to variant
            catalogue). Indexed for fast patient-variant lookups.
        old_class: Previous ClinVar clinical significance.
        new_class: Updated ClinVar clinical significance.
        clinvar_accession: ClinVar accession number (e.g. RCV000123456).
        clinvar_variation_id: ClinVar variation ID (integer, e.g. 12345).
        clinvar_date: Date the reclassification was published in ClinVar.
        detected_at: Timestamp when the GenomeForge daemon first detected
            this reclassification.
        fhir_task_id: FHIR Task resource ID created for recontact workflow,
            populated after Task creation in fhir_task.py.
        recontact_required: Whether clinical recontact is required per
            ACGS 2024 §9 criteria (P/LP↔VUS or P/LP↔benign transitions).
        notes: Free-text notes from clinical scientist review.
        patient_variants: Related PatientVariant records for affected patients.
    """

    __tablename__ = "reclassification_events"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    variant_id: Mapped[str] = mapped_column(String(64), nullable=False, index=True)
    old_class: Mapped[str] = mapped_column(
        Enum(ClinicalSignificance), nullable=False
    )
    new_class: Mapped[str] = mapped_column(
        Enum(ClinicalSignificance), nullable=False
    )
    clinvar_accession: Mapped[Optional[str]] = mapped_column(String(32), nullable=True)
    clinvar_variation_id: Mapped[Optional[int]] = mapped_column(nullable=True)
    clinvar_date: Mapped[date] = mapped_column(Date, nullable=False)
    detected_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        nullable=False,
        server_default=func.now(),
    )
    fhir_task_id: Mapped[Optional[str]] = mapped_column(
        String(128), nullable=True, index=True
    )
    recontact_required: Mapped[bool] = mapped_column(
        nullable=False, default=False
    )
    notes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    # Relationships
    patient_variants: Mapped[list["PatientVariant"]] = relationship(
        "PatientVariant",
        back_populates="reclassification_event",
        lazy="select",
    )

    __table_args__ = (
        # Ensure we don't double-record the same reclassification
        UniqueConstraint(
            "variant_id", "clinvar_date", "old_class", "new_class",
            name="uq_reclassification_event"
        ),
        Index("ix_reclassification_detected_at", "detected_at"),
    )

    def __repr__(self) -> str:
        return (
            f"<ReclassificationEvent variant={self.variant_id!r} "
            f"{self.old_class} -> {self.new_class} on {self.clinvar_date}>"
        )


class ClinVarSubmissionQueue(Base):
    """Queue of variants awaiting NHS-mandated ClinVar submission.

    Per the ACGS 2024 Introduction: "The NHS Genomic Medicine Service (GMS)
    requires that all laboratories participating in the NHS WGS programme
    submit clinically interpreted variants to ClinVar. Submission of
    pathogenic and likely pathogenic variants is mandatory within 3 months
    of clinical report issue."

    This model tracks the full submission lifecycle, including MANE Select
    HGVSc nomenclature (required by NCBI for ClinVar submissions), the
    BayesACMG posterior probability supporting the classification, and
    NCBI API responses.

    Attributes:
        id: Auto-incremented primary key.
        variant_id: Internal GenomeForge variant identifier.
        gene_symbol: HGNC gene symbol (e.g. BRCA1).
        chromosome: Chromosome (e.g. '17', 'X').
        position_grch38: GRCh38 genomic position (1-based).
        ref_allele: Reference allele sequence.
        alt_allele: Alternate allele sequence.
        mane_select_hgvsc: MANE Select transcript HGVSc notation, required
            for ClinVar submission. MANE Select transcripts are the single
            representative transcript per gene agreed between Ensembl and
            RefSeq (Morales et al. 2022, Nat Methods).
        clinical_significance: ACMG/AMP classification tier.
        condition_name: MedGen/OMIM condition name.
        condition_id: MedGen CUI or OMIM ID for structured condition lookup.
        bayesacmg_probability: Posterior probability from BayesACMG model;
            provides quantitative support for the clinical significance tier
            per Tavtigian et al. 2020 (PMID:31479589).
        evidence_codes: JSON array of ACMG/AMP criteria codes applied
            (e.g. ["PVS1", "PS1", "PM2"]).
        submission_status: Current lifecycle status.
        submitted_at: Timestamp of NCBI API submission.
        ncbi_response: Raw NCBI API response JSON for audit trail.
        ncbi_submission_id: NCBI-assigned submission batch ID.
        error_message: Error details if submission failed.
        created_at: Record creation timestamp.
        updated_at: Last modification timestamp.
    """

    __tablename__ = "clinvar_submission_queue"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    variant_id: Mapped[str] = mapped_column(String(64), nullable=False, index=True)

    # Variant description fields
    gene_symbol: Mapped[str] = mapped_column(String(32), nullable=False)
    chromosome: Mapped[str] = mapped_column(String(4), nullable=False)
    position_grch38: Mapped[int] = mapped_column(nullable=False)
    ref_allele: Mapped[str] = mapped_column(String(1000), nullable=False)
    alt_allele: Mapped[str] = mapped_column(String(1000), nullable=False)

    # MANE Select HGVSc — mandatory for NCBI ClinVar submissions
    mane_select_hgvsc: Mapped[Optional[str]] = mapped_column(
        String(512), nullable=True,
        comment=(
            "MANE Select transcript HGVSc notation required by NCBI. "
            "Must use RefSeq NM_ accession matching MANE Select v1.3+."
        )
    )

    # Classification fields
    clinical_significance: Mapped[str] = mapped_column(
        Enum(ClinicalSignificance), nullable=False
    )
    condition_name: Mapped[str] = mapped_column(String(512), nullable=False)
    condition_id: Mapped[Optional[str]] = mapped_column(String(32), nullable=True)

    # BayesACMG posterior probability (0.0–1.0)
    bayesacmg_probability: Mapped[Optional[float]] = mapped_column(nullable=True)

    # JSON-encoded ACMG criteria codes list
    evidence_codes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    # Submission lifecycle
    submission_status: Mapped[str] = mapped_column(
        Enum(SubmissionStatus),
        nullable=False,
        default=SubmissionStatus.PENDING,
    )
    submitted_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    ncbi_response: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    ncbi_submission_id: Mapped[Optional[str]] = mapped_column(
        String(128), nullable=True
    )
    error_message: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    # Audit timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        nullable=False,
        server_default=func.now(),
        onupdate=func.now(),
    )

    __table_args__ = (
        Index("ix_submission_status", "submission_status"),
        Index("ix_submission_created_at", "created_at"),
    )

    def __repr__(self) -> str:
        return (
            f"<ClinVarSubmissionQueue variant={self.variant_id!r} "
            f"status={self.submission_status} gene={self.gene_symbol}>"
        )


class PatientVariant(Base):
    """Links patients to variant classifications for reclassification tracking.

    This junction table enables efficient identification of all patients
    affected by a reclassification event. Patient identifiers are stored
    as pseudonymous GMS IDs to comply with UK GDPR Article 89 research
    pseudonymisation requirements.

    Attributes:
        id: Auto-incremented primary key.
        patient_gms_id: NHS GMS pseudonymous patient identifier.
        variant_id: Internal GenomeForge variant identifier.
        reclassification_event_id: FK to ReclassificationEvent if this
            patient-variant pair has been affected by a reclassification.
        current_classification: Classification at time of patient report.
        report_date: Date of original diagnostic report.
        recontact_sent_at: Timestamp when recontact letter was generated.
        recontact_fhir_task_id: FHIR Task ID for tracking recontact workflow.
        lab_sample_id: Laboratory sample identifier for traceability.
        consent_research: Whether patient consented to research use of data.
    """

    __tablename__ = "patient_variants"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)

    # Pseudonymous NHS GMS patient identifier (not NHS number)
    patient_gms_id: Mapped[str] = mapped_column(
        String(64), nullable=False, index=True
    )
    variant_id: Mapped[str] = mapped_column(String(64), nullable=False, index=True)

    # FK to reclassification event (nullable — not all patient-variants
    # will have been reclassified)
    reclassification_event_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("reclassification_events.id", ondelete="SET NULL"),
        nullable=True,
        index=True,
    )
    reclassification_event: Mapped[Optional["ReclassificationEvent"]] = relationship(
        "ReclassificationEvent",
        back_populates="patient_variants",
    )

    # Classification and report details
    current_classification: Mapped[str] = mapped_column(
        Enum(ClinicalSignificance), nullable=False
    )
    report_date: Mapped[date] = mapped_column(Date, nullable=False)
    lab_sample_id: Mapped[Optional[str]] = mapped_column(String(64), nullable=True)
    consent_research: Mapped[bool] = mapped_column(nullable=False, default=False)

    # Recontact workflow tracking
    recontact_sent_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    recontact_fhir_task_id: Mapped[Optional[str]] = mapped_column(
        String(128), nullable=True
    )

    __table_args__ = (
        UniqueConstraint(
            "patient_gms_id", "variant_id",
            name="uq_patient_variant"
        ),
    )

    def __repr__(self) -> str:
        return (
            f"<PatientVariant patient={self.patient_gms_id!r} "
            f"variant={self.variant_id!r} class={self.current_classification}>"
        )


class VUSReviewSchedule(Base):
    """Schedules periodic re-review of variants of uncertain significance.

    Per ACGS 2024 §9: "Variants of uncertain significance (VUS) should be
    re-evaluated every 2 years, or sooner if significant new evidence
    becomes available (e.g. gnomAD v4.1 population frequency data or new
    functional studies)."

    Attributes:
        id: Auto-incremented primary key.
        variant_id: Internal GenomeForge variant identifier.
        patient_gms_id: Pseudonymous patient identifier (variant may be
            present in multiple patients; one schedule row per patient-variant
            pair to allow individual tracking).
        initial_classification_date: Date of original VUS classification.
        review_due_date: Date by which re-review must be completed.
            Set to initial_classification_date + 2 years per ACGS 2024 §9.
        review_completed_at: Timestamp when re-review was completed.
        post_review_classification: Classification following review,
            if changed from VUS.
        review_notes: Clinical scientist notes from review.
        reminder_sent_at: Timestamp when reminder notification was sent.
        is_overdue: Computed property; True if review_due_date < today
            and review_completed_at is None.
    """

    __tablename__ = "vus_review_schedule"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    variant_id: Mapped[str] = mapped_column(String(64), nullable=False, index=True)
    patient_gms_id: Mapped[str] = mapped_column(String(64), nullable=False, index=True)

    initial_classification_date: Mapped[date] = mapped_column(Date, nullable=False)

    # review_due_date = initial_classification_date + 2 years (ACGS 2024 §9)
    review_due_date: Mapped[date] = mapped_column(Date, nullable=False, index=True)

    review_completed_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    post_review_classification: Mapped[Optional[str]] = mapped_column(
        Enum(ClinicalSignificance), nullable=True
    )
    review_notes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    reminder_sent_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )

    __table_args__ = (
        UniqueConstraint(
            "variant_id", "patient_gms_id",
            name="uq_vus_review_schedule"
        ),
        Index("ix_vus_review_due_date", "review_due_date"),
    )

    @property
    def is_overdue(self) -> bool:
        """Return True if VUS review is past due and not yet completed."""
        if self.review_completed_at is not None:
            return False
        return date.today() > self.review_due_date

    def __repr__(self) -> str:
        return (
            f"<VUSReviewSchedule variant={self.variant_id!r} "
            f"patient={self.patient_gms_id!r} due={self.review_due_date}>"
        )
