"""
beacon_api.vrs_utils
====================
GA4GH VRS v2.0 variant identifier utilities.

Wagner et al. 2021 Cell Genomics PMID:35072137
VRS spec: https://vrs.ga4gh.org/

VRS v2.0 identifiers are deterministic, globally unique identifiers for
genomic variants.  They are computed from the canonical allele representation
using a SHA-512 digest with URL-safe base64 encoding.

The canonical serialisation used here follows the VRS v2.0 Allele schema:
    ga4gh:VA.<chrom>:<pos>:<ref>:<alt>

Production implementations should use the official ga4gh.vrs Python SDK:
    pip install ga4gh.vrs

References:
    Wagner et al. 2021 Cell Genomics PMID:35072137 — VRS specification.
    VRS v2.0 spec: https://vrs.ga4gh.org/en/stable/
    ga4gh.vrs Python SDK: https://github.com/ga4gh/vrs-python
"""

from __future__ import annotations

import base64
import hashlib
from dataclasses import dataclass


@dataclass
class VRSAllele:
    """GA4GH VRS v2.0 Allele representation.

    Attributes:
        vrs_id: Computed GA4GH VRS v2.0 identifier (24-char digest prefixed
            with ``"ga4gh:VA."``).
        chrom: Chromosome (GRCh38 notation, e.g. ``"chr17"``).
        pos: 1-based genomic position (GRCh38).
        ref: Reference allele (VCF notation, uppercase).
        alt: Alternate allele (VCF notation, uppercase).
        digest: Raw 24-character base64url digest.
    """

    vrs_id: str      # e.g. "ga4gh:VA.AbCdEfGhIjKlMnOpQrStUvWx"
    chrom: str
    pos: int
    ref: str
    alt: str
    digest: str      # 24-char base64url


def compute_vrs_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Compute GA4GH VRS v2.0 computed identifier (24-char digest).

    VRS identifiers are deterministic from the canonical allele tuple.
    The identifier uses SHA-512 with URL-safe base64 encoding per VRS v2.0 spec.

    This is a simplified implementation for demonstration purposes.
    Production code should use the ga4gh.vrs Python SDK for full compliance
    with the VRS Digest algorithm (VRS Digest Algorithm Spec v2.0).

    Args:
        chrom: Chromosome in GRCh38 notation (e.g. ``"chr17"``).
        pos: 1-based genomic position.
        ref: Reference allele (uppercase, VCF notation).
        alt: Alternate allele (uppercase, VCF notation).

    Returns:
        Full GA4GH VRS v2.0 identifier string in the format
        ``"ga4gh:VA.<24-char-digest>"``.

    Examples:
        >>> compute_vrs_id("chr17", 43044295, "G", "A")
        'ga4gh:VA.XXXXXXXXXXXXXXXXXXXXXXXX'

    References:
        Wagner et al. 2021 Cell Genomics PMID:35072137 — VRS spec.
        VRS Digest Algorithm: https://vrs.ga4gh.org/en/stable/impl-guide/computed_identifiers.html
    """
    # VRS canonical serialisation
    # Full VRS v2.0 uses the VRS Digest algorithm with JSON canonicalisation
    # and ga4gh-digest-multibase encoding.  Here we use a simplified form.
    canonical = f"ga4gh:VA.{chrom}:{pos}:{ref}:{alt}"
    digest_bytes = hashlib.sha512(canonical.encode("utf-8")).digest()[:18]
    # URL-safe base64 encoding, no padding, truncated to 24 chars
    digest_str = base64.urlsafe_b64encode(digest_bytes).decode("ascii").rstrip("=")[:24]
    return f"ga4gh:VA.{digest_str}"


def make_vrs_allele(chrom: str, pos: int, ref: str, alt: str) -> VRSAllele:
    """Create a VRSAllele object with a computed VRS v2.0 identifier.

    Args:
        chrom: Chromosome in GRCh38 notation (e.g. ``"chr17"``).
        pos: 1-based genomic position.
        ref: Reference allele (uppercase, VCF notation).
        alt: Alternate allele (uppercase, VCF notation).

    Returns:
        VRSAllele with a computed GA4GH VRS v2.0 identifier.

    References:
        Wagner et al. 2021 Cell Genomics PMID:35072137.
    """
    full_id = compute_vrs_id(chrom, pos, ref, alt)
    digest = full_id.replace("ga4gh:VA.", "")
    return VRSAllele(
        vrs_id=full_id,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        digest=digest,
    )


def vrs_allele_to_dict(allele: VRSAllele) -> dict[str, object]:
    """Serialise a VRSAllele to a GA4GH VRS v2.0 JSON-compatible dict.

    Args:
        allele: VRSAllele to serialise.

    Returns:
        Dict with ``"id"``, ``"type"``, ``"location"``, and ``"state"``
        keys per VRS v2.0 Allele schema.

    References:
        VRS v2.0 Allele schema: https://vrs.ga4gh.org/en/stable/terms_and_model.html#allele
    """
    return {
        "id": allele.vrs_id,
        "type": "Allele",
        "digest": allele.digest,
        "location": {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": f"SQ.{allele.chrom}",  # simplified; real uses GA4GH refget
            },
            "start": allele.pos - 1,   # VRS uses 0-based interbase coordinates
            "end": allele.pos - 1 + len(allele.ref),
        },
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": allele.alt,
        },
    }
