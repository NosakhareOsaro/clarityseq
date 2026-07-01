"""GA4GH Beacon v2.1.1 API package for ClaritySeq WGS platform.

Implements the GA4GH Beacon v2.1.1 specification (released December 13, 2024)
for federated genomic data discovery across NHS GMS sequencing centres.

Specification references:
    - GA4GH Beacon v2.1.1 (December 13, 2024):
      https://github.com/ga4gh-beacon/beacon-v2
    - GA4GH VRS v2.0: Variant representation in Beacon responses uses
      GA4GH Variation Representation Specification v2.0 identifiers.
    - GA4GH Passports v1.2: Authentication and authorisation via GA4GH
      Passport JWT tokens for controlled access to variant data.

Endpoints implemented:
    GET /info            — Beacon metadata and capabilities.
    POST /g_variants     — Genomic variant queries with VRS identifiers.
    POST /individuals    — Individual-level queries (GA4GH Passport protected).
"""

from beacon_api.main import app

__all__ = ["app"]
__version__ = "2.1.1"
