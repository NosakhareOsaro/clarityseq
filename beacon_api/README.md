# beacon_api/

GA4GH Beacon v2.1.1 API with GA4GH VRS v2.0 variant identifiers and GA4GH Passport authentication.

## Specification versions

- **Beacon**: v2.1.1 (released December 13, 2024) — bug fixes over v2.0; no breaking changes; improved VRS alignment
- **VRS**: v2.0 — canonical variant identifiers (24-char computed digests)
- **GA4GH Passports**: JWT-based access control

## Endpoints

| Endpoint | Description |
|----------|-------------|
| `GET /info` | Beacon metadata and capabilities |
| `GET /g_variants` | Genomic variant queries (VRS identifiers in response) |
| `GET /individuals` | Individual-level queries (Passport-gated) |

## VRS identifiers

All variant responses include GA4GH VRS v2.0 computed identifiers.
Reference: Wagner et al. 2021 Cell Genomics PMID:35072137

## Running locally

```bash
# With Docker Compose
docker compose up beacon

# Or directly
uvicorn beacon_api.main:app --host 0.0.0.0 --port 8080 --reload
```

## Authentication

JWT tokens validated as GA4GH Passports. See `auth/passports.py`.
Register with EGA: `BEACON_EGA_REGISTRATION_URL` in `.env`.
