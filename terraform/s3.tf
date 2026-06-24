# GenomeForge — S3 Buckets (AWS provider v5)
#
# IMPORTANT: AWS provider v5 breaking change —
#   aws_s3_bucket_acl REMOVED. Use aws_s3_bucket_ownership_controls instead.
#   Set object_ownership = "BucketOwnerEnforced" to disable ACLs entirely.
#   All access controlled via bucket policies and IAM only.
#
# Buckets:
#   1. genomeforge-results    — pipeline output VCFs, reports, BAMs (clinical data)
#   2. genomeforge-reference  — DRAGMAP index, gnomAD v4.1, VEP cache v111
#   3. genomeforge-clinvar    — Weekly ClinVar diff cache (Celery daemon)
#
# All buckets: encryption at rest (SSE-S3 or SSE-KMS), versioning on,
#              public access blocked, lifecycle rules to reduce cost.

# ── KMS key for clinical data encryption ─────────────────────────────────────
resource "aws_kms_key" "s3" {
  description             = "GenomeForge S3 encryption key (clinical genomic data)"
  deletion_window_in_days = 10
  enable_key_rotation     = true  # Annual rotation

  tags = { Name = "genomeforge-s3-kms" }
}

resource "aws_kms_alias" "s3" {
  name          = "alias/genomeforge-s3"
  target_key_id = aws_kms_key.s3.key_id
}

# ── Results bucket (clinical pipeline output) ─────────────────────────────────
resource "aws_s3_bucket" "results" {
  bucket = var.s3_results_bucket
  tags   = { Name = var.s3_results_bucket, DataClass = "Clinical" }
}

resource "aws_s3_bucket_versioning" "results" {
  bucket = aws_s3_bucket.results.id
  versioning_configuration { status = "Enabled" }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "results" {
  bucket = aws_s3_bucket.results.id
  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm     = "aws:kms"
      kms_master_key_id = aws_kms_key.s3.arn
    }
    bucket_key_enabled = true  # Reduces KMS API call cost
  }
}

# v5 replacement for aws_s3_bucket_acl: disable ACLs, enforce bucket owner
resource "aws_s3_bucket_ownership_controls" "results" {
  bucket = aws_s3_bucket.results.id
  rule { object_ownership = "BucketOwnerEnforced" }
}

resource "aws_s3_bucket_public_access_block" "results" {
  bucket                  = aws_s3_bucket.results.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# Lifecycle: transition to Glacier after 90 days; delete after 7 years (clinical retention)
resource "aws_s3_bucket_lifecycle_configuration" "results" {
  bucket = aws_s3_bucket.results.id
  rule {
    id     = "archive-old-results"
    status = "Enabled"
    transition {
      days          = 90
      storage_class = "GLACIER"
    }
    expiration { days = 2555 }  # ~7 years clinical data retention
  }
}

# ── Reference data bucket ─────────────────────────────────────────────────────
resource "aws_s3_bucket" "reference" {
  bucket = var.s3_reference_bucket
  tags   = { Name = var.s3_reference_bucket, DataClass = "Reference" }
}

resource "aws_s3_bucket_versioning" "reference" {
  bucket = aws_s3_bucket.reference.id
  versioning_configuration { status = "Enabled" }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "reference" {
  bucket = aws_s3_bucket.reference.id
  rule {
    apply_server_side_encryption_by_default { sse_algorithm = "AES256" }
  }
}

resource "aws_s3_bucket_ownership_controls" "reference" {
  bucket = aws_s3_bucket.reference.id
  rule { object_ownership = "BucketOwnerEnforced" }
}

resource "aws_s3_bucket_public_access_block" "reference" {
  bucket                  = aws_s3_bucket.reference.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# Intelligent tiering: reference data rarely changes after upload
resource "aws_s3_bucket_intelligent_tiering_configuration" "reference" {
  bucket = aws_s3_bucket.reference.id
  name   = "reference-intelligent-tiering"
  tiering {
    access_tier = "DEEP_ARCHIVE_ACCESS"
    days        = 180
  }
  tiering {
    access_tier = "ARCHIVE_ACCESS"
    days        = 90
  }
}

# ── ClinVar diff cache bucket ─────────────────────────────────────────────────
resource "aws_s3_bucket" "clinvar_cache" {
  bucket = var.s3_clinvar_cache_bucket
  tags   = { Name = var.s3_clinvar_cache_bucket, DataClass = "Reference" }
}

resource "aws_s3_bucket_versioning" "clinvar_cache" {
  bucket = aws_s3_bucket.clinvar_cache.id
  versioning_configuration { status = "Enabled" }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "clinvar_cache" {
  bucket = aws_s3_bucket.clinvar_cache.id
  rule {
    apply_server_side_encryption_by_default { sse_algorithm = "AES256" }
  }
}

resource "aws_s3_bucket_ownership_controls" "clinvar_cache" {
  bucket = aws_s3_bucket.clinvar_cache.id
  rule { object_ownership = "BucketOwnerEnforced" }
}

resource "aws_s3_bucket_public_access_block" "clinvar_cache" {
  bucket                  = aws_s3_bucket.clinvar_cache.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# Keep only 52 weekly ClinVar diffs (1 year)
resource "aws_s3_bucket_lifecycle_configuration" "clinvar_cache" {
  bucket = aws_s3_bucket.clinvar_cache.id
  rule {
    id     = "expire-weekly-diffs"
    status = "Enabled"
    expiration { days = 365 }
    noncurrent_version_expiration { noncurrent_days = 7 }
  }
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "s3_results_bucket_name" {
  description = "Results S3 bucket name"
  value       = aws_s3_bucket.results.id
}

output "s3_reference_bucket_name" {
  description = "Reference data S3 bucket name"
  value       = aws_s3_bucket.reference.id
}

output "s3_results_bucket_arn" {
  description = "Results S3 bucket ARN (for IAM policies)"
  value       = aws_s3_bucket.results.arn
}
