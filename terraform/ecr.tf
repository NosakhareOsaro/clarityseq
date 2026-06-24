# GenomeForge — ECR Repositories
# Three images: pipeline (GATK/DRAGMAP/VEP), beacon (FastAPI), daemon (Celery)
# All repositories: IMMUTABLE tags (audit trail); vulnerability scanning on push

resource "aws_ecr_repository" "pipeline" {
  name                 = "genomeforge/pipeline"
  image_tag_mutability = var.ecr_image_tag_mutability

  image_scanning_configuration {
    scan_on_push = var.ecr_scan_on_push
  }

  encryption_configuration {
    encryption_type = "KMS"
    kms_key         = aws_kms_key.s3.arn
  }

  tags = { Name = "genomeforge-pipeline-ecr" }
}

resource "aws_ecr_repository" "beacon" {
  name                 = "genomeforge/beacon"
  image_tag_mutability = var.ecr_image_tag_mutability

  image_scanning_configuration {
    scan_on_push = var.ecr_scan_on_push
  }

  encryption_configuration {
    encryption_type = "KMS"
    kms_key         = aws_kms_key.s3.arn
  }

  tags = { Name = "genomeforge-beacon-ecr" }
}

resource "aws_ecr_repository" "daemon" {
  name                 = "genomeforge/daemon"
  image_tag_mutability = var.ecr_image_tag_mutability

  image_scanning_configuration {
    scan_on_push = var.ecr_scan_on_push
  }

  encryption_configuration {
    encryption_type = "KMS"
    kms_key         = aws_kms_key.s3.arn
  }

  tags = { Name = "genomeforge-daemon-ecr" }
}

# ── Lifecycle policies: keep last 10 tagged images + remove untagged after 1 day ──
locals {
  ecr_lifecycle_policy = jsonencode({
    rules = [
      {
        rulePriority = 1
        description  = "Remove untagged images after 1 day"
        selection = {
          tagStatus   = "untagged"
          countType   = "sinceImagePushed"
          countUnit   = "days"
          countNumber = 1
        }
        action = { type = "expire" }
      },
      {
        rulePriority = 2
        description  = "Keep last 10 tagged images"
        selection = {
          tagStatus     = "tagged"
          tagPrefixList = ["v", "sha-"]
          countType     = "imageCountMoreThan"
          countNumber   = 10
        }
        action = { type = "expire" }
      }
    ]
  })
}

resource "aws_ecr_lifecycle_policy" "pipeline" {
  repository = aws_ecr_repository.pipeline.name
  policy     = local.ecr_lifecycle_policy
}

resource "aws_ecr_lifecycle_policy" "beacon" {
  repository = aws_ecr_repository.beacon.name
  policy     = local.ecr_lifecycle_policy
}

resource "aws_ecr_lifecycle_policy" "daemon" {
  repository = aws_ecr_repository.daemon.name
  policy     = local.ecr_lifecycle_policy
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "ecr_pipeline_url" {
  description = "ECR URL for pipeline image (use in Nextflow config)"
  value       = aws_ecr_repository.pipeline.repository_url
}

output "ecr_beacon_url" {
  description = "ECR URL for Beacon API image"
  value       = aws_ecr_repository.beacon.repository_url
}

output "ecr_daemon_url" {
  description = "ECR URL for reclassification daemon image"
  value       = aws_ecr_repository.daemon.repository_url
}
