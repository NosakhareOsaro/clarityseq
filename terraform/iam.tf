# GenomeForge — IAM Roles and Policies
# Principle of least privilege: each service gets only what it needs.
#
# Roles:
#   genomeforge-batch-job     — Nextflow pipeline jobs (S3 + ECR + Secrets Manager)
#   genomeforge-batch-instance — EC2 instance profile for Batch compute
#   genomeforge-ecs-task      — ECS Fargate tasks (Beacon API + Celery daemon)
#   genomeforge-github-actions — OIDC role for GitHub Actions CI/CD

# ── Batch job role ────────────────────────────────────────────────────────────
resource "aws_iam_role" "batch_job" {
  name        = "genomeforge-batch-job"
  description = "IAM role for GenomeForge Nextflow pipeline Batch jobs"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ecs-tasks.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "batch_job_s3" {
  name = "genomeforge-batch-s3"
  role = aws_iam_role.batch_job.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "ReadReference"
        Effect = "Allow"
        Action = ["s3:GetObject", "s3:ListBucket"]
        Resource = [
          aws_s3_bucket.reference.arn,
          "${aws_s3_bucket.reference.arn}/*",
        ]
      },
      {
        Sid    = "WriteResults"
        Effect = "Allow"
        Action = ["s3:PutObject", "s3:GetObject", "s3:DeleteObject", "s3:ListBucket"]
        Resource = [
          aws_s3_bucket.results.arn,
          "${aws_s3_bucket.results.arn}/*",
        ]
      },
      {
        Sid    = "KMSDecrypt"
        Effect = "Allow"
        Action = ["kms:Decrypt", "kms:GenerateDataKey"]
        Resource = [aws_kms_key.s3.arn]
      },
    ]
  })
}

resource "aws_iam_role_policy_attachment" "batch_job_ecr" {
  role       = aws_iam_role.batch_job.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly"
}

# ── Batch instance profile ────────────────────────────────────────────────────
resource "aws_iam_role" "batch_instance" {
  name = "genomeforge-batch-instance"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ec2.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "batch_instance_ecs" {
  role       = aws_iam_role.batch_instance.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_instance_profile" "batch" {
  name = "genomeforge-batch"
  role = aws_iam_role.batch_instance.name
}

# ── ECS task role (Beacon API + Celery daemon) ────────────────────────────────
resource "aws_iam_role" "ecs_task" {
  name        = "genomeforge-ecs-task"
  description = "IAM role for ECS Fargate tasks (Beacon API + Celery daemon)"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ecs-tasks.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "ecs_task_secrets" {
  name = "genomeforge-ecs-secrets"
  role = aws_iam_role.ecs_task.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "ReadSecrets"
        Effect = "Allow"
        Action = ["secretsmanager:GetSecretValue"]
        Resource = [
          aws_secretsmanager_secret.db_password.arn,
        ]
      },
      {
        Sid    = "ReadClinVarCache"
        Effect = "Allow"
        Action = ["s3:GetObject", "s3:PutObject", "s3:ListBucket"]
        Resource = [
          aws_s3_bucket.clinvar_cache.arn,
          "${aws_s3_bucket.clinvar_cache.arn}/*",
        ]
      },
    ]
  })
}

# ECS execution role (pulls images from ECR, writes logs to CloudWatch)
resource "aws_iam_role" "ecs_execution" {
  name = "genomeforge-ecs-execution"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ecs-tasks.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "ecs_execution" {
  role       = aws_iam_role.ecs_execution.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

# ── GitHub Actions OIDC role ──────────────────────────────────────────────────
# Enables GitHub Actions to push to ECR without long-lived access keys
data "aws_iam_openid_connect_provider" "github" {
  url = "https://token.actions.githubusercontent.com"
}

resource "aws_iam_role" "github_actions" {
  name        = "genomeforge-github-actions"
  description = "GitHub Actions OIDC role for GenomeForge CI/CD"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Principal = {
        Federated = data.aws_iam_openid_connect_provider.github.arn
      }
      Action = "sts:AssumeRoleWithWebIdentity"
      Condition = {
        StringEquals = {
          "token.actions.githubusercontent.com:aud" = "sts.amazonaws.com"
        }
        StringLike = {
          # Restrict to the GenomeForge repository only
          "token.actions.githubusercontent.com:sub" = "repo:genomeforge/genomeforge:*"
        }
      }
    }]
  })
}

resource "aws_iam_role_policy" "github_actions_ecr" {
  name = "genomeforge-github-ecr-push"
  role = aws_iam_role.github_actions.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "ECRAuth"
        Effect = "Allow"
        Action = ["ecr:GetAuthorizationToken"]
        Resource = ["*"]
      },
      {
        Sid    = "ECRPush"
        Effect = "Allow"
        Action = [
          "ecr:BatchCheckLayerAvailability",
          "ecr:GetDownloadUrlForLayer",
          "ecr:BatchGetImage",
          "ecr:InitiateLayerUpload",
          "ecr:UploadLayerPart",
          "ecr:CompleteLayerUpload",
          "ecr:PutImage",
        ]
        Resource = [
          aws_ecr_repository.pipeline.arn,
          aws_ecr_repository.beacon.arn,
          aws_ecr_repository.daemon.arn,
        ]
      },
    ]
  })
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "batch_job_role_arn" {
  description = "ARN of the Batch job IAM role (use in Nextflow config)"
  value       = aws_iam_role.batch_job.arn
}

output "github_actions_role_arn" {
  description = "ARN for GitHub Actions OIDC role (add to repo Actions secrets)"
  value       = aws_iam_role.github_actions.arn
}
