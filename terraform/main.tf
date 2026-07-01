# ClaritySeq — Terraform Root Module
# AWS provider v5.x (REQUIRED: v4→v5 breaking changes)
# Breaking change: aws_s3_bucket_acl removed — use aws_s3_bucket_ownership_controls
# Breaking change: AWS provider ≥5.0 requires Terraform ≥1.5.0
#
# Architecture:
#   - AWS Batch (Spot) for Nextflow pipeline execution (60–80% cost saving)
#   - RDS PostgreSQL 16 on db.t4g.micro (Graviton2; 30% cheaper than t3.micro)
#   - ECR for Docker image storage (pipeline, beacon, daemon)
#   - S3 for FASTQs, results, reference data, ClinVar cache
#   - ECS Fargate for Beacon API + Celery daemon
#   - VPC with private subnets (no public internet for clinical data)
#
# Estimated monthly cost (us-east-1, light usage):
#   RDS db.t4g.micro (PostgreSQL 16):  ~$13/month
#   S3 (1 TB reference + results):     ~$23/month
#   ECR (3 images, ~10 GB):           ~$1/month
#   AWS Batch Spot (per-run):          ~$8–15/WGS run
#   ECS Fargate (Beacon API):          ~$5/month (0.25 vCPU, 0.5 GB)
#   Total idle (no runs):              ~$42/month

terraform {
  required_version = ">= 1.5.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"  # v5 required; v4 has aws_s3_bucket_acl which is removed in v5
    }
  }

  # Remote state — use S3 backend in production
  # backend "s3" {
  #   bucket         = "clarityseq-tfstate"
  #   key            = "clarityseq/terraform.tfstate"
  #   region         = "us-east-1"
  #   encrypt        = true
  #   dynamodb_table = "clarityseq-tfstate-lock"
  # }
}

provider "aws" {
  region = var.aws_region

  default_tags {
    tags = {
      Project     = "ClaritySeq"
      Environment = var.environment
      ManagedBy   = "Terraform"
    }
  }
}

# ── Data sources ──────────────────────────────────────────────────────────────
data "aws_caller_identity" "current" {}
data "aws_region" "current" {}

# ── VPC ───────────────────────────────────────────────────────────────────────
resource "aws_vpc" "main" {
  cidr_block           = var.vpc_cidr
  enable_dns_support   = true
  enable_dns_hostnames = true

  tags = { Name = "clarityseq-vpc" }
}

resource "aws_internet_gateway" "main" {
  vpc_id = aws_vpc.main.id
  tags   = { Name = "clarityseq-igw" }
}

resource "aws_subnet" "private" {
  count             = length(var.private_subnet_cidrs)
  vpc_id            = aws_vpc.main.id
  cidr_block        = var.private_subnet_cidrs[count.index]
  availability_zone = var.availability_zones[count.index]

  tags = { Name = "clarityseq-private-${count.index + 1}" }
}

resource "aws_subnet" "public" {
  count                   = length(var.public_subnet_cidrs)
  vpc_id                  = aws_vpc.main.id
  cidr_block              = var.public_subnet_cidrs[count.index]
  availability_zone       = var.availability_zones[count.index]
  map_public_ip_on_launch = false  # Never auto-assign public IPs

  tags = { Name = "clarityseq-public-${count.index + 1}" }
}

# NAT Gateway for private subnet outbound (reference data downloads)
resource "aws_eip" "nat" {
  domain = "vpc"
  tags   = { Name = "clarityseq-nat-eip" }
}

resource "aws_nat_gateway" "main" {
  allocation_id = aws_eip.nat.id
  subnet_id     = aws_subnet.public[0].id
  tags          = { Name = "clarityseq-nat" }
  depends_on    = [aws_internet_gateway.main]
}

resource "aws_route_table" "private" {
  vpc_id = aws_vpc.main.id
  route {
    cidr_block     = "0.0.0.0/0"
    nat_gateway_id = aws_nat_gateway.main.id
  }
  tags = { Name = "clarityseq-rt-private" }
}

resource "aws_route_table_association" "private" {
  count          = length(aws_subnet.private)
  subnet_id      = aws_subnet.private[count.index].id
  route_table_id = aws_route_table.private.id
}

# ── Security Groups ───────────────────────────────────────────────────────────
resource "aws_security_group" "batch" {
  name        = "clarityseq-batch"
  description = "AWS Batch compute environment (WGS pipeline workers)"
  vpc_id      = aws_vpc.main.id

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound (reference data + S3 + ECR)"
  }

  tags = { Name = "clarityseq-batch-sg" }
}

resource "aws_security_group" "rds" {
  name        = "clarityseq-rds"
  description = "PostgreSQL 16 (BayesACMG results + VUS tracking)"
  vpc_id      = aws_vpc.main.id

  ingress {
    from_port       = 5432
    to_port         = 5432
    protocol        = "tcp"
    security_groups = [aws_security_group.batch.id, aws_security_group.ecs.id]
    description     = "PostgreSQL from Batch workers and ECS services"
  }

  tags = { Name = "clarityseq-rds-sg" }
}

resource "aws_security_group" "ecs" {
  name        = "clarityseq-ecs"
  description = "ECS Fargate (Beacon API + Celery daemon)"
  vpc_id      = aws_vpc.main.id

  ingress {
    from_port   = 8080
    to_port     = 8080
    protocol    = "tcp"
    cidr_blocks = var.allowed_beacon_cidrs
    description = "GA4GH Beacon v2.1.1 API"
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
    description = "Allow all outbound (ClinVar API + NCBI)"
  }

  tags = { Name = "clarityseq-ecs-sg" }
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "vpc_id" {
  description = "VPC ID"
  value       = aws_vpc.main.id
}

output "private_subnet_ids" {
  description = "Private subnet IDs for Batch and RDS"
  value       = aws_subnet.private[*].id
}

output "aws_account_id" {
  description = "AWS account ID"
  value       = data.aws_caller_identity.current.account_id
}
