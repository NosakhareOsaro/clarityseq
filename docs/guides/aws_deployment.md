# AWS Deployment Guide

Deploy GenomeForge on AWS using Terraform (provider v5) + AWS Batch + RDS PostgreSQL 16.

## Architecture overview

```
┌─────────────────────────────────────────────────────────────────┐
│ AWS (eu-west-2)                                                 │
│                                                                 │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────┐  │
│  │ S3 (inputs)  │    │ S3 (results) │    │ S3 (references)  │  │
│  └──────┬───────┘    └──────────────┘    └──────────────────┘  │
│         │                                                       │
│  ┌──────▼───────────────────────────────┐                      │
│  │ AWS Batch (Spot fleet)               │                      │
│  │ - ECS tasks from ECR images          │                      │
│  │ - GPU instances for DeepVariant      │                      │
│  └──────────────────────────────────────┘                      │
│                                                                 │
│  ┌──────────────┐    ┌──────────────┐                          │
│  │ RDS           │    │ ElastiCache  │                          │
│  │ PostgreSQL 16 │    │ Redis 7.x    │                          │
│  │ db.t4g.micro  │    │ t4g.micro    │                          │
│  └──────────────┘    └──────────────┘                          │
│                                                                 │
│  ┌──────────────────────────────────────┐                      │
│  │ ECS Fargate (always-on services)     │                      │
│  │ - Beacon API (beacon:latest)         │                      │
│  │ - Reclassification daemon            │                      │
│  └──────────────────────────────────────┘                      │
└─────────────────────────────────────────────────────────────────┘
```

## Prerequisites

```bash
# AWS CLI v2
aws --version  # ≥ 2.x

# Terraform ≥ 1.7
terraform --version

# GitHub CLI (for CI/CD setup)
gh --version
```

## Initial deployment

```bash
# 1. Configure AWS credentials
aws configure --profile genomeforge

# 2. Create Terraform state bucket
aws s3 mb s3://genomeforge-terraform-state --region eu-west-2

# 3. Deploy infrastructure
cd terraform
terraform init
terraform plan -var="aws_region=eu-west-2" -out=tfplan
terraform apply tfplan

# 4. Push Docker images to ECR
make docker-build
make ecr-push

# 5. Run a test sample
aws batch submit-job \
  --job-name genomeforge-test-HG001 \
  --job-queue genomeforge-spot \
  --job-definition genomeforge-pipeline \
  --container-overrides '{"command": ["nextflow", "run", "pipelines/wgs_grch38.nf", "-profile", "test"]}'
```

## Cost optimisation

- **Batch Spot fleet**: 60–90% cost reduction vs On-Demand
- **db.t4g.micro**: ~30% cheaper than db.t3.micro for same workload
- **Budget alert**: $100/month CloudWatch alarm (configured in `terraform/main.tf`)

## AWS Provider v5 note

This configuration uses `hashicorp/aws ~> 5.0`. **Breaking change from v4:**
- `aws_s3_bucket_acl` removed → use `aws_s3_bucket_ownership_controls`
- Migration guide: https://registry.terraform.io/providers/hashicorp/aws/latest/docs/guides/version-5-upgrade
