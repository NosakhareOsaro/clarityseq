# ClaritySeq Terraform Variables
# Override defaults in terraform.tfvars or via -var flags

variable "aws_region" {
  type        = string
  description = "AWS region for all resources"
  default     = "us-east-1"
}

variable "environment" {
  type        = string
  description = "Deployment environment tag (dev / staging / prod)"
  default     = "dev"

  validation {
    condition     = contains(["dev", "staging", "prod"], var.environment)
    error_message = "environment must be dev, staging, or prod"
  }
}

variable "vpc_cidr" {
  type        = string
  description = "CIDR block for the ClaritySeq VPC"
  default     = "10.0.0.0/16"
}

variable "private_subnet_cidrs" {
  type        = list(string)
  description = "CIDR blocks for private subnets (Batch + RDS + ECS)"
  default     = ["10.0.1.0/24", "10.0.2.0/24"]
}

variable "public_subnet_cidrs" {
  type        = list(string)
  description = "CIDR blocks for public subnets (NAT Gateway, Load Balancer)"
  default     = ["10.0.101.0/24", "10.0.102.0/24"]
}

variable "availability_zones" {
  type        = list(string)
  description = "Availability zones (must match subnet count)"
  default     = ["us-east-1a", "us-east-1b"]
}

variable "allowed_beacon_cidrs" {
  type        = list(string)
  description = "CIDR blocks allowed to reach the Beacon API (restrict to institutional networks)"
  default     = ["0.0.0.0/0"]  # Restrict in production to NHS/institutional IP ranges
}

# ── RDS (PostgreSQL 16) ───────────────────────────────────────────────────────
variable "db_instance_class" {
  type        = string
  description = "RDS instance class. db.t4g.micro saves 30% vs db.t3.micro (Graviton2)"
  default     = "db.t4g.micro"
  # Cost comparison (us-east-1, on-demand):
  #   db.t4g.micro:  ~$0.016/hr = ~$11.50/month
  #   db.t3.micro:   ~$0.017/hr = ~$12.50/month  (x86; slightly more expensive)
  #   db.t4g.small:  ~$0.032/hr = ~$23/month      (use for production/staging)
}

variable "db_name" {
  type        = string
  description = "PostgreSQL database name"
  default     = "clarityseq"
}

variable "db_username" {
  type        = string
  description = "PostgreSQL admin username"
  default     = "clarityseq"
  sensitive   = true
}

variable "db_password" {
  type        = string
  description = "PostgreSQL admin password (use AWS Secrets Manager in production)"
  sensitive   = true
  # No default — must be provided via TF_VAR_db_password or -var flag
}

variable "db_storage_gb" {
  type        = number
  description = "RDS allocated storage in GB"
  default     = 20
}

variable "db_backup_retention_days" {
  type        = number
  description = "RDS automated backup retention (days). Minimum 7 for clinical data"
  default     = 7
}

# ── S3 ────────────────────────────────────────────────────────────────────────
variable "s3_results_bucket" {
  type        = string
  description = "S3 bucket name for pipeline results and reports"
  default     = "clarityseq-results"
}

variable "s3_reference_bucket" {
  type        = string
  description = "S3 bucket for reference data (DRAGMAP index, gnomAD v4.1, VEP cache)"
  default     = "clarityseq-reference"
}

variable "s3_clinvar_cache_bucket" {
  type        = string
  description = "S3 bucket for ClinVar weekly diff cache (Celery daemon)"
  default     = "clarityseq-clinvar-cache"
}

# ── AWS Batch ────────────────────────────────────────────────────────────────
variable "batch_spot_bid_percentage" {
  type        = number
  description = "Max Spot bid as % of on-demand price (60% saves ~40% vs on-demand)"
  default     = 60
}

variable "batch_max_vcpus" {
  type        = number
  description = "Maximum vCPUs across all Batch jobs (cost guard)"
  default     = 256
}

variable "batch_instance_types" {
  type        = list(string)
  description = "EC2 instance types for Batch Spot fleet (memory-optimised for WGS)"
  default     = ["r6i.8xlarge", "r6i.4xlarge", "r5.8xlarge", "r5.4xlarge"]
  # r6i.8xlarge: 32 vCPU / 256 GB RAM — ideal for GATK HaplotypeCaller
  # Minimum for HPRC pangenome: r6i.8xlarge (needs 256 GB RAM for GBZ graph)
}

# ── ECR ───────────────────────────────────────────────────────────────────────
variable "ecr_image_tag_mutability" {
  type        = string
  description = "ECR image tag mutability (IMMUTABLE for production audit trail)"
  default     = "IMMUTABLE"
}

variable "ecr_scan_on_push" {
  type        = bool
  description = "Enable ECR image scanning on push (CVE detection)"
  default     = true
}
