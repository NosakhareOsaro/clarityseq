# GenomeForge — AWS Batch (Spot fleet for Nextflow pipeline)
#
# Cost model: Spot instances at 60% of on-demand price = ~40% cost saving.
# A typical 30× WGS run (DRAGMAP + GATK4 + DeepVariant) costs ~$8–15 on Spot.
# vs. ~$20–40 on-demand.
#
# Instance selection: r6i.8xlarge preferred (32 vCPU / 256 GB RAM).
# The HPRC pangenome GBZ graph requires ≥256 GB RAM; r6i.4xlarge is insufficient.
# For standard GRCh38 WGS (no pangenome), r6i.4xlarge (16 vCPU / 128 GB) works.
#
# Nextflow config:
#   process.executor = 'awsbatch'
#   process.queue    = 'genomeforge-spot'
#   aws.region       = 'us-east-1'
#   aws.batch.jobRole = '<batch_job_role_arn>'

# ── Compute environment (Spot) ────────────────────────────────────────────────
resource "aws_batch_compute_environment" "spot" {
  compute_environment_name = "genomeforge-spot"
  type                     = "MANAGED"
  state                    = "ENABLED"

  compute_resources {
    type                = "SPOT"
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"  # Prefer pools with most capacity
    bid_percentage      = var.batch_spot_bid_percentage  # Default: 60% of on-demand

    min_vcpus     = 0   # Scale to zero when idle (cost saving)
    max_vcpus     = var.batch_max_vcpus
    desired_vcpus = 0

    instance_type = var.batch_instance_types

    subnets            = aws_subnet.private[*].id
    security_group_ids = [aws_security_group.batch.id]

    instance_role = aws_iam_instance_profile.batch.arn

    # ECS-optimised AMI (latest; managed by AWS)
    # Ubuntu 24.04 ECS-optimised AMI would require custom AMI — use default ECS AMI
    # (GATK tools run inside Docker containers, so host OS is less critical)

    tags = { Name = "genomeforge-batch-spot" }
  }

  service_role = aws_iam_role.batch_service.arn

  depends_on = [
    aws_iam_role_policy_attachment.batch_service,
    aws_iam_instance_profile.batch,
  ]
}

# On-demand fallback (for jobs requiring guaranteed capacity: clinical urgency)
resource "aws_batch_compute_environment" "ondemand" {
  compute_environment_name = "genomeforge-ondemand"
  type                     = "MANAGED"
  state                    = "ENABLED"

  compute_resources {
    type                = "EC2"
    allocation_strategy = "BEST_FIT_PROGRESSIVE"

    min_vcpus     = 0
    max_vcpus     = 32  # Limited: on-demand is expensive; use for clinical urgency only
    desired_vcpus = 0

    instance_type = ["r6i.4xlarge"]

    subnets            = aws_subnet.private[*].id
    security_group_ids = [aws_security_group.batch.id]

    instance_role = aws_iam_instance_profile.batch.arn

    tags = { Name = "genomeforge-batch-ondemand" }
  }

  service_role = aws_iam_role.batch_service.arn

  depends_on = [
    aws_iam_role_policy_attachment.batch_service,
    aws_iam_instance_profile.batch,
  ]
}

# ── Job queues ────────────────────────────────────────────────────────────────
# Spot queue with on-demand fallback (Nextflow uses this for all pipeline steps)
resource "aws_batch_job_queue" "main" {
  name     = "genomeforge-spot"
  state    = "ENABLED"
  priority = 1

  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.spot.arn
  }

  compute_environment_order {
    order               = 2
    compute_environment = aws_batch_compute_environment.ondemand.arn  # Fallback
  }
}

# High-priority queue for clinical-urgent runs (on-demand only)
resource "aws_batch_job_queue" "urgent" {
  name     = "genomeforge-urgent"
  state    = "ENABLED"
  priority = 100  # Higher priority than main queue

  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.ondemand.arn
  }
}

# ── Batch service role ────────────────────────────────────────────────────────
resource "aws_iam_role" "batch_service" {
  name = "genomeforge-batch-service"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "batch.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "batch_service" {
  role       = aws_iam_role.batch_service.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "batch_job_queue_arn" {
  description = "Spot job queue ARN (use as process.queue in nextflow.config aws profile)"
  value       = aws_batch_job_queue.main.arn
}

output "batch_job_queue_name" {
  description = "Spot job queue name"
  value       = aws_batch_job_queue.main.name
}

output "batch_urgent_queue_arn" {
  description = "On-demand urgent queue ARN (clinical-urgent runs)"
  value       = aws_batch_job_queue.urgent.arn
}
