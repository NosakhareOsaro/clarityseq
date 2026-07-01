# ClaritySeq — Terraform Outputs
# Reference these in nextflow.config aws profile and CI/CD secrets

output "region" {
  description = "AWS region"
  value       = var.aws_region
}

output "environment" {
  description = "Deployment environment"
  value       = var.environment
}

# ── Summary output block for CI/CD ───────────────────────────────────────────
output "deployment_summary" {
  description = "Key deployment values for nextflow.config and GitHub Actions"
  value = {
    region              = var.aws_region
    results_bucket      = aws_s3_bucket.results.id
    reference_bucket    = aws_s3_bucket.reference.id
    batch_queue         = aws_batch_job_queue.main.name
    batch_urgent_queue  = aws_batch_job_queue.urgent.name
    db_endpoint         = "${aws_db_instance.main.address}:${aws_db_instance.main.port}"
    ecr_pipeline        = aws_ecr_repository.pipeline.repository_url
    ecr_beacon          = aws_ecr_repository.beacon.repository_url
    ecr_daemon          = aws_ecr_repository.daemon.repository_url
    github_actions_role = aws_iam_role.github_actions.arn
    batch_job_role      = aws_iam_role.batch_job.arn
  }
  sensitive = false
}
