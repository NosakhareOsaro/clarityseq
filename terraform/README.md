# terraform/

AWS infrastructure as code using Terraform with AWS Provider v5.

## Architecture

```
VPC (eu-west-2)
├── Public subnets (Beacon API Load Balancer)
├── Private subnets
│   ├── AWS Batch (Spot fleet) — pipeline compute
│   ├── RDS PostgreSQL 16 (db.t4g.micro) — variant/audit database
│   └── ElastiCache Redis 7.x — Celery broker + BayesACMG cache
└── ECR repositories (pipeline, beacon, daemon images)
```

## AWS Provider v5 note

This Terraform uses `hashicorp/aws ~> 5.0`. **Breaking change from v4:**
`aws_s3_bucket_acl` is removed — use `aws_s3_bucket_ownership_controls` instead.
See: https://registry.terraform.io/providers/hashicorp/aws/latest/docs/guides/version-5-upgrade

## Database

- **Engine**: PostgreSQL 16 (`db.t4g.micro` — ARM Graviton2)
- **vs db.t3.micro**: ~30% cheaper (~$12/month vs ~$16/month)
- **PostgreSQL 16 benefits**: Better JSONB query performance (ACMG evidence stored as JSONB); improved logical replication for Beacon read replicas

## Deploy

```bash
cd terraform
terraform init
terraform plan -out=tfplan
terraform apply tfplan
```

## Cost estimate

- Batch (Spot): ~$0.10–$0.50/sample (variable)
- RDS db.t4g.micro: ~$12/month
- ElastiCache cache.t4g.micro: ~$12/month
- ECR: ~$0.10/GB/month
- Budget alert: $100/month (aws_budgets_budget)
