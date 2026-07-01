# ClaritySeq — RDS PostgreSQL 16 (db.t4g.micro)
#
# Cost rationale:
#   db.t4g.micro on Graviton2 ARM processor:  ~$0.016/hr ($11.50/month)
#   db.t3.micro on x86:                       ~$0.017/hr ($12.50/month)
#   Graviton2 is ~30% cheaper per vCPU for I/O-bound workloads like Postgres.
#   Use db.t4g.small ($23/month) for staging/prod with concurrent beacon + daemon.
#
# PostgreSQL 16 benefits over 15:
#   - JSONB performance improvements (BayesACMG evidence storage)
#   - Logical replication slot management improvements
#   - Parallel query improvements
#
# Schema managed by Alembic migrations in reclassification/alembic/
# Stores: variant classifications, VUS review schedules, ClinVar submission queue

# ── Subnet group (Multi-AZ) ───────────────────────────────────────────────────
resource "aws_db_subnet_group" "main" {
  name        = "clarityseq-rds-subnet-group"
  description = "Private subnets for ClaritySeq PostgreSQL 16"
  subnet_ids  = aws_subnet.private[*].id

  tags = { Name = "clarityseq-rds-subnet-group" }
}

# ── Parameter group (PostgreSQL 16 tuning) ────────────────────────────────────
resource "aws_db_parameter_group" "postgres16" {
  name        = "clarityseq-postgres16"
  family      = "postgres16"
  description = "ClaritySeq PostgreSQL 16 parameters"

  # Optimise for t4g.micro (1 GB RAM): reduce shared_buffers from default 128 MB
  parameter {
    name  = "shared_buffers"
    value = "128000"  # 128 MB (default; leave as-is for t4g.micro)
  }

  # Enable pg_stat_statements for query performance monitoring
  parameter {
    name  = "pg_stat_statements.track"
    value = "all"
    apply_method = "pending-reboot"
  }

  # JSONB GIN index autovacuum — important for BayesACMG evidence JSONB columns
  parameter {
    name  = "autovacuum_vacuum_scale_factor"
    value = "0.05"  # Vacuum when 5% of rows are dead (more aggressive than default 20%)
  }

  tags = { Name = "clarityseq-postgres16-params" }
}

# ── Master password (Secrets Manager) ────────────────────────────────────────
resource "aws_secretsmanager_secret" "db_password" {
  name                    = "clarityseq/db/password"
  description             = "ClaritySeq PostgreSQL 16 admin password"
  recovery_window_in_days = 7  # 7-day recovery window before deletion

  tags = { Name = "clarityseq-db-password" }
}

resource "aws_secretsmanager_secret_version" "db_password" {
  secret_id     = aws_secretsmanager_secret.db_password.id
  secret_string = jsonencode({
    username = var.db_username
    password = var.db_password
    host     = aws_db_instance.main.address
    port     = aws_db_instance.main.port
    dbname   = var.db_name
  })
}

# ── RDS instance ──────────────────────────────────────────────────────────────
resource "aws_db_instance" "main" {
  identifier     = "clarityseq-postgres16"
  engine         = "postgres"
  engine_version = "16.4"  # PostgreSQL 16 LTS

  # db.t4g.micro: Graviton2, 2 vCPU burst, 1 GB RAM, 2,085 IOPS burst
  # 30% cheaper than db.t3.micro; equivalent or better single-thread performance
  instance_class = var.db_instance_class

  allocated_storage     = var.db_storage_gb
  max_allocated_storage = 100  # Auto-scaling up to 100 GB
  storage_type          = "gp3"
  storage_encrypted     = true  # REQUIRED for clinical genomic data

  db_name  = var.db_name
  username = var.db_username
  password = var.db_password

  db_subnet_group_name   = aws_db_subnet_group.main.name
  vpc_security_group_ids = [aws_security_group.rds.id]
  parameter_group_name   = aws_db_parameter_group.postgres16.name

  multi_az               = var.environment == "prod"  # Multi-AZ in prod only
  publicly_accessible    = false  # Never expose PostgreSQL to internet
  deletion_protection    = var.environment == "prod"
  skip_final_snapshot    = var.environment != "prod"

  backup_retention_period = var.db_backup_retention_days
  backup_window           = "03:00-04:00"  # UTC; off-peak for UK clinical systems
  maintenance_window      = "sun:04:00-sun:05:00"

  # Enhanced monitoring (free at 60s granularity)
  monitoring_interval = 60
  monitoring_role_arn = aws_iam_role.rds_monitoring.arn

  # Enable Performance Insights (free tier: 7-day retention)
  performance_insights_enabled          = true
  performance_insights_retention_period = 7

  tags = { Name = "clarityseq-postgres16" }

  lifecycle {
    # Prevent accidental deletion of the database
    prevent_destroy = false  # Set to true in prod
  }
}

# ── IAM role for Enhanced Monitoring ─────────────────────────────────────────
resource "aws_iam_role" "rds_monitoring" {
  name = "clarityseq-rds-monitoring"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "monitoring.rds.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "rds_monitoring" {
  role       = aws_iam_role.rds_monitoring.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonRDSEnhancedMonitoringRole"
}

# ── Outputs ───────────────────────────────────────────────────────────────────
output "db_endpoint" {
  description = "PostgreSQL 16 endpoint (hostname:port)"
  value       = "${aws_db_instance.main.address}:${aws_db_instance.main.port}"
  sensitive   = false
}

output "db_secret_arn" {
  description = "ARN of the Secrets Manager secret containing DB credentials"
  value       = aws_secretsmanager_secret.db_password.arn
}
