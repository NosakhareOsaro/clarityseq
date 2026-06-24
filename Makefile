# =============================================================================
# GenomeForge — Makefile
# =============================================================================
# Convenience targets for development, testing, and deployment.
# All targets are documented. Run `make help` to see available targets.
#
# Prerequisites: see docs/setup/prerequisites.md
# Python 3.12 required exactly.
# =============================================================================

SHELL := /bin/bash
.ONESHELL:
.DEFAULT_GOAL := help

# Python and tool paths
PYTHON := python3.12
PIP    := $(PYTHON) -m pip
PYTEST := $(PYTHON) -m pytest
BLACK  := $(PYTHON) -m black
RUFF   := $(PYTHON) -m ruff
MYPY   := $(PYTHON) -m mypy

# Project directories
VENV_DIR   := .venv
SRC_DIRS   := bayesacmg/src beacon_api annotation reclassification reporting pgx prioritisation
TEST_DIRS  := bayesacmg/tests beacon_api/tests annotation/tests reclassification/tests

# Docker image names (ECR prefix set in .env)
PIPELINE_IMAGE := genomeforge/pipeline
BEACON_IMAGE   := genomeforge/beacon
DAEMON_IMAGE   := genomeforge/daemon

# Colours for output
RESET  := \033[0m
BOLD   := \033[1m
GREEN  := \033[0;32m
YELLOW := \033[0;33m
CYAN   := \033[0;36m

# =============================================================================
.PHONY: help
help: ## Show this help message
	@printf "$(BOLD)GenomeForge — Available Make Targets$(RESET)\n"
	@printf "$(CYAN)═══════════════════════════════════════════════════════════$(RESET)\n"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  $(GREEN)%-20s$(RESET) %s\n", $$1, $$2}'
	@printf "\n"

# =============================================================================
.PHONY: install
install: ## Install all Python dependencies in a virtual environment
	@printf "$(BOLD)Installing GenomeForge dependencies...$(RESET)\n"
	$(PYTHON) -m venv $(VENV_DIR)
	source $(VENV_DIR)/bin/activate && \
		$(PIP) install --upgrade pip && \
		$(PIP) install -e ".[dev,docs]" && \
		$(PIP) install -e bayesacmg/ && \
		pre-commit install && \
		pre-commit install --hook-type commit-msg
	@printf "$(GREEN)Installation complete. Activate with: source $(VENV_DIR)/bin/activate$(RESET)\n"

# =============================================================================
.PHONY: test
test: ## Run the full test suite with coverage (requires ≥ 90% coverage)
	@printf "$(BOLD)Running tests with coverage...$(RESET)\n"
	$(PYTEST) \
		--cov=. \
		--cov-report=html:htmlcov \
		--cov-report=xml:coverage.xml \
		--cov-report=term-missing \
		--cov-fail-under=90 \
		-v \
		--tb=short \
		$(TEST_DIRS)
	@printf "$(GREEN)Tests complete. Coverage report: htmlcov/index.html$(RESET)\n"

# =============================================================================
.PHONY: test-unit
test-unit: ## Run only unit tests (fast; no external services required)
	@printf "$(BOLD)Running unit tests...$(RESET)\n"
	$(PYTEST) \
		-m "not integration and not slow" \
		-v \
		--tb=short \
		$(TEST_DIRS)

# =============================================================================
.PHONY: test-acmg
test-acmg: ## Run ACMG rule unit tests (uses real ClinVar RCV accessions)
	@printf "$(BOLD)Running ACMG classification tests...$(RESET)\n"
	$(PYTEST) \
		-m "acmg" \
		-v \
		--tb=long \
		bayesacmg/tests/

# =============================================================================
.PHONY: lint
lint: ## Run all linters (ruff, black --check, mypy) without modifying files
	@printf "$(BOLD)Running linters...$(RESET)\n"
	@printf "$(CYAN)→ black (check mode)$(RESET)\n"
	$(BLACK) --check --diff $(SRC_DIRS)
	@printf "$(CYAN)→ ruff$(RESET)\n"
	$(RUFF) check $(SRC_DIRS)
	@printf "$(CYAN)→ mypy (strict, Python 3.12)$(RESET)\n"
	$(MYPY) $(SRC_DIRS)
	@printf "$(GREEN)All linters passed.$(RESET)\n"

# =============================================================================
.PHONY: format
format: ## Auto-format code with black and ruff --fix
	@printf "$(BOLD)Formatting code...$(RESET)\n"
	$(BLACK) $(SRC_DIRS)
	$(RUFF) check --fix $(SRC_DIRS)
	@printf "$(GREEN)Formatting complete.$(RESET)\n"

# =============================================================================
.PHONY: docs
docs: ## Build Sphinx documentation (HTML output in docs/_build/html/)
	@printf "$(BOLD)Building documentation...$(RESET)\n"
	$(PYTHON) -m sphinx \
		-b html \
		-d docs/_build/.doctrees \
		docs/ \
		docs/_build/html/
	@printf "$(GREEN)Docs built: docs/_build/html/index.html$(RESET)\n"

# =============================================================================
.PHONY: figures
figures: ## Generate all figures for JOSS paper and documentation
	@printf "$(BOLD)Generating figures...$(RESET)\n"
	@# Architecture diagram (Mermaid → PNG via mmdc)
	@if command -v mmdc >/dev/null 2>&1; then \
		mmdc -i joss_paper/figures/architecture.mmd -o joss_paper/figures/architecture.png -w 1200; \
		printf "$(GREEN)Architecture diagram generated.$(RESET)\n"; \
	else \
		printf "$(YELLOW)mmdc not found. Install with: npm install -g @mermaid-js/mermaid-cli$(RESET)\n"; \
	fi
	@# Calibration plots — requires calibration module
	$(PYTHON) calibration/run_calibration.py --figures-only || \
		printf "$(YELLOW)Calibration figures skipped (run calibration first).$(RESET)\n"
	@printf "$(GREEN)Figure generation complete.$(RESET)\n"

# =============================================================================
.PHONY: docker-build
docker-build: ## Build all Docker images (pipeline, beacon, daemon)
	@printf "$(BOLD)Building Docker images...$(RESET)\n"
	@printf "$(CYAN)→ Pipeline image (Ubuntu 24.04, GATK 4.6.0.0)$(RESET)\n"
	docker build \
		--file docker/Dockerfile.pipeline \
		--tag $(PIPELINE_IMAGE):0.1.0 \
		--tag $(PIPELINE_IMAGE):latest \
		--label "org.opencontainers.image.version=0.1.0" \
		--label "org.opencontainers.image.source=https://github.com/genomeforge/genome-forge" \
		.
	@printf "$(CYAN)→ Beacon API image$(RESET)\n"
	docker build \
		--file docker/Dockerfile.beacon \
		--tag $(BEACON_IMAGE):0.1.0 \
		--tag $(BEACON_IMAGE):latest \
		.
	@printf "$(CYAN)→ Reclassification daemon image$(RESET)\n"
	docker build \
		--file docker/Dockerfile.daemon \
		--tag $(DAEMON_IMAGE):0.1.0 \
		--tag $(DAEMON_IMAGE):latest \
		.
	@printf "$(GREEN)All Docker images built successfully.$(RESET)\n"

# =============================================================================
.PHONY: ci-test
ci-test: ## Run the chr22 CI test profile (matches GitHub Actions ci.yml)
	@printf "$(BOLD)Running chr22 CI test profile...$(RESET)\n"
	@if ! command -v nextflow >/dev/null 2>&1; then \
		printf "$(YELLOW)ERROR: nextflow not found. Install from https://nextflow.io$(RESET)\n"; \
		exit 1; \
	fi
	nextflow run pipelines/wgs_grch38.nf \
		-profile test,docker \
		-resume \
		--outdir results/ci_test_$$(date +%Y%m%d_%H%M%S)
	@printf "$(GREEN)chr22 CI test complete.$(RESET)\n"

# =============================================================================
.PHONY: beacon-start
beacon-start: ## Start the GA4GH Beacon v2.1.1 API server (development mode)
	@printf "$(BOLD)Starting Beacon v2.1.1 API (development mode)...$(RESET)\n"
	@if [ ! -f .env ]; then \
		printf "$(YELLOW)ERROR: .env file not found. Copy .env.example to .env and fill in values.$(RESET)\n"; \
		exit 1; \
	fi
	@# Validate required Beacon environment variables
	@source .env && \
		[ -n "$$BEACON_JWT_SECRET" ] || { \
			printf "$(YELLOW)ERROR: BEACON_JWT_SECRET not set in .env$(RESET)\n"; exit 1; \
		}
	uvicorn beacon_api.main:app \
		--reload \
		--host 0.0.0.0 \
		--port 5050 \
		--log-level info \
		--workers 1
	# For production: use gunicorn with uvicorn workers
	# gunicorn beacon_api.main:app -k uvicorn.workers.UvicornWorker -w 4

# =============================================================================
.PHONY: daemon-start
daemon-start: ## Start the ClinVar reclassification Celery daemon
	@printf "$(BOLD)Starting ClinVar reclassification daemon...$(RESET)\n"
	@if [ ! -f .env ]; then \
		printf "$(YELLOW)ERROR: .env file not found. Copy .env.example to .env and fill in values.$(RESET)\n"; \
		exit 1; \
	fi
	@source .env && \
		[ -n "$$REDIS_URL" ] || { \
			printf "$(YELLOW)ERROR: REDIS_URL not set in .env$(RESET)\n"; exit 1; \
		}
	celery \
		--app=reclassification.daemon:app \
		worker \
		--loglevel=INFO \
		--concurrency=2 \
		--queues=clinvar_diff,fhir_tasks,notifications \
		--beat \
		--scheduler=celery.beat:PersistentScheduler &
	@printf "$(GREEN)Daemon started. Check logs with: celery --app=reclassification.daemon:app inspect active$(RESET)\n"

# =============================================================================
.PHONY: clean
clean: ## Remove build artifacts, caches, and temporary files
	@printf "$(BOLD)Cleaning build artifacts...$(RESET)\n"
	@# Python build artifacts
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	find . -type f -name "*.pyo" -delete 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	@# Coverage artifacts
	rm -rf htmlcov/ coverage.xml .coverage 2>/dev/null || true
	@# Documentation build
	rm -rf docs/_build/ 2>/dev/null || true
	@# Nextflow work directory (can be very large — ask before deleting in CI)
	@if [ -d work/ ]; then \
		printf "$(YELLOW)Removing Nextflow work/ directory (can be large)...$(RESET)\n"; \
		rm -rf work/; \
	fi
	@printf "$(GREEN)Clean complete.$(RESET)\n"

# =============================================================================
.PHONY: secrets-baseline
secrets-baseline: ## Create/update the detect-secrets baseline file
	@printf "$(BOLD)Updating detect-secrets baseline...$(RESET)\n"
	detect-secrets scan \
		--exclude-files '\.env\.example$$' \
		--exclude-files '\.secrets\.baseline$$' \
		> .secrets.baseline
	@printf "$(GREEN)Baseline updated: .secrets.baseline$(RESET)\n"

# =============================================================================
.PHONY: check-prereqs
check-prereqs: ## Check all prerequisite tools are installed and at correct versions
	@printf "$(BOLD)Checking prerequisites...$(RESET)\n"
	@$(PYTHON) -c "import sys; assert sys.version_info[:2] == (3, 12), f'Python 3.12 required, got {sys.version}'" && \
		printf "$(GREEN)✓ Python $$($(PYTHON) --version)$(RESET)\n" || \
		printf "$(YELLOW)✗ Python 3.12 required$(RESET)\n"
	@git --version | grep -qE "git version ([2-9]\.[4-9]|[3-9])" && \
		printf "$(GREEN)✓ $$(git --version)$(RESET)\n" || \
		printf "$(YELLOW)✗ git ≥ 2.40 required$(RESET)\n"
	@nextflow -version 2>/dev/null | head -1 | grep -q "nextflow" && \
		printf "$(GREEN)✓ Nextflow $$(nextflow -version 2>&1 | head -1)$(RESET)\n" || \
		printf "$(YELLOW)✗ Nextflow not found$(RESET)\n"
	@docker --version 2>/dev/null && \
		printf "$(GREEN)✓ $$(docker --version)$(RESET)\n" || \
		printf "$(YELLOW)✗ Docker not found$(RESET)\n"
	@terraform --version 2>/dev/null | head -1 && \
		printf "$(GREEN)✓ Terraform found$(RESET)\n" || \
		printf "$(YELLOW)✗ Terraform not found$(RESET)\n"
	@java -version 2>&1 | head -1 && \
		printf "$(GREEN)✓ Java found$(RESET)\n" || \
		printf "$(YELLOW)✗ Java ≥ 17 required (for GATK 4.6.0.0)$(RESET)\n"
