# ADR-003: Why Celery over APScheduler for the ClinVar reclassification daemon

**Status:** Accepted  
**Date:** 2026-06-22  
**Deciders:** ClaritySeq core team  
**Category:** Task scheduling / background processing

---

## Context

ClaritySeq includes a **ClinVar reclassification daemon** (`reclassification/daemon.py`)
that:

1. Downloads weekly ClinVar XML diffs every Monday
2. Identifies variants that have changed classification since last checked
3. Matches changed variants against reported patient variants in the PostgreSQL database
4. Creates FHIR R4 Task resources (Genomics Reporting IG v3.0.0) for labs to action
5. Sends Slack/email notifications to clinical teams

Two task scheduling frameworks were evaluated:

1. **Celery v5.3.x with Redis 7** — distributed task queue and beat scheduler
2. **APScheduler v3.x** — lightweight Python in-process job scheduler

---

## Decision

**Celery v5.3.x with Redis 7 as broker is the reclassification daemon framework.**

---

## Rationale

### 1. Celery provides distributed, fault-tolerant task execution

The reclassification workflow involves multiple discrete tasks:

```
clinvar_diff (download + parse) → match_variants → create_fhir_tasks → notify
```

Celery executes these as a chain, with:
- **Retry semantics**: Failed tasks are retried with exponential backoff
- **Result tracking**: Task results stored in Redis for audit logging
- **Priority queues**: Urgent reclassifications (P/LP → B/LB) prioritised over
  informational updates

APScheduler is a single-process scheduler. If the process crashes mid-task (e.g.,
during a large ClinVar XML download), the task state is lost.

### 2. Redis 7 is already required for caching

ClaritySeq already uses Redis 7 for:
- API response caching (Beacon v2.1.1 frequency queries)
- ClinVar annotation caching (avoid repeated API calls)
- FHIR server response caching

Using Redis as the Celery broker adds zero new infrastructure dependencies.

### 3. Celery Beat provides production-grade cron scheduling

`celery beat` runs scheduled tasks on a cron-like schedule. The weekly ClinVar
diff is scheduled as:

```python
CELERY_BEAT_SCHEDULE = {
    "clinvar-weekly-diff": {
        "task": "reclassification.daemon.run_clinvar_diff",
        "schedule": crontab(hour=6, minute=0, day_of_week=1),
        # Every Monday at 06:00 UTC — ClinVar releases weekly updates on Mondays
    },
}
```

Celery beat persists schedule state to Redis, surviving process restarts.

### 4. Auditable task history

ClinVar reclassification is a clinical safety function. ACGS 2024 v1.2 requires
an audit trail of:
- When reclassification events were detected
- Which patients were affected
- When FHIR Tasks were created
- Whether notifications were sent

Celery's task result backend (Redis) stores the full task execution history,
which `reclassification/audit_logger.py` reads for compliance reporting.

APScheduler has no built-in result backend or audit trail.

### 5. NHS England operational requirements

NHS GMS labs must demonstrate continuous monitoring of ClinVar reclassifications
(ACGS 2024 §8.3). The daemon must:
- Run continuously (not just when the main API server is up)
- Recover automatically from crashes
- Alert on-call teams if the daemon has not run for > 72 hours

Celery workers are standard production services with health-check endpoints:
```bash
celery --app=reclassification.daemon:app inspect ping
```

APScheduler provides no health-check endpoint.

---

## Architecture

```
Redis 7 (broker + result backend)
    ↑
Celery Beat (scheduler) → Celery Worker × 2
                                ↓
    clinvar_diff_queue: ClinVarDiffTask
    fhir_tasks_queue:   FHIRTaskCreatorTask
    notifications_queue: SlackNotifierTask, EmailNotifierTask
```

The daemon is started with:
```bash
celery --app=reclassification.daemon:app worker \
    --beat --loglevel=INFO \
    --queues=clinvar_diff,fhir_tasks,notifications \
    --concurrency=2
```

---

## Consequences

**Positive:**
- Fault-tolerant: tasks survive worker crashes
- Distributed: scale workers independently
- Auditable: full task history in Redis
- Production-grade health checks
- No additional infrastructure (Redis already required)

**Negative:**
- More complex than APScheduler (worker process, beat process, Redis)
- Celery can be difficult to debug locally (use `celery --app=... events`)

**Mitigation:**
- `make daemon-start` wraps the Celery invocation
- Flower dashboard (optional): `pip install flower && celery flower`

---

## Alternatives Considered

| Framework | Reason Not Chosen |
|-----------|------------------|
| APScheduler v3.x | No result backend; single-process; no fault tolerance |
| cron (system) | Not Python-native; no retry logic; no audit trail |
| AWS EventBridge + Lambda | Vendor lock-in; adds complexity; not needed for weekly job |
| Airflow | Massive overhead for a single weekly task |
| Prefect / Dagster | Overkill; adds new infrastructure |
| asyncio tasks | No persistence across restarts; no audit trail |

---

## References

- Celery documentation: <https://docs.celeryq.dev/>
- ACGS 2024 v1.2 §8.3: ClinVar reclassification monitoring requirements
- FHIR Genomics Reporting IG v3.0.0: Recontact Task profile
- Redis 7 documentation: <https://redis.io/docs/>
