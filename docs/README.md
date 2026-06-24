# docs/

GenomeForge documentation (Sphinx + autodoc + napoleon + mermaid + furo theme).

## Structure

```
docs/
├── conf.py              Sphinx configuration
├── index.rst            Documentation root
├── setup/
│   └── prerequisites.md Tool version requirements
├── adr/                 Architecture Decision Records
│   ├── 001-why-dragmap-over-bwa-mem2.md
│   ├── 002-why-pymc-over-stan.md
│   ├── 003-why-celery-over-apscheduler.md
│   ├── 004-why-vg-giraffe-over-vg-map.md
│   ├── 005-why-alphamissense-primary-pp3-bp4.md
│   └── 006-why-pm2-supporting-not-moderate.md
├── pipeline/            Pipeline usage docs
├── api/                 Beacon API docs
└── guides/
    ├── quickstart.md
    ├── data_setup.md
    ├── aws_deployment.md
    ├── acmg_classification.md  Full ACGS 2024 v1.2 implementation guide
    └── clinvar_submission.md   NHS ClinVar submission workflow
```

## Build

```bash
make docs
# Opens docs/_build/html/index.html
```
