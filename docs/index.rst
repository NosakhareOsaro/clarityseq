GenomeForge Documentation
=========================

**GenomeForge** is a production-grade, research-novel whole-genome sequencing (WGS)
clinical variant interpretation platform.

.. image:: https://img.shields.io/badge/ACGS-2024%20v1.2-blue
   :alt: ACGS 2024 v1.2 compliant

.. image:: https://img.shields.io/badge/Python-3.12-blue
   :alt: Python 3.12

.. image:: https://img.shields.io/badge/GATK-4.6.0.0-orange
   :alt: GATK 4.6.0.0

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   setup/prerequisites
   guides/quickstart
   guides/data_setup
   guides/aws_deployment

.. toctree::
   :maxdepth: 2
   :caption: Pipeline

   guides/acmg_classification
   guides/clinvar_submission

.. toctree::
   :maxdepth: 2
   :caption: Architecture decisions

   adr/001-why-dragmap-over-bwa-mem2
   adr/002-why-pymc-over-stan
   adr/003-why-celery-over-apscheduler
   adr/004-why-vg-giraffe-over-vg-map
   adr/005-why-alphamissense-primary-pp3-bp4
   adr/006-why-pm2-supporting-not-moderate

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/bayesacmg
   api/beacon
   api/annotation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
