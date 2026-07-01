"""
bayesacmg.cli
=============

Click command-line interface for BayesACMG variant classification.

Commands:
    classify  — Classify a single variant from command-line arguments.
    classify-vcf — (future) Classify all variants in a VCF file.
    show-spec — Display the VCEP specification for a gene.

References:
    Richards et al. 2015 PMID:25741868
    ClinGen SVI Working Group 2024
    ACGS 2024 v1.2 §5
"""

from __future__ import annotations

import asyncio
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from bayesacmg.combinations import classify_variant
from bayesacmg.models import GeneData, TranscriptData, VariantInput
from bayesacmg.rules import benign as benign_rules
from bayesacmg.rules import pathogenic as path_rules
from bayesacmg.vcep_client import get_vcep_spec

console = Console()

# ---------------------------------------------------------------------------
# Classification colour mapping
# ---------------------------------------------------------------------------

_CLASSIFICATION_COLOURS: dict[str, str] = {
    "Pathogenic": "bold red",
    "Likely_Pathogenic": "red",
    "VUS": "yellow",
    "Likely_Benign": "green",
    "Benign": "bold green",
}


def _colour_classification(classification: str) -> str:
    """Wrap a classification string in Rich markup for coloured output.

    Args:
        classification: One of the five ACMG/AMP classification strings.

    Returns:
        Rich markup string, e.g. ``"[bold red]Pathogenic[/bold red]"``.
    """
    colour = _CLASSIFICATION_COLOURS.get(classification, "white")
    return f"[{colour}]{classification}[/{colour}]"


def _quick_classify(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str,
    gnomad_af: float | None,
    variant_type: str,
    alphamissense: float | None,
    is_de_novo: bool,
    lof_mechanism: bool,
    is_mane_select: bool,
) -> Any:
    """Run a quick classification with minimal inputs using default rules.

    Args:
        chrom: Chromosome string.
        pos: Genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        gene: Gene symbol.
        gnomad_af: gnomAD v4.1 allele frequency or None.
        variant_type: Variant type string.
        alphamissense: AlphaMissense score or None.
        is_de_novo: True if confirmed de novo.
        lof_mechanism: True if LoF is the disease mechanism.
        is_mane_select: True if the variant is on the MANE Select transcript.

    Returns:
        ClassificationResult from combine_and_classify.
    """
    variant = VariantInput(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_type=variant_type,
        gnomad_af=gnomad_af,
        alphamissense_score=alphamissense,
        is_de_novo=is_de_novo,
        gene_symbol=gene,
    )
    transcript = TranscriptData(
        transcript_id="MANE_SELECT",
        is_mane_select=is_mane_select,
        gene_symbol=gene,
        consequence=variant_type,
    )
    gene_data = GeneData(
        gene_symbol=gene,
        lof_is_disease_mechanism=lof_mechanism,
    )

    rules = []

    # Evaluate core pathogenic rules
    rules.append(path_rules.rule_pvs1(variant, transcript, gene_data))
    rules.append(path_rules.rule_pm2(variant))
    rules.append(path_rules.rule_ps2(variant))
    rules.append(path_rules.rule_pp3(variant, alphamissense, None))

    # Evaluate core benign rules
    rules.append(benign_rules.rule_ba1(variant))
    rules.append(benign_rules.rule_bp4(variant, alphamissense, None))

    return classify_variant(rules, variant)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


@click.group()
@click.version_option(package_name="bayesacmg")
def main() -> None:
    """BayesACMG — Bayesian ACMG/AMP variant classifier.

    Implements ACGS 2024 v1.2, ClinGen SVI 2024, and all 28 ACMG/AMP criteria.
    PM2 is applied at Supporting weight (1 pt) per ClinGen SVI 2024.
    """


@main.command("classify")
@click.option("--chrom", required=True, help="Chromosome (e.g. chr17)")
@click.option(
    "--pos", required=True, type=int, help="Genomic position (GRCh38, 1-based)"
)
@click.option("--ref", required=True, help="Reference allele")
@click.option("--alt", required=True, help="Alternate allele")
@click.option("--gene", required=True, help="HGNC gene symbol")
@click.option(
    "--gnomad-af", default=None, type=float, help="gnomAD v4.1 allele frequency"
)
@click.option(
    "--variant-type",
    default="snv",
    type=click.Choice(
        [
            "snv",
            "indel",
            "frameshift",
            "nonsense",
            "splice_canonical",
            "splice_region",
            "start_loss",
            "stop_loss",
            "inframe_insertion",
            "inframe_deletion",
            "synonymous",
        ]
    ),
    help="Variant type",
)
@click.option(
    "--alphamissense", default=None, type=float, help="AlphaMissense score (0-1)"
)
@click.option(
    "--de-novo", is_flag=True, default=False, help="Confirmed de novo variant"
)
@click.option(
    "--lof-mechanism", is_flag=True, default=False, help="LoF is disease mechanism"
)
@click.option(
    "--mane-select",
    is_flag=True,
    default=True,
    help="Variant on MANE Select transcript",
)
@click.option(
    "--json-output", is_flag=True, default=False, help="Output result as JSON"
)
def classify_cmd(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str,
    gnomad_af: float | None,
    variant_type: str,
    alphamissense: float | None,
    de_novo: bool,
    lof_mechanism: bool,
    mane_select: bool,
    json_output: bool,
) -> None:
    """Classify a single variant using ACMG/AMP 2024 criteria.

    \b
    Example:
        bayesacmg classify --chrom chr17 --pos 43057103 --ref G --alt A \\
            --gene BRCA1 --gnomad-af 0 --variant-type nonsense \\
            --lof-mechanism --mane-select

    \b
    References:
        Richards et al. 2015 PMID:25741868
        ClinGen SVI 2024 (PM2 at Supporting weight)
        ACGS 2024 v1.2
    """
    result = _quick_classify(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene=gene,
        gnomad_af=gnomad_af,
        variant_type=variant_type,
        alphamissense=alphamissense,
        is_de_novo=de_novo,
        lof_mechanism=lof_mechanism,
        is_mane_select=mane_select,
    )

    if json_output:
        import json as json_mod

        out: dict[str, Any] = {
            "variant": f"{chrom}:{pos}:{ref}>{alt}",
            "gene": gene,
            "classification": result.classification,
            "total_points": result.total_points,
            "novel_combination": result.novel_combination,
            "rules_applied": [
                {"rule_id": r.rule_id, "strength": r.strength.value, "points": r.points}
                for r in result.rules_applied
            ],
        }
        click.echo(json_mod.dumps(out, indent=2))
        return

    # Rich formatted output
    title = f"BayesACMG Classification: " f"{chrom}:{pos} {ref}>{alt} ({gene})"
    table = Table(title="Applied Rules", show_header=True, header_style="bold cyan")
    table.add_column("Rule", style="cyan", width=12)
    table.add_column("Strength", width=16)
    table.add_column("Points", justify="right", width=8)
    table.add_column("Evidence", overflow="fold")

    for rule in result.rules_applied:
        table.add_row(
            rule.rule_id,
            rule.strength.value,
            f"{rule.points:+d}",
            "; ".join(rule.evidence_items[:1]),  # first evidence item
        )

    console.print(Panel(table, title=title, border_style="blue"))

    classification_str = _colour_classification(result.classification)
    console.print(
        f"\nTotal points: [bold]{result.total_points:+d}[/bold] → {classification_str}"
    )

    if result.novel_combination:
        console.print(
            f"[yellow]Novel combination applied: {result.novel_combination}[/yellow]"
        )


@main.command("show-spec")
@click.argument("gene_symbol")
def show_spec_cmd(gene_symbol: str) -> None:
    """Display the ClinGen VCEP specification for a gene.

    \b
    Example:
        bayesacmg show-spec BRCA1
    """
    spec = asyncio.run(get_vcep_spec(gene_symbol))

    table = Table(title=f"VCEP Specification: {gene_symbol}", show_header=True)
    table.add_column("Field", style="cyan")
    table.add_column("Value")
    table.add_row("VCEP Name", spec.vcep_name or "(none)")
    table.add_row("VCEP ID", spec.vcep_id or "(none)")
    table.add_row("PM2 Weight", spec.pm2_weight)
    table.add_row(
        "PP3 AlphaMissense threshold",
        (
            str(spec.pp3_threshold_alphamissense)
            if spec.pp3_threshold_alphamissense
            else "default (0.564)"
        ),
    )
    table.add_row(
        "BP4 AlphaMissense threshold",
        (
            str(spec.bp4_threshold_alphamissense)
            if spec.bp4_threshold_alphamissense
            else "default (0.340)"
        ),
    )
    console.print(table)

    if spec.pm2_weight == "moderate":
        console.print(
            "[yellow]Note: this VCEP overrides PM2 to Moderate (2 pts). "
            "Default ClinGen SVI 2024 is Supporting (1 pt).[/yellow]"
        )
