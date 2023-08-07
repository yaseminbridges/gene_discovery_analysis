from pathlib import Path

import click

from gene_discovery_analysis.filter.filter_exomiser_variants_output import filter_exomiser_results


@click.command("filter-results")
@click.option(
    "--results-dir",
    "-r",
    required=True,
    help="Path to results directory containing raw results directory.",
    type=Path,
)
@click.option(
    "--configuration",
    "-c",
    required=True,
    help="Path to analysis configurations file.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
def filter_results_command(results_dir: Path, configuration: Path, output_dir: Path):
    """Filter Exomiser variant tsv results."""
    filter_exomiser_results(results_dir, configuration, output_dir)
