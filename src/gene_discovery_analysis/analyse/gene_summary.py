from pathlib import Path
from statistics import mean

import click
import pandas as pd

from gene_discovery_analysis.filter.filter_exomiser_variants_output import (
    FilterExomiserVariantsOutput,
    read_variant_tsv,
    write_filtered_output,
)
from gene_discovery_analysis.utils.file_utils import get_variant_tsv_files, parse_analysis_config


class GeneSummaryCreator:
    gene_summary = {}
    number_of_patients = 0

    def add_reciprocal_rank_for_gene(self, row: pd.Series):
        """Add reciprocal rank for a gene."""
        gene_symbol = row["GENE_SYMBOL"]
        if gene_symbol not in self.gene_summary:
            self.gene_summary[gene_symbol] = []
        self.gene_summary[gene_symbol].append(1 / row["NEW_RANK"])

    def add_ranks_from_filtered_df(self, filtered_variants: pd.DataFrame):
        unique_genes = filtered_variants.drop_duplicates(subset=["GENE_SYMBOL"])
        for _index, row in unique_genes.iterrows():
            self.add_reciprocal_rank_for_gene(row)

    @staticmethod
    def calculate_mrr(reciprocal_ranks: [int]) -> float:
        return mean(reciprocal_ranks)

    def write_gene_summary(self, output_dir: Path) -> None:
        with open(output_dir.joinpath("filtered_gene_summary.tsv"), "w") as output:
            output.write("GENE_SYMBOL\tMRR\tNUMBER_OF_SAMPLES\n")
            for gene, ranks in self.gene_summary.items():
                output.write(f"{gene}\t{self.calculate_mrr(ranks)}\t{len(ranks)}\n")
        output.close()


def analyse_filtered_variants(results_dir: Path, analysis_configuration: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    configurations = parse_analysis_config(analysis_configuration)
    gene_summary = GeneSummaryCreator()
    for result in get_variant_tsv_files(results_dir.joinpath("raw_results")):
        variants_df = read_variant_tsv(result)
        filtered_variants_df = FilterExomiserVariantsOutput().create_filtered_output(
            variants_df=variants_df,
            score_name=configurations.score_name,
            max_rank=configurations.max_rank,
            whitelist_filter=configurations.whitelist_filter,
            contributing_variant_filter=configurations.contributing_variant_filter,
        )
        gene_summary.add_ranks_from_filtered_df(filtered_variants_df)
        write_filtered_output(result, output_dir, filtered_variants_df)
    gene_summary.write_gene_summary(output_dir)


@click.command("gene-summary")
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
def gene_summary_command(results_dir: Path, configuration: Path, output_dir: Path):
    """Filter Exomiser variant tsv results and get gene summary."""
    analyse_filtered_variants(results_dir, configuration, output_dir)
