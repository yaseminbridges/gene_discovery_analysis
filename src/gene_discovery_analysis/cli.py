import click

from gene_discovery_analysis.analyse.gene_summary import gene_summary_command
from gene_discovery_analysis.filter.filter import filter_results_command


@click.group()
def main():
    """Exomiser results filter."""


main.add_command(filter_results_command)
main.add_command(gene_summary_command)
if __name__ == "__main__":
    main()
