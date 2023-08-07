import click

from gene_discovery_analysis.filter.filter import filter_results_command


@click.group()
def main():
    """Exomiser results filter."""


main.add_command(filter_results_command)

if __name__ == "__main__":
    main()
