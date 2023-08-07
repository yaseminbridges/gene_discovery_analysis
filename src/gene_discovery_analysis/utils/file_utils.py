from pathlib import Path

import yaml

from gene_discovery_analysis.analysis_config.parse_configuration_settings import Configurations


def parse_analysis_config(analysis_config: Path) -> Configurations:
    """Read the config file."""
    with open(analysis_config, "r") as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    return Configurations.model_validate(config)


def get_variant_tsv_files(results_dir: Path) -> [Path]:
    """Get a list of paths for Exomiser variant tsv output files."""
    files = [path for path in results_dir.iterdir() if str(path).endswith("variants.tsv")]
    files.sort()
    return files
