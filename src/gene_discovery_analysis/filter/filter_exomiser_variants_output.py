from pathlib import Path

import pandas as pd

from gene_discovery_analysis.utils.file_utils import get_variant_tsv_files, parse_analysis_config


def read_variant_tsv(variant_tsv: Path) -> pd.DataFrame:
    """Read variant tsv output from Exomiser."""
    return pd.read_csv(variant_tsv, delimiter="\t")


class FilterExomiserVariantsOutput:
    rank: int = 0
    current_score: float = float("inf")
    count: int = 0

    @staticmethod
    def sort_results(variants_df: pd.DataFrame, score_name: str) -> pd.DataFrame:
        """Sort results."""
        return variants_df.sort_values(by=score_name, ascending=False)

    @staticmethod
    def rank_scores(variants_df: pd.DataFrame) -> int:
        """Add ranks to a result, equal scores are given the same rank e.g., 1,1,3."""
        current_rank = 0
        current_score = None
        num_changes = 0
        ranks = []
        for index, row in variants_df.iterrows():
            if current_score is None or row['EXOMISER_GENE_COMBINED_SCORE'] != current_score:
                current_rank += 1
                current_rank += num_changes
                num_changes = 0
                current_score = row['EXOMISER_GENE_COMBINED_SCORE']
            if index > 0 and row['EXOMISER_GENE_COMBINED_SCORE'] == variants_df.loc[
                index - 1, 'EXOMISER_GENE_COMBINED_SCORE']:
                if row['GENE_SYMBOL'] != variants_df.loc[index - 1, 'GENE_SYMBOL']:
                    num_changes += 1
                elif row['GENE_SYMBOL'] == variants_df.loc[index - 1, 'GENE_SYMBOL'] and row['MOI'] != variants_df.loc[
                    index - 1, 'MOI']:
                    num_changes += 1
            ranks.append(current_rank)
        variants_df["NEW_RANK"] = ranks
        return variants_df

    @staticmethod
    def filter_max_rank(variants_df: pd.DataFrame, max_rank: int) -> pd.DataFrame:
        """Filter variant results for maximum rank."""
        if max_rank != 0:
            return variants_df[variants_df["NEW_RANK"] <= max_rank]
        return variants_df

    @staticmethod
    def get_transcript_id(variants_df: pd.DataFrame) -> pd.DataFrame:
        """Get transcript ID for variant."""
        copied_df = variants_df.copy()
        copied_df["TRANSCRIPT_ID"] = copied_df["HGVS"].str.split(":").str[1]
        return copied_df

    @staticmethod
    def filter_whitelist(variants_df: pd.DataFrame) -> pd.DataFrame:
        """Filter variants that are defined as a whitelist variant."""
        return variants_df[variants_df["WHITELIST_VARIANT"] != 0]

    @staticmethod
    def filter_contributing_variant(variants_df: pd.DataFrame) -> pd.DataFrame:
        """Filter variants that are defined as a contributing variant."""
        return variants_df[variants_df["CONTRIBUTING_VARIANT"] != 0]

    def filter_output(
            self,
            variants_df: pd.DataFrame,
            max_rank: int,
            whitelist_filter: bool,
            contributing_variant_filter: bool,
    ) -> pd.DataFrame:
        """Return filtered prioritised variant output."""
        filtered_variants = self.filter_max_rank(variants_df, max_rank)
        if whitelist_filter is True:
            filtered_variants = self.filter_whitelist(filtered_variants)
        if contributing_variant_filter is True:
            filtered_variants = self.filter_contributing_variant(filtered_variants)
        return filtered_variants

    def create_filtered_output(
            self,
            variants_df: pd.DataFrame,
            score_name: str,
            max_rank: int,
            whitelist_filter: bool,
            contributing_variant_filter: bool,
    ) -> pd.DataFrame:
        """Create the filtered prioritised variant output with annotated transcript ID."""
        sorted_variants_df = self.rank_scores(
            self.sort_results(variants_df, score_name)
        )
        filtered_variants = self.filter_output(
            sorted_variants_df, max_rank, whitelist_filter, contributing_variant_filter
        )
        return self.get_transcript_id(filtered_variants)


def write_filtered_output(
        variant_tsv_file: Path, output_dir: Path, filtered_variants_df: pd.DataFrame
) -> None:
    """Write filtered Exomiser variant tsv output to a tsv file."""
    filtered_variants_df.to_csv(
        output_dir.joinpath(f"{variant_tsv_file.stem}-filtered.variants.tsv"), sep="\t"
    )


def filter_exomiser_results(results_dir: Path, analysis_configuration: Path, output_dir: Path):
    """Filter a directory of Exomiser variant tsv results."""
    output_dir.mkdir(parents=True, exist_ok=True)
    configurations = parse_analysis_config(analysis_configuration)
    for result in get_variant_tsv_files(results_dir.joinpath("raw_results")):
        variants_df = read_variant_tsv(result)
        filtered_variants_df = FilterExomiserVariantsOutput().create_filtered_output(
            variants_df=variants_df,
            score_name=configurations.score_name,
            max_rank=configurations.max_rank,
            whitelist_filter=configurations.whitelist_filter,
            contributing_variant_filter=configurations.contributing_variant_filter,
        )
        write_filtered_output(result, output_dir, filtered_variants_df)
