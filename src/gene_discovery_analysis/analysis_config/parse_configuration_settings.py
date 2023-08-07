from pydantic import BaseModel, Field


class Configurations(BaseModel):
    """
    Class for defining the application.properties configurations.
    Args:
        score_name (str): Name of the score to sort prioritised variant results.
        max_rank (int): Maximum rank to be included in results.
        whitelist_filter (bool): Filter for whitelist variants
        contributing_variant_filter (bool): Filter for contributing variants.

    """

    score_name: str = Field(None)
    max_rank: int = Field(None)
    whitelist_filter: bool = Field(None)
    contributing_variant_filter: bool = Field(None)
