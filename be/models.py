from pydantic import BaseModel, Field
from typing import Literal


class MaveDBData(BaseModel):
    urn: str
    title: str
    shortDescription: str
    sequenceType: str
    geneName: str
    geneCategory: str
    publicationUrl: str
    publicationYear: int
    numVariants: int


class MaveDBResponse(BaseModel):
    total: int
    start: int
    size: int
    results: list[MaveDBData]
    aggregations: dict[str, dict]


class MaveDBDetailsResponse(BaseModel):
    results: list[MaveDBData]


class SearchParams(BaseModel):
    model_config = {
        "populate_by_name": True,
        "extra": "forbid"
    }

    q: str | None = Field(None, description="Search query string")
    publication_year: str | None = Field(None,
                                         description="PublicationYear query",
                                         alias="publicationYear")
    gene_category: str | None = Field(None,
                                      description="GeneCategory query",
                                      alias="geneCategory")
    sequence_type: str | None = Field(None,
                                      description="SequenceType query",
                                      alias="sequenceType")
    start: int = Field(0, description="Starting point of the results")
    size: int = Field(10, gt=0, description="Number of results per page")
    sort_field: str | None = Field("publicationYear", description="Sort field")
    sort_order: Literal["desc", "asc"] = "desc"
