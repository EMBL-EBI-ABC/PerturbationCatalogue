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
    # TODO: undef values are strings, convert it to None
    publicationYear: int|str
    numVariants: int|str

class AggregationBucket(BaseModel):
    key: int|str
    doc_count: int

class Aggregation(BaseModel):
    doc_count_error_upper_bound: int
    sum_other_doc_count: int
    buckets: list[AggregationBucket]

class AggregationResponse(BaseModel):
    publicationYear: Aggregation
    sequenceType: Aggregation
    geneCategory: Aggregation


class MaveDBResponse(BaseModel):
    total: int
    start: int
    size: int
    results: list[MaveDBData]
    aggregations: AggregationResponse


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
