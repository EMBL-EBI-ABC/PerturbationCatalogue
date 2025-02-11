from pydantic import BaseModel, Field
from typing import Generic, List, Literal, TypeVar

T = TypeVar("T")  # Datasource data type
A = TypeVar("A")  # Datasource aggregation type


# Generic aggregation classes.

class AggregationBucket(BaseModel):
    key: int|str
    doc_count: int

class Aggregation(BaseModel):
    doc_count_error_upper_bound: int
    sum_other_doc_count: int
    buckets: list[AggregationBucket]

def get_list_of_aggregations(aggregation_class):
    return sorted(aggregation_class.schema()["properties"].keys())

# Generic Elastic response classes.

class ElasticResponse(BaseModel, Generic[T, A]):
    total: int
    start: int
    size: int
    results: List[T]
    aggregations: A

class ElasticDetailsResponse(BaseModel, Generic[T]):
    results: List[T]


# Base Elastic query class.

class SearchParams(BaseModel):
    model_config = {
        "populate_by_name": True,
        "extra": "forbid",
    }
    # Basic query parameters.
    q: str | None = Field(None, description="Search query string")
    start: int = Field(0, description="Starting point of the results")
    size: int = Field(10, gt=0, description="Number of results per page")
    # No sorting by default, child classes can override this.
    sort_field: str | None = None
    sort_order: Literal["desc", "asc"] = "asc"


# MaveDB.

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

class MaveDBAggregationResponse(BaseModel):
    publicationYear: Aggregation
    sequenceType: Aggregation
    geneCategory: Aggregation

class MaveDBSearchParams(SearchParams):
    sort_field: str | None = Field("publicationYear", description="Sort field")
    sort_order: Literal["desc", "asc"] = "desc"
    publicationYear: str | None = Field(
        None, description="PublicationYear query", alias="publicationYear"
    )
    geneCategory: str | None = Field(
        None, description="GeneCategory query", alias="geneCategory"
    )
    sequenceType: str | None = Field(
        None, description="SequenceType query", alias="sequenceType"
    )
