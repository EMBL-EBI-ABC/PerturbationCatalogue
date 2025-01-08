from enum import Enum
from pydantic import BaseModel

class SortDirections(str, Enum):
    desc = "desc"
    asc = "asc"

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
