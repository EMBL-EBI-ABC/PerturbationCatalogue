"""Pydantic models used by the backend API."""

from typing import Optional, List, Dict, Any

from pydantic import BaseModel, Field


class SearchRequest(BaseModel):
    query: Optional[str] = Field(
        None, description="Search query for perturbed_target_symbol"
    )
    filters: Optional[Dict[str, List[str]]] = Field(
        None,
        description=(
            "Filter by facet values. Keys: licenses_tested, data_modalities, tissues_tested, "
            "cell_types_tested, cell_lines_tested, sex_tested, "
            "developmental_stages_tested, diseases_tested"
        ),
    )
    page: int = Field(1, ge=1, description="Page number (1-indexed)")
    size: int = Field(6, ge=1, le=100, description="Number of results per page")


class FacetValue(BaseModel):
    value: str
    count: int


class Facets(BaseModel):
    licenses_tested: List[FacetValue]
    data_modalities: List[FacetValue]
    tissues_tested: List[FacetValue]
    cell_types_tested: List[FacetValue]
    cell_lines_tested: List[FacetValue]
    sex_tested: List[FacetValue]
    developmental_stages_tested: List[FacetValue]
    diseases_tested: List[FacetValue]


class SearchResponse(BaseModel):
    total: int
    page: int
    size: int
    total_pages: int
    results: List[Dict[str, Any]]
    facets: Facets


class SummaryTopEntry(BaseModel):
    value: str
    n_datasets: int


class LandingPageSummary(BaseModel):
    n_datasets: int
    n_experiments: int
    min_year: int
    max_year: int
    n_targets: int
    n_tissues: int
    n_cell_types: int
    n_cell_lines: int
    n_diseases: int
    top_modalities: List[SummaryTopEntry]
    top_tissues: List[SummaryTopEntry]
    top_cell_types: List[SummaryTopEntry]
    top_cell_lines: List[SummaryTopEntry]
    top_perturbation_types: List[SummaryTopEntry]
    top_diseases: List[SummaryTopEntry]
    top_sexes: List[SummaryTopEntry]
    top_dev_stages: List[SummaryTopEntry]
