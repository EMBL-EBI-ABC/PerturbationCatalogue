import dash

PAGE_ORDER = ("Data Portal", "Dashboards", "API documentation", "About")


def get_all_pages_dict():
    return {x["name"]: x for x in dash.page_registry.values()}


def get_pages(
    include_home=True, include_other=True, require_button=False, require_icon=False
):
    all_pages_dict = get_all_pages_dict()
    pages = []
    if include_home:
        pages.append(all_pages_dict["Home"])
    if include_other:
        pages.extend([all_pages_dict[p] for p in PAGE_ORDER])
    if require_button:
        pages = [p for p in pages if p.get("button")]
    if require_icon:
        pages = [p for p in pages if p.get("icon")]
    return pages


def get_home_page():
    return get_all_pages_dict()["Home"]


REQUEST_DATASET_PAGE = {
    "name": "Request dataset",
    "description": "Use this link to request a dataset to be ingested into the Perturbation Catalogue.",
    "href": "https://pollunit.com/polls/zl6y1cje-smx6jikfvfwaw",
    "icon": "bi-box-arrow-up-right",
    "button": "Request dataset",
}
