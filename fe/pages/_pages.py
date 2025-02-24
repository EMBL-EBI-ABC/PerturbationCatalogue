from collections import namedtuple

# Define app pages.
Page = namedtuple("Page", ["name", "selector", "button", "description", "icon"])
pages = [
    Page(
        name="Home",
        selector="home",
        button=None,
        description=None,
        icon="bi-house",
    ),
    Page(
        name="Data Portal",
        selector="data-portal",
        button="Open Data Portal",
        description="The Data Portal allows users to sort and filter metadata using a set of predefined filters, and it also has free-text search capabilities.",
        icon="bi-table",
    ),
    Page(
        name="Dashboards",
        selector="dashboards",
        button="Explore Dashboards",
        description="The Dashboards tab provides a visual overview of the existing data.",
        icon="bi-bar-chart-line",
    ),
    Page(
        name="Data Analytics",
        selector="data-analytics",
        button="Data Analytics",
        description="The Data Analytics tab allows users to search the data using the Data Warehouse.",
        icon="bi-graph-up",
    ),
    Page(
        name="About",
        selector="about",
        button="Contact Us",
        description="Here you can find more information about how to get help or propose new features.",
        icon="bi-question-circle",
    ),
]
