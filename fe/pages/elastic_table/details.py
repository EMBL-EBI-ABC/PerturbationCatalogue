import requests

from dash import html


class ElasticTableDetails:
    def __init__(self, api_endpoint):
        self.api_endpoint = api_endpoint

    def get_detail(self, urn):
        response = requests.get(f"{self.api_endpoint}/{urn}")

        if response.status_code != 200:
            return html.Div(
                "Error: Unable to fetch data from the API.",
                className="alert alert-danger",
            )

        return response.json()
