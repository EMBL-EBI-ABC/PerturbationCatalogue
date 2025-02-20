from dash import html
import requests


def layout(urn):
    api_url = f"https://perturbation-catalogue-be-959149465821.europe-west2.run.app/mavedb/search/{urn}"
    response = requests.get(api_url)

    if response.status_code != 200:
        return html.Div(
            "Error: Unable to fetch data from the API.", className="alert alert-danger"
        )

    data = response.json().get("results", [{}])[0]

    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H4(
                                        data.get("title", "N/A"), className="card-title"
                                    ),
                                    html.P(
                                        data.get("shortDescription", "N/A"),
                                        className="card-text",
                                    ),
                                    html.Div(
                                        [
                                            html.A(
                                                "View on MaveDB",
                                                href=f"https://www.mavedb.org/score-sets/{urn}/",
                                                className="btn btn-primary mb-3",
                                            )
                                        ]
                                    ),
                                    html.P(
                                        [html.Strong("URN: "), data.get("urn", "N/A")],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Sequence Type: "),
                                            data.get("sequenceType", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Gene Name: "),
                                            data.get("geneName", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Gene Category: "),
                                            data.get("geneCategory", "N/A"),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Publication: "),
                                            html.A(
                                                data.get("publicationUrl", "N/A"),
                                                href=data.get("publicationUrl", "#"),
                                                target="_blank",
                                            ),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Publication Year: "),
                                            str(data.get("publicationYear", "N/A")),
                                        ],
                                        className="card-text",
                                    ),
                                    html.P(
                                        [
                                            html.Strong("Number of Variants: "),
                                            str(data.get("numVariants", "N/A")),
                                        ],
                                        className="card-text",
                                    ),
                                ],
                                className="card-body",
                            )
                        ],
                        className="card",
                    )
                ],
                className="col-md-8 mx-auto",
            )
        ],
        className="container mt-4",
    )
