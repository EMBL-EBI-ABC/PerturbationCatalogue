import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

# Register the page
dash.register_page(__name__, path="/privacy-notice")

layout = dbc.Container(
    [
        html.H1("Privacy notice for Perturbation Catalogue", className="mt-4 mb-4"),
        html.P(
            "This Privacy notice explains what personal data is collected by the specific EMBL-EBI service you are requesting, for what purposes, how it is processed, and how we keep it secure."
        ),
        html.H2(
            "Who controls your personal data and how to contact us?", className="mt-4"
        ),
        html.P("The EMBL data controller's contact details are:"),
        html.Ul(
            [
                html.Li("Rolf Apweiler and Ewan Birney, EMBL-EBI Directors"),
                html.Li("Email: data-controller@ebi.ac.uk"),
                html.Li(
                    "EMBL-EBI, Wellcome Genome Campus, CB10 1SD Hinxton, Cambridgeshire, UK."
                ),
            ]
        ),
        html.P("The EMBL Data Protection Officer's contact details are:"),
        html.Ul(
            [
                html.Li("EMBL Data Protection Officer"),
                html.Li("Email: dpo@embl.org"),
                html.Li("EMBL Heidelberg, Meyerhofstraße 1, 69117 Heidelberg, Germany"),
            ]
        ),
        html.H2(
            "Which is the lawful basis for processing your personal data?",
            className="mt-4",
        ),
        html.P(
            "Processing your personal data is necessary for our legitimate interest of delivering a service that provides broader societal benefits."
        ),
        html.H2(
            "What personal data do we collect? How do we use this personal data?",
            className="mt-4",
        ),
        html.P("We collect the following personal data:"),
        html.Ul(
            [
                html.Li("IP addresses"),
                html.Li("Date and time of a visit to the service website"),
                html.Li("Operating system"),
                html.Li("Browser"),
            ]
        ),
        html.P("We will use the personal data:"),
        html.Ul(
            [
                html.Li("To provide the user access to the service"),
                html.Li(
                    "To better understand the needs of the users and guide future improvements of the service"
                ),
                html.Li("To create anonymous usage statistics"),
                html.Li("To conduct and monitor security activities"),
            ]
        ),
        html.H2("Who will have access to your personal data?", className="mt-4"),
        html.P("The personal data will be disclosed to:"),
        html.Ul([html.Li("Authorised EMBL-EBI staff.")]),
        html.H2(
            "Will your personal data be transferred to third countries (i.e. countries not part of EU/EAA) and/or international organisations?",
            className="mt-4",
        ),
        html.P(
            "There are no personal data transfers to third countries or to international organisations."
        ),
        html.H2("How long do we keep your personal data?", className="mt-4"),
        html.P(
            "Any personal data directly obtained from you will be retained as long as the service is live, even if you stop using the service. We will keep the personal data for the minimum amount of time possible to ensure legal compliance and to facilitate internal and external audits if they arise."
        ),
        html.H2("Your rights regarding your personal data", className="mt-4"),
        html.P("You have the right to:"),
        html.Ol(
            [
                html.Li(
                    "Not be subject to decisions based solely on an automated processing of data (i.e. without human intervention) without you having your views taken into consideration."
                ),
                html.Li(
                    "Request at reasonable intervals and without excessive delay or expense, information about the personal data processed about you. Under your request we will inform you in writing about, for example, the origin of the personal data or the preservation period."
                ),
                html.Li(
                    "Request information to understand data processing activities when the results of these activities are applied to you."
                ),
                html.Li(
                    "Object at any time to the processing of your personal data unless we can demonstrate that we have legitimate reasons to process your personal data."
                ),
                html.Li(
                    "Request free of charge and without excessive delay rectification or erasure of your personal data if we have not been processing it respecting the EMBL Internal Policy for Data Protection."
                ),
            ]
        ),
        html.P(
            "It must be clarified that rights 4 and 5 are only available whenever the processing of your personal data is not necessary to:"
        ),
        html.Ol(
            [
                html.Li("Comply with a legal obligation."),
                html.Li("Perform a task carried out in the public interest."),
                html.Li("Exercise authority as a data controller."),
                html.Li(
                    "Archive for purposes in the public interest, or for historical research purposes, or for statistical purposes."
                ),
                html.Li("Establish, exercise or defend legal claims."),
            ]
        ),
        html.P("Published at: Tuesday, October 29, 2024 - 14:40", className="mt-4"),
    ],
    className="p-4",
)
