from dash import dcc, html
import dash_bootstrap_components as dbc

from .elastic_table import ElasticTable


def layout(filter_fields):
    from ..data_portal import elastic_table

    return elastic_table.layout()
