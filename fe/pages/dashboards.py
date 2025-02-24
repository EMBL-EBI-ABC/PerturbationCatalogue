import dash

from . import _iframe

dash.register_page(__name__)

layout = _iframe.layout(
    "https://lookerstudio.google.com/embed/reporting/86ab32f9-151b-4f91-87eb-060a22f2f890/page/QJFZE"
)
