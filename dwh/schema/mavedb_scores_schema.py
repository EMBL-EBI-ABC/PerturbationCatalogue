mavedb_scores_schema = {
    "fields": [
        {
            "name": "accession",
            "type": "string",
            "mode": "REQUIRED"
        },
        {
            "name": "hgvs_nt",
            "type": "string",
            "mode": "NULLABLE"
        },
        {
            "name": "hgvs_splice",
            "type": "string",
            "mode": "NULLABLE"
        },
        {
            "name": "hgvs_pro",
            "type": "string",
            "mode": "NULLABLE"
        },
        {
            "name": "score",
            "type": "FLOAT",
            "mode": "NULLABLE"
        }
    ]
}