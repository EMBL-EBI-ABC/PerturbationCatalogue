import asyncio
from unittest.mock import Mock

import data_query
from fastapi import HTTPException
import pytest


class _Acquire:
    def __init__(self, connection):
        self.connection = connection

    async def __aenter__(self):
        return self.connection

    async def __aexit__(self, *_):
        return False


class _Connection:
    async def fetchval(self, query, dataset_id):
        assert "COUNT(*)" in query
        assert dataset_id == "demo"
        return 18

    async def fetch(self, query, *params):
        self.query = query
        self.params = params
        return [
            {
                "perturbed_target_ensg": "ENSG1",
                "term": "Pathway A",
                "es": 0.5,
                "nes": 1.2,
                "pval": 0.01,
                "sidak": 0.03,
                "fdr": 0.02,
                "geneset_size": 12,
                "leading_edge": ["A", "B"],
                "cell_type": "T cell",
            }
        ]


class _Pool:
    def __init__(self, connection):
        self.connection = connection

    def acquire(self):
        return _Acquire(self.connection)


def test_dataset_gsea_page_uses_stable_order_and_offset(monkeypatch):
    connection = _Connection()
    monkeypatch.setitem(data_query.db_pools, "pg", _Pool(connection))

    response = asyncio.run(
        data_query.get_dataset_perturb_seq_gsea("demo", limit=5, offset=10)
    )

    assert response["total_rows_count"] == 18
    assert response["results"][0]["leading_edge"] == ["A", "B"]
    assert "ORDER BY sidak ASC NULLS LAST" in connection.query
    assert "ctid ASC" in connection.query
    assert "LIMIT $2 OFFSET $3" in connection.query
    assert connection.params == ("demo", 5, 10)


def test_target_modal_gsea_fetch_keeps_its_50_row_cap():
    connection = _Connection()
    rows = asyncio.run(
        data_query._fetch_perturb_seq_gsea(
            connection, "demo", {"perturbation_gene_name": "ENSG1"}
        )
    )

    assert rows
    assert "LIMIT 50" in connection.query


def test_gsea_routes_do_not_replace_target_modal_route():
    paths = {route.path for route in data_query.router.routes}

    assert "/v1/perturb-seq-gsea" in paths
    assert "/v1/perturb-seq/{dataset_id}/gsea" in paths
    assert "/v1/perturb-seq/{dataset_id}/gsea/download" in paths


def test_gsea_download_rejects_unknown_format_and_redirects(monkeypatch):
    monkeypatch.setattr(
        data_query,
        "_release_signed_url",
        lambda modality, dataset_id, download_format, artifact=None: (
            f"https://storage.example/{modality}/{dataset_id}.{artifact}.{download_format}"
        ),
    )

    with pytest.raises(HTTPException) as error:
        asyncio.run(data_query.download_dataset_perturb_seq_gsea("demo", "json"))
    assert error.value.status_code == 400

    response = asyncio.run(
        data_query.download_dataset_perturb_seq_gsea("demo", "parquet")
    )
    assert response.status_code == 307
    assert response.headers["location"].endswith("/perturb-seq/demo.gsea.parquet")


def test_gsea_release_signing_uses_sibling_artifact(monkeypatch):
    import google.auth
    from google.cloud import storage

    credentials = Mock(valid=True, token="token", service_account_email="service")
    blob = Mock(exists=Mock(return_value=True))
    blob.generate_signed_url.return_value = "https://storage.example/signed"
    client = Mock()
    client.bucket.return_value.blob.return_value = blob
    monkeypatch.setenv("RELEASE_BUCKET", "release-bucket")
    monkeypatch.setattr(google.auth, "default", lambda: (credentials, "project"))
    monkeypatch.setattr(storage, "Client", lambda **_: client)

    result = data_query._release_signed_url(
        "perturb-seq", "demo", "parquet", artifact="gsea"
    )

    assert result == "https://storage.example/signed"
    client.bucket.return_value.blob.assert_called_once_with(
        "perturb-seq/demo.gsea.parquet"
    )
