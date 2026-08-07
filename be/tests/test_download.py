"""Local integration tests for streamed dataset downloads."""

import asyncio

import data_query
import pytest
from starlette.requests import Request


DATASET_ID = "replogle_2022_rpe1_essential_normalized"


def test_download_streams_csv_from_postgres(run_with_dev_db):
    async def download():
        response = await data_query.download_dataset_data(
            "perturb-seq",
            DATASET_ID,
            Request({"type": "http", "query_string": b""}),
            limit=3,
            offset=0,
        )
        chunks = []
        async for chunk in response.body_iterator:
            chunks.append(chunk)
        return b"".join(chunks)

    content = run_with_dev_db(download)
    lines = content.decode().splitlines()

    assert lines[0] == (
        "Perturbed Target ENSG,Perturbed Target Name,Effect Gene ENSG,"
        "Effect Gene Name,Log2FC,Padj,Score Name,Score Value,Cell Type"
    )
    assert {tuple(line.split(",")[:4]) for line in lines[1:]} == {
        ("ENSG00000171421", "MRPL36", "ENSG00000198804", "MT-CO1"),
        ("ENSG00000075624", "ACTB", "ENSG00000180914", "OXTR"),
        ("ENSG00000108064", "TFAM", "ENSG00000198804", "MT-CO1"),
    }
    assert len(lines) == 4


def test_download_without_limit_builds_unbounded_query(run_with_dev_db):
    async def prepare():
        return await data_query._search_dataset_impl(
            "perturb-seq",
            DATASET_ID,
            {"limit": None, "offset": 0},
            return_query=True,
        )

    query, _, _ = run_with_dev_db(prepare)
    assert "LIMIT" not in query


def test_stream_failure_propagates_without_successful_completion(monkeypatch):
    class Context:
        async def __aenter__(self):
            return self

        async def __aexit__(self, *_):
            return False

    class BrokenConnection:
        def transaction(self):
            return Context()

        async def cursor(self, *_args, **_kwargs):
            yield {"perturbed_target_ensg": "ENSG1"}
            raise RuntimeError("database connection lost")

    class BrokenPool:
        def acquire(self):
            class Acquire(Context):
                async def __aenter__(self):
                    return BrokenConnection()

            return Acquire()

    monkeypatch.setitem(data_query.db_pools, "pg", BrokenPool())

    async def no_gene_symbols():
        return {}

    monkeypatch.setattr(data_query, "_fetch_all_gene_symbols", no_gene_symbols)

    async def consume():
        stream = data_query._stream_dataset_csv(
            "SELECT 1",
            [],
            data_query.get_api_to_db_mapping("perturb-seq"),
            "perturb-seq",
        )
        assert (await anext(stream)).startswith(b"Perturbed Target ENSG")
        with pytest.raises(RuntimeError, match="database connection lost"):
            await anext(stream)

    asyncio.run(consume())


def test_full_release_download_redirects_to_signed_url(monkeypatch):
    monkeypatch.setattr(
        data_query,
        "_release_signed_url",
        lambda modality, dataset_id, download_format: (
            f"https://storage.example/{modality}/{dataset_id}.{download_format}"
        ),
    )

    async def download():
        return await data_query.download_dataset_data(
            "mave",
            "dataset-1",
            Request({"type": "http", "query_string": b"format=parquet"}),
            download_format="parquet",
        )

    response = asyncio.run(download())
    assert response.status_code == 307
    assert response.headers["location"].endswith("/mave/dataset-1.parquet")
