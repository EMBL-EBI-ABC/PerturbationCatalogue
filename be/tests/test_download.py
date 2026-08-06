"""Local integration tests for streamed dataset downloads."""

import asyncio

import data_query
import pytest


DATASET_ID = "replogle_2022_rpe1_essential_normalized"


def test_download_streams_csv_from_postgres(run_with_dev_db):
    async def download():
        response = await data_query.download_dataset_data(
            "perturb-seq",
            DATASET_ID,
            limit=3,
            offset=0,
            sort=None,
            perturbation_gene_name=None,
            effect_gene_name=None,
            effect_log2fc=None,
            effect_padj=None,
            effect_score_name=None,
            effect_score_value=None,
            effect_cell_type=None,
            effect_significant=None,
            effect_significance_criteria=None,
            perturbation_name=None,
            perturbation_position=None,
            perturbation_aa_wt=None,
            perturbation_aa_change=None,
        )
        chunks = []
        async for chunk in response.body_iterator:
            chunks.append(chunk)
        return b"".join(chunks)

    content = run_with_dev_db(download)
    lines = content.decode().splitlines()

    assert lines[0].startswith("Perturbed Target ENSG,Effect Gene ENSG")
    assert len(lines) == 4


def test_download_without_limit_builds_unbounded_query(run_with_dev_db):
    async def prepare():
        return await data_query._prepare_dataset_query(
            "perturb-seq", DATASET_ID, {"limit": None, "offset": 0}
        )

    query = run_with_dev_db(prepare)
    assert "LIMIT" not in query["data_query"]


def test_stream_failure_propagates_without_successful_completion(monkeypatch):
    class BrokenCursor:
        def __aiter__(self):
            return self

        async def __anext__(self):
            raise RuntimeError("database connection lost")

    class Context:
        async def __aenter__(self):
            return self

        async def __aexit__(self, *_):
            return False

    class BrokenConnection:
        def transaction(self):
            return Context()

        def cursor(self, *_args, **_kwargs):
            return BrokenCursor()

    class BrokenPool:
        def acquire(self):
            class Acquire(Context):
                async def __aenter__(self):
                    return BrokenConnection()

            return Acquire()

    monkeypatch.setitem(data_query.db_pools, "pg", BrokenPool())
    query = {
        "data_query": "SELECT 1",
        "pg_params": [],
        "api_to_db": data_query.get_api_to_db_mapping("perturb-seq"),
    }

    async def consume():
        stream = data_query._stream_dataset_csv(query, "perturb-seq")
        assert (await anext(stream)).startswith(b"Perturbed Target ENSG")
        with pytest.raises(RuntimeError, match="database connection lost"):
            await anext(stream)

    asyncio.run(consume())
