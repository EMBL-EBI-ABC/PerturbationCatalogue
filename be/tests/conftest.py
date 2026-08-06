"""Shared local BE test helpers.

These tests intentionally use the configured development PostgreSQL database.
Source pc_secrets before running pytest; the fixture fails early otherwise.
"""

import asyncio
import os
import sys
from pathlib import Path

import asyncpg
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
import data_query  # noqa: E402


@pytest.fixture
def run_with_dev_db():
    """Run one async operation with a short-lived dev database pool."""
    required = ("PG_HOST", "PG_PORT", "PG_USER", "PG_PASSWORD", "PG_DB")
    missing = [name for name in required if not os.getenv(name)]
    if missing:
        pytest.fail(
            "Source pc_secrets dev before running BE tests; missing: "
            + ", ".join(missing)
        )

    def run(operation):
        async def execute():
            pool = await asyncpg.create_pool(
                host=os.environ["PG_HOST"],
                port=int(os.environ["PG_PORT"]),
                user=os.environ["PG_USER"],
                password=os.environ["PG_PASSWORD"],
                database=os.environ["PG_DB"],
                min_size=1,
                max_size=1,
            )
            previous_pool = data_query.db_pools.get("pg")
            data_query.db_pools["pg"] = pool
            try:
                return await operation()
            finally:
                if previous_pool is None:
                    data_query.db_pools.pop("pg", None)
                else:
                    data_query.db_pools["pg"] = previous_pool
                await pool.close()

        return asyncio.run(execute())

    return run
