# BE tests

These are local PostgreSQL and Elasticsearch integration tests. They are
intentionally not part of CI and use the development services configured by
`pc_secrets`.

From the repository root:

```bash
pc_secrets -v dev
be/fastapi-env/bin/pytest -q be/tests
```

The external development PostgreSQL and Elasticsearch connections must be
reachable from the current machine. The suite uses the existing
`replogle_2022_rpe1_essential_normalized` dataset and does not modify data.

The tests cover the streamed CSV response, the no-limit query shape, and
propagation of a mid-stream database failure. Keep BE tests in this directory
and run them manually; do not add them to deployment or CI actions yet.
