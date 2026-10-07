# Study curation source caches

These files are the live inputs to the
[study curation workflow](../curation_tools/study_curation/README.md).

```text
curation_cache/
├── publications/
│   ├── raw/              # Downloaded PDF and XML articles
│   ├── markdown/         # Converted publication text
│   └── download.log
└── mavedb/
    ├── metadata/         # Cached MaveDB JSON records
    ├── urn_to_dois.json
    └── doi_to_fulltext.json
```

Publication filenames use normalized DOI identifiers. The URN-to-DOI map links
MaveDB datasets to publications; the DOI-to-fulltext map points to local raw
article files. The mapping JSON files are versioned; downloaded records, article
files, and logs are local artifacts ignored by Git.

Default paths are defined in `curation_tools/study_curation/paths.py`. Refresh
these caches from the repository root with the project environment active:

```bash
export PYTHONPATH=data_exploration
.venv/bin/python -m curation_tools.study_curation.sources.mavedb
```

Use the source collector's `--help` to override locations or selection settings.
Creating a curation run snapshots selected cache contents into its database;
the live caches can subsequently change without altering that run's inputs.
