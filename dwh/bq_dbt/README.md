# dbt project for data marts

### Installing dependencies
- ```python -m venv .venv```
- ```source .venv/bin/activate```
- ```pip install -r requirements.txt```

### Initialize dbt:
- ```dbt init```

### Update location in ~/.dbt/profiles.yml to europe-west2


### Run dbt:
- ```dbt run --exclude target_summary gene_summary```
- ```dbt run --select target_summary gene_summary```
