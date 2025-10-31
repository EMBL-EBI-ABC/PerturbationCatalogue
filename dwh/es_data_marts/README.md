# dbt project for data marts

### Installing dependencies

Run following commands:
- ```python -m venv .venv```
- ```source .venv/bin/activate```
- ```pip install -r requirements.txt```


### Run dbt:
- ```dbt run --exclude target_summary gene_summary```
- ```dbt run --select target_summary gene_summary```
