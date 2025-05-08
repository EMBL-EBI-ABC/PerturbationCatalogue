# Data exploration notebooks

Data exploration notebooks use the common environment. To set it up on your machine, run from this directory:

```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

Install the kernel:
```bash
python -m ipykernel install --user --name=data_exploration_env
```

Verify that the kernel has been added:

```bash
jupyter kernelspec list | grep data_exploration_env
```

Run Jupyter notebook:
```bash
jupyter notebook
```

Open any notebook and change the kernel: Kernel → Change kernel → data_exploration_env.
