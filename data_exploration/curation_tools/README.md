# Curation Tools

This package provides utilities for data curation in the PerturbationCatalogue project.

## Installation

You can install the package locally in editable mode using either `pip` or `uv`. This allows you to import `curation_tools` anywhere in your project without using relative imports.

### Using pip

Open a terminal in the `data_exploration` directory and run:

```bash
python -m pip install -e .
```

### Using uv

If you have [uv](https://github.com/astral-sh/uv) installed, run:

```bash
uv pip install -e .
```

## Usage

After installation, you can import modules and classes from `curation_tools` in your scripts or notebooks:

```python
from curation_tools.curation_tools import CuratedDataset, ObsSchema, VarSchema, Experiment
```

Any changes to the source code will be reflected immediately thanks to the editable install.


## Development

For development, keep using the editable install (`-e`) so your changes are instantly available.