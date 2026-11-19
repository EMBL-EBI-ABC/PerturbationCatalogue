"""Utilities for loading Pydantic extraction schemas from CLI references."""

import importlib
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Type

from pydantic import BaseModel


def _load_module_from_file(module_path: Path) -> ModuleType:
    """Import a Python module from an explicit filesystem path."""
    module_path = module_path.resolve()
    if not module_path.is_file():
        raise FileNotFoundError(f"Schema module file not found: {module_path}")

    module_name = f"_metadata_extraction_schema_{module_path.stem}"
    module_spec = importlib.util.spec_from_file_location(module_name, module_path)
    if module_spec is None or module_spec.loader is None:
        raise ImportError(f"Unable to load schema module from file: {module_path}")

    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


def load_extraction_schema(schema_reference: str) -> Type[BaseModel]:
    """Load a user-provided Pydantic extraction schema.

    Accepted formats:
    - package.module:SchemaClass
    - /path/to/schema.py:SchemaClass
    - relative/path/to/schema.py:SchemaClass
    """
    if ":" not in schema_reference:
        raise ValueError(
            "Schema must be provided as 'module:ClassName' or 'path/to/schema.py:ClassName'."
        )

    module_reference, class_name = schema_reference.rsplit(":", maxsplit=1)
    if not module_reference or not class_name:
        raise ValueError(
            "Schema must include both a module/path and class name, for example "
            "'my_schema:MySchema'."
        )

    module_path = Path(module_reference)
    if module_reference.endswith(".py") or module_path.exists():
        module = _load_module_from_file(module_path)
    else:
        module = importlib.import_module(module_reference)

    schema_class = getattr(module, class_name, None)
    if schema_class is None:
        raise AttributeError(
            f"Schema class '{class_name}' was not found in '{module_reference}'."
        )
    if not isinstance(schema_class, type) or not issubclass(schema_class, BaseModel):
        raise TypeError(
            f"Schema '{schema_reference}' must reference a pydantic.BaseModel subclass."
        )
    return schema_class
