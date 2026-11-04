"""Shared file publication, provenance and CLI validation for reference pipelines."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd


def read_tsv(path: Path, required: set[str] | None = None) -> pd.DataFrame:
    """Read reference TSVs with explicit missing-value and column rules."""
    if not path.is_file():
        raise FileNotFoundError(f"Required input table not found: {path}")
    frame = pd.read_csv(
        path, sep="\t", keep_default_na=False, na_values=["", "NA"], low_memory=False
    )
    if required and (missing := required - set(frame.columns)):
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    return frame


def file_record(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def verify_file(record: dict[str, Any]) -> Path:
    """Reject source bytes that no longer match the published bundle."""
    path = Path(record["path"])
    if not path.is_file():
        raise FileNotFoundError(f"Downloaded source file not found: {path}")
    if (
        path.stat().st_size != record["bytes"]
        or file_record(path)["sha256"] != record["sha256"]
    ):
        raise ValueError(
            f"Source file changed since the bundle was downloaded: {path}; "
            "rerun download_reference_data.py before processing"
        )
    return path


def write_json(path: Path, value: Any) -> None:
    """Replace a JSON file only after the complete write succeeds."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.close()
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)


def publish_tables(
    staged: Path,
    output: Path,
    *,
    manifest: str,
    managed: Iterable[Path] | None = None,
) -> None:
    """Install a complete result, removing obsolete managed files with rollback.

    The manifest is installed last. Unrelated destination files are preserved.
    Callers must wait for successful completion before reading the result.
    """
    manifest_path = Path(manifest)
    if not (staged / manifest_path).is_file():
        raise FileNotFoundError(f"Staged manifest is missing: {staged / manifest_path}")
    names = (
        set(managed)
        if managed is not None
        else {path.relative_to(staged) for path in staged.rglob("*") if path.is_file()}
    )
    names.add(manifest_path)
    ordered = sorted(names - {manifest_path}) + [manifest_path]
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".reference-backup-", dir=output.parent
    ) as directory:
        backup = Path(directory)
        saved: list[Path] = []
        installed: list[Path] = []
        try:
            for relative in ordered:
                destination = output / relative
                if destination.exists():
                    if not destination.is_file():
                        raise OSError(f"Output is not a regular file: {destination}")
                    previous = backup / relative
                    previous.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(destination, previous)
                    saved.append(relative)
                source = staged / relative
                if source.is_file():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(source, destination)
                    installed.append(relative)
        except (OSError, KeyboardInterrupt):
            for relative in reversed(installed):
                (output / relative).unlink(missing_ok=True)
            for relative in reversed(saved):
                os.replace(backup / relative, output / relative)
            raise


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return number


def nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Must be nonnegative")
    return number


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("Must be finite and positive")
    return number


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("Must be finite and nonnegative")
    return number
