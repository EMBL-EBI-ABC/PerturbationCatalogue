"""Source catalogue, bundle paths, and canonical ID mappings for MaveDB."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import pandas as pd

from data_sources.mavedb.analysis._reference_io import read_tsv

MAVEDB_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = MAVEDB_ROOT / "tmp" / "reference_data"
DEFAULT_SOURCE_ROOT = MAVEDB_ROOT / "references"
SOURCE_MANIFEST = "manifest.json"
SOURCE_SCHEMA_VERSION = 3
DEFAULT_SOURCE_MANIFEST = DEFAULT_SOURCE_ROOT / SOURCE_MANIFEST
DEFAULT_ID_MAP = Path(__file__).parent / "data" / "mavedb_uniprot_to_ensembl.tsv"
DEFAULT_METADATA = MAVEDB_ROOT / "output" / "processed" / "metadata.tsv"
DEFAULT_ENSEMBL = "https://rest.ensembl.org"
DEFAULT_EVIDENCE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/evidence_eva/"
DEFAULT_CODING = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.06/view/coding_variant_table/"
OT_ROOT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output"
SOURCE_URLS = {
    "coding_variants": DEFAULT_CODING,
    "evidence": DEFAULT_EVIDENCE,
    "credible_sets": f"{OT_ROOT}/credible_set/",
    "pharmacogenomics": f"{OT_ROOT}/pharmacogenomics/",
    "drugs": f"{OT_ROOT}/drug_molecule/",
}
POCKET_FILENAME = "2024.05.28_pockets.tsv"
POCKET_COLUMNS = {
    "struct_id",
    "pocket_id",
    "pocket_resid",
    "pocket_score_combined_scaled",
    "pocket_rad_gyration",
    "pocket_plddt_mean",
}
DEFAULT_POCKETS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/pockets/"
    f"{POCKET_FILENAME}.gz"
)


def reference_directory(root: Path, name: str, url: str) -> Path:
    """Name stored datasets by provider, source release, and dataset role."""
    custom = "custom-" + hashlib.sha256(url.rstrip("/").encode()).hexdigest()[:12]
    if name in SOURCE_URLS:
        release = re.search(r"/platform/(\d{2}\.\d{2})/", urlsplit(url).path)
        dataset = "evidence_eva" if name == "evidence" else name
        return root / "opentargets" / (release[1] if release else custom) / dataset
    if name == "pockets":
        date = re.search(r"\d{4}\.\d{2}\.\d{2}", Path(urlsplit(url).path).name)
        return root / "protvar" / (date[0] if date else custom)
    raise ValueError(f"Unknown downloaded reference source: {name}")


def load_id_map(path: Path) -> pd.DataFrame:
    """Load canonical UniProt/Ensembl pairs, preserving duplicate rows."""
    frame = read_tsv(path).rename(columns={"uniprot_id": "From", "ensembl_id": "To"})
    if not {"From", "To"}.issubset(frame.columns):
        raise ValueError("ID mapping requires From/To or uniprot_id/ensembl_id columns")
    frame = frame[["From", "To"]].dropna().copy()
    frame["From"] = frame["From"].astype(str).str.strip()
    frame["To"] = (
        frame["To"].astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)
    )
    frame = frame.loc[frame["From"].ne("") & frame["To"].ne("")]
    if frame.empty:
        raise ValueError("ID mapping contains no usable UniProt/Ensembl pairs")
    return frame.reset_index(drop=True)


def _resolve_path(record: dict[str, Any], source: str, manifest: Path) -> None:
    relative = record.get("path")
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError(
            f"Source manifest entry {source} requires a relative file path: {manifest}"
        )
    root = manifest.resolve().parent
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(
            f"Source manifest entry {source} escapes its bundle: {manifest}"
        )
    record["path"] = str(resolved)


def _validate_file_identity(record: Any, source: str, manifest: Path) -> None:
    if not isinstance(record, dict):
        raise ValueError(
            f"Source manifest entry {source} requires file identity: {manifest}"
        )
    size = record.get("bytes")
    digest = record.get("sha256")
    _resolve_path(record, source, manifest)
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise ValueError(
            f"Source manifest entry {source} requires a nonnegative byte count: {manifest}"
        )
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise ValueError(
            f"Source manifest entry {source} requires a SHA-256 digest: {manifest}"
        )


def load_sources(path: Path) -> dict[str, Any]:
    """Validate portable bundle descriptors and resolve paths for consumers."""
    value = json.loads(path.read_text())
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != SOURCE_SCHEMA_VERSION
    ):
        raise ValueError(
            f"Unsupported source manifest: {path}; regenerate it with "
            "scripts/download_reference_data.py (existing downloads can be reused offline)"
        )
    if not isinstance(value.get("sources"), dict):
        raise ValueError(f"Source manifest is missing its sources: {path}")
    for name, record in value["sources"].items():
        if not isinstance(record, dict):
            raise ValueError(f"Source manifest entry {name} must be an object: {path}")
        if name in SOURCE_URLS:
            files = record.get("files")
            _resolve_path(record, name, path)
            if not isinstance(files, list) or not files:
                raise ValueError(
                    f"Source manifest entry {name} requires nonempty file descriptors: {path}"
                )
            seen: set[str] = set()
            for index, descriptor in enumerate(files):
                _validate_file_identity(descriptor, f"{name}[{index}]", path)
                file_path = descriptor["path"]
                if not Path(file_path).is_relative_to(Path(record["path"])):
                    raise ValueError(
                        f"Source {name} file is outside its dataset: {path}"
                    )
                if file_path in seen:
                    raise ValueError(
                        f"Source manifest entry {name} contains duplicate file paths: {path}"
                    )
                seen.add(file_path)
            if "download_manifest" in record:
                raise ValueError(
                    f"Source manifest entry {name} cannot reference a child manifest: {path}"
                )
        else:
            _validate_file_identity(record, name, path)
    return value
