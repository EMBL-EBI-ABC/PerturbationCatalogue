"""Validate and publish an unmodified MaveDB ZIP as a raw ingestion bundle."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import zipfile

LOGGER = logging.getLogger(__name__)


def sha256_file(path: Path) -> str:
    """Hash source bytes without loading the archive into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    """Calculate a file's MD5 for comparison with Zenodo metadata."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def validate_members(archive: zipfile.ZipFile) -> None:
    """Reject paths and file types that cannot form an ordinary dump tree."""
    seen: set[str] = set()
    for member in archive.infolist():
        path = PurePosixPath(member.filename)
        if (
            path.is_absolute()
            or ".." in path.parts
            or "\\" in member.filename
            or not path.parts
            or ":" in path.parts[0]
        ):
            raise ValueError(f"Unsafe archive member: {member.filename}")
        if path.as_posix() in seen:
            raise ValueError(f"Duplicate archive member: {member.filename}")
        seen.add(path.as_posix())
        mode = member.external_attr >> 16
        file_type = stat.S_IFMT(mode)
        if file_type not in (0, stat.S_IFREG, stat.S_IFDIR):
            raise ValueError(f"Unsupported archive member type: {member.filename}")


def validate_dump(root: Path) -> tuple[Path, dict[str, int]]:
    """Locate one dump and validate the metadata shape and score CSV headers."""
    candidates = [path for path in root.rglob("main.json") if path.is_file()]
    if len(candidates) != 1:
        raise ValueError("Archive must contain exactly one main.json")
    dump = candidates[0].parent
    with candidates[0].open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    experiment_sets = (
        metadata.get("experimentSets") if isinstance(metadata, dict) else metadata
    )
    if not isinstance(experiment_sets, list):
        raise ValueError("main.json must contain a list or an experimentSets list")
    score_set_ids: set[str] = set()
    for experiment_set in experiment_sets:
        if not isinstance(experiment_set, dict) or not isinstance(
            experiment_set.get("experiments"), list
        ):
            raise ValueError("Each experiment set must contain an experiments list")
        for experiment in experiment_set["experiments"]:
            if not isinstance(experiment, dict) or not isinstance(
                experiment.get("scoreSets"), list
            ):
                raise ValueError("Each experiment must contain a scoreSets list")
            for score_set in experiment["scoreSets"]:
                if (
                    not isinstance(score_set, dict)
                    or not isinstance(score_set.get("urn"), str)
                    or not score_set["urn"]
                ):
                    raise ValueError("Each score set must contain a urn")
                if score_set["urn"] in score_set_ids:
                    raise ValueError(f"Duplicate score-set urn: {score_set['urn']}")
                score_set_ids.add(score_set["urn"])

    csv_dir = dump / "csv"
    score_files = sorted(csv_dir.glob("*.scores.csv"))
    if not score_files:
        raise ValueError("Dump must contain csv/*.scores.csv")
    for path in score_files:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            header = next(csv.reader(handle), [])
        if "score" not in header:
            raise ValueError(f"Score CSV has no score column: {path.name}")
    return dump, {
        "experiment_sets": len(experiment_sets),
        "metadata_score_sets": len(score_set_ids),
        "score_files": len(score_files),
        "count_files": len(list(csv_dir.glob("*.counts.csv"))),
    }


def ingest_archive(
    archive_path: Path,
    output_dir: Path,
    source: str,
    expected_sha256: str = "",
    zenodo_metadata: dict[str, str | int | bool] | None = None,
) -> Path:
    """Stage, check and atomically rename a content-addressed raw bundle."""
    digest = sha256_file(archive_path)
    if expected_sha256 and digest.lower() != expected_sha256.lower():
        raise ValueError(
            f"Archive SHA-256 mismatch: expected {expected_sha256}, got {digest}"
        )
    if zenodo_metadata:
        if not zenodo_metadata.get("verified"):
            raise ValueError("Zenodo checksum metadata was not marked verified")
        actual_md5 = md5_file(archive_path)
        expected_md5 = str(zenodo_metadata["md5"])
        if actual_md5 != expected_md5:
            raise ValueError(
                f"Archive MD5 mismatch: expected {expected_md5}, got {actual_md5}"
            )
        expected_size = int(zenodo_metadata["size"])
        actual_size = archive_path.stat().st_size
        if actual_size != expected_size:
            raise ValueError(
                f"Archive size mismatch: expected {expected_size}, got {actual_size}"
            )
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"sha256-{digest}"
    if destination.exists():
        raise FileExistsError(f"Bundle already exists: {destination}")
    with tempfile.TemporaryDirectory(prefix=".ingest-", dir=output_dir) as temporary:
        stage = Path(temporary)
        extracted = stage / "extracted"
        with zipfile.ZipFile(archive_path) as archive:
            validate_members(archive)
            # Extraction checks ZIP CRCs while streaming every member to disk.
            archive.extractall(extracted)
        dump, counts = validate_dump(extracted)
        bundle = stage / "bundle"
        bundle.mkdir()
        shutil.move(str(dump), str(bundle / "mavedb_dump"))
        manifest = {
            "schema_version": 1,
            "source": source,
            "archive_sha256": digest,
            "archive_bytes": archive_path.stat().st_size,
            "expected_sha256": expected_sha256.lower() or None,
            "zenodo_file": zenodo_metadata,
            "ingested_at": datetime.now(timezone.utc).isoformat(),
            "dump_path": "mavedb_dump",
            "metadata_path": "mavedb_dump/main.json",
            "scores_path": "mavedb_dump/csv",
            "counts": counts,
        }
        with (bundle / "manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        bundle.rename(destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--expected-sha256", default="")
    parser.add_argument("--zenodo-metadata", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        zenodo_metadata = (
            json.loads(args.zenodo_metadata.read_text(encoding="utf-8"))
            if args.zenodo_metadata
            else None
        )
        destination = ingest_archive(
            args.archive,
            args.output,
            args.source,
            args.expected_sha256,
            zenodo_metadata=zenodo_metadata,
        )
    except (
        OSError,
        ValueError,
        RuntimeError,
        EOFError,
        csv.Error,
        zipfile.BadZipFile,
    ) as error:
        parser.exit(1, f"MaveDB ingest failed: {error}\n")
    LOGGER.info("Published raw ingestion bundle: %s", destination)


if __name__ == "__main__":
    main()
