"""Download an NCBI taxonomy categories snapshot and its provenance manifest."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import shutil
import tarfile
import tempfile
from urllib.request import Request, urlopen

DEFAULT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz"
LOGGER = logging.getLogger(__name__)


def file_identity(path: Path) -> dict[str, int | str]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def download_taxonomy(output: Path, url: str = DEFAULT_URL) -> None:
    """Stage a verified tar extraction before exposing the taxonomy bundle."""
    if output.exists():
        raise FileExistsError(f"Taxonomy output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".taxonomy-", dir=output.parent
    ) as temporary:
        stage = Path(temporary)
        archive_path = stage / "taxcat.tar.gz"
        request = Request(
            url, headers={"User-Agent": "PerturbationCatalogue taxonomy downloader"}
        )
        with urlopen(request, timeout=60) as response, archive_path.open(
            "wb"
        ) as handle:
            shutil.copyfileobj(response, handle)
        bundle = stage / "taxonomy"
        bundle.mkdir()
        taxonomy = bundle / "taxcat.dmp"
        try:
            with tarfile.open(archive_path, mode="r:gz") as archive:
                members = [
                    member
                    for member in archive.getmembers()
                    if Path(member.name).name == "categories.dmp" and member.isfile()
                ]
                if len(members) != 1:
                    raise ValueError(
                        "NCBI archive must contain one categories.dmp file"
                    )
                source = archive.extractfile(members[0])
                if source is None:
                    raise ValueError("Could not read categories.dmp from NCBI archive")
                # Only copy the selected regular file; archive paths are never extracted.
                with source, taxonomy.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
        except tarfile.TarError as error:
            raise ValueError(f"Invalid NCBI taxonomy archive: {error}") from error
        manifest = {
            "schema_version": 1,
            "source": url,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "archive": file_identity(archive_path),
            "path": taxonomy.name,
            **file_identity(taxonomy),
        }
        with (bundle / "manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        bundle.rename(output)
    LOGGER.info("Published taxonomy snapshot: %s", output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    download_taxonomy(args.output, args.url)


if __name__ == "__main__":
    main()
