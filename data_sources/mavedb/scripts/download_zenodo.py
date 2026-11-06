"""Download a MaveDB dump from Zenodo and optionally unpack it."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from urllib.error import URLError
from urllib.parse import unquote, urlsplit
from urllib.request import Request, urlopen
import zipfile

from tqdm import tqdm

DEFAULT_URL = (
    "https://zenodo.org/records/18511521/files/"
    "mavedb-dump.20260206153444.zip?download=1"
)
_RECORD_FILE_PATH = re.compile(r"/records/(\d+)/files/([^/]+)$")


def zenodo_file_metadata(url: str) -> dict[str, str | int]:
    """Read the published checksum and size for a direct Zenodo file URL."""
    parsed_url = urlsplit(url)
    if parsed_url.scheme != "https" or not parsed_url.hostname:
        raise ValueError("Zenodo file URL must use HTTPS")
    if parsed_url.hostname != "zenodo.org" and not parsed_url.hostname.endswith(
        ".zenodo.org"
    ):
        raise ValueError("Expected a file URL hosted on zenodo.org")
    match = _RECORD_FILE_PATH.search(parsed_url.path)
    if not match:
        raise ValueError(
            "Expected a versioned Zenodo file URL like "
            "https://zenodo.org/records/<id>/files/<filename>"
        )

    record_id, encoded_filename = match.groups()
    filename = unquote(encoded_filename)
    api_url = f"{parsed_url.scheme}://{parsed_url.netloc}/api/records/{record_id}"
    request = Request(api_url, headers={"Accept": "application/json"})
    with urlopen(request, timeout=60) as response:
        record = json.load(response)

    files = record.get("files") if isinstance(record, dict) else None
    if not isinstance(files, list):
        raise ValueError(f"Zenodo record {record_id} returned no files list")
    matches = [
        item
        for item in files
        if isinstance(item, dict) and item.get("key") == filename
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {filename!r} entry in Zenodo record {record_id}; "
            f"found {len(matches)}"
        )
    file_info = matches[0]
    checksum = file_info.get("checksum", "")
    algorithm, separator, digest = checksum.partition(":")
    if not separator or algorithm.lower() != "md5" or not re.fullmatch(
        r"[a-fA-F0-9]{32}", digest
    ):
        raise ValueError(f"Zenodo returned an unsupported file checksum: {checksum!r}")
    if "size" not in file_info:
        raise ValueError(f"Zenodo returned no file size for {filename!r}")
    try:
        size = int(file_info["size"])
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Zenodo returned an invalid file size for {filename!r}"
        ) from error
    return {
        "record_id": record_id,
        "filename": filename,
        "md5": digest.lower(),
        "size": size,
    }


def md5_file(path: Path) -> str:
    """Calculate the MD5 of a file using bounded memory."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_archive(path: Path, metadata: dict[str, str | int]) -> None:
    """Reject a cached or downloaded file that differs from Zenodo metadata."""
    actual_size = path.stat().st_size
    expected_size = int(metadata["size"])
    if actual_size != expected_size:
        raise ValueError(
            f"Archive size mismatch for {path.name}: "
            f"expected {expected_size} bytes, got {actual_size}"
        )
    actual_md5 = md5_file(path)
    expected_md5 = str(metadata["md5"])
    if actual_md5 != expected_md5:
        raise ValueError(
            f"Archive MD5 mismatch for {path.name}: "
            f"expected {expected_md5}, got {actual_md5}"
        )


def archive_name_from_url(url: str) -> str:
    """Return the archive filename from a URL, with a useful fallback."""
    name = Path(unquote(urlsplit(url).path)).name
    return name or "mavedb-dump.zip"


def download_archive(url: str, output_path: Path) -> None:
    """Stream a URL to disk, replacing the destination only after success."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Starting download: {url}")
    request = Request(
        url, headers={"User-Agent": "PerturbationCatalogue MaveDB downloader"}
    )
    temporary_path = None

    try:
        with urlopen(request, timeout=60) as response:
            try:
                total_bytes = int(response.headers.get("Content-Length", ""))
            except ValueError:
                total_bytes = 0

            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=output_path.parent,
                prefix=f".{output_path.name}.",
                suffix=".part",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                with tqdm(
                    total=total_bytes or None,
                    desc="Downloading",
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                ) as progress:
                    while chunk := response.read(1024 * 1024):
                        temporary_file.write(chunk)
                        progress.update(len(chunk))

        temporary_path.replace(output_path)
    except (OSError, URLError):
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def unpack_archive(archive_path: Path, output_dir: Path) -> None:
    """Unpack a ZIP archive after checking that its members stay in output_dir."""
    resolved_output_dir = output_dir.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target_path = (output_dir / member.filename).resolve()
            if not target_path.is_relative_to(resolved_output_dir):
                raise ValueError(
                    f"Archive member would be written outside the output directory: "
                    f"{member.filename}"
                )
        archive.extractall(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download and optionally unpack a MaveDB dump from Zenodo."
    )
    parser.add_argument(
        "--url",
        default=DEFAULT_URL,
        help=f"Zenodo archive URL (default: {DEFAULT_URL})",
    )
    parser.add_argument(
        "--output",
        "--output-dir",
        dest="output_dir",
        type=Path,
        required=True,
        help="Directory where the ZIP archive is saved and, by default, unpacked.",
    )
    parser.add_argument(
        "--no-unpack",
        dest="no_unpack",
        action="store_true",
        default=False,
        help=(
            "Skip unpacking the archive; by default, the archive is unpacked "
            "into the output directory."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for persistent copies of downloaded archives.",
    )
    parser.add_argument(
        "--metadata-output",
        type=Path,
        help="Write the verified Zenodo file metadata as JSON to this path.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    archive_path = output_dir / archive_name_from_url(args.url)

    try:
        metadata = zenodo_file_metadata(args.url)
        archive_path = output_dir / str(metadata["filename"])
        if args.cache_dir:
            checksum_dir = args.cache_dir / f"md5-{metadata['md5']}"
            cached_archive = checksum_dir / archive_path.name
            if cached_archive.is_file():
                try:
                    verify_archive(cached_archive, metadata)
                except (OSError, ValueError):
                    print(f"Discarding invalid cached archive: {cached_archive}")
                    cached_archive.unlink(missing_ok=True)
            if not cached_archive.is_file():
                download_archive(args.url, cached_archive)
                try:
                    verify_archive(cached_archive, metadata)
                except (OSError, ValueError):
                    cached_archive.unlink(missing_ok=True)
                    raise
                print(f"Downloaded and verified archive in cache: {cached_archive}")
            else:
                print(f"Using verified cached archive: {cached_archive}")
            output_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cached_archive, archive_path)
        else:
            download_archive(args.url, archive_path)
            try:
                verify_archive(archive_path, metadata)
            except (OSError, ValueError):
                archive_path.unlink(missing_ok=True)
                raise
            print(f"Downloaded and verified archive to {archive_path}")

        if args.metadata_output:
            args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
            args.metadata_output.write_text(
                json.dumps({**metadata, "verified": True}, indent=2) + "\n",
                encoding="utf-8",
            )

        if not args.no_unpack:
            unpack_archive(archive_path, output_dir)
            print(f"Unpacked archive into {output_dir}")
    except (OSError, URLError, ValueError, zipfile.BadZipFile) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
