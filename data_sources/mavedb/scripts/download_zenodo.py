"""Download a MaveDB dump from Zenodo and optionally unpack it."""

from __future__ import annotations

import argparse
from pathlib import Path
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
    args = parser.parse_args()

    output_dir = args.output_dir
    archive_path = output_dir / archive_name_from_url(args.url)

    try:
        download_archive(args.url, archive_path)
        print(f"Downloaded archive to {archive_path}")

        if not args.no_unpack:
            unpack_archive(archive_path, output_dir)
            print(f"Unpacked archive into {output_dir}")
    except (OSError, URLError, ValueError, zipfile.BadZipFile) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
