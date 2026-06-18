"""Load Open Targets Platform targets into the configured BigQuery table.

This loader is intentionally scoped to the Open Targets reference table used by
the ENSG dev data stack. It downloads the release target parquet parts,
assembles a single parquet file, creates the configured reference dataset if
needed, and loads the parquet into BigQuery.
"""

from __future__ import annotations

import argparse
import os
import re
import time
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable
from urllib.error import URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import pyarrow.parquet as pq
from google.cloud import bigquery


DEFAULT_RELEASE = "26.03"
DEFAULT_BASE_TEMPLATE = (
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/"
    "{release}/output/target/"
)
DEFAULT_WORK_DIR = "/tmp/perturbation_catalogue/opentargets"
CHUNK_SIZE = 1024 * 1024
EXPECTED_COLUMNS = {
    "id",
    "approvedSymbol",
    "approvedName",
    "symbolSynonyms",
    "nameSynonyms",
    "obsoleteSymbols",
    "obsoleteNames",
    "synonyms",
}
BIGQUERY_IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class LinkParser(HTMLParser):
    """Extract href values from an Apache-style directory listing."""

    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for name, value in attrs:
            if name.lower() == "href" and value:
                self.links.append(value)


def format_bytes(size: int | None) -> str:
    """Format byte counts for progress output."""
    if size is None:
        return "unknown"
    units = ["B", "KB", "MB", "GB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GB"


def release_slug(release: str) -> str:
    """Return a filesystem-friendly release label."""
    return release.replace(".", "_").replace("-", "_")


def validate_bigquery_identifier(value: str, label: str) -> None:
    """Reject identifiers that cannot be used as BigQuery dataset/table IDs."""
    if not BIGQUERY_IDENTIFIER_RE.fullmatch(value):
        raise ValueError(f"{label} must be a BigQuery identifier, got: {value!r}")


def validate_dev_destination(dataset: str, table: str, allow_non_dev: bool) -> None:
    """Avoid accidental writes to production reference tables."""
    if allow_non_dev:
        return
    if "ensg_dev" in dataset or "ensg_dev" in table:
        return
    raise ValueError(
        "Refusing to load Open Targets into a non-dev destination. "
        "Use a BQ_REFERENCE_DATASET/table containing 'ensg_dev', or pass "
        "--allow-non-dev-destination explicitly."
    )


def list_parquet_parts(base_url: str) -> list[str]:
    """List Parquet part URLs in an Open Targets target directory."""
    request = Request(base_url, headers={"User-Agent": "PerturbationCatalogue/1.0"})
    with urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")

    parser = LinkParser()
    parser.feed(html)

    part_names = sorted(
        {
            link
            for link in parser.links
            if link.endswith(".parquet") and not link.startswith("/")
        }
    )
    if not part_names:
        raise RuntimeError(f"No Parquet parts found at {base_url}")

    return [urljoin(base_url, part_name) for part_name in part_names]


def parquet_part_path(download_dir: Path, part_url: str) -> Path:
    """Return the local path for a remote part URL."""
    return download_dir / part_url.rstrip("/").split("/")[-1]


def is_valid_parquet(path: Path) -> bool:
    """Return True when a local file can be opened as Parquet."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        pq.ParquetFile(path)
    except Exception:
        return False
    return True


def parquet_row_count(path: Path) -> int:
    """Return the number of rows in a parquet file."""
    return pq.ParquetFile(path).metadata.num_rows


def download_part(part_url: str, destination: Path, overwrite: bool) -> None:
    """Download one Parquet part atomically."""
    if destination.exists() and is_valid_parquet(destination) and not overwrite:
        print(
            f"    already downloaded: {format_bytes(destination.stat().st_size)}",
            flush=True,
        )
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    tmp_path.unlink(missing_ok=True)

    request = Request(part_url, headers={"User-Agent": "PerturbationCatalogue/1.0"})
    with urlopen(request, timeout=120) as response, tmp_path.open("wb") as handle:
        total_header = response.headers.get("Content-Length")
        total_bytes = int(total_header) if total_header else None
        downloaded_bytes = 0
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            handle.write(chunk)
            downloaded_bytes += len(chunk)

            if total_bytes:
                percent = downloaded_bytes / total_bytes * 100
                print(
                    (
                        f"    downloaded {format_bytes(downloaded_bytes)} / "
                        f"{format_bytes(total_bytes)} ({percent:.1f}%)"
                    ),
                    flush=True,
                )
            else:
                print(f"    downloaded {format_bytes(downloaded_bytes)}", flush=True)

    if not is_valid_parquet(tmp_path):
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file is not valid Parquet: {part_url}")

    tmp_path.replace(destination)
    print(f"    saved: {format_bytes(destination.stat().st_size)}", flush=True)


def download_parts(
    part_urls: Iterable[str],
    download_dir: Path,
    overwrite: bool,
) -> list[Path]:
    """Download all Parquet parts and return local paths."""
    part_urls = list(part_urls)
    local_paths: list[Path] = []
    total = len(part_urls)

    for index, part_url in enumerate(part_urls, start=1):
        path = parquet_part_path(download_dir, part_url)
        print(f"[{index}/{total}] {path.name}", flush=True)
        try:
            download_part(part_url, path, overwrite=overwrite)
        except URLError as exc:
            raise RuntimeError(f"Failed to download {part_url}: {exc}") from exc
        local_paths.append(path)

    return local_paths


def assemble_parquet(part_paths: list[Path], output_path: Path, overwrite: bool) -> int:
    """Write one Parquet file by streaming downloaded parts."""
    if output_path.exists() and is_valid_parquet(output_path) and not overwrite:
        print(f"Using existing assembled parquet: {output_path}", flush=True)
        return parquet_row_count(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output_path.with_suffix(output_path.suffix + ".tmp")
    tmp_output.unlink(missing_ok=True)

    writer: pq.ParquetWriter | None = None
    total_rows = 0
    try:
        for path in part_paths:
            print(f"Assembling {path.name}", flush=True)
            table = pq.read_table(path)
            if writer is None:
                writer = pq.ParquetWriter(
                    tmp_output,
                    table.schema,
                    compression="snappy",
                    use_dictionary=True,
                )
            elif table.schema != writer.schema:
                raise RuntimeError(f"Schema mismatch in {path}")

            writer.write_table(table)
            total_rows += table.num_rows
    finally:
        if writer is not None:
            writer.close()

    if not is_valid_parquet(tmp_output):
        tmp_output.unlink(missing_ok=True)
        raise RuntimeError(f"Assembled output is not valid Parquet: {tmp_output}")

    tmp_output.replace(output_path)
    return total_rows


def prepare_parquet(args: argparse.Namespace) -> tuple[Path, int]:
    """Return a local Open Targets target parquet, downloading it if needed."""
    if args.parquet_path:
        path = Path(args.parquet_path)
        if not is_valid_parquet(path):
            raise ValueError(f"Input parquet is missing or invalid: {path}")
        return path, parquet_row_count(path)

    slug = release_slug(args.release)
    work_dir = Path(args.work_dir)
    download_dir = Path(args.download_dir or work_dir / f"target_parts_{slug}")
    output_path = Path(args.output or work_dir / f"opentargets_targets_{slug}.parquet")
    base_url = args.base_url or DEFAULT_BASE_TEMPLATE.format(release=args.release)

    print(f"Listing Open Targets target parts: {base_url}", flush=True)
    part_urls = list_parquet_parts(base_url)
    print(f"Found {len(part_urls)} Parquet parts", flush=True)

    part_paths = download_parts(part_urls, download_dir, overwrite=args.overwrite)
    print(f"Assembling {len(part_paths)} parts into {output_path}", flush=True)
    total_rows = assemble_parquet(part_paths, output_path, overwrite=args.overwrite)
    return output_path, total_rows


def wait_for_job(job: bigquery.LoadJob, label: str) -> None:
    """Poll a BigQuery job and print progress while it runs."""
    while not job.done():
        print(f"{label}: {job.state}", flush=True)
        time.sleep(15)
    job.result()


def load_reference_table(
    args: argparse.Namespace, parquet_path: Path
) -> bigquery.Table:
    """Create the configured dataset if needed and load the parquet into BigQuery."""
    client = bigquery.Client(project=args.project)
    dataset_ref = bigquery.Dataset(f"{args.project}.{args.dataset}")
    dataset_ref.location = args.location
    client.create_dataset(dataset_ref, exists_ok=True)

    destination = f"{args.project}.{args.dataset}.{args.table}"
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        write_disposition=args.write_disposition,
    )

    print(f"Loading {parquet_path} into {destination}", flush=True)
    with parquet_path.open("rb") as parquet_file:
        load_job = client.load_table_from_file(
            parquet_file,
            destination,
            job_config=job_config,
            location=args.location,
            rewind=True,
        )
    wait_for_job(load_job, f"BigQuery load job {load_job.job_id}")

    table = client.get_table(destination)
    schema_columns = {field.name for field in table.schema}
    missing_columns = sorted(EXPECTED_COLUMNS - schema_columns)
    if missing_columns:
        raise RuntimeError(
            f"{destination} is missing expected columns: {', '.join(missing_columns)}"
        )

    return table


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Load Open Targets target parquet into the configured BQ table."
    )
    parser.add_argument("--project", default=os.getenv("GCLOUD_PROJECT"))
    parser.add_argument("--location", default=os.getenv("BQ_LOCATION", "EU"))
    parser.add_argument("--dataset", default=os.getenv("BQ_REFERENCE_DATASET"))
    parser.add_argument("--table", default=os.getenv("BQ_OPENTARGETS_TARGETS_TABLE"))
    parser.add_argument(
        "--release",
        default=os.getenv("OPENTARGETS_RELEASE", DEFAULT_RELEASE),
        help=f"Open Targets Platform release. Default: {DEFAULT_RELEASE}",
    )
    parser.add_argument("--base-url", default=None)
    parser.add_argument(
        "--work-dir",
        default=os.getenv("OPENTARGETS_WORK_DIR", DEFAULT_WORK_DIR),
    )
    parser.add_argument("--download-dir", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument(
        "--parquet-path",
        default=None,
        help=(
            "Use an existing assembled Open Targets target parquet instead of "
            "downloading."
        ),
    )
    parser.add_argument(
        "--write-disposition",
        default=bigquery.WriteDisposition.WRITE_TRUNCATE,
        choices=[
            bigquery.WriteDisposition.WRITE_APPEND,
            bigquery.WriteDisposition.WRITE_EMPTY,
            bigquery.WriteDisposition.WRITE_TRUNCATE,
        ],
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download existing parts and replace an existing assembled parquet.",
    )
    parser.add_argument(
        "--allow-non-dev-destination",
        action="store_true",
        help="Allow loading into a destination that does not contain 'ensg_dev'.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Validate configuration and print the planned destination without "
            "loading."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the Open Targets reference load."""
    args = parse_args()
    if not args.project:
        raise ValueError("GCLOUD_PROJECT or --project is required")
    if not args.dataset:
        raise ValueError("BQ_REFERENCE_DATASET or --dataset is required")
    if not args.table:
        raise ValueError("BQ_OPENTARGETS_TARGETS_TABLE or --table is required")

    validate_bigquery_identifier(args.dataset, "dataset")
    validate_bigquery_identifier(args.table, "table")
    validate_dev_destination(args.dataset, args.table, args.allow_non_dev_destination)

    destination = f"{args.project}.{args.dataset}.{args.table}"
    print("Open Targets reference load", flush=True)
    print(f"  release:     {args.release}", flush=True)
    print(f"  location:    {args.location}", flush=True)
    print(f"  destination: {destination}", flush=True)
    print(f"  disposition: {args.write_disposition}", flush=True)

    if args.dry_run:
        print(
            "Dry run only; no files downloaded and no BigQuery table updated.",
            flush=True,
        )
        return

    parquet_path, parquet_rows = prepare_parquet(args)
    print(f"Prepared parquet: {parquet_path}", flush=True)
    print(f"Parquet rows: {parquet_rows}", flush=True)

    table = load_reference_table(args, parquet_path)
    print(f"Loaded rows: {table.num_rows}", flush=True)
    print(
        f"Loaded table: {table.project}.{table.dataset_id}.{table.table_id}",
        flush=True,
    )


if __name__ == "__main__":
    main()
