#!/usr/bin/env python3
"""
Import the CRISPRGeneEffect CSV for a DepMap release into required_data/.

The DepMap download API is behind a CAPTCHA, so download CRISPRGeneEffect.csv
from https://depmap.org/portal/data_page/?tab=allData and save it in
required_data/ before running this script.

Usage:
    python step1_fetch_depmap.py --release NAME [--input-file FILE] [--output-dir DIR]
    python step1_fetch_depmap.py --from-api [--output-dir DIR]   (blocked by CAPTCHA)

Outputs (written to output_dir):
    CRISPRGeneEffect_<version>.csv  — CRISPR gene effect scores
    depmap_version.txt              — DepMap release name, e.g. "DepMap Public 26Q2"
"""

import argparse
import hashlib
import io
import os
import re
import shutil
import sys

import pandas as pd
import requests

FILES_API    = "https://depmap.org/portal/api/download/files"
TARGET_FILE  = "CRISPRGeneEffect.csv"

_HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; depmap-pipeline/1.0)"
}

# Gene columns in CRISPRGeneEffect.csv are named "SYMBOL (EntrezID)", e.g. "A1BG (1)".
_GENE_COLUMN_RE = re.compile(r"^\S+ \(\d+\)$")

_CAPTCHA_MESSAGE = (
    "The DepMap portal returned its CAPTCHA verification page instead of the file "
    "listing, so the file cannot be downloaded by a script.\n"
    "  Download CRISPRGeneEffect.csv manually from "
    "https://depmap.org/portal/data_page/?tab=allData, save it in required_data/ and run:\n"
    "    python step1_fetch_depmap.py --release \"DepMap Public <version>\""
)


class CaptchaBlockedError(RuntimeError):
    pass


def _is_verification_page(resp):
    content_type = resp.headers.get("content-type", "")
    head = resp.text[:2000].lower()
    return "text/html" in content_type or "<html" in head or "turnstile" in head


def fetch_file_listing():
    resp = requests.get(FILES_API, headers=_HEADERS, timeout=30)
    resp.raise_for_status()
    if _is_verification_page(resp):
        raise CaptchaBlockedError(_CAPTCHA_MESSAGE)
    df = pd.read_csv(io.StringIO(resp.text))

    required = {"release", "release_date", "filename", "url", "md5_hash"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"DepMap files API returned unexpected columns.\n"
            f"  Expected: {sorted(required)}\n"
            f"  Got:      {sorted(df.columns)}\n"
            f"  The API schema may have changed — update this script accordingly."
        )
    return df


def find_latest_entry(df, filename=TARGET_FILE):
    matches = df[df["filename"] == filename].copy()
    if matches.empty:
        raise ValueError(
            f"'{filename}' not found in the DepMap file listing.\n"
            f"  Available filenames (sample): {sorted(df['filename'].unique())[:10]}"
        )
    matches["release_date"] = pd.to_datetime(matches["release_date"])
    return matches.sort_values("release_date", ascending=False).iloc[0]


def extract_version(release_name):
    """Extract short version tag from release name, e.g. 'DepMap Public 26Q1' → '26Q1'."""
    match = re.search(r"(\d+Q\d+)", release_name, re.IGNORECASE)
    if not match:
        raise ValueError(
            f"Could not parse version from release name: {release_name!r}\n"
            f"  Expected a pattern like '26Q1' or '25Q3'."
        )
    return match.group(1)


def file_md5(path):
    hasher = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def validate_gene_effect_csv(path):
    """Fail fast on a wrong file rather than hours into step 2."""
    with open(path, "rb") as fh:
        start = fh.read(1024).lower()
    if b"<html" in start or b"<!doctype" in start:
        raise ValueError(
            f"{path} is an HTML page, not a CSV — the browser probably saved the "
            f"portal's verification page. Download CRISPRGeneEffect.csv again."
        )

    try:
        head = pd.read_csv(path, index_col=0, nrows=5)
    except Exception as exc:
        raise ValueError(f"Could not read {path} as a CSV: {exc}") from exc

    gene_cols = [c for c in head.columns if _GENE_COLUMN_RE.match(str(c))]
    if len(gene_cols) < 0.9 * len(head.columns) or len(gene_cols) < 1000:
        raise ValueError(
            f"{path} does not look like CRISPRGeneEffect.csv.\n"
            f"  Expected ~17,000 gene columns named 'SYMBOL (EntrezID)'; found "
            f"{len(gene_cols):,} of {len(head.columns):,} columns.\n"
            f"  First columns: {list(head.columns[:5])}"
        )
    if not head.index.astype(str).str.startswith("ACH-").all():
        raise ValueError(
            f"{path} does not look like CRISPRGeneEffect.csv: expected cell-line "
            f"IDs starting with 'ACH-' in the first column, got {list(head.index[:3])}."
        )
    return len(gene_cols)


def write_version(output_dir, release):
    version_path = os.path.join(output_dir, "depmap_version.txt")
    with open(version_path, "w") as fh:
        fh.write(f"{release}\n")
    print(f"      Version recorded -> {version_path}")


def download_with_md5(url, dest_path, expected_md5):
    resp = requests.get(url, headers=_HEADERS, stream=True, timeout=600)
    resp.raise_for_status()

    hasher = hashlib.md5()
    bytes_written = 0
    with open(dest_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=4 * 1024 * 1024):
            fh.write(chunk)
            hasher.update(chunk)
            bytes_written += len(chunk)

    size_mb = bytes_written / (1024 * 1024)
    print(f"      {size_mb:.1f} MB written")

    actual_md5 = hasher.hexdigest()
    if actual_md5 != expected_md5:
        os.remove(dest_path)
        raise ValueError(
            f"MD5 mismatch — file removed to avoid corrupt data.\n"
            f"  Expected: {expected_md5}\n"
            f"  Got:      {actual_md5}"
        )
    print(f"      MD5 verified: {actual_md5}")


def run_from_file(input_file, release, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    version = extract_version(release)

    print(f"[1/3] Validating {input_file} ...")
    if not os.path.isfile(input_file):
        raise FileNotFoundError(
            f"Input file not found: {input_file}\n"
            f"  Download CRISPRGeneEffect.csv from "
            f"https://depmap.org/portal/data_page/?tab=allData and save it in {output_dir}"
        )
    n_genes = validate_gene_effect_csv(input_file)
    print(f"      Looks like CRISPRGeneEffect.csv ({n_genes:,} gene columns)")

    dest = os.path.join(output_dir, f"CRISPRGeneEffect_{version}.csv")
    same_dir = os.path.samefile(os.path.dirname(os.path.abspath(input_file)), output_dir)
    print(f"[2/3] {'Renaming' if same_dir else 'Copying'} to {dest} ...")
    if os.path.abspath(input_file) != os.path.abspath(dest):
        if same_dir:
            os.replace(input_file, dest)
        else:
            shutil.copyfile(input_file, dest)
    size_mb = os.path.getsize(dest) / (1024 * 1024)
    print(f"      {size_mb:.1f} MB written")
    print(f"      MD5: {file_md5(dest)}")  # for the record; nothing to verify against

    print(f"[3/3] Recording release {release!r} (version {version}) ...")
    write_version(output_dir, release)

    print("Done.")
    return version


def run_from_api(output_dir):
    os.makedirs(output_dir, exist_ok=True)

    print(f"[1/3] Fetching DepMap file listing from {FILES_API} ...")
    df = fetch_file_listing()
    print(f"      {len(df):,} files listed across {df['release'].nunique()} release(s)")

    print(f"[2/3] Identifying latest {TARGET_FILE} release ...")
    entry   = find_latest_entry(df, TARGET_FILE)
    release = entry["release"]
    version = extract_version(release)
    print(f"      Release : {release}")
    print(f"      Date    : {entry['release_date'].date()}")
    print(f"      Version : {version}")

    dest = os.path.join(output_dir, f"CRISPRGeneEffect_{version}.csv")
    print(f"[3/3] Downloading {TARGET_FILE} to {dest} ...")
    download_with_md5(entry["url"], dest, expected_md5=entry["md5_hash"])

    write_version(output_dir, release)

    print("Done.")
    return version


def main():
    parser = argparse.ArgumentParser(
        description="Import the CRISPRGeneEffect CSV for a DepMap release",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input-file",
        help="CRISPRGeneEffect.csv downloaded manually from the DepMap portal "
             "(default: <output-dir>/CRISPRGeneEffect.csv)",
    )
    parser.add_argument(
        "--release",
        help='Release the file belongs to, e.g. "DepMap Public 26Q2" (required unless --from-api)',
    )
    parser.add_argument(
        "--from-api",
        action="store_true",
        help="Fetch the latest release from the DepMap download API "
             "(currently blocked by the portal's CAPTCHA)",
    )
    parser.add_argument(
        "--output-dir",
        default="./required_data",
        help="Directory to save the file to (default: ./required_data)",
    )
    args = parser.parse_args()

    if args.from_api and (args.input_file or args.release):
        parser.error("--from-api finds the release itself; don't combine it with "
                     "--input-file or --release")
    if not args.from_api and not args.release:
        parser.error('--release is required, e.g. --release "DepMap Public 26Q2"')

    try:
        if args.from_api:
            run_from_api(args.output_dir)
        else:
            input_file = args.input_file or os.path.join(args.output_dir, TARGET_FILE)
            run_from_file(input_file, args.release, args.output_dir)
    except (CaptchaBlockedError, ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
