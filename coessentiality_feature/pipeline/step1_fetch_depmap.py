#!/usr/bin/env python3
"""
Fetch the latest CRISPRGeneEffect CSV from the DepMap portal API.

Usage:
    python step1_fetch_depmap.py [--output-dir DIR]

Outputs (written to output_dir):
    CRISPRGeneEffect_<version>.csv  — CRISPR gene effect scores (latest release)
    depmap_version.txt              — DepMap release name, e.g. "DepMap Public 26Q2"
"""

import argparse
import hashlib
import io
import os
import re
import sys

import pandas as pd
import requests

FILES_API    = "https://depmap.org/portal/api/download/files"
TARGET_FILE  = "CRISPRGeneEffect.csv"
MODEL_FILE   = "Model.csv"

_HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; depmap-pipeline/1.0)"
}


def fetch_file_listing():
    resp = requests.get(FILES_API, headers=_HEADERS, timeout=30)
    resp.raise_for_status()
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


def run(output_dir):
    os.makedirs(output_dir, exist_ok=True)

    print(f"[1/4] Fetching DepMap file listing from {FILES_API} ...")
    df = fetch_file_listing()
    print(f"      {len(df):,} files listed across {df['release'].nunique()} release(s)")

    print(f"[2/4] Identifying latest {TARGET_FILE} release ...")
    entry   = find_latest_entry(df, TARGET_FILE)
    release = entry["release"]
    version = extract_version(release)
    print(f"      Release : {release}")
    print(f"      Date    : {entry['release_date'].date()}")
    print(f"      Version : {version}")

    dest = os.path.join(output_dir, f"CRISPRGeneEffect_{version}.csv")
    print(f"[3/4] Downloading {TARGET_FILE} to {dest} ...")
    download_with_md5(entry["url"], dest, expected_md5=entry["md5_hash"])

    model_dest = os.path.join(output_dir, "Model.csv")
    print(f"[4/4] Downloading {MODEL_FILE} (cancer model metadata) to {model_dest} ...")
    try:
        model_entry = find_latest_entry(df, MODEL_FILE)
        download_with_md5(model_entry["url"], model_dest, expected_md5=model_entry["md5_hash"])
    except ValueError as exc:
        print(f"      WARNING: {exc}\n      Skipping Model.csv — cancer subtype info will be unavailable.")

    version_path = os.path.join(output_dir, "depmap_version.txt")
    with open(version_path, "w") as fh:
        fh.write(f"{release}\n")
    print(f"      Version recorded -> {version_path}")

    print("Done.")
    return version


def main():
    parser = argparse.ArgumentParser(
        description="Fetch the latest CRISPRGeneEffect CSV from the DepMap portal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output-dir",
        default="./data",
        help="Directory to save the downloaded file (default: ./data)",
    )
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()