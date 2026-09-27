#!/usr/bin/env python3
"""Download the official Zhu Cell Ranger matrices with bounded parallelism."""

import argparse
import concurrent.futures
from pathlib import Path
import time
import urllib.request


def https_url(url):
    if not url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        raise ValueError(f"Unexpected GEO URL: {url}")
    return "https://ftp.ncbi.nlm.nih.gov/" + url.split("ftp.ncbi.nlm.nih.gov/", 1)[1]


def download(url, outdir):
    url = https_url(url)
    output = outdir / url.rsplit("/", 1)[1]
    temporary = output.with_suffix(output.suffix + ".part")
    for attempt in range(5):
        try:
            request = urllib.request.Request(url)
            with urllib.request.urlopen(request, timeout=180) as response:
                expected = int(response.headers["Content-Length"])
                if output.is_file() and output.stat().st_size == expected:
                    return output.name, expected, "existing"
                written = 0
                with temporary.open("wb") as handle:
                    while block := response.read(4 * 1024 * 1024):
                        handle.write(block)
                        written += len(block)
            if written != expected:
                raise RuntimeError(
                    f"Short download for {output.name}: {written}/{expected}"
                )
            temporary.replace(output)
            return output.name, written, "downloaded"
        except Exception:
            temporary.unlink(missing_ok=True)
            if attempt == 4:
                raise
            time.sleep(2**attempt)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert (
            https_url("ftp://ftp.ncbi.nlm.nih.gov/a")
            == "https://ftp.ncbi.nlm.nih.gov/a"
        )
        return
    if not 1 <= args.workers <= 32:
        parser.error("--workers must be between 1 and 32")
    urls = [
        line.strip() for line in args.manifest.read_text().splitlines() if line.strip()
    ]
    if len(urls) != 284 or len(urls) != len(set(urls)):
        raise ValueError(f"Expected 284 unique matrix URLs, found {len(urls)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    total = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download, url, args.outdir) for url in urls]
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            name, size, state = future.result()
            total += size
            print(f"files={done}/{len(urls)} bytes={total} {state}={name}", flush=True)


if __name__ == "__main__":
    main()
