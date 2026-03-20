#!/usr/bin/env python3
import argparse
import json
import urllib.request
import urllib.error
import subprocess
import os
import sys
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


def get_srr_accessions(sample_id):
    """Fetch SRR run accessions for a given sample ID from the ENA portal API."""
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sample_id}&result=read_run&fields=run_accession&format=json"
    print(f"Fetching SRR accessions for {sample_id} from ENA API...")
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode("utf-8"))
            srr_ids = [
                item["run_accession"] for item in data if "run_accession" in item
            ]
            return srr_ids
    except urllib.error.URLError as e:
        print(f"Error fetching data from ENA: {e}")
        sys.exit(1)


def process_srr(srr_id, out_dir, temp_base):
    """Download, extract, and compress FASTQ for a single SRR accession using SRA toolkit."""
    srr_temp = os.path.join(temp_base, srr_id)
    os.makedirs(srr_temp, exist_ok=True)

    try:
        # 1. Prefetch the SRA file
        prefetch_cmd = ["prefetch", srr_id, "-O", srr_temp, "--max-size", "100G"]
        subprocess.run(
            prefetch_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        sra_path = None
        for root, _, files in os.walk(srr_temp):
            for file in files:
                if file.endswith(".sra"):
                    sra_path = os.path.join(root, file)
                    break
            if sra_path:
                break

        if not sra_path:
            sra_path = srr_id

        # 2. Extract FASTQ files
        dump_cmd = [
            "fasterq-dump",
            sra_path,
            "--split-files",
            "--include-technical",
            "-O",
            out_dir,
            "-t",
            srr_temp,
            "-e",
            "1",
        ]
        subprocess.run(
            dump_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # 3. Compress the extracted FASTQ files
        fastq_files = [
            f
            for f in os.listdir(out_dir)
            if f.startswith(srr_id) and f.endswith(".fastq")
        ]
        for fq in fastq_files:
            fq_path = os.path.join(out_dir, fq)
            subprocess.run(["gzip", fq_path], check=True)

        return srr_id, True, None
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode("utf-8").strip() if e.stderr else str(e)
        return srr_id, False, error_msg
    except Exception as e:
        return srr_id, False, str(e)
    finally:
        # Cleanup temporary files for this SRR
        if os.path.exists(srr_temp):
            shutil.rmtree(srr_temp)


def format_time(seconds):
    """Helper to format seconds into HH:MM:SS."""
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def main():
    parser = argparse.ArgumentParser(
        description="Download FASTQ data (including technical reads) from SRA directly."
    )
    parser.add_argument(
        "--sample-id", required=True, help="Sample accession ID (e.g., SAMN40972597)"
    )
    parser.add_argument("--out-dir", required=True, help="Base output directory")
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Number of concurrent downloads (default: 8)",
    )
    args = parser.parse_args()

    target_dir = os.path.join(args.out_dir, args.sample_id)
    os.makedirs(target_dir, exist_ok=True)

    srr_ids = get_srr_accessions(args.sample_id)
    if not srr_ids:
        print(f"No run accessions found for sample {args.sample_id}.")
        sys.exit(1)

    print(f"Found {len(srr_ids)} runs for {args.sample_id}.")
    print(f"Output directory: {target_dir}")

    temp_base = os.path.join(target_dir, "_temp")
    os.makedirs(temp_base, exist_ok=True)

    success_count = 0
    fail_count = 0
    total = len(srr_ids)

    print(f"Starting download with {args.jobs} workers...")
    start_time = time.time()

    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(process_srr, srr, target_dir, temp_base): srr
                for srr in srr_ids
            }

            for i, future in enumerate(as_completed(futures), 1):
                srr_id, success, error_msg = future.result()
                if success:
                    success_count += 1
                else:
                    fail_count += 1
                    sys.stdout.write("\r\033[K")
                    print(f"[ERROR] Failed {srr_id}: {error_msg}")

                elapsed = time.time() - start_time
                rate = i / elapsed
                remaining = total - i
                eta_seconds = remaining / rate if rate > 0 else 0

                # Dynamic progress update with Elapsed and ETA
                sys.stdout.write(
                    f"\r\033[KProgress: [{i}/{total}] | Success: {success_count} | Failed: {fail_count} | Elapsed: {format_time(elapsed)} | ETA: {format_time(eta_seconds)} | Last: {srr_id}"
                )
                sys.stdout.flush()

    except KeyboardInterrupt:
        print("\n\nDownload interrupted by user.")
    finally:
        print("\n\nCleaning up temporary files...")
        if os.path.exists(temp_base):
            shutil.rmtree(temp_base)

        elapsed = time.time() - start_time
        print(f"Finished in {format_time(elapsed)}.")
        print(f"Total successful: {success_count}/{total}")
        if fail_count > 0:
            print(f"Total failed: {fail_count}/{total}")
            sys.exit(1)


if __name__ == "__main__":
    main()
