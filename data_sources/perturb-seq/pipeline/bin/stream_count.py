#!/usr/bin/env python3
"""Overlap bounded SRA downloading, disk extraction and one sample's kb counting."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import queue
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from urllib.parse import urlsplit


def download_sra(accession, directory, logdir, execute):
    """Fetch one full-quality archive using 32 checked HTTPS ranges."""
    if not re.fullmatch(r"(SRR|ERR|DRR)\d+", accession):
        raise ValueError("Invalid SRA accession")
    target = directory / accession
    target.mkdir()
    options = [
        "--fail",
        "--silent",
        "--show-error",
        "--location",
        "--proto",
        "=https",
        "--proto-redir",
        "=https",
        "--connect-timeout",
        "30",
        "--max-time",
        "900",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "--retry-max-time",
        "1800",
    ]
    metadata = target / "locator.json"
    execute(
        [
            "curl",
            *options,
            "--max-filesize",
            "1048576",
            "--output",
            str(metadata),
            "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc="
            + accession
            + "&accept-proto=https",
        ],
        logdir / (accession + ".locator.log"),
    )
    bundles = json.loads(metadata.read_text())["result"]
    files = [
        f
        for b in bundles
        if b.get("bundle") == accession and b.get("status") == 200
        for f in b.get("files", [])
        if f.get("type") == "sra" and f.get("accession") == accession
    ]
    if len(files) != 1:
        raise RuntimeError("Expected one full-quality SRA archive: " + accession)
    source = files[0]
    size, digest = source.get("size"), source.get("md5", "")
    urls = [
        loc["link"]
        for loc in source.get("locations", [])
        if urlsplit(loc.get("link", "")).scheme == "https"
    ]
    if (
        type(size) is not int
        or not 0 < size <= 100 * 1024**3
        or not re.fullmatch(r"[0-9a-fA-F]{32}", digest)
        or not urls
    ):
        raise RuntimeError("Missing or invalid archive size, checksum or HTTPS URL")
    url = urls[0]
    if not urlsplit(url).hostname or urlsplit(url).username or urlsplit(url).password:
        raise RuntimeError("Invalid archive URL")
    count = min(32, size)
    ranges = [(size * i // count, size * (i + 1) // count - 1) for i in range(count)]
    command = [
        "curl",
        "--parallel",
        "--parallel-immediate",
        "--parallel-max",
        str(count),
    ]
    for i, (start, end) in enumerate(ranges):
        if i:
            command += ["--next"]
        command += [
            *options,
            "--range",
            f"{start}-{end}",
            "--max-filesize",
            str(end - start + 1),
            "--dump-header",
            str(target / f"headers-{i}"),
            "--output",
            str(target / f"part-{i}"),
            url,
        ]
    execute(command, logdir / (accession + ".curl.log"))
    checksum = hashlib.md5()
    temporary = target / (accession + ".sra.partial")
    with temporary.open("xb") as output:
        for i, (start, end) in enumerate(ranges):
            part, header = target / f"part-{i}", target / f"headers-{i}"
            returned = re.findall(
                r"(?im)^content-range:\s*bytes (\d+)-(\d+)/(\d+)\s*$",
                header.read_text(),
            )
            if (
                not returned
                or tuple(map(int, returned[-1])) != (start, end, size)
                or part.stat().st_size != end - start + 1
            ):
                raise RuntimeError(
                    f"Invalid or incomplete HTTP range for {accession}: {i}"
                )
            with part.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    output.write(block)
                    checksum.update(block)
            # Reclaim each part as it is assembled, bounding the extra space to one range.
            part.unlink()
            header.unlink()
    if temporary.stat().st_size != size or checksum.hexdigest() != digest.lower():
        raise RuntimeError("Archive checksum mismatch: " + accession)
    temporary.replace(target / (accession + ".sra"))
    metadata.unlink()
    return dict(
        archive_bytes=size, archive_md5=digest.lower(), download_connections=count
    )


def extraction_summary(path):
    text = path.read_text()
    result = {}
    for label in ("spots read", "reads read", "reads written"):
        found = re.search(label + r"\s*:\s*([\d,]+)", text)
        if not found:
            raise RuntimeError("Missing extraction summary: " + str(path))
        result[label] = int(found[1].replace(",", ""))
    if result["spots read"] <= 0 or result["reads read"] != result["reads written"]:
        raise RuntimeError("Incomplete extraction: " + str(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accessions", nargs="+", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--t2g", required=True)
    parser.add_argument("--chemistry", default="10xv3")
    parser.add_argument("--workflow", choices=("standard", "kite"), required=True)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--sra-bin", default="")
    args = parser.parse_args()
    if (
        args.cpus < 4
        or len(set(args.accessions)) != len(args.accessions)
        or any(not re.fullmatch(r"(SRR|ERR|DRR)\d+", acc) for acc in args.accessions)
    ):
        parser.error("Use at least four CPUs and unique SRA run accessions")
    extract_threads = min(4, args.cpus - 3)
    count_threads = args.cpus - extract_threads - 2
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    logdir = Path("stream_logs")
    logdir.mkdir()
    router = Path("read_router").resolve()
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++17",
            "-Wall",
            "-Wextra",
            str(Path(__file__).with_name("read_router.cpp")),
            "-o",
            str(router),
        ],
        check=True,
    )
    archives, fastqs = queue.Queue(maxsize=1), queue.Queue(maxsize=1)
    archive_slot, fastq_slot = threading.Semaphore(1), threading.Semaphore(1)
    stopped = threading.Event()
    lock = threading.RLock()
    processes, errors = set(), []
    metrics = {acc: {} for acc in args.accessions}

    def acquire(semaphore):
        while not stopped.is_set():
            if semaphore.acquire(timeout=0.2):
                return
        raise RuntimeError("Pipeline stopped")

    def receive(channel):
        while not stopped.is_set():
            try:
                return channel.get(timeout=0.2)
            except queue.Empty:
                pass
        raise RuntimeError("; ".join(errors) or "Pipeline stopped")

    def execute(command, log, stdout=None, timeout=3600):
        with log.open("wb") as stream:
            with lock:
                if stopped.is_set():
                    raise RuntimeError("Pipeline stopped")
                child = subprocess.Popen(
                    command,
                    stdout=stdout if stdout is not None else stream,
                    stderr=stream,
                    start_new_session=True,
                )
                processes.add(child)
            try:
                code = child.wait(timeout=timeout)
                if code:
                    with lock:
                        concurrent = "; ".join(errors)
                    detail = (
                        "; concurrent producer error: " + concurrent
                        if concurrent
                        else ""
                    )
                    raise RuntimeError(f"{command[0]} exited {code}; see {log}{detail}")
            finally:
                # Keep live children registered until shutdown can kill them.
                if child.poll() is not None:
                    with lock:
                        processes.discard(child)

    def tool(name):
        return str(Path(args.sra_bin) / name) if args.sra_bin else name

    def stage(name, accession, action):
        begin = time.monotonic()
        action()
        with lock:
            metrics[accession][name + "_seconds"] = time.monotonic() - begin

    def guard(action):
        try:
            action()
        except BaseException as error:
            with lock:
                errors.append(repr(error))
            stopped.set()
            with lock:
                for child in list(processes):
                    if child.poll() is None:
                        try:
                            os.killpg(child.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    counter = None
    workers = []
    with tempfile.TemporaryDirectory(prefix="sra-", dir=Path.cwd()) as scratch:
        scratch = Path(scratch)
        downloads, extracted = scratch / "archives", scratch / "fastqs"
        downloads.mkdir()
        extracted.mkdir()

        def download():
            for accession in args.accessions:
                acquire(archive_slot)  # Reserve space before producing another archive.

                def action():
                    metrics[accession].update(
                        download_sra(accession, downloads, logdir, execute)
                    )

                stage("download", accession, action)
                archives.put(accession)
            acquire(archive_slot)
            archives.put(None)

        def extract():
            while True:
                acquire(
                    fastq_slot
                )  # Includes in-progress extraction, so completed FASTQs cannot accumulate.
                accession = receive(archives)
                if accession is None:
                    fastqs.put(None)
                    return
                archive_slot.release()

                def action():
                    target = extracted / accession
                    target.mkdir()
                    logfile = logdir / (accession + ".fasterq.log")
                    execute(
                        [
                            tool("fasterq-dump"),
                            str(downloads / accession),
                            "--split-spot",
                            "--include-technical",
                            "--seq-defline",
                            "@$ac.$si/$ri",
                            "--qual-defline",
                            "+",
                            "-e",
                            str(extract_threads),
                            "-t",
                            str(target),
                            "-o",
                            str(target / "reads.fastq"),
                        ],
                        logfile,
                    )
                    summary = extraction_summary(logfile)
                    metrics[accession].update(
                        spots=summary["spots read"],
                        input_records=summary["reads written"],
                        fastq_bytes=(target / "reads.fastq").stat().st_size,
                    )
                    shutil.rmtree(downloads / accession)

                stage("extract", accession, action)
                fastqs.put(accession)

        command = [
            "kb",
            "count",
            "-i",
            args.index,
            "-g",
            args.t2g,
            "-x",
            args.chemistry,
            "-o",
            "out",
            "--h5ad",
            "--inleaved",
            "-t",
            str(count_threads),
        ]
        command += (
            ["--filter", "bustools"]
            if args.workflow == "standard"
            else ["--workflow", "kite"]
        )
        command += ["-"]
        try:
            counter = subprocess.Popen(
                command, stdin=subprocess.PIPE, start_new_session=True, bufsize=0
            )
            processes.add(counter)
            for action in (download, extract):
                worker = threading.Thread(target=guard, args=(action,))
                workers.append(worker)
                worker.start()
            while True:
                accession = receive(fastqs)
                if accession is None:
                    break
                fastq_slot.release()

                def feed():
                    logfile = logdir / (accession + ".router.log")
                    target = extracted / accession
                    execute(
                        [str(router), accession, str(target / "reads.fastq")],
                        logfile,
                        stdout=counter.stdin,
                        timeout=12 * 3600,
                    )
                    actual = json.loads(logfile.read_text())
                    expected = {
                        key: metrics[accession][key]
                        for key in ("spots", "input_records")
                    }
                    if actual != expected:
                        raise RuntimeError(
                            f"{accession}: router/extraction counts differ: {actual} != {expected}"
                        )
                    shutil.rmtree(target)

                stage("feed", accession, feed)
            counter.stdin.close()
            if counter.wait(timeout=12 * 3600):
                raise RuntimeError("kb count failed")
            for worker in workers:
                worker.join(timeout=5)
                if worker.is_alive():
                    raise RuntimeError("Producer did not finish")
            if errors:
                raise RuntimeError("; ".join(errors))
            info = json.loads(Path("out/run_info.json").read_text())
            if info["n_processed"] != sum(row["spots"] for row in metrics.values()):
                raise RuntimeError(
                    "Counter processed-pair total does not match all extracted spots"
                )
            output = (
                "counts_filtered"
                if args.workflow == "standard"
                else "counts_unfiltered"
            )
            if not (Path("out") / output / "adata.h5ad").is_file():
                raise RuntimeError("Missing count matrix")
        finally:
            stopped.set()
            with lock:
                children = list(processes)
                for child in children:
                    if child.poll() is None:
                        try:
                            os.killpg(child.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
            for child in children:
                child.wait()
            for worker in workers:
                worker.join()
            if counter and counter.stdin and not counter.stdin.closed:
                counter.stdin.close()
    result = dict(
        started_at=started_at,
        finished_at=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic() - started,
        workflow=args.workflow,
        total_spots=sum(row["spots"] for row in metrics.values()),
        runs=metrics,
        spots_per_accession={acc: row["spots"] for acc, row in metrics.items()},
    )
    Path("stream_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":

    def interrupted(signum, frame):
        raise KeyboardInterrupt("Received signal " + str(signum))

    signal.signal(signal.SIGTERM, interrupted)
    main()
