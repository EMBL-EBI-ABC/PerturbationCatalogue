#!/usr/bin/env python3
"""Overlap bounded SRA downloading, disk extraction and one sample's kb counting."""
import argparse
from datetime import datetime, timezone
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


def disk_bytes(root):
    """Conservative file bytes; shared filesystems can delay allocated-block accounting."""
    total = 0
    for directory, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = [d for d in dirs if not (Path(directory) / d).is_symlink()]
        for name in files:
            try:
                stat = (Path(directory) / name).lstat()
                total += max(stat.st_size, stat.st_blocks * 512)
            except FileNotFoundError:
                pass
    return total


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
    state = {"download": None, "extract": None, "feed": None}
    completed = {key: 0 for key in state}
    metrics = {acc: {} for acc in args.accessions}
    peaks = {"task_bytes": 0, "buffer_bytes": 0}
    events = Path("stream_events.jsonl").open("w", buffering=1)

    def event(kind, **values):
        with lock:
            events.write(
                json.dumps(dict(event=kind, timestamp=time.time(), **values)) + "\n"
            )

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
                    raise RuntimeError(f"{command[0]} exited {code}; see {log}")
            finally:
                # Keep live children registered until shutdown can kill them.
                if child.poll() is not None:
                    with lock:
                        processes.discard(child)

    def tool(name):
        return str(Path(args.sra_bin) / name) if args.sra_bin else name

    def stage(name, accession, action):
        with lock:
            state[name] = accession
        event(name + "_start", accession=accession)
        begin = time.monotonic()
        action()
        with lock:
            metrics[accession][name + "_seconds"] = time.monotonic() - begin
            completed[name] += 1
            state[name] = None
        event(name + "_complete", accession=accession, **metrics[accession])

    def guard(action):
        try:
            action()
        except BaseException as error:
            with lock:
                errors.append(repr(error))
            stopped.set()
            event("failure", error=repr(error))
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
                    execute(
                        [
                            tool("prefetch"),
                            accession,
                            "-O",
                            str(downloads),
                            "--max-size",
                            "100G",
                        ],
                        logdir / (accession + ".prefetch.log"),
                    )
                    archive = downloads / accession
                    if not (archive / (accession + ".sra")).is_file():
                        raise RuntimeError("Missing downloaded archive: " + accession)
                    metrics[accession]["archive_bytes"] = disk_bytes(archive)

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
                with lock:
                    state["extract"] = accession
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

        def sample():
            while not stopped.is_set():
                task_size, buffer_size = disk_bytes(Path.cwd()), disk_bytes(scratch)
                with lock:
                    peaks["task_bytes"] = max(peaks["task_bytes"], task_size)
                    peaks["buffer_bytes"] = max(peaks["buffer_bytes"], buffer_size)
                    snapshot = dict(
                        timestamp=time.time(),
                        elapsed_seconds=time.monotonic() - started,
                        workflow=args.workflow,
                        state=dict(state),
                        completed=dict(completed),
                        runs_total=len(args.accessions),
                        task_bytes=task_size,
                        buffer_bytes=buffer_size,
                        peaks=dict(peaks),
                        allocated_cpus=args.cpus,
                        extraction_threads=extract_threads,
                        counting_threads=count_threads,
                        runs={acc: dict(values) for acc, values in metrics.items()},
                    )
                temporary = Path("stream_status.json.tmp")
                temporary.write_text(json.dumps(snapshot) + "\n")
                temporary.replace("stream_status.json")
                event(
                    "sample",
                    task_bytes=task_size,
                    buffer_bytes=buffer_size,
                    state=snapshot["state"],
                    completed=snapshot["completed"],
                )
                stopped.wait(5)

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
            for action in (download, extract, sample):
                worker = threading.Thread(target=guard, args=(action,))
                workers.append(worker)
                worker.start()
            while True:
                accession = receive(fastqs)
                if accession is None:
                    break
                with lock:
                    state["feed"] = accession
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
            event("count_finalize_start")
            if counter.wait(timeout=12 * 3600):
                raise RuntimeError("kb count failed")
            for worker in workers[:2]:
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
            event("count_finalize_complete")
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
            events.close()
    result = dict(
        started_at=started_at,
        finished_at=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic() - started,
        workflow=args.workflow,
        total_spots=sum(row["spots"] for row in metrics.values()),
        runs=metrics,
        spots_per_accession={acc: row["spots"] for acc, row in metrics.items()},
        sampled_peaks=peaks,
        sample_interval_seconds=5,
    )
    Path("stream_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":

    def interrupted(signum, frame):
        raise KeyboardInterrupt("Received signal " + str(signum))

    signal.signal(signal.SIGTERM, interrupted)
    main()
