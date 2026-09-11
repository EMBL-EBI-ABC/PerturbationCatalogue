#!/usr/bin/env python3
"""Download SRA runs and stream validated barcode/biological pairs into kb count."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time


def pairs(stream, accession):
    """Select the two read roles used by this pipeline, preserving whole spots."""
    spot, reads, layout = None, {}, None

    def select():
        nonlocal layout
        current = {rid: len(rec[1]) for rid, rec in reads.items()}
        if layout is None:
            layout = current
        if current != layout:
            raise ValueError(
                f"{accession}: read layout changed at spot {spot}: {current} != {layout}"
            )
        barcode = [rec for rec in reads.values() if 20 <= len(rec[1]) <= 40]
        biological = [rec for rec in reads.values() if len(rec[1]) > 40]
        if len(barcode) != 1 or len(biological) != 1:
            raise ValueError(
                f"{accession}: expected one barcode and one biological read at spot {spot}: {current}"
            )
        return tuple(barcode + biological)

    while True:
        header = stream.readline()
        if not header:
            break
        sequence, plus, quality = [stream.readline() for _ in range(3)]
        match = re.fullmatch(rb"@" + accession.encode() + rb"\.(\d+)/(\d+)\n", header)
        if (
            not match
            or not plus.startswith(b"+")
            or not sequence.endswith(b"\n")
            or not quality.endswith(b"\n")
            or len(sequence) != len(quality)
            or len(sequence) < 2
        ):
            raise ValueError(
                f"{accession}: malformed or truncated FASTQ record {header!r}"
            )
        next_spot, read_id = map(int, match.groups())
        if next_spot != spot:
            if spot is not None:
                if next_spot <= spot:
                    raise ValueError(f"{accession}: spot order is not increasing")
                yield select()
            spot, reads = next_spot, {}
        if read_id in reads:
            raise ValueError(f"{accession}: duplicate read {read_id} in spot {spot}")
        reads[read_id] = (
            header.rstrip(b"\n"),
            sequence.rstrip(b"\n"),
            quality.rstrip(b"\n"),
        )
    if reads:
        yield select()
    else:
        raise ValueError(f"{accession}: empty stream")


def stop(process):
    if process is not None and process.poll() is None:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


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
    prefetch = str(Path(args.sra_bin) / "prefetch") if args.sra_bin else "prefetch"
    fasterq = (
        str(Path(args.sra_bin) / "fasterq-dump") if args.sra_bin else "fasterq-dump"
    )
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
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
        str(args.cpus - 3),
    ]
    command += (
        ["--filter", "bustools"]
        if args.workflow == "standard"
        else ["--workflow", "kite"]
    )
    command += ["-"]
    counts = {}
    # One archive is being extracted and at most one more is being downloaded.
    with tempfile.TemporaryDirectory(prefix="sra-", dir=Path.cwd()) as scratch:
        scratch = Path(scratch)
        pending = counter = extractor = None
        logdir = Path("stream_logs")
        logdir.mkdir()

        def download(acc):
            with (logdir / (acc + ".prefetch.log")).open("wb") as log:
                return subprocess.Popen(
                    [prefetch, acc, "-O", str(scratch), "--max-size", "100G"],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )

        try:
            pending = download(args.accessions[0])
            counter = subprocess.Popen(
                command, stdin=subprocess.PIPE, start_new_session=True
            )
            for index, accession in enumerate(args.accessions):
                if pending.wait(timeout=3600) != 0:
                    raise RuntimeError(f"{accession}: prefetch failed; see stream_logs")
                sra = scratch / accession / (accession + ".sra")
                if not sra.is_file():
                    raise FileNotFoundError(sra)
                pending = (
                    download(args.accessions[index + 1])
                    if index + 1 < len(args.accessions)
                    else None
                )
                with (logdir / (accession + ".fasterq.log")).open("wb") as log:
                    extractor = subprocess.Popen(
                        [
                            fasterq,
                            str(sra),
                            "--split-spot",
                            "--stdout",
                            "--include-technical",
                            "--seq-defline",
                            "@$ac.$si/$ri",
                            "--qual-defline",
                            "+",
                            "-e",
                            "2",
                            "-t",
                            str(scratch),
                        ],
                        stdout=subprocess.PIPE,
                        stderr=log,
                        start_new_session=True,
                    )
                    n_spots = 0
                    for pair in pairs(extractor.stdout, accession):
                        counter.stdin.write(
                            b"".join(
                                h + b"\n" + seq + b"\n+\n" + qual + b"\n"
                                for h, seq, qual in pair
                            )
                        )
                        n_spots += 1
                    extractor.stdout.close()
                    if extractor.wait(timeout=120) != 0:
                        raise RuntimeError(
                            f"{accession}: fasterq-dump failed; see stream_logs"
                        )
                counts[accession] = n_spots
                print(
                    f"{accession}: streamed {n_spots} barcode/biological pairs",
                    flush=True,
                )
                shutil.rmtree(sra.parent)
            counter.stdin.close()
            if counter.wait() != 0:
                raise RuntimeError("kb count failed")
            info = json.loads(Path("out/run_info.json").read_text())
            if info["n_processed"] != sum(counts.values()):
                raise RuntimeError(
                    f'Counter consumed {info["n_processed"]} pairs; producer sent {sum(counts.values())}'
                )
            output = (
                "counts_filtered"
                if args.workflow == "standard"
                else "counts_unfiltered"
            )
            if not (Path("out") / output / "adata.h5ad").is_file():
                raise RuntimeError("Missing count matrix")
        finally:
            for child in (pending, extractor, counter):
                stop(child)
        result = {
            "started_at": started_at,
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - started,
            "spots_per_accession": counts,
            "total_spots": sum(counts.values()),
            "workflow": args.workflow,
            "fastq_files_written": 0,
        }
        Path("stream_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
