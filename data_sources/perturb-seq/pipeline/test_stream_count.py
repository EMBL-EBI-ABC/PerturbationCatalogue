"""Run with python3 -B test_stream_count.py; compile router and exercise bounded handoffs."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

BIN = Path(__file__).parent.resolve() / "bin"


def fastq(spot, read_id, length, accession="SRR1"):
    return (
        f"@{accession}.{spot}/{read_id}\n"
        + "ACGT"[read_id % 4] * length
        + "\n+\n"
        + "I" * length
        + "\n"
    )


def test():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        router = root / "router"
        subprocess.run(
            [
                "g++",
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                str(BIN / "read_router.cpp"),
                "-o",
                str(router),
            ],
            check=True,
        )
        valid = "".join(
            fastq(spot, rid, length)
            for spot in (1, 2)
            for rid, length in ((1, 8), (2, 30), (3, 98))
        )
        expected = "".join(
            fastq(spot, rid, length)
            for spot in (1, 2)
            for rid, length in ((2, 30), (3, 98))
        )
        source = root / "reads.fastq"
        source.write_text(valid)
        result = subprocess.run(
            [str(router), "SRR1", str(source)],
            capture_output=True,
            check=True,
            text=True,
        )
        assert result.stdout == expected
        assert json.loads(result.stderr) == dict(spots=2, input_records=6)
        for malformed in (
            "",
            valid[:-1],
            valid[:-2],
            fastq(1, 2, 30),
            fastq(1, 2, 30) * 2,
            fastq(1, 2, 30) + fastq(1, 3, 98) + fastq(1, 4, 98),
            valid + valid,
            valid + fastq(3, 2, 30) + fastq(3, 3, 98),
            valid.replace(".2/", ".3/"),
        ):
            source.write_text(malformed)
            assert (
                subprocess.run(
                    [str(router), "SRR1", str(source)], capture_output=True
                ).returncode
                != 0
            )

        # Real subprocesses, queues and compiled routing; only SRA/kb commands are simulated.
        fake = root / "tools"
        fake.mkdir()
        program = """#!/usr/bin/env python3
import json, os, pathlib, sys, time
p = pathlib.Path
args = sys.argv[1:]
name = p(sys.argv[0]).name
if name == 'prefetch':
    base = p(args[args.index('-O')+1]); acc = args[0]
    assert len(list(base.glob('*/*.sra'))) <= 1
    target = base / acc; target.mkdir()
    time.sleep(0.03)
    (target / (acc+'.sra')).write_text(acc)
elif name == 'fasterq-dump':
    target = p(args[args.index('-o')+1]); acc = p(args[0]).name
    assert len(list(target.parent.parent.glob('*/reads.fastq'))) <= 1
    if os.environ.get('FAIL_EXTRACT') == acc: sys.exit(7)
    time.sleep(0.06)
    with target.open('w') as f:
        for spot in range(1, 4001):
            for rid, length in ((1,8),(2,30),(3,98)):
                f.write('@'+acc+'.'+str(spot)+'/'+str(rid)+'\\n'+'A'*length+'\\n+\\n'+'I'*length+'\\n')
    print('spots read : 4,000\\nreads read : 12,000\\nreads written : 12,000', file=sys.stderr)
else:
    nlines = 0
    while True:
        chunk = sys.stdin.buffer.read(65536)
        if not chunk: break
        nlines += chunk.count(b'\\n'); time.sleep(0.006)
    assert nlines % 8 == 0
    p('out/counts_unfiltered').mkdir(parents=True)
    p('out/run_info.json').write_text(json.dumps({'n_processed': nlines//8}))
    p('out/counts_unfiltered/adata.h5ad').write_text('fixture')
"""
        for name in ("prefetch", "fasterq-dump", "kb"):
            path = fake / name
            path.write_text(program)
            path.chmod(0o755)
        env = dict(os.environ, PATH=str(fake) + os.pathsep + os.environ["PATH"])
        command = [
            sys.executable,
            "-B",
            str(BIN / "stream_count.py"),
            "--workflow",
            "kite",
            "--index",
            "index",
            "--t2g",
            "t2g",
            "--cpus",
            "4",
            "--accessions",
            "SRR1",
            "SRR2",
            "SRR3",
            "SRR4",
        ]
        for fail in (False, True):
            work = root / ("failure" if fail else "success")
            work.mkdir()
            result = subprocess.run(
                command,
                cwd=work,
                env=dict(env, FAIL_EXTRACT="SRR2" if fail else ""),
                capture_output=True,
                text=True,
                timeout=30,
            )
            assert bool(result.returncode) == fail, result.stderr
            assert not list(work.glob("sra-*"))
            if not fail:
                m = json.loads((work / "stream_metrics.json").read_text())
                assert m["total_spots"] == 16000 and len(m["runs"]) == 4
                events = [
                    json.loads(line)
                    for line in (work / "stream_events.jsonl").read_text().splitlines()
                ]
                starts = {
                    (e["event"], e.get("accession")): e["timestamp"] for e in events
                }
                assert (
                    starts["download_start", "SRR3"] < starts["feed_complete", "SRR2"]
                )
                assert starts["extract_start", "SRR2"] < starts["feed_complete", "SRR1"]
    print(
        "Native routing, malformed-input rejection, bounded overlap, totals and producer-failure cleanup passed"
    )


if __name__ == "__main__":
    test()
