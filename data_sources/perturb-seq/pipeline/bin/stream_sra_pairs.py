#!/usr/bin/env python3
"""Extract two equal SRA row blocks as validated, paired gzip FASTQs."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from collections import Counter


def sequence_rows(vdb_dump, archive):
    result = subprocess.run(
        [str(vdb_dump), "--info", str(archive)],
        check=True,
        capture_output=True,
        text=True,
    )
    matches = re.findall(
        r"^\s*SEQ\s*:\s*([\d,]+)(?:\s+rows)?\s*$", result.stdout, re.M | re.I
    )
    if len(matches) != 1:
        raise RuntimeError("vdb-dump --info did not report one SEQ row count")
    return int(matches[0].replace(",", ""))


def fastq_record(stream, label):
    header = stream.readline()
    if not header:
        return None
    sequence, plus, quality = (stream.readline() for _ in range(3))
    if not sequence or not plus or not quality:
        raise RuntimeError(f"Truncated {label} FASTQ record")
    lines = tuple(line.rstrip("\r\n") for line in (header, sequence, plus, quality))
    if not lines[0].startswith("@") or not lines[2].startswith("+"):
        raise RuntimeError(f"Malformed {label} FASTQ record")
    if len(lines[1]) != len(lines[3]):
        raise RuntimeError(f"{label} FASTQ sequence/quality lengths differ")
    return lines


def spot_name(record, label):
    name = record[0][1:].split(None, 1)[0]
    if not name:
        raise RuntimeError(f"Empty {label} FASTQ read name")
    return re.sub(r"/[12]$", "", name)


def extract_pairs(
    archive, accession, sra_bin, r1_path, r2_path, metrics_path, max_pairs=0
):
    archive = Path(archive).resolve(strict=True)
    sra_bin = Path(sra_bin).resolve(strict=True)
    vdb_dump, fastq_dump = sra_bin / "vdb-dump", sra_bin / "fastq-dump"
    if not vdb_dump.is_file() or not fastq_dump.is_file():
        raise FileNotFoundError("sra-bin must contain vdb-dump and fastq-dump")
    if not re.fullmatch(r"(?:SRR|ERR|DRR)\d+", accession):
        raise ValueError("Invalid SRA accession")

    rows = sequence_rows(vdb_dump, archive)
    if rows < 2 or rows % 2:
        raise RuntimeError(f"Expected an even paired-row count, found {rows}")
    if max_pairs < 0:
        raise ValueError("max-pairs must be nonnegative")
    half = rows // 2
    count_limit = min(max_pairs, half) if max_pairs else half
    outputs = [Path(value).resolve() for value in (r1_path, r2_path, metrics_path)]
    if len(set(outputs)) != 3:
        raise ValueError("FASTQ and metrics paths must be distinct")
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            raise FileExistsError(path)
    partials = [Path(str(path) + ".partial") for path in outputs]
    if any(path.exists() for path in partials):
        raise FileExistsError("A partial output already exists")

    commands = [
        [
            str(fastq_dump),
            "--origfmt",
            "--stdout",
            "-N",
            "1",
            "-X",
            str(count_limit),
            str(archive),
        ],
        [
            str(fastq_dump),
            "--origfmt",
            "--stdout",
            "-N",
            str(half + 1),
            "-X",
            str(half + count_limit),
            str(archive),
        ],
    ]
    processes = []
    count = 0
    r1_lengths, r2_lengths = Counter(), Counter()
    r1_names, r2_names, name_pairs = (
        hashlib.sha256(),
        hashlib.sha256(),
        hashlib.sha256(),
    )
    try:
        for command in commands:
            processes.append(
                subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    text=True,
                    encoding="ascii",
                    errors="strict",
                    bufsize=1024 * 1024,
                )
            )
        with gzip.open(
            partials[0], "wt", encoding="ascii", newline="\n", compresslevel=1
        ) as out1, gzip.open(
            partials[1], "wt", encoding="ascii", newline="\n", compresslevel=1
        ) as out2:
            while True:
                rec1 = fastq_record(processes[0].stdout, "R1")
                rec2 = fastq_record(processes[1].stdout, "R2")
                if rec1 is None or rec2 is None:
                    if rec1 is not None or rec2 is not None:
                        raise RuntimeError(
                            "R1 and R2 FASTQ blocks have different record counts"
                        )
                    break
                name1, name2 = spot_name(rec1, "R1"), spot_name(rec2, "R2")
                if name1 != name2:
                    raise RuntimeError(
                        f"SRA row blocks are not matching pairs: {name1!r} != {name2!r}"
                    )
                for record, output in ((rec1, out1), (rec2, out2)):
                    output.write("\n".join(record) + "\n")
                encoded = name1.encode("ascii")
                r1_names.update(encoded + b"\n")
                r2_names.update(name2.encode("ascii") + b"\n")
                name_pairs.update(encoded + b"\0" + name2.encode("ascii") + b"\n")
                r1_lengths[len(rec1[1])] += 1
                r2_lengths[len(rec2[1])] += 1
                count += 1
        return_codes = [process.wait() for process in processes]
        if return_codes != [0, 0]:
            raise RuntimeError(f"fastq-dump failed: exit codes {return_codes}")
        if count != count_limit:
            raise RuntimeError(
                f"Expected {count_limit} paired records, extracted {count}"
            )
        metrics = {
            "accession": accession,
            "sra_sequence_rows": rows,
            "split_row": half,
            "is_partial": count_limit < half,
            "requested_pairs": max_pairs,
            "paired_records": count,
            "r1_length_counts": {str(k): v for k, v in sorted(r1_lengths.items())},
            "r2_length_counts": {str(k): v for k, v in sorted(r2_lengths.items())},
            "r1_names_sha256": r1_names.hexdigest(),
            "r2_names_sha256": r2_names.hexdigest(),
            "name_pairs_sha256": name_pairs.hexdigest(),
        }
        partials[2].write_text(json.dumps(metrics, sort_keys=True) + "\n")
        for partial, final in zip(partials, outputs):
            partial.replace(final)
        return metrics
    except BaseException:
        for process in processes:
            if process.poll() is None:
                process.kill()
        for process in processes:
            if process.stdout:
                process.stdout.close()
            process.wait()
        for path in (*partials, *outputs):
            path.unlink(missing_ok=True)
        raise
    finally:
        for process in processes:
            if process.stdout and not process.stdout.closed:
                process.stdout.close()


def self_test():
    with tempfile.TemporaryDirectory(prefix="stream-sra-pairs-") as temp:
        root = Path(temp)
        bindir = root / "bin"
        bindir.mkdir()
        archive = root / "SRR1.sra"
        archive.touch()
        vdb_dump = bindir / "vdb-dump"
        vdb_dump.write_text("#!/usr/bin/env python3\nprint('SEQ    : 4')\n")
        vdb_dump.chmod(0o755)
        records = [
            ("spot-a/1", "A" * 28),
            ("spot-b/1", "C" * 28),
            ("spot-a/2", "G" * 47),
            ("spot-b/2", "T" * 58),
        ]
        fastq_dump = bindir / "fastq-dump"
        fastq_dump.write_text(
            "#!/usr/bin/env python3\n"
            "import os, sys\n"
            "args = sys.argv\n"
            "start, end = int(args[args.index('-N')+1]), int(args[args.index('-X')+1])\n"
            "rows = " + repr(records) + "\n"
            "if os.environ.get('MISMATCH') and start == 3: rows[2] = ('wrong-name', rows[2][1])\n"
            "for name, seq in rows[start-1:end]: print('@'+name+'\\n'+seq+'\\n+\\n'+'I'*len(seq))\n"
        )
        fastq_dump.chmod(0o755)
        r1, r2, metrics = (
            root / name for name in ("r1.fastq.gz", "r2.fastq.gz", "metrics.json")
        )
        result = extract_pairs(archive, "SRR1", bindir, r1, r2, metrics)
        assert result["paired_records"] == 2
        assert result["r1_length_counts"] == {"28": 2}
        assert result["r2_length_counts"] == {"47": 1, "58": 1}
        with gzip.open(r1, "rt", encoding="ascii") as in1, gzip.open(
            r2, "rt", encoding="ascii"
        ) as in2:
            pairs = []
            for _ in range(2):
                pairs.append(
                    (fastq_record(in1, "test R1"), fastq_record(in2, "test R2"))
                )
            assert fastq_record(in1, "test R1") is None
            assert fastq_record(in2, "test R2") is None
        assert [spot_name(pair[0], "R1") for pair in pairs] == ["spot-a", "spot-b"]
        assert [spot_name(pair[1], "R2") for pair in pairs] == ["spot-a", "spot-b"]
        assert pairs[0][0][0] == "@spot-a/1" and pairs[0][1][0] == "@spot-a/2"
        assert result["r1_names_sha256"] == result["r2_names_sha256"]
        assert json.loads(metrics.read_text()) == result
        for path in (r1, r2, metrics):
            path.unlink()
        partial = extract_pairs(archive, "SRR1", bindir, r1, r2, metrics, max_pairs=1)
        assert partial["is_partial"] and partial["requested_pairs"] == 1
        assert partial["paired_records"] == 1
        for path in (r1, r2, metrics):
            path.unlink()
        os.environ["MISMATCH"] = "1"
        try:
            try:
                extract_pairs(archive, "SRR1", bindir, r1, r2, metrics)
            except RuntimeError as error:
                assert "not matching pairs" in str(error)
            else:
                raise AssertionError("A mismatched pair was accepted")
            assert not any(path.exists() for path in (r1, r2, metrics))
        finally:
            os.environ.pop("MISMATCH", None)
    print("stream_sra_pairs self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sra")
    parser.add_argument("--accession")
    parser.add_argument("--sra-bin")
    parser.add_argument("--r1")
    parser.add_argument("--r2")
    parser.add_argument("--metrics")
    parser.add_argument("--max-pairs", type=int, default=0)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not all(
        (args.sra, args.accession, args.sra_bin, args.r1, args.r2, args.metrics)
    ):
        parser.error(
            "--sra, --accession, --sra-bin, --r1, --r2 and --metrics are required"
        )
    print(
        json.dumps(
            extract_pairs(
                args.sra,
                args.accession,
                args.sra_bin,
                args.r1,
                args.r2,
                args.metrics,
                args.max_pairs,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
