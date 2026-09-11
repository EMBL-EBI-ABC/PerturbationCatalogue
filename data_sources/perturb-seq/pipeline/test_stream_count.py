"""Run with python3 -B test_stream_count.py; no SRA connection required."""

import io
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent / "bin"))
from stream_count import pairs


def fastq(spot, read_id, length):
    return f"@SRR1.{spot}/{read_id}\n" + "A" * length + "\n+\n" + "I" * length + "\n"


def test():
    valid = "".join(
        fastq(spot, rid, length)
        for spot in (1, 2)
        for rid, length in ((1, 8), (2, 30), (3, 98))
    )
    records = list(pairs(io.BytesIO(valid.encode()), "SRR1"))
    assert len(records) == 2
    assert [[len(rec[1]) for rec in pair] for pair in records] == [[30, 98], [30, 98]]
    for malformed in (
        "",
        valid[:-2],
        fastq(1, 2, 30),
        fastq(1, 2, 30) * 2,
        fastq(1, 2, 30) + fastq(1, 3, 98) + fastq(1, 4, 98),
        valid + fastq(1, 2, 30),
        valid + fastq(3, 2, 30) + fastq(3, 3, 98),
    ):
        try:
            list(pairs(io.BytesIO(malformed.encode()), "SRR1"))
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid stream accepted")
    print(
        "Stream pairing, technical-read selection and truncated/missing/duplicate input checks passed"
    )


if __name__ == "__main__":
    test()
