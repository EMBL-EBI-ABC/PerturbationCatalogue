#!/usr/bin/env python3
"""Stream original 10x BAM tags as interleaved barcode/feature FASTQ."""

import argparse
from contextlib import contextmanager
import json
import re
import subprocess
import sys
from urllib.parse import urlsplit
from urllib.request import urlopen


LOCATOR = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc={}&accept-proto=https"
ACCESSION = re.compile(r"(?:SRR|ERR|DRR)\d+$")


def locator(accession):
    with urlopen(LOCATOR.format(accession), timeout=60) as response:
        result = json.load(response)["result"]
    files = [
        file
        for bundle in result
        if bundle.get("status") == 200
        for file in bundle.get("files", [])
        if file.get("accession") == accession
        and file.get("type", "").lower() in {"bam", "tenx"}
    ]
    urls = [
        location["link"]
        for file in files
        for location in file.get("locations", [])
        if urlsplit(location.get("link", "")).scheme == "https"
    ]
    if len(files) != 1 or not urls:
        raise RuntimeError(f"No unique HTTPS BAM locator result for {accession}")
    return urls[0]


def source_parts(source):
    if source.startswith("BAMFILE:"):
        value = source[len("BAMFILE:") :]
        path, _, group = value.partition("#")
        return path, int(group) if group else None
    if not source.startswith("BAM:"):
        raise ValueError(f"Invalid BAM source: {source}")
    value = source[len("BAM:") :]
    value, _, group = value.partition("#")
    if ACCESSION.fullmatch(value):
        value = locator(value)
    if not urlsplit(value).scheme and not value.startswith("/"):
        raise ValueError(f"BAM source is not an accession, URL or path: {source}")
    if urlsplit(value).scheme != "https" and not value.startswith("/"):
        raise ValueError("BAM source must use HTTPS")
    return value, int(group) if group else None


@contextmanager
def stream_sam(source):
    value, group = source_parts(source)
    curl = None
    handle = None
    if value.startswith("/"):
        handle = open(value, "rb")
    else:
        curl = subprocess.Popen(
            [
                "curl",
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
                "172800",
                "--retry",
                "8",
                "--retry-all-errors",
                "--retry-delay",
                "2",
                value,
            ],
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
        )
        handle = curl.stdout
    sam = subprocess.Popen(
        ["samtools", "view", "-h", "-"],
        stdin=handle,
        stdout=subprocess.PIPE,
        text=True,
        errors="replace",
    )
    if curl:
        handle.close()
    try:
        yield sam.stdout, group
    finally:
        sam.stdout.close()
        sam_code = sam.wait()
        if curl:
            curl_code = curl.wait()
            if curl_code:
                raise RuntimeError(f"curl exited with status {curl_code}")
        if sam_code:
            raise RuntimeError(f"samtools view exited with status {sam_code}")
        if handle and not handle.closed:
            handle.close()


def parse_specs(headers):
    specs = {}
    pattern = re.compile(r"^@CO\t10x_bam_to_fastq:(\S+)\((\S+)\)$")
    for line in headers:
        match = pattern.match(line.rstrip("\n"))
        if match:
            specs[match[1]] = [tuple(part.split(":")) for part in match[2].split(",")]
    return specs or {
        "R1": [("SEQ", "QUAL")],
        "R2": [("UR", "UQ")],
        "I1": [("CR", "CQ")],
        "I2": [("BC", "QT")],
    }


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def build_read(spec, tags, sequence, quality, reverse):
    reads, qualities = [], []
    for index, (read_tag, quality_tag) in enumerate(spec):
        last = index == len(spec) - 1
        if read_tag == "SEQ":
            read, qual = sequence, quality
            if reverse:
                read, qual = reverse_complement(read), qual[::-1]
        else:
            read = tags.get(read_tag, "")
            qual = tags.get(quality_tag, "")
            if not read and not last:
                raise RuntimeError(f"BAM record is missing required tag {read_tag}")
            if not qual and not last:
                raise RuntimeError(f"BAM record is missing required tag {quality_tag}")
        reads.append(read)
        qualities.append(qual)
    read, qual = "".join(reads), "".join(qualities)
    if len(read) != len(qual):
        raise RuntimeError("BAM-derived FASTQ sequence and quality lengths differ")
    return read, qual


def tags_from_fields(fields):
    tags = {}
    for item in fields[11:]:
        parts = item.split(":", 2)
        if len(parts) == 3 and parts[1] in {"A", "Z"}:
            tags[parts[0]] = parts[2]
    return tags


def gem_group(tags):
    for key in ("CB", "BX"):
        match = re.search(r"-(\d+)$", tags.get(key, ""))
        if match:
            return int(match[1])
    return None


def choose_reads(reads, workflow, cr11):
    def nonempty(name):
        return reads.get(name, ("", ""))[0]

    r1, r2 = nonempty("R1"), nonempty("R2")
    i1, i2 = nonempty("I1"), nonempty("I2")
    if cr11:
        barcode = i1 + r2
        if not barcode:
            raise RuntimeError(
                "Cell barcode/UMI tags are missing from a Cell Ranger 1.x BAM"
            )
        feature = i2 if workflow == "kite" and i2 else r1
        if not feature:
            raise RuntimeError(
                "No feature read was reconstructed from a Cell Ranger 1.x BAM"
            )
        return (barcode, reads["I1"][1] + reads["R2"][1]), (
            feature,
            reads["I2"][1] if workflow == "kite" and i2 else reads["R1"][1],
        )

    candidates = [(name, value) for name, value in reads.items() if value[0]]
    long_reads = [item for item in candidates if len(item[1][0]) > 40]
    barcode_reads = [item for item in candidates if 20 <= len(item[1][0]) <= 40]
    if len(long_reads) == 1 and barcode_reads:
        barcode = barcode_reads[0][1]
        return barcode, long_reads[0][1]
    raise RuntimeError(
        "Could not identify one 20-40 bp barcode/UMI read and one >40 bp feature read; "
        + repr({name: len(value[0]) for name, value in candidates})
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--workflow", choices=("standard", "kite"), required=True)
    parser.add_argument("--max-spots", type=int, default=0)
    args = parser.parse_args()
    if args.max_spots < 0:
        parser.error("--max-spots must be nonnegative")

    output = bytearray()
    spots = records = 0
    headers = []
    specs = None
    group = None
    try:
        with stream_sam(args.source) as (lines, requested_group):
            group = requested_group
            for line in lines:
                if line.startswith("@"):
                    headers.append(line)
                    continue
                if specs is None:
                    specs = parse_specs(headers)
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                tags = tags_from_fields(fields)
                if group is not None and gem_group(tags) != group:
                    continue
                sequence, quality = fields[9], fields[10]
                if sequence == "*" or quality == "*":
                    continue
                records += 1
                reads = {
                    name: build_read(
                        spec, tags, sequence, quality, bool(int(fields[1]) & 16)
                    )
                    for name, spec in specs.items()
                }
                (barcode, barcode_quality), (feature, feature_quality) = choose_reads(
                    reads,
                    args.workflow,
                    "I1" in specs and "I2" in specs and specs["R2"] == [("UR", "UQ")],
                )
                if not barcode or not feature:
                    continue
                spot = f"bam.{spots + 1}"
                output.extend(f"@{spot}/1\n{barcode}\n+\n{barcode_quality}\n".encode())
                output.extend(f"@{spot}/2\n{feature}\n+\n{feature_quality}\n".encode())
                spots += 1
                if len(output) >= 1 << 20:
                    sys.stdout.buffer.write(output)
                    output.clear()
                if args.max_spots and spots >= args.max_spots:
                    break
    finally:
        if output:
            sys.stdout.buffer.write(output)
        sys.stdout.buffer.flush()
    print(
        json.dumps({"spots": spots, "input_records": records}),
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
