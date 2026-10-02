#!/usr/bin/env python3
"""Stream original 10x BAM tags as interleaved barcode/feature FASTQ."""

import argparse
import concurrent.futures
from contextlib import contextmanager
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


LOCATOR = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc={}&accept-proto=https"
ACCESSION = re.compile(r"(?:SRR|ERR|DRR)\d+$")
GEM_GROUP = re.compile(r"-(\d+)$")
BAM_RANGE_SIZE = 128 * 1024**2
BAM_RANGE_WORKERS = 16


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
    return {
        "url": urls[0],
        "size": files[0].get("size"),
        "md5": files[0].get("md5", ""),
    }


def download_ranges(
    url, size, digest, output, chunk_size=BAM_RANGE_SIZE, workers=BAM_RANGE_WORKERS
):
    """Download one locator-pinned object via checked parallel HTTP ranges."""
    parsed = urlsplit(url)
    if (
        parsed.scheme not in {"http", "https"}
        or not parsed.hostname
        or parsed.username
        or parsed.password
        or type(size) is not int
        or size <= 0
        or not re.fullmatch(r"[0-9a-fA-F]{32}", digest or "")
        or chunk_size <= 0
        or workers <= 0
    ):
        raise ValueError("Invalid ranged-download URL, size, checksum or settings")

    output = Path(output)
    partial = output.with_name(output.name + ".partial")
    if any(path.is_symlink() or path.exists() for path in (output, partial)):
        raise FileExistsError(f"Download target already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fd = None
    created_partial = False
    created_output = False
    try:
        fd = os.open(
            partial,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
        )
        created_partial = True
        os.ftruncate(fd, size)
        ranges = [
            (start, min(size - 1, start + chunk_size - 1))
            for start in range(0, size, chunk_size)
        ]

        def fetch(start, end):
            expected = end - start + 1
            for attempt in range(5):
                try:
                    request = Request(url, headers={"Range": f"bytes={start}-{end}"})
                    with urlopen(request, timeout=180) as response:
                        try:
                            content_length = int(
                                response.headers.get("Content-Length", "")
                            )
                        except ValueError as error:
                            raise RuntimeError("missing HTTP range length") from error
                        returned = response.headers.get("Content-Range", "")
                        final_url = urlsplit(response.geturl())
                        if (
                            response.status != 206
                            or returned != f"bytes {start}-{end}/{size}"
                            or content_length != expected
                            or final_url.scheme != parsed.scheme
                        ):
                            raise RuntimeError(
                                f"unexpected HTTP range response: {response.status} {returned}"
                            )
                        offset = start
                        while True:
                            block = response.read(4 * 1024 * 1024)
                            if not block:
                                break
                            view = memoryview(block)
                            while view:
                                written = os.pwrite(fd, view, offset)
                                if not written:
                                    raise OSError("short ranged write")
                                offset += written
                                view = view[written:]
                        if offset != end + 1:
                            raise RuntimeError(
                                f"short HTTP range: {offset - start}/{expected}"
                            )
                    return expected
                except Exception:
                    if attempt == 4:
                        raise
                    time.sleep(2**attempt)

        started = time.monotonic()
        downloaded = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(fetch, start, end) for start, end in ranges]
            for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
                try:
                    downloaded += future.result()
                except BaseException:
                    for pending in futures:
                        pending.cancel()
                    raise
                if size >= 1024**3 and (index % 32 == 0 or index == len(ranges)):
                    elapsed = time.monotonic() - started
                    rate = downloaded / elapsed / 1024**2
                    print(
                        f"source_download ranges={index}/{len(ranges)} "
                        f"bytes={downloaded} MiB_per_second={rate:.1f}",
                        file=sys.stderr,
                        flush=True,
                    )
        if downloaded != size:
            raise RuntimeError(f"Downloaded {downloaded} of {size} bytes")
        os.fsync(fd)
        os.close(fd)
        fd = None

        checksum = hashlib.md5()
        with partial.open("rb") as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                checksum.update(block)
        if partial.stat().st_size != size or checksum.hexdigest() != digest.lower():
            raise RuntimeError("BAM download size or MD5 mismatch")
        os.replace(partial, output)
        created_partial = False
        created_output = True
        return {
            "bytes": size,
            "md5": digest.lower(),
            "ranges": len(ranges),
            "workers": min(workers, len(ranges)),
            "seconds": round(time.monotonic() - started, 1),
        }
    except BaseException:
        if fd is not None:
            os.close(fd)
        if created_partial:
            partial.unlink(missing_ok=True)
        if created_output:
            output.unlink(missing_ok=True)
        raise


def source_parts(source):
    if source.startswith("BAMFILE:"):
        value = source[len("BAMFILE:") :]
        path, _, group = value.partition("#")
        if not path or urlsplit(path).scheme:
            raise ValueError(f"BAMFILE source must be a local path: {source}")
        return str(Path(path).resolve(strict=True)), int(group) if group else None
    if not source.startswith("BAM:"):
        raise ValueError(f"Invalid BAM source: {source}")
    value = source[len("BAM:") :]
    value, _, group = value.partition("#")
    if ACCESSION.fullmatch(value):
        value = locator(value)["url"]
    if not urlsplit(value).scheme and not value.startswith("/"):
        raise ValueError(f"BAM source is not an accession, URL or path: {source}")
    if urlsplit(value).scheme != "https" and not value.startswith("/"):
        raise ValueError("BAM source must use HTTPS")
    return value, int(group) if group else None


@contextmanager
def stream_sam(source, tolerate_sigpipe=False, threads=0):
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
    samtools = ["samtools", "view"]
    if threads:
        samtools += ["-@", str(threads)]
    sam = subprocess.Popen(
        samtools + ["-h", "-"],
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
        if sam_code and not (tolerate_sigpipe and sam_code == -13):
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
        match = GEM_GROUP.search(tags.get(key, ""))
        if match:
            return int(match[1])
    return None


def sam_record_group(fields, groups):
    if len(fields) < 12:
        return None
    tags = {}
    for field in fields[11:]:
        if field.startswith(("CB:", "BX:")):
            parts = field.split(":", 2)
            if len(parts) == 3 and parts[1] in {"A", "Z"}:
                tags[parts[0]] = parts[2]
    group = gem_group(tags)
    return group if group in groups else None


def download_source_bam(source, output_dir):
    if source.startswith("BAMFILE:"):
        return None, None
    if not source.startswith("BAM:"):
        raise ValueError("--download-source is only valid for BAM inputs")
    value, _, group = source[len("BAM:") :].partition("#")
    if group:
        raise ValueError("--download-source cannot select one GEM group")
    if value.startswith("/"):
        return None, None
    if not ACCESSION.fullmatch(value):
        raise ValueError(
            "--download-source needs a BAM accession with locator metadata"
        )
    pinned = locator(value)
    if (
        type(pinned["size"]) is not int
        or pinned["size"] <= 0
        or not re.fullmatch(r"[0-9a-fA-F]{32}", pinned["md5"] or "")
        or urlsplit(pinned["url"]).scheme != "https"
    ):
        raise RuntimeError(f"BAM locator lacks a valid size, MD5 or HTTPS URL: {value}")
    output_dir = Path(output_dir).resolve()
    staged = output_dir.parent / f"{value}.source.bam"
    metrics = download_ranges(pinned["url"], pinned["size"], pinned["md5"], staged)
    metrics["accession"] = value
    return staged, metrics


def split_bam_by_gem_group(
    source, groups, output_dir, threads=0, download_source=False
):
    groups = tuple(groups)
    if (
        not groups
        or len(set(groups)) != len(groups)
        or any(group < 1 for group in groups)
    ):
        raise ValueError("Gem groups must be unique positive integers")
    selected_groups = set(groups)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = {group: output_dir / f"group_{group}.bam" for group in groups}
    temporary = {group: output_dir / f"group_{group}.partial.bam" for group in groups}
    metrics_path = output_dir / "split_metrics.json"
    metrics_temporary = output_dir / "split_metrics.json.partial"
    if any(
        path.exists()
        for path in (
            *targets.values(),
            *temporary.values(),
            metrics_path,
            metrics_temporary,
        )
    ):
        raise FileExistsError("Gem-group output already exists")

    writers = {}
    staged_source = None
    records = dict.fromkeys(groups, 0)
    unassigned_records = 0
    try:
        if download_source:
            staged_source, download_metrics = download_source_bam(source, output_dir)
        else:
            download_metrics = None
        split_source = f"BAMFILE:{staged_source}" if staged_source else source
        with stream_sam(split_source, threads=threads) as (lines, requested_group):
            if requested_group is not None:
                raise ValueError("The split source must not select one gem group")
            headers = []
            for line in lines:
                if line.startswith("@"):
                    headers.append(line)
                    continue
                if not writers:
                    for group in groups:
                        writers[group] = subprocess.Popen(
                            [
                                "samtools",
                                "view",
                                "-b",
                                "-1",
                                "-o",
                                str(temporary[group]),
                                "-",
                            ],
                            stdin=subprocess.PIPE,
                            start_new_session=True,
                        )
                        for header in headers:
                            writers[group].stdin.write(header.encode())
                group = sam_record_group(line.rstrip("\n").split("\t"), selected_groups)
                if group is not None:
                    writers[group].stdin.write(line.encode())
                    records[group] += 1
                else:
                    unassigned_records += 1
        if not writers:
            raise RuntimeError(f"No alignment records found in {source}")
        for writer in writers.values():
            writer.stdin.close()
        for group, writer in writers.items():
            if writer.wait():
                raise RuntimeError(f"samtools failed writing gem group {group}")
        empty = [group for group, count in records.items() if not count]
        if empty:
            raise RuntimeError(f"No alignments for gem groups: {empty}")
        subprocess.run(
            [
                "samtools",
                "quickcheck",
                "-v",
                *(str(temporary[group]) for group in groups),
            ],
            check=True,
        )
        for group in groups:
            temporary[group].replace(targets[group])
        result = {
            "source": source,
            "groups": {
                str(group): {
                    "bam": str(targets[group]),
                    "input_records": records[group],
                }
                for group in groups
            },
            "unassigned_records": unassigned_records,
        }
        if download_metrics:
            result["source_download"] = download_metrics
        metrics_temporary.write_text(json.dumps(result, sort_keys=True) + "\n")
        metrics_temporary.replace(metrics_path)
        if staged_source:
            staged_source.unlink()
    except BaseException:
        for writer in writers.values():
            if writer.poll() is None:
                try:
                    os.killpg(writer.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        for writer in writers.values():
            if writer.stdin and not writer.stdin.closed:
                try:
                    writer.stdin.close()
                except BrokenPipeError:
                    pass
            writer.wait()
        for path in (
            *temporary.values(),
            *targets.values(),
            metrics_temporary,
            metrics_path,
        ):
            path.unlink(missing_ok=True)
        if staged_source:
            staged_source.unlink(missing_ok=True)
        raise
    return result


def choose_reads(reads, workflow, cr11):
    def nonempty(name):
        return reads.get(name, ("", ""))[0]

    r1, r2 = nonempty("R1"), nonempty("R2")
    i1, i2 = nonempty("I1"), nonempty("I2")
    if cr11:
        if not i1 or not r2:
            raise RuntimeError(
                "Cell barcode/UMI tags are missing from a Cell Ranger 1.x BAM"
            )
        # In the legacy 10x v1 layout I2/BC is the sample index.  The
        # captured guide or cDNA is the long SEQ/R1 read in both workflows.
        feature = r1
        if not feature:
            raise RuntimeError(
                "No feature read was reconstructed from a Cell Ranger 1.x BAM"
            )
        return (
            reads["I1"],
            reads["R2"],
            (
                feature,
                reads["R1"][1],
            ),
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


def test_range_download():
    data = bytes(range(256)) * 4096
    seen = []
    lock = threading.Lock()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            match = re.fullmatch(r"bytes=(\d+)-(\d+)", self.headers.get("Range", ""))
            if not match:
                self.send_error(400)
                return
            start, end = map(int, match.groups())
            with lock:
                seen.append((start, end))
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
            self.send_header("Content-Length", str(end - start + 1))
            self.end_headers()
            self.wfile.write(data[start : end + 1])

        def log_message(self, *_):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}/synthetic.bam"
    with tempfile.TemporaryDirectory(prefix="bam-range-test-", dir=Path.cwd()) as temp:
        target = Path(temp) / "synthetic.bam"
        digest = hashlib.md5(data).hexdigest()
        try:
            metrics = download_ranges(
                url, len(data), digest, target, chunk_size=64 * 1024, workers=4
            )
            assert target.read_bytes() == data
            assert metrics["ranges"] == 16 and len(seen) == 16

            rejected = Path(temp) / "bad-checksum.bam"
            try:
                download_ranges(
                    url,
                    len(data),
                    "0" * 32,
                    rejected,
                    chunk_size=64 * 1024,
                    workers=4,
                )
            except RuntimeError as error:
                assert "MD5 mismatch" in str(error)
            else:
                raise AssertionError("incorrect checksum was accepted")
            assert not rejected.exists()
            assert not rejected.with_name(rejected.name + ".partial").exists()

            sentinel = Path(temp) / "preserve.bam"
            sentinel.write_bytes(b"keep")
            try:
                download_ranges(
                    url,
                    len(data),
                    digest,
                    sentinel,
                    chunk_size=64 * 1024,
                    workers=4,
                )
            except FileExistsError:
                pass
            else:
                raise AssertionError("existing output target was overwritten")
            assert sentinel.read_bytes() == b"keep"

            partial = Path(temp) / "occupied.bam.partial"
            partial.write_bytes(b"protect")
            try:
                download_ranges(
                    url,
                    len(data),
                    digest,
                    partial.with_suffix(""),
                    chunk_size=64 * 1024,
                    workers=4,
                )
            except FileExistsError:
                pass
            else:
                raise AssertionError("existing partial target was removed or reused")
            assert partial.read_bytes() == b"protect"
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source")
    parser.add_argument("--workflow", choices=("standard", "kite"))
    parser.add_argument("--feature-offset", type=int, default=0)
    parser.add_argument("--feature-length", type=int, default=0)
    parser.add_argument("--max-spots", type=int, default=0)
    parser.add_argument("--samtools-threads", type=int, default=0)
    parser.add_argument("--split-groups", type=int, default=0)
    parser.add_argument("--output-dir")
    parser.add_argument("--download-source", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.source:
        parser.error("--source is required")
    if (
        args.max_spots < 0
        or args.feature_offset < 0
        or args.feature_length < 0
        or args.samtools_threads < 0
    ):
        parser.error("spot and feature trim values must be nonnegative")
    if args.split_groups:
        if args.split_groups < 1 or not args.output_dir or args.workflow:
            parser.error(
                "split mode needs --split-groups, --output-dir, and no --workflow"
            )
        if args.feature_offset or args.feature_length or args.max_spots:
            parser.error("read trimming and spot limits are not valid in split mode")
        if args.download_source and not args.source.startswith(("BAM:", "BAMFILE:")):
            parser.error("--download-source requires a BAM source")
        print(
            json.dumps(
                split_bam_by_gem_group(
                    args.source,
                    range(1, args.split_groups + 1),
                    args.output_dir,
                    args.samtools_threads,
                    args.download_source,
                )
            ),
            file=sys.stderr,
            flush=True,
        )
        return
    if args.download_source:
        parser.error("--download-source is only valid in split mode")
    if args.output_dir or not args.workflow:
        parser.error("count mode needs --workflow and does not use --output-dir")
    if args.feature_length == 1:
        parser.error("--feature-length must be zero or at least two")
    if args.workflow != "kite" and (args.feature_offset or args.feature_length):
        parser.error("feature trimming is only valid for the kite workflow")

    output = bytearray()
    spots = records = 0
    headers = []
    specs = None
    group = None
    try:
        with stream_sam(
            args.source,
            tolerate_sigpipe=bool(args.max_spots),
            threads=args.samtools_threads,
        ) as (
            lines,
            requested_group,
        ):
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
                cr11 = "I1" in specs and "I2" in specs and specs["R2"] == [("UR", "UQ")]
                selected = choose_reads(reads, args.workflow, cr11)
                if cr11:
                    (
                        (cell, cell_quality),
                        (umi, umi_quality),
                        (
                            feature,
                            feature_quality,
                        ),
                    ) = selected
                    output_reads = [
                        (cell, cell_quality),
                        (umi, umi_quality),
                        (feature, feature_quality),
                    ]
                else:
                    (barcode, barcode_quality), (feature, feature_quality) = selected
                    output_reads = [
                        (barcode, barcode_quality),
                        (feature, feature_quality),
                    ]
                if args.feature_offset or args.feature_length:
                    end = (
                        args.feature_offset + args.feature_length
                        if args.feature_length
                        else len(feature)
                    )
                    if end > len(feature):
                        raise RuntimeError(
                            "Feature trim exceeds reconstructed feature read length"
                        )
                    feature = feature[args.feature_offset : end]
                    feature_quality = feature_quality[args.feature_offset : end]
                    output_reads[-1] = (feature, feature_quality)
                if any(not read for read, _ in output_reads):
                    continue
                spot = f"bam.{spots + 1}"
                for index, (read, quality) in enumerate(output_reads, 1):
                    output.extend(f"@{spot}/{index}\n{read}\n+\n{quality}\n".encode())
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


def self_test():
    if not shutil.which("samtools"):
        raise RuntimeError("--self-test requires samtools")
    test_range_download()

    groups = set(range(1, 11))
    fields = [
        "spot",
        "0",
        "*",
        "0",
        "0",
        "*",
        "*",
        "0",
        "0",
        "ACGT",
        "IIII",
        "CB:Z:ACGT-7",
        "UR:Z:UMI1",
        "CR:Z:ACGT",
        "CQ:Z:IIII",
        "UQ:Z:IIII",
        "BC:Z:ACGT",
        "QT:Z:IIII",
    ]
    assert sam_record_group(fields, groups) == 7
    assert sam_record_group(fields, {1, 2}) is None
    fallback = fields[:11] + ["BX:Z:ACGT-3"]
    assert sam_record_group(fallback, groups) == 3

    specs = {
        "R1": [("SEQ", "QUAL")],
        "R2": [("UR", "UQ")],
        "I1": [("CR", "CQ")],
        "I2": [("BC", "QT")],
    }
    tags = tags_from_fields(fields)
    reads = {
        name: build_read(spec, tags, "A" * 80, "I" * 80, False)
        for name, spec in specs.items()
    }
    cell, umi, (feature, _) = choose_reads(reads, "standard", True)
    assert (cell[0], umi[0], feature) == ("ACGT", "UMI1", "A" * 80)

    def record(name, cell_group=None, bx_group=None):
        fields = [
            name,
            "0",
            "chr1",
            "1",
            "60",
            "80M",
            "*",
            "0",
            "0",
            "A" * 80,
            "I" * 80,
            "UR:Z:ACGTACGTAC",
            "UQ:Z:IIIIIIIIII",
            "CR:Z:ACGTACGTACGTACGT",
            "CQ:Z:IIIIIIIIIIIIIIII",
            "BC:Z:ACGT",
            "QT:Z:IIII",
        ]
        if cell_group is not None:
            fields.append(f"CB:Z:ACGTACGTACGTACGT{cell_group}")
        if bx_group is not None:
            fields.append(f"BX:Z:ACGTACGTACGTACGT-{bx_group}")
        return "\t".join(fields) + "\n"

    with tempfile.TemporaryDirectory(
        prefix="bam-split-self-test-", dir=Path.cwd()
    ) as temp:
        temp = Path(temp)
        sam_path, bam_path = temp / "input.sam", temp / "input.bam"
        sam_path.write_text(
            "@HD\tVN:1.6\tSO:unsorted\n"
            "@SQ\tSN:chr1\tLN:100000\n"
            "@CO\t10x_bam_to_fastq:R1(SEQ:QUAL)\n"
            "@CO\t10x_bam_to_fastq:R2(UR:UQ)\n"
            "@CO\t10x_bam_to_fastq:I1(CR:CQ)\n"
            "@CO\t10x_bam_to_fastq:I2(BC:QT)\n"
            + record("group1a", "-1")
            + record("group2", "-2", 9)
            + record("group1b", "-1")
            + record("bx_fallback", None, 2)
            + record("unassigned", "-3")
        )
        subprocess.run(
            ["samtools", "view", "-b", "-o", str(bam_path), str(sam_path)],
            check=True,
        )
        relative_bam_path = bam_path.relative_to(Path.cwd())
        result = split_bam_by_gem_group(
            f"BAMFILE:{relative_bam_path}", (1, 2), temp / "groups", threads=1
        )
        flagged_result = split_bam_by_gem_group(
            f"BAMFILE:{relative_bam_path}",
            (1, 2),
            temp / "groups-with-download-flag",
            threads=1,
            download_source=True,
        )
        assert result["groups"]["1"]["input_records"] == 2
        assert result["groups"]["2"]["input_records"] == 2
        assert {
            key: value["input_records"]
            for key, value in flagged_result["groups"].items()
        } == {key: value["input_records"] for key, value in result["groups"].items()}
        assert result["unassigned_records"] == 1
        metrics = json.loads((temp / "groups" / "split_metrics.json").read_text())
        assert metrics == result

        script = str(Path(__file__).resolve())
        for group in (1, 2):
            original = subprocess.run(
                [
                    sys.executable,
                    script,
                    "--source",
                    f"BAMFILE:{relative_bam_path}#{group}",
                    "--workflow",
                    "standard",
                ],
                capture_output=True,
                check=True,
            )
            split = subprocess.run(
                [
                    sys.executable,
                    script,
                    "--source",
                    f"BAMFILE:{result['groups'][str(group)]['bam']}",
                    "--workflow",
                    "standard",
                ],
                capture_output=True,
                check=True,
            )
            assert split.stdout == original.stdout
            flagged = subprocess.run(
                [
                    sys.executable,
                    script,
                    "--source",
                    f"BAMFILE:{flagged_result['groups'][str(group)]['bam']}",
                    "--workflow",
                    "standard",
                ],
                capture_output=True,
                check=True,
            )
            assert flagged.stdout == original.stdout


if __name__ == "__main__":
    main()
