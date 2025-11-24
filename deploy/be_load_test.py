"""
Perturbation Catalogue Backend Load Testing Script

This script performs a load test on the Perturbation Catalogue backend API. It simulates real-world usage patterns
to verify the system's stability and autoscaling capabilities under load.

Key Features:
1.  **Concurrency**: Asyncio workers simulating users.
2.  **Simulation**: Realistic queries (genes, modalities) with user think-time.
3.  **Live Dashboard**:
    -   Full terminal width utilization.
    -   Consistent time-binning for graphs.
    -   Live graphs: Error Rate, RPS, Latency Percentiles.
    -   Live Latency Histogram (Capped at 30s).
"""

import asyncio
import aiohttp
import argparse
import random
import time
import sys
import math
import os

# --- Configuration & Constants ---

# ANSI Colors
C_RESET = "\033[0m"
C_RED = "\033[91m"
C_GREEN = "\033[92m"
C_YELLOW = "\033[93m"
C_BLUE = "\033[94m"
C_MAGENTA = "\033[95m"
C_CYAN = "\033[96m"
C_WHITE = "\033[97m"
C_BOLD = "\033[1m"
C_DIM = "\033[2m"

GENES = [
    "TP53",
    "MYC",
    "EGFR",
    "KRAS",
    "BRCA1",
    "BRCA2",
    "AKT1",
    "PIK3CA",
    "PTEN",
    "APC",
    "ESR1",
    "VEGFA",
    "GAPDH",
    "MKI67",
    "CCND1",
    "CDKN2A",
    "MTOR",
    "TNF",
    "IL6",
    "IL1B",
    "CD44",
    "CD8A",
    "CD4",
    "INS",
    "ALB",
    "ACTB",
    "B2M",
    "RPLP0",
    "HPRT1",
    "TUBB",
    "PPIA",
    "TBP",
    "GUSB",
    "UBC",
    "YWHAZ",
    "AR",
    "ATM",
    "BRAF",
    "CCNE1",
    "CDK4",
    "CDK6",
    "CTNNB1",
    "ERBB2",
    "FGFR1",
    "HIF1A",
    "HRAS",
    "IDH1",
    "JAK2",
    "KIT",
    "MAP2K1",
    "MET",
    "NRAS",
    "PDGFRA",
    "RB1",
    "RET",
    "ROS1",
    "SMAD4",
    "SMO",
    "SRC",
    "STK11",
    "TERT",
    "VHL",
    "NOTCH1",
    "FLT3",
    "GNAS",
    "KDR",
    "MPL",
    "NPM1",
    "PTPN11",
    "ABL1",
    "ALK",
    "BCL2",
    "BCL6",
    "CCND2",
    "CD79B",
    "CSF1R",
    "EP300",
    "FGF1",
    "FGF2",
    "FGF3",
    "FGF4",
    "FLT1",
    "FLT4",
    "GATA1",
    "GATA2",
    "HNF1A",
    "IGF1R",
    "JAK1",
    "JAK3",
    "KMT2A",
    "MAPK1",
    "MDM2",
    "MDM4",
    "MLH1",
    "MSH2",
    "MSH6",
    "NF1",
    "NF2",
    "PDGFRB",
    "PIK3R1",
]

MODALITIES = ["perturb-seq", "crispr-screen", "mave"]

# --- Network Logic ---


async def make_request(session, url, params):
    start_time = time.time()
    try:
        async with session.get(url, params=params) as response:
            await response.read()
            duration = time.time() - start_time
            return response.status, duration
    except Exception:
        return 0, time.time() - start_time


async def worker(worker_id, session, base_url, stop_event, stats, wait_time):
    while not stop_event.is_set():
        jitter = random.uniform(-0.1 * wait_time, 0.1 * wait_time)
        actual_wait = max(0, wait_time + jitter)

        if wait_time > 0:
            try:
                await asyncio.sleep(actual_wait)
            except asyncio.CancelledError:
                break

        if stop_event.is_set():
            break

        modality = random.choice(MODALITIES)
        search_url = f"{base_url}/v1/{modality}/search"
        params = {}

        if modality == "perturb-seq":
            query_type = random.choice(["perturbation", "effect", "both"])
            if query_type == "perturbation":
                params["perturbation_gene_name"] = random.choice(GENES)
            elif query_type == "effect":
                params["effect_gene_name"] = random.choice(GENES)
            else:
                params["perturbation_gene_name"] = random.choice(GENES)
                params["effect_gene_name"] = random.choice(GENES)
        else:
            params["perturbation_gene_name"] = random.choice(GENES)

        params["dataset_limit"] = 10
        params["rows_per_dataset_limit"] = 10

        # Record absolute time of request completion
        req_start = time.time()
        status, duration = await make_request(session, search_url, params)
        req_end = time.time()

        is_error = not (200 <= status < 300)

        # Append to history: (timestamp, duration, is_error)
        stats["history"].append((req_end, duration, is_error))

        stats["requests"] += 1
        if not is_error:
            stats["success"] += 1
        else:
            stats["errors"] += 1


# --- Visualization Helpers ---


def get_terminal_width():
    try:
        return os.get_terminal_size().columns
    except OSError:
        return 80


def calculate_bin_config(duration, available_width):
    """Calculates a sensible bin size to fit the duration within available width."""
    # Margins for Y-axis labels (approx 8 chars)
    graph_width = max(10, available_width - 10)

    if duration <= 0:
        return 1.0, graph_width

    ideal_bin_size = duration / graph_width

    # Nice round numbers for bin sizes (seconds)
    candidates = [0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 60.0]

    bin_size = candidates[-1]
    for c in candidates:
        if c >= ideal_bin_size:
            bin_size = c
            break

    return bin_size, graph_width


def draw_binned_graph(
    title, buckets, series_defs, bin_size, num_bins, width=80, height=6
):
    """
    Draws a graph based on pre-calculated buckets.
    series_defs: list of dict {'key': func(bucket), 'color': str, 'char': str, 'label': str}
    """
    # Extract values for Y-scaling
    all_values = []
    for s in series_defs:
        for b in buckets:
            val = s["extractor"](b)
            if val is not None:
                all_values.append(val)

    if not all_values:
        return [f"{title} (Waiting for data...)"]

    min_y = 0
    max_y = max(all_values)
    if max_y == 0:
        max_y = 1

    grid = [[" " for _ in range(num_bins)] for _ in range(height)]

    for s in series_defs:
        color = s.get("color", C_WHITE)
        char = s.get("char", "*")
        extractor = s["extractor"]

        for x, b in enumerate(buckets):
            if x >= num_bins:
                break
            val = extractor(b)
            if val is None:
                continue

            ratio = (val - min_y) / (max_y - min_y)
            y = int(ratio * (height - 1))
            y = max(0, min(height - 1, y))
            row = height - 1 - y
            grid[row][x] = f"{color}{char}{C_RESET}"

    lines = [f"{C_BOLD}{title}{C_RESET}"]
    for r in range(height):
        if r == 0:
            label = f"{max_y:>6.2f} ┤"
        elif r == height - 1:
            label = f"{min_y:>6.2f} ┤"
        else:
            label = "       │"
        row_str = "".join(grid[r])
        lines.append(f"{C_DIM}{label}{C_RESET}{row_str}")

    lines.append(f"{C_DIM}{' '*7}└{'─' * num_bins}{C_RESET}")
    # X-Axis Label
    lines.append(
        f"{C_DIM}{' '*7}0s (1 char = {bin_size}s){' '*(max(0, num_bins-20))}{num_bins*bin_size:.0f}s{C_RESET}"
    )

    legend = []
    for s in series_defs:
        # Get last non-None value
        last_val = 0
        for b in reversed(buckets):
            v = s["extractor"](b)
            if v is not None:
                last_val = v
                break
        legend.append(
            f"{s.get('color')}{s['char']} {s['label']}: {last_val:.2f}{C_RESET}"
        )

    lines.append(" " * 9 + "   ".join(legend))
    return lines


def render_histogram(latencies, height=6, num_bins=60, max_x_val=30.0):
    """Generates a vertical histogram of latencies (Capped at max_x_val)."""
    if not latencies:
        return ["Waiting for data..."]

    min_val = 0.0
    # Cap the max range logic, but actual data might be clamped
    capped_data = [min(l, max_x_val) for l in latencies]

    bins = [0] * num_bins
    bin_step = (max_x_val - min_val) / num_bins

    for l in capped_data:
        # Logic to find bin index
        if l >= max_x_val:
            idx = num_bins - 1
        else:
            idx = int((l - min_val) / bin_step)

        if 0 <= idx < num_bins:
            bins[idx] += 1

    max_count = max(bins)
    if max_count == 0:
        return ["No data"]

    lines = [f"{C_BOLD}LATENCY DISTRIBUTION (Max 30s){C_RESET}"]

    for r in range(height - 1, -1, -1):
        row_chars = ""
        threshold = r / (height - 1) if height > 1 else 0

        if r == height - 1:
            label = f"{max_count:>5} ┤"
        elif r == 0:
            label = f"{0:>5} ┤"
        else:
            label = "      │"

        for count in bins:
            norm_count = count / max_count if max_count else 0
            if norm_count > threshold:
                row_chars += "█"
            else:
                row_chars += " "
        lines.append(f"{C_DIM}{label}{C_RESET}{C_BLUE}{row_chars}{C_RESET}")

    lines.append(f"{C_DIM}      └{'─' * num_bins}{C_RESET}")
    x_label = f"      0s{' ' * (num_bins - 8)}{max_x_val:.1f}s"
    lines.append(f"{C_DIM}{x_label}{C_RESET}")

    return lines


async def reporter(stop_event, stats, start_time, args):
    """Main loop for the live dashboard."""

    # Calculate Layout Config once
    term_width = get_terminal_width()
    bin_size, graph_width = calculate_bin_config(args.duration, term_width)

    # Initialize Buckets
    # Each bucket: { 'reqs': 0, 'errs': 0, 'lats': [] }
    # We might need more bins if duration isn't perfectly divisible, but max_bins is purely for width calculation
    total_bins = math.ceil(args.duration / bin_size)
    buckets = [{"reqs": 0, "errs": 0, "lats": []} for _ in range(total_bins + 1)]

    last_processed_idx = 0

    # Initial clear
    sys.stdout.write("\033[2J\033[3J\033[H")

    while not stop_event.is_set():
        await asyncio.sleep(0.5)

        # Dynamic Terminal Resize Support (simple check)
        curr_width = get_terminal_width()
        if curr_width != term_width:
            term_width = curr_width
            bin_size, graph_width = calculate_bin_config(args.duration, term_width)

        # Process new events
        history = stats["history"]
        curr_len = len(history)

        # Aggregate into buckets
        for i in range(last_processed_idx, curr_len):
            ts, duration, is_error = history[i]
            rel_time = ts - start_time

            if rel_time < 0:
                rel_time = 0
            b_idx = int(rel_time / bin_size)

            if b_idx < len(buckets):
                buckets[b_idx]["reqs"] += 1
                if is_error:
                    buckets[b_idx]["errs"] += 1
                buckets[b_idx]["lats"].append(duration)

        last_processed_idx = curr_len

        # --- Render ---

        sys.stdout.write("\033[2J\033[3J\033[H")

        lines = []
        lines.append(f"{C_BOLD}{C_CYAN}PERTURBATION CATALOGUE LOAD TEST{C_RESET}")
        lines.append(f"{C_DIM}{'-'*term_width}{C_RESET}")

        curr_time = time.time() - start_time

        # Settings (Two Lines)
        lines.append(f"Target: {C_BOLD}{args.url}{C_RESET}")
        lines.append(
            f"Workers: {args.concurrency} | Wait: {args.wait}s | Duration: {args.duration}s"
        )

        lines.append("")

        # Metrics
        rps_overall = stats["requests"] / curr_time if curr_time > 0 else 0
        lines.append(
            f"Total Reqs: {C_BOLD}{C_BLUE}{stats['requests']}{C_RESET} | "
            f"Total Succ: {C_BOLD}{C_GREEN}{stats['success']}{C_RESET} | "
            f"Total Errs: {C_BOLD}{C_RED}{stats['errors']}{C_RESET} | "
            f"Overall RPS: {C_BOLD}{C_CYAN}{rps_overall:.1f}{C_RESET}"
        )
        lines.append(f"{C_DIM}{'-'*term_width}{C_RESET}")
        lines.append("")

        # --- Extractor Functions ---
        def get_error_rate(b):
            return b["errs"] / b["reqs"] if b["reqs"] > 0 else 0.0

        def get_rps(b):
            return b["reqs"] / bin_size

        def get_lat_percentile(b, p):
            if not b["lats"]:
                return None
            b["lats"].sort()  # Timsort is efficient on sorted data
            idx = int(len(b["lats"]) * p)
            return b["lats"][min(idx, len(b["lats"]) - 1)]

        # Graph 1: Error Rate
        lines.extend(
            draw_binned_graph(
                "ERROR RATE (Errors/Total per bin)",
                buckets,
                [
                    {
                        "extractor": get_error_rate,
                        "color": C_RED,
                        "char": "x",
                        "label": "ErrRate",
                    }
                ],
                bin_size,
                len(buckets),
                width=term_width,
                height=6,
            )
        )
        lines.append("")

        # Graph 2: RPS
        lines.extend(
            draw_binned_graph(
                "RPS (Requests/sec)",
                buckets,
                [{"extractor": get_rps, "color": C_CYAN, "char": "≈", "label": "RPS"}],
                bin_size,
                len(buckets),
                width=term_width,
                height=6,
            )
        )
        lines.append("")

        # Graph 3: Latencies
        lines.extend(
            draw_binned_graph(
                "LATENCY PERCENTILES",
                buckets,
                [
                    {
                        "extractor": lambda b: get_lat_percentile(b, 1.0),
                        "color": C_RED,
                        "char": "^",
                        "label": "Max",
                    },
                    {
                        "extractor": lambda b: get_lat_percentile(b, 0.99),
                        "color": C_MAGENTA,
                        "char": "~",
                        "label": "P99",
                    },
                    {
                        "extractor": lambda b: get_lat_percentile(b, 0.95),
                        "color": C_YELLOW,
                        "char": "-",
                        "label": "P95",
                    },
                    {
                        "extractor": lambda b: get_lat_percentile(b, 0.50),
                        "color": C_GREEN,
                        "char": ".",
                        "label": "P50",
                    },
                ],
                bin_size,
                len(buckets),
                width=term_width,
                height=6,
            )
        )
        lines.append("")

        # Histogram
        # Get all latencies from stats for the global histogram
        all_lats = [x[1] for x in history]
        # Use graph_width for num_bins to align with other graphs
        lines.extend(
            render_histogram(all_lats, height=6, num_bins=graph_width, max_x_val=30.0)
        )

        sys.stdout.write("\n".join(lines))
        sys.stdout.flush()


async def main():
    parser = argparse.ArgumentParser(
        description="Load test the Perturbation Catalogue BE."
    )
    parser.add_argument("url", help="Deployment URL")
    parser.add_argument("concurrency", type=int, help="Number of concurrent requests")
    parser.add_argument(
        "--duration", type=int, default=30, help="Duration in seconds (default: 30)"
    )
    parser.add_argument(
        "--wait",
        type=float,
        default=5.0,
        help="Wait time between requests (default: 5.0)",
    )

    args = parser.parse_args()
    base_url = args.url.rstrip("/")

    # stats['history'] stores list of (timestamp, duration, is_error)
    stats = {"requests": 0, "success": 0, "errors": 0, "history": []}

    stop_event = asyncio.Event()
    start_time = time.time()

    connector = aiohttp.TCPConnector(limit=args.concurrency + 10)

    async with aiohttp.ClientSession(connector=connector) as session:
        workers = [
            asyncio.create_task(
                worker(i, session, base_url, stop_event, stats, args.wait)
            )
            for i in range(args.concurrency)
        ]

        reporter_task = asyncio.create_task(
            reporter(stop_event, stats, start_time, args)
        )

        try:
            await asyncio.sleep(args.duration)
        except KeyboardInterrupt:
            pass

        stop_event.set()

        reporter_task.cancel()
        try:
            await reporter_task
        except asyncio.CancelledError:
            pass

        await asyncio.gather(*workers)

    print(f"\n{C_BOLD}TEST COMPLETE{C_RESET}")
    print(f"Total Requests: {stats['requests']}")
    print(f"Errors: {stats['errors']}")
    if stats["history"]:
        lats = [x[1] for x in stats["history"]]
        avg = sum(lats) / len(lats)
        print(f"Avg Latency: {avg:.4f}s")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        sys.exit(0)
