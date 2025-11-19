"""
Perturbation Catalogue Backend Load Testing Script

This script performs a load test on the Perturbation Catalogue backend API. It simulates real-world usage patterns
to verify the system's stability and autoscaling capabilities under load.

Key Features and Logic:
1.  **Concurrency & Asynchrony**: Utilizes `aiohttp` and `asyncio` to spawn multiple concurrent 'workers',
    representing individual users or client connections interacting with the API simultaneously.

2.  **Realistic Query Simulation**:
    -   Workers randomly select one of the three available modalities: `perturb-seq`, `crispr-screen`, or `mave`.
    -   Queries are constructed using a pre-defined list of 100 common biologically relevant genes (e.g., TP53, KRAS, MYC).
    -   For `perturb-seq` specifically, the script intelligently varies the query structure to test different database indices:
        filtering by perturbation gene, effect gene, or both simultaneously.

3.  **User Pacing**: Includes a configurable 'wait time' (default: 5s) between requests for each worker.
    This mimics human 'think time' or the natural delay between interactions in a frontend application,
    preventing an unrealistic flood of back-to-back requests from a single connection.

4.  **Live Monitoring**: Provides real-time feedback on request rates, latencies, and error rates, allowing
    immediate observation of how the backend responds to the applied load (e.g., seeing latency spikes during scaling).
"""

import asyncio
import aiohttp
import argparse
import random
import time
import sys

# List of 100 likely genes (common cancer drivers, essential genes, etc.)
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


async def make_request(session, url, params):
    start_time = time.time()
    try:
        async with session.get(url, params=params) as response:
            # Read the full response to ensure the request completes
            await response.read()
            duration = time.time() - start_time
            return response.status, duration
    except Exception as e:
        return 0, time.time() - start_time


async def worker(worker_id, session, base_url, stop_event, stats, wait_time):
    while not stop_event.is_set():
        # Sleep first (except maybe first run, but consistent spacing is better)
        # We add a small random jitter to avoid "thundering herd" where all workers wake up exactly at the same time
        jitter = random.uniform(-0.1 * wait_time, 0.1 * wait_time)
        actual_wait = max(0, wait_time + jitter)

        # Check stop event before and after sleep
        if stop_event.is_set():
            break
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
            else:  # both
                params["perturbation_gene_name"] = random.choice(GENES)
                params["effect_gene_name"] = random.choice(GENES)
        else:
            # crispr-screen or mave
            # Focus on perturbation_gene_name as the primary query filter
            params["perturbation_gene_name"] = random.choice(GENES)

        # Add limits to keep response sizes reasonable for a load test
        params["dataset_limit"] = 10
        params["rows_per_dataset_limit"] = 10

        status, duration = await make_request(session, search_url, params)

        # Update stats
        stats["requests"] += 1
        if 200 <= status < 300:
            stats["success"] += 1
        else:
            stats["errors"] += 1

        stats["latencies"].append(duration)


async def reporter(stop_event, stats, start_time):
    """Prints live statistics every second."""
    print(
        f"{ 'Time':<10} | { 'Reqs':<8} | { 'Succ':<8} | { 'Err':<5} | { 'RPS':<6} | { 'AvgLat':<7} | { 'LastLat':<7}"
    )
    print("-" * 70)

    last_requests = 0

    while not stop_event.is_set():
        await asyncio.sleep(1)

        current_time = time.time() - start_time
        total_requests = stats["requests"]
        success = stats["success"]
        errors = stats["errors"]
        latencies = stats["latencies"]

        # Calculate RPS for the last second (approx)
        # For simpler logic, we'll just do overall RPS or instant if we tracked last count
        # Let's do instantaneous RPS based on delta
        delta_req = total_requests - last_requests
        last_requests = total_requests

        avg_latency = sum(latencies) / len(latencies) if latencies else 0.0
        last_latency = latencies[-1] if latencies else 0.0

        print(
            f"{current_time:<10.1f} | {total_requests:<8} | {success:<8} | {errors:<5} | {delta_req:<6} | {avg_latency:<7.3f} | {last_latency:<7.3f}"
        )


async def main():
    parser = argparse.ArgumentParser(
        description="Load test the Perturbation Catalogue BE."
    )
    parser.add_argument("url", help="Deployment URL (e.g. https://my-be.run.app)")
    parser.add_argument("concurrency", type=int, help="Number of concurrent requests")
    parser.add_argument(
        "--duration",
        type=int,
        default=30,
        help="Duration of test in seconds (default: 30)",
    )
    parser.add_argument(
        "--wait",
        type=float,
        default=5.0,
        help="Wait time (seconds) between requests per worker (default: 5.0)",
    )

    args = parser.parse_args()

    base_url = args.url.rstrip("/")

    print(f"Starting load test on {base_url}")
    print(f"Concurrency: {args.concurrency}")
    print(f"Wait Time: {args.wait}s")
    print(f"Duration: {args.duration} seconds")
    print(f"Targeting modalities: {MODALITIES}")
    print("-" * 40)

    stats = {"requests": 0, "success": 0, "errors": 0, "latencies": []}

    stop_event = asyncio.Event()
    start_time = time.time()

    # Use a connector with a higher limit just in case
    connector = aiohttp.TCPConnector(limit=args.concurrency + 10)

    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        # Start workers
        for i in range(args.concurrency):
            tasks.append(
                asyncio.create_task(
                    worker(i, session, base_url, stop_event, stats, args.wait)
                )
            )

        # Start reporter
        reporter_task = asyncio.create_task(reporter(stop_event, stats, start_time))

        # Run for the specified duration
        try:
            await asyncio.sleep(args.duration)
        except KeyboardInterrupt:
            pass  # Handle Ctrl+C gracefully by just setting stop event

        stop_event.set()

        # Wait for all workers to finish
        # We cancel the reporter immediately so it doesn't print one last line that might be confusing or redundant
        reporter_task.cancel()
        try:
            await reporter_task
        except asyncio.CancelledError:
            pass

        await asyncio.gather(*tasks)

    print("\n" + "=" * 40)
    print("FINAL RESULTS")
    print("=" * 40)
    print(f"Total Requests: {stats['requests']}")
    print(f"Successful: {stats['success']}")
    print(f"Errors: {stats['errors']}")

    if stats["latencies"]:
        avg_latency = sum(stats["latencies"]) / len(stats["latencies"])
        max_latency = max(stats["latencies"])
        min_latency = min(stats["latencies"])
        # Calculate percentiles
        sorted_latencies = sorted(stats["latencies"])
        p95_index = int(len(sorted_latencies) * 0.95)
        p99_index = int(len(sorted_latencies) * 0.99)
        p95_latency = sorted_latencies[p95_index] if sorted_latencies else 0
        p99_latency = sorted_latencies[p99_index] if sorted_latencies else 0

        print(f"Min Latency: {min_latency:.4f}s")
        print(f"Average Latency: {avg_latency:.4f}s")
        print(f"P95 Latency: {p95_latency:.4f}s")
        print(f"P99 Latency: {p99_latency:.4f}s")
        print(f"Max Latency: {max_latency:.4f}s")

        total_duration = time.time() - start_time
        rps = stats["requests"] / total_duration
        print(f"Overall Requests Per Second (RPS): {rps:.2f}")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nTest interrupted.")
        sys.exit(0)
