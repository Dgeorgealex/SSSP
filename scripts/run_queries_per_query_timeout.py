#!/usr/bin/env python3
"""
Run query files one query at a time with a hard per-query timeout.

Behavior:
- Reads a query file line-by-line (skips blank lines and comments).
- Executes each query through the existing C++ binary by writing the query
  to a temporary single-query file.
- Enforces timeout per query (default: 300s).
- On timeout, appends a marker entry to:
    time_<basename(graph_filename)>
  in the same folder as the graph file from that query.
"""

from __future__ import annotations

import argparse
import os
import shlex
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


DEFAULT_TIMEOUT_SECONDS = 3000
DEFAULT_BINARY = "../build/Main"


@dataclass
class QueryRunResult:
    index: int
    query_line: str
    returncode: int
    elapsed_seconds: float
    timed_out: bool
    stderr: str


def parse_query_graph_path(query_line: str) -> Optional[str]:
    parts = query_line.split()
    if len(parts) < 4:
        return None

    algo_type = parts[0]
    query_type = parts[1]
    if algo_type == "SSSP":
        # SSSP <check|time> <alg> <graph_file> <source> [iters]
        return parts[3] if query_type in {"check", "time"} else None
    if algo_type == "NegCycle":
        # NegCycle <check|time> <alg> <graph_file> [iters]
        return parts[3] if query_type in {"check", "time"} else None
    return None


def build_timeout_file_path(graph_path: str) -> Path:
    graph = Path(graph_path)
    return graph.parent / f"time_{graph.stem}.txt"


def append_timeout_marker(
    timeout_file: Path,
    query_index: int,
    elapsed_seconds: float,
    timeout_seconds: int,
    query_line: str,
) -> None:
    timeout_file.parent.mkdir(parents=True, exist_ok=True)
    with timeout_file.open("a", encoding="utf-8") as f:
        f.write("\n\n")
        f.write(
            f"TIMEOUT: query={query_index}, elapsed={elapsed_seconds:.2f}s, "
            f"limit={timeout_seconds}s\n"
        )
        f.write(f"QUERY: {query_line}\n")


def run_single_query(
    binary: str,
    config_args: List[str],
    query_line: str,
    query_index: int,
    timeout_seconds: int,
) -> QueryRunResult:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="utf-8") as tmp:
        tmp.write(query_line + "\n")
        single_query_file = tmp.name

    cmd = [binary, *config_args, single_query_file]
    start = time.monotonic()
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            _stdout, stderr = proc.communicate(timeout=timeout_seconds)
            elapsed = time.monotonic() - start
            return QueryRunResult(
                index=query_index,
                query_line=query_line,
                returncode=proc.returncode,
                elapsed_seconds=elapsed,
                timed_out=False,
                stderr=stderr,
            )
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            _stdout, stderr = proc.communicate()
            elapsed = time.monotonic() - start
            return QueryRunResult(
                index=query_index,
                query_line=query_line,
                returncode=-1,
                elapsed_seconds=elapsed,
                timed_out=True,
                stderr=stderr,
            )
    finally:
        try:
            os.unlink(single_query_file)
        except FileNotFoundError:
            pass


def load_queries(queries_file: Path) -> List[str]:
    queries: List[str] = []
    with queries_file.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            queries.append(line)
    return queries


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run each query from a query file separately with a hard timeout."
    )
    parser.add_argument("queries_file", type=Path, help="Path to query file.")
    parser.add_argument(
        "--binary",
        default=DEFAULT_BINARY,
        help=f"Path to C++ binary (default: {DEFAULT_BINARY}).",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=DEFAULT_TIMEOUT_SECONDS,
        help=f"Per-query timeout in seconds (default: {DEFAULT_TIMEOUT_SECONDS}).",
    )
    parser.add_argument(
        "--config-arg",
        action="append",
        default=[],
        help="Config argument to forward to Main (example: --config-arg cycle_detection=1).",
    )
    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    if args.timeout_seconds <= 0:
        print("Error: --timeout-seconds must be > 0", file=sys.stderr)
        return 2

    if not args.queries_file.exists():
        print(f"Error: queries file does not exist: {args.queries_file}", file=sys.stderr)
        return 2

    if not Path(args.binary).exists():
        print(
            f"Error: binary does not exist at {args.binary}. "
            "Build first or pass --binary.",
            file=sys.stderr,
        )
        return 2

    queries = load_queries(args.queries_file)
    print(f"Loaded {len(queries)} queries from {args.queries_file}")
    print(f"Binary: {args.binary}")
    print(f"Timeout per query: {args.timeout_seconds}s")
    if args.config_arg:
        print(f"Config args: {' '.join(shlex.quote(v) for v in args.config_arg)}")

    timed_out = 0
    failed = 0

    for idx, query in enumerate(queries, start=1):
        print(f"[{idx}/{len(queries)}] Running...")
        result = run_single_query(
            binary=args.binary,
            config_args=args.config_arg,
            query_line=query,
            query_index=idx,
            timeout_seconds=args.timeout_seconds,
        )

        if result.timed_out:
            timed_out += 1
            graph_path = parse_query_graph_path(query)
            if graph_path is not None:
                timeout_file = build_timeout_file_path(graph_path)
                append_timeout_marker(
                    timeout_file=timeout_file,
                    query_index=idx,
                    elapsed_seconds=result.elapsed_seconds,
                    timeout_seconds=args.timeout_seconds,
                    query_line=query,
                )
                print(f"  TIMEOUT after {result.elapsed_seconds:.2f}s -> {timeout_file}")
            else:
                print(f"  TIMEOUT after {result.elapsed_seconds:.2f}s (no graph file parsed)")
            continue

        if result.returncode != 0:
            failed += 1
            print(f"  FAILED (exit={result.returncode}, {result.elapsed_seconds:.2f}s)")
            if result.stderr.strip():
                print(f"  stderr: {result.stderr.strip()}")
            continue

        print(f"  OK ({result.elapsed_seconds:.2f}s)")

    print("\nSummary:")
    print(f"  total:    {len(queries)}")
    print(f"  timeout:  {timed_out}")
    print(f"  failed:   {failed}")
    print(f"  success:  {len(queries) - timed_out - failed}")

    return 1 if (timed_out > 0 or failed > 0) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
