#!/usr/bin/env python3
"""Auditable, resumable HTTP range downloader for large public snapshots."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import sys
import time
from datetime import datetime, timezone

import requests


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_chunk(args):
    url, folder, index, start, end, retries = args
    target = pathlib.Path(folder) / ("part_%06d" % index)
    expected = end - start + 1
    if target.exists() and target.stat().st_size == expected:
        return str(target), expected
    headers = {"Range": "bytes=%d-%d" % (start, end), "User-Agent": "phase0-feasibility-audit/1.0"}
    for attempt in range(retries):
        try:
            with requests.get(url, headers=headers, stream=True, timeout=(30, 180)) as r:
                if r.status_code != 206:
                    raise RuntimeError("expected HTTP 206, received %s" % r.status_code)
                content_range = r.headers.get("Content-Range", "")
                if not content_range.startswith("bytes %d-%d/" % (start, end)):
                    raise RuntimeError("unexpected Content-Range: %s" % content_range)
                tmp = target.with_suffix(".partial")
                n = 0
                with open(tmp, "wb") as f:
                    for block in r.iter_content(1024 * 1024):
                        if block:
                            f.write(block)
                            n += len(block)
                if n != expected:
                    raise RuntimeError("short range: %d != %d" % (n, expected))
                os.replace(tmp, target)
                return str(target), n
        except Exception:
            if attempt + 1 >= retries:
                raise
            time.sleep(min(2 ** attempt, 30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--size", required=True, type=int)
    ap.add_argument("--chunk-mib", type=int, default=16)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--retries", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    out = pathlib.Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not args.force:
        sys.exit("Refusing to overwrite existing file: %s (use --force explicitly)" % out)
    chunk_bytes = args.chunk_mib * 1024 * 1024
    folder = out.parent.parent / "interim" / "download_chunks" / out.name
    folder.mkdir(parents=True, exist_ok=True)
    ranges = []
    for i, start in enumerate(range(0, args.size, chunk_bytes)):
        ranges.append((args.url, str(folder), i, start, min(start + chunk_bytes, args.size) - 1, args.retries))
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(fetch_chunk, job) for job in ranges]
        for fut in concurrent.futures.as_completed(futures):
            _, n = fut.result()
            completed += n
            print("completed_bytes=%d/%d" % (completed, args.size), flush=True)
    tmp = out.with_suffix(out.suffix + ".assembling")
    with open(tmp, "wb") as dest:
        for i in range(len(ranges)):
            part = folder / ("part_%06d" % i)
            with open(part, "rb") as src:
                for block in iter(lambda: src.read(8 * 1024 * 1024), b""):
                    dest.write(block)
    if tmp.stat().st_size != args.size:
        sys.exit("assembled size mismatch: %d != %d" % (tmp.stat().st_size, args.size))
    os.replace(tmp, out)
    meta = {
        "url": args.url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "bytes": out.stat().st_size,
        "sha256": sha256(out),
        "chunk_mib": args.chunk_mib,
        "workers": args.workers,
    }
    sidecar = out.with_suffix(out.suffix + ".download.json")
    with open(sidecar, "x" if not sidecar.exists() else "w") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
