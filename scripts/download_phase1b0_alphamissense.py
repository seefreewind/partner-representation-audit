#!/usr/bin/env python3
"""Fetch the frozen official AlphaMissense amino-acid substitution table."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/raw/phase1B0/AlphaMissense_aa_substitutions.tsv.gz'
URL = 'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz'
EXPECTED_SIZE = 1207278510
EXPECTED_MD5 = 'b9ccb339e0de6cb0a8d1973ad2026576'
OUT.parent.mkdir(parents=True, exist_ok=True)
if OUT.exists():
    h = hashlib.md5(OUT.read_bytes()).hexdigest()
    if OUT.stat().st_size == EXPECTED_SIZE and h == EXPECTED_MD5:
        print(f'reused {OUT}; size={OUT.stat().st_size}; md5={h}')
        raise SystemExit(0)
    raise SystemExit(f'existing file checksum mismatch; preserving it: {OUT}')
part = OUT.with_suffix(OUT.suffix + '.part')
session = requests.Session()
last_error = None
for attempt in range(1, 16):
    offset = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={offset}-", "Accept-Encoding": "identity"} if offset else {"Accept-Encoding": "identity"}
    try:
        with session.get(URL, headers=headers, stream=True, timeout=(30, 180)) as r:
            if offset and r.status_code != 206:
                raise RuntimeError(f"server did not honor Range request at offset {offset}: HTTP {r.status_code}")
            if not offset:
                r.raise_for_status()
            h = hashlib.md5(); n_chunk = 0
            mode = 'ab' if offset else 'wb'
            with part.open(mode) as f:
                for chunk in r.iter_content(chunk_size=4 * 1024 * 1024):
                    if not chunk: continue
                    f.write(chunk); n_chunk += len(chunk)
            print(f"attempt={attempt} resumed_from={offset} received={n_chunk} total={part.stat().st_size}", flush=True)
        if part.stat().st_size == EXPECTED_SIZE:
            break
        if part.stat().st_size > EXPECTED_SIZE:
            raise RuntimeError(f"partial file exceeds expected size: {part.stat().st_size}")
        last_error = RuntimeError("connection ended before the expected object size")
    except Exception as exc:
        last_error = exc
        current = part.stat().st_size if part.exists() else 0
        print(f"attempt={attempt} interrupted_at={current}: {type(exc).__name__}: {exc}", flush=True)
    import time
    time.sleep(min(2 * attempt, 12))
else:
    raise SystemExit(f"download retries exhausted at {part.stat().st_size if part.exists() else 0} bytes: {last_error}")

md5 = hashlib.md5(); sha = hashlib.sha256()
with part.open('rb') as f:
    for chunk in iter(lambda: f.read(4 * 1024 * 1024), b''):
        md5.update(chunk); sha.update(chunk)
if part.stat().st_size != EXPECTED_SIZE or md5.hexdigest() != EXPECTED_MD5:
    raise SystemExit(f'checksum/size mismatch: size={part.stat().st_size}, md5={md5.hexdigest()}')
part.replace(OUT)
print(json.dumps({'path': str(OUT.relative_to(ROOT)), 'bytes': OUT.stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest(), 'zenodo_record': '10.5281/zenodo.10813168'}, indent=2))
