#!/usr/bin/env python3
"""Fetch frozen AlphaFold DB monomer CIFs for the 42 Phase 1A genes."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/phase1B0/alphafold_db"
META = ROOT / "metadata/phase1B0_afdb_freeze.json"
COV = pd.read_csv(ROOT / "results/phase1A/coverage_by_gene.tsv", sep="\t", low_memory=False)
GENES = set(COV.loc[(COV.primary_cohort_eligible == "YES") &
                    (COV.passes_ge70_representation_threshold == "YES"), "gene"])
MAP = pd.read_csv(ROOT / "results/phase1A/protein_mapping.tsv", sep="\t", low_memory=False)
ACC = sorted(set(MAP.loc[MAP.gene.isin(GENES) &
                         MAP.mapping_status.isin(["EXACT_CANONICAL", "SAFE_ISOFORM_REMAP"]),
                         "canonical_UniProt_accession"].dropna().astype(str)) - {""})
OUT.mkdir(parents=True, exist_ok=True)
META.parent.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(acc: str) -> dict:
    sess = requests.Session()
    api = sess.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=(15, 45))
    api.raise_for_status()
    hits = api.json()
    if not hits:
        return {"accession": acc, "status": "NO_AFDB_ENTRY"}
    hit = next((x for x in hits if x.get("uniprotAccession") == acc), hits[0])
    url = hit.get("cifUrl") or ""
    version = hit.get("latestVersion")
    if not url:
        return {"accession": acc, "status": "NO_CIF_URL", "api_entry": hit}
    out = OUT / f"AF-{acc}-F1-model_v{version}.cif"
    if out.exists():
        return {"accession": acc, "status": "REUSED", "version": version,
                "url": url, "path": str(out.relative_to(ROOT)),
                "bytes": out.stat().st_size, "sha256": sha256(out)}
    part = out.with_suffix(out.suffix + ".part")
    r = sess.get(url, stream=True, timeout=(20, 120))
    r.raise_for_status()
    h = hashlib.sha256(); n = 0
    with part.open("wb") as f:
        for chunk in r.iter_content(1024 * 1024):
            if chunk:
                f.write(chunk); h.update(chunk); n += len(chunk)
    part.replace(out)
    return {"accession": acc, "status": "DOWNLOADED", "version": version,
            "url": url, "path": str(out.relative_to(ROOT)), "bytes": n,
            "sha256": h.hexdigest()}


results = []
with ThreadPoolExecutor(max_workers=1) as pool:
    futures = {pool.submit(fetch, acc): acc for acc in ACC}
    for fut in as_completed(futures):
        try:
            results.append(fut.result())
        except Exception as exc:
            results.append({"accession": futures[fut], "status": "FETCH_FAILED",
                            "error": f"{type(exc).__name__}: {exc}"[:300]})
results.sort(key=lambda x: x["accession"])
META.write_text(json.dumps({"phase": "PHASE1B0", "source": "AlphaFold Protein Structure Database",
                            "accessions_requested": ACC, "model_results": results},
                           indent=2) + "\n", encoding="utf-8")
print(json.dumps({"requested": len(ACC),
                  "downloaded_or_reused": sum(x["status"] in {"DOWNLOADED", "REUSED"} for x in results),
                  "failed": [x for x in results if x["status"] not in {"DOWNLOADED", "REUSED"}],
                  "metadata": str(META.relative_to(ROOT))}, indent=2))
