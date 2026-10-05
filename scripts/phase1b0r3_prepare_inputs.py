#!/usr/bin/env python3
"""Map frozen variant-partner rows to available experimental assembly chains.

This is a data-preparation wrapper for the frozen MutPred-PPI v1.0.0 graph
builder. It makes no disease-label comparisons and does not edit the release
inference code.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import itertools
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import Atom as BioAtom
from Bio.PDB import Chain as BioChain
from Bio.PDB import MMCIFIO, Polypeptide
from Bio.PDB import Residue as BioResidue
from Bio.Data import PDBData
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure as BioStructure
import gemmi
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase1B0R3"
MODEL_INPUT = OUT / "structure_inputs"
CHECKPOINT_DIR = OUT / "structure_mapping_checkpoints"
UNIVERSE = ROOT / "results/phase1B0/variant_partner_scoring_universe.tsv"
SEQ_DIR = ROOT / "data/interim/phase1A_uniprot_sequences"
CIF_DIR = ROOT / "data/interim/phase1A_parsed_cif"
EDGE_CUTOFF_A = 4.5
MAPPING_VERSION = "gemmi-mmCIF-v2-null-label-seq"
AA = dict(Polypeptide.protein_letters_3to1_extended)
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2.0
ALIGNER.mismatch_score = -1.0
ALIGNER.open_gap_score = -8.0
ALIGNER.extend_gap_score = -0.5
ALIGNMENT_CACHE: dict[tuple[str, str], dict] = {}


def clean_seq(value: object) -> str:
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]", "", str(value or "").upper())


def read_fasta(path: Path) -> str:
    return clean_seq("".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">")))


def map_entity_to_reference(source: str, target: str) -> dict:
    """Map 1-based entity sequence positions to UniProt positions."""
    source, target = clean_seq(source), clean_seq(target)
    cache_key = (source, target)
    if cache_key in ALIGNMENT_CACHE:
        return ALIGNMENT_CACHE[cache_key]
    if not source or not target:
        answer = {"mapping": {}, "identity": 0.0, "coverage": 0.0, "ambiguous_targets": set()}
        ALIGNMENT_CACHE[cache_key] = answer
        return answer
    if source == target:
        answer = {"mapping": {i: i for i in range(1, len(source) + 1)},
                  "identity": 1.0, "coverage": 1.0, "ambiguous_targets": set()}
        ALIGNMENT_CACHE[cache_key] = answer
        return answer
    try:
        alignments = ALIGNER.align(source, target)
        best = alignments[0]
    except Exception:
        answer = {"mapping": {}, "identity": 0.0, "coverage": 0.0, "ambiguous_targets": set()}
        ALIGNMENT_CACHE[cache_key] = answer
        return answer
    mapping = {}
    aligned = identical = 0
    for si, ti in zip(best.indices[0], best.indices[1]):
        if si >= 0 and ti >= 0:
            mapping[int(si) + 1] = int(ti) + 1
            aligned += 1
            identical += source[si] == target[ti]
    ambiguous = set()
    top = best.score
    # PairwiseAlignments can have astronomically many equally optimal paths;
    # len(alignments) may overflow even before iteration starts. Inspect only
    # the first eleven alternatives. Any observed coordinate ambiguity is
    # conservatively excluded downstream.
    for other in itertools.islice(iter(alignments), 1, 12):
        if abs(other.score - top) > 1e-8:
            break
        alt_map = {}
        for si, ti in zip(other.indices[0], other.indices[1]):
            if si >= 0 and ti >= 0:
                alt_map[int(si) + 1] = int(ti) + 1
        for spos, tpos in mapping.items():
            if alt_map.get(spos) != tpos:
                ambiguous.add(tpos)
    answer = {"mapping": mapping, "identity": identical / aligned if aligned else 0.0,
              "coverage": aligned / len(source) if source else 0.0,
              "ambiguous_targets": ambiguous}
    ALIGNMENT_CACHE[cache_key] = answer
    return answer


def acc_sequence(accession: str, seqs: dict[str, str]) -> str:
    accession = str(accession).split("-")[0]
    if accession not in seqs:
        p = SEQ_DIR / f"{accession}.fasta"
        if not p.is_file():
            return ""
        seqs[accession] = read_fasta(p)
    return seqs[accession]


def candidate_ids(value: str) -> list[tuple[str, str]]:
    return [(pid.upper(), f"assembly{assembly}") for pid, assembly in
            re.findall(r"([A-Za-z0-9]{4,})\s*:\s*assembly\s*(\d+)", str(value or ""))]


@lru_cache(maxsize=12)
def parse_cif(path: Path):
    # Gemmi's C++ parser avoids the repeated Python tokenization and object
    # construction cost of Biopython on large biological assemblies.
    structure = gemmi.read_structure(str(path))
    sequences = {
        str(entity.name): clean_seq("".join(AA.get(str(resname).upper(), "X")
                                             for resname in entity.full_sequence))
        for entity in structure.entities
        if entity.entity_type == gemmi.EntityType.Polymer
    }
    asym_entity: dict[str, str] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.subchain and residue.entity_id:
                    asym_entity.setdefault(str(residue.subchain), str(residue.entity_id))
    # Label-asym groups are the coordinate chains consumed by the official
    # graph builder. Assembly CIFs may store several such groups per auth chain.
    chains: dict[str, list] = {}
    for chain in structure[0]:
        for residue in chain:
            label = str(residue.subchain or chain.name)
            try:
                label_seq = int(residue.label_seq)
            except (TypeError, ValueError):
                continue
            if label_seq > 0 and AA.get(str(residue.name).upper()):
                chains.setdefault(label, []).append(residue)
    for label in chains:
        chains[label].sort(key=lambda residue: int(residue.label_seq))
    geometries = {}
    for label, residues in chains.items():
        xyz = np.asarray([[atom.pos.x, atom.pos.y, atom.pos.z]
                          for residue in residues for atom in residue if not atom.is_hydrogen()],
                         dtype=np.float32)
        geometries[label] = (xyz, cKDTree(xyz) if len(xyz) else None)
    return structure, sequences, asym_entity, chains, geometries


def observed_residues(residues: list) -> tuple[list, dict[int, int], dict[int, str]]:
    label_to_index: dict[int, int] = {}
    label_to_aa: dict[int, str] = {}
    for i, residue in enumerate(residues):
        try: label = int(residue.label_seq)
        except (TypeError, ValueError): continue
        label_to_index.setdefault(label, i)
        label_to_aa.setdefault(label, AA.get(str(residue.name).upper(), "X"))
    return residues, label_to_index, label_to_aa


def min_chain_distance(geometry_a: tuple, geometry_b: tuple) -> float:
    xyz_a, tree_a = geometry_a
    xyz_b, _ = geometry_b
    if not len(xyz_a) or not len(xyz_b) or tree_a is None: return float("inf")
    distances, _ = tree_a.query(xyz_b, k=1)
    return float(np.min(distances))


def eligible_entities(entity_seqs: dict[str, str], target_seq: str, partner_seq: str,
                      same_protein: bool) -> tuple[dict, dict]:
    tmap, pmap = {}, {}
    for eid, seq in entity_seqs.items():
        tm = map_entity_to_reference(seq, target_seq)
        tq = tm["identity"] * tm["coverage"]
        tgood = tm["identity"] >= 0.85 and tm["coverage"] >= 0.70
        if same_protein:
            if tgood:
                tmap[eid] = tm
                pmap[eid] = tm
            continue
        pm = map_entity_to_reference(seq, partner_seq)
        pq = pm["identity"] * pm["coverage"]
        pgood = pm["identity"] >= 0.85 and pm["coverage"] >= 0.70
        if tgood and (not pgood or tq >= pq + 0.03):
            tmap[eid] = tm
        elif pgood and (not tgood or pq >= tq + 0.03):
            pmap[eid] = pm
    return tmap, pmap


def variant_coordinate_map(variants: pd.DataFrame, tmap: dict, target_chain: list,
                           target_entity: str) -> dict[str, dict]:
    _, label_index, label_aa = observed_residues(target_chain)
    entity_map = tmap.get(target_entity, {}).get("mapping", {})
    inverse: dict[int, list[int]] = {}
    for entity_pos, canonical_pos in entity_map.items():
        inverse.setdefault(canonical_pos, []).append(entity_pos)
    result = {}
    for row in variants.itertuples(index=False):
        key = str(row.canonical_variant)
        try: canonical_pos = int(row.canonical_residue)
        except Exception:
            result[key] = {"status": "VARIANT_REFERENCE_MISMATCH"}; continue
        variant = str(row.canonical_variant)
        wt, alt = variant[0], variant[-1]
        positions = inverse.get(canonical_pos, [])
        if len(positions) != 1 or canonical_pos in tmap.get(target_entity, {}).get("ambiguous_targets", set()):
            result[key] = {"status": "STRUCTURE_RESIDUE_FAILURE"}; continue
        entity_pos = positions[0]
        if entity_pos not in label_index:
            result[key] = {"status": "STRUCTURE_RESIDUE_FAILURE"}; continue
        observed = label_aa.get(entity_pos, "X")
        if observed != wt:
            result[key] = {"status": "VARIANT_REFERENCE_MISMATCH", "observed_aa": observed}; continue
        result[key] = {
            "status": "READY",
            "local_index": label_index[entity_pos],
            "local_variant": f"{wt}{label_index[entity_pos] + 1}{alt}",
            "entity_position": entity_pos,
            "observed_aa": observed,
        }
    return result


def save_two_chain_cif(chains: dict[str, list], target_chain_id: str,
                       partner_chain_id: str, path: Path) -> None:
    new_structure = BioStructure("selected_pair")
    new_model = Model(0)
    new_structure.add(new_model)
    for new_chain_id, source_chain_id in (("A", target_chain_id), ("B", partner_chain_id)):
        new_chain = BioChain.Chain(new_chain_id)
        for local_index, residue in enumerate(chains[source_chain_id], start=1):
            one_letter = AA.get(str(residue.name).upper())
            resname = PDBData.protein_letters_1to3.get(one_letter, "")
            if not resname:
                continue
            bio_residue = BioResidue.Residue((" ", local_index, " "), resname, "")
            # Resolve alternate locations by highest occupancy, with blank and
            # then A as deterministic tie breakers, and retain heavy atoms only.
            selected_atoms = {}
            for atom in residue:
                if atom.is_hydrogen():
                    continue
                atom_name = str(atom.name).strip()
                altloc = " " if atom.altloc in ("\x00", " ") else str(atom.altloc)
                rank = (float(atom.occ), altloc == " ", altloc == "A")
                if atom_name not in selected_atoms or rank > selected_atoms[atom_name][0]:
                    selected_atoms[atom_name] = (rank, atom, altloc)
            for serial, (atom_name, (_, atom, altloc)) in enumerate(selected_atoms.items(), start=1):
                coord = np.asarray([atom.pos.x, atom.pos.y, atom.pos.z], dtype=np.float32)
                bio_atom = BioAtom.Atom(atom_name, coord, float(atom.b_iso), float(atom.occ),
                                        altloc, atom_name, serial,
                                        element=atom.element.name.upper())
                bio_residue.add(bio_atom)
            if len(bio_residue):
                new_chain.add(bio_residue)
        new_model.add(new_chain)
    writer = MMCIFIO()
    writer.set_structure(new_structure)
    writer.save(str(path))


def write_pair_checkpoint(path: Path, pair_key: str, fingerprint: str,
                          rows: list[dict], edge: dict) -> None:
    record = {"pair_key": pair_key, "fingerprint": fingerprint,
              "mapping_version": MAPPING_VERSION,
              "rows": rows, "edge": edge}
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(record, ensure_ascii=False, default=str))
    tmp.replace(path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    MODEL_INPUT.mkdir(parents=True, exist_ok=True)
    STRUCT_DIR = MODEL_INPUT / "structures"
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    existing_checkpoints = [p for p in CHECKPOINT_DIR.glob("*.json")
                            if not p.name.startswith("._")]
    has_current_checkpoint = False
    checkpoint_records = {}
    for checkpoint in existing_checkpoints:
        try:
            record = json.loads(checkpoint.read_text())
            checkpoint_records[record["pair_key"]] = record
            if record.get("mapping_version") == MAPPING_VERSION:
                has_current_checkpoint = True
        except Exception:
            continue
    # If a previous run has no group checkpoints, its CIFs are unassociated
    # partial artifacts. Clear only those script-generated files. Later runs
    # preserve outputs referenced by completed checkpoints.
    if not has_current_checkpoint:
        for stale in STRUCT_DIR.glob("*.cif"):
            try:
                stale.unlink()
            except FileNotFoundError:
                pass
    universe = pd.read_csv(UNIVERSE, sep="\t", dtype={"VariationID": str})
    if len(universe) != 11150 or universe[["gene", "VariationID"]].drop_duplicates().shape[0] != 5357:
        raise SystemExit("INPUT_FREEZE_VIOLATION: frozen scoring universe dimensions changed")

    seqs: dict[str, str] = {}
    group_keys = ["canonical_UniProt_accession", "partner_accession"]
    edge_rows, row_maps = [], []
    unique_pairs = universe[group_keys].drop_duplicates().reset_index(drop=True)
    # Read the small local sequence cache in one pass before opening large CIFs.
    all_accessions = sorted(set(unique_pairs.iloc[:, 0].astype(str)) |
                            set(unique_pairs.iloc[:, 1].astype(str)))
    for seq_index, accession in enumerate(all_accessions, start=1):
        acc_sequence(accession, seqs)
        if seq_index == 1 or seq_index % 50 == 0 or seq_index == len(all_accessions):
            print(f"sequence_cache={seq_index}/{len(all_accessions)}", flush=True)
    for pair_index, pair in enumerate(unique_pairs.itertuples(index=False), start=1):
        target_acc, partner_acc = str(pair[0]), str(pair[1])
        if pair_index == 1 or pair_index % 10 == 0 or pair_index == len(unique_pairs):
            print(f"preflight_pair={pair_index}/{len(unique_pairs)} target={target_acc} partner={partner_acc}", flush=True)
        rows = universe[(universe.canonical_UniProt_accession.astype(str) == target_acc) &
                        (universe.partner_accession.astype(str) == partner_acc)].copy()
        pair_key = f"{target_acc}|{partner_acc}"
        fingerprint = hashlib.sha256(rows.to_csv(sep="\t", index=False).encode()).hexdigest()
        checkpoint_path = CHECKPOINT_DIR / (hashlib.sha256(pair_key.encode()).hexdigest()[:20] + ".json")
        cached = checkpoint_records.get(pair_key)
        cache_is_current = bool(cached and cached.get("mapping_version") == MAPPING_VERSION)
        if cached and cached.get("fingerprint") == fingerprint and len(cached.get("rows", [])) == len(rows) and cache_is_current:
            saved_structure = str(cached.get("edge", {}).get("structure_file", ""))
            if not saved_structure or (ROOT / saved_structure).is_file():
                row_maps.extend(cached["rows"])
                edge_rows.append(cached["edge"])
                continue
        row_start = len(row_maps)
        target_seq, partner_seq = acc_sequence(target_acc, seqs), acc_sequence(partner_acc, seqs)
        if not target_seq or not partner_seq:
            for r in rows.itertuples(index=False):
                row_maps.append({"gene": r.gene, "VariationID": r.VariationID, "variant": r.canonical_variant,
                                 "partner": r.partner, "target_accession": target_acc,
                                 "partner_accession": partner_acc, "status": "SEQUENCE_MISMATCH",
                                 "failure_reason": "SEQUENCE_MISMATCH"})
            edge_rows.append({"target_accession": target_acc, "partner_accession": partner_acc,
                              "status": "NO_LOCAL_REFERENCE_SEQUENCE", "scoreable_unique_variants": 0,
                              "unique_variants": rows.canonical_variant.nunique()})
            write_pair_checkpoint(checkpoint_path, pair_key, fingerprint,
                                  row_maps[row_start:], edge_rows[-1])
            continue

        # Frozen candidate order is preserved. If every row agrees on its candidate list,
        # this is simply the row-level list; otherwise use the ordered intersection.
        row_candidate_lists = [candidate_ids(x) for x in rows.experimental_PDB_IDs.fillna("").unique()]
        candidates = row_candidate_lists[0] if row_candidate_lists else []
        if any(lst != candidates for lst in row_candidate_lists[1:]):
            common = set(candidates)
            for lst in row_candidate_lists[1:]: common.intersection_update(lst)
            candidates = [x for x in candidates if x in common]
        unique_variants = rows[["canonical_variant", "canonical_residue"]].drop_duplicates()
        same_protein = target_acc.split("-")[0] == partner_acc.split("-")[0]
        candidate_evaluations = []
        edge_reason = "STRUCTURE_CHAIN_FAILURE"
        for candidate_i, (pid, assembly) in enumerate(candidates):
            cif_path = CIF_DIR / f"{pid}-{assembly}.cif"
            if not cif_path.is_file():
                edge_reason = "STRUCTURE_CHAIN_FAILURE"
                continue
            try:
                structure, entity_seqs, asym_entity, chains, geometries = parse_cif(cif_path)
            except Exception:
                edge_reason = "STRUCTURE_CHAIN_FAILURE"
                continue
            target_maps, partner_maps = eligible_entities(entity_seqs, target_seq, partner_seq, same_protein)
            target_chains = [cid for cid, eid in asym_entity.items() if eid in target_maps and cid in chains]
            partner_chains = [cid for cid, eid in asym_entity.items() if eid in partner_maps and cid in chains]
            if not target_chains or not partner_chains:
                edge_reason = "STRUCTURE_CHAIN_FAILURE"
                continue
            for tc in target_chains:
                target_entity = asym_entity[tc]
                variant_map = variant_coordinate_map(unique_variants, target_maps, chains[tc], target_entity)
                for pc in partner_chains:
                    if tc == pc: continue
                    distance = min_chain_distance(geometries[tc], geometries[pc])
                    if not np.isfinite(distance) or distance > EDGE_CUTOFF_A:
                        continue
                    covered = sum(x["status"] == "READY" for x in variant_map.values())
                    candidate_evaluations.append({
                        "candidate_index": candidate_i, "pdb_id": pid, "assembly": assembly,
                        "cif_path": cif_path, "structure": structure, "target_chain": tc,
                        "partner_chain": pc, "target_entity": target_entity,
                        "partner_entity": asym_entity[pc], "target_maps": target_maps,
                        "partner_maps": partner_maps, "variant_map": variant_map,
                        "chains": chains,
                        "covered_variants": covered, "unique_variants": len(unique_variants),
                        "min_interchain_distance_A": distance,
                    })
            if candidate_evaluations:
                edge_reason = "STRUCTURE_RESIDUE_FAILURE"

        if not candidate_evaluations:
            for r in rows.itertuples(index=False):
                row_maps.append({"gene": r.gene, "VariationID": r.VariationID, "variant": r.canonical_variant,
                                 "partner": r.partner, "target_accession": target_acc,
                                 "partner_accession": partner_acc, "status": edge_reason,
                                 "failure_reason": edge_reason, "structure_candidates": r.experimental_PDB_IDs})
            edge_rows.append({"target_accession": target_acc, "partner_accession": partner_acc,
                              "status": edge_reason, "scoreable_unique_variants": 0,
                              "unique_variants": len(unique_variants), "candidate_n": len(candidates)})
            write_pair_checkpoint(checkpoint_path, pair_key, fingerprint,
                                  row_maps[row_start:], edge_rows[-1])
            continue

        # Pick the structure/chain pair that covers the most unique target variants;
        # ties follow frozen candidate order, then shortest physical chain distance,
        # then stable chain IDs. No disease labels or scores enter this choice.
        chosen = sorted(candidate_evaluations, key=lambda x: (
            -x["covered_variants"], x["candidate_index"],
            x["min_interchain_distance_A"], x["target_chain"], x["partner_chain"],
        ))[0]
        token = hashlib.sha256(
            f"{target_acc}|{partner_acc}|{chosen['pdb_id']}|{chosen['assembly']}|"
            f"{chosen['target_chain']}|{chosen['partner_chain']}".encode()
        ).hexdigest()[:12]
        target_id, partner_id = "T" + token, "P" + token
        cif_out = STRUCT_DIR / f"{target_id}_{partner_id}_{chosen['pdb_id']}-{chosen['assembly']}.cif"
        save_two_chain_cif(chosen["chains"], chosen["target_chain"], chosen["partner_chain"], cif_out)
        variant_map = chosen["variant_map"]
        for r in rows.itertuples(index=False):
            info = variant_map.get(str(r.canonical_variant), {"status": "STRUCTURE_RESIDUE_FAILURE"})
            status = info["status"]
            row_maps.append({
                "gene": r.gene, "VariationID": r.VariationID, "variant": r.canonical_variant,
                "partner": r.partner, "target_accession": target_acc,
                "partner_accession": partner_acc,
                "target_id": target_id, "partner_id": partner_id,
                "complex_id": f"{target_id}_{partner_id}",
                "selected_structure": f"{chosen['pdb_id']}:{chosen['assembly']}",
                "target_chain_in_source": chosen["target_chain"],
                "partner_chain_in_source": chosen["partner_chain"],
                "local_variant": info.get("local_variant", "NA"),
                "local_index_zero_based": info.get("local_index", "NA"),
                "observed_reference_aa": info.get("observed_aa", "NA"),
                "minimum_interchain_distance_A": chosen["min_interchain_distance_A"],
                "status": status,
                "failure_reason": "NA" if status == "READY" else status,
                "structure_candidates": r.experimental_PDB_IDs,
                "tier": r.tier,
            })
        edge_rows.append({
            "target_accession": target_acc, "partner_accession": partner_acc,
            "target_id": target_id, "partner_id": partner_id,
            "selected_structure": f"{chosen['pdb_id']}:{chosen['assembly']}",
            "target_chain_in_source": chosen["target_chain"],
            "partner_chain_in_source": chosen["partner_chain"],
            "min_interchain_distance_A": chosen["min_interchain_distance_A"],
            "unique_variants": chosen["unique_variants"],
            "scoreable_unique_variants": chosen["covered_variants"],
            "variant_residue_coverage_pct": round(100*chosen["covered_variants"]/max(1,chosen["unique_variants"]), 3),
            "candidate_n": len(candidates),
            "evaluated_chain_pairs": len(candidate_evaluations),
            "status": "READY" if chosen["covered_variants"] else "STRUCTURE_RESIDUE_FAILURE",
            "structure_file": str(cif_out.relative_to(ROOT)),
            "selection_rule": "max unique variant residues with coordinates; tie: frozen candidate order, min chain distance, chain IDs",
        })
        write_pair_checkpoint(checkpoint_path, pair_key, fingerprint,
                              row_maps[row_start:], edge_rows[-1])

    rows_df = pd.DataFrame(row_maps)
    edges_df = pd.DataFrame(edge_rows)
    rows_df.to_csv(OUT / "row_structure_mapping.tsv", sep="\t", index=False, na_rep="NA")
    edges_df.to_csv(OUT / "structure_pair_manifest.tsv", sep="\t", index=False, na_rep="NA")
    ready = rows_df[rows_df.status.eq("READY")].copy()
    variants_out = ready[["target_id", "local_variant", "partner_id"]].drop_duplicates()
    variants_out.to_csv(MODEL_INPUT / "variants_ready.tsv", sep="\t", index=False, header=False)

    print(f"frozen_rows={len(universe)} mapped_rows={len(rows_df)} ready_rows={len(ready)}")
    print(f"unique_accession_pairs={len(unique_pairs)} prepared_pair_structures={edges_df.status.eq('READY').sum()}")
    print("row_status_counts=" + repr(rows_df.status.value_counts().to_dict()))
    print(f"structure_pair_manifest={OUT/'structure_pair_manifest.tsv'}")


if __name__ == "__main__":
    main()
