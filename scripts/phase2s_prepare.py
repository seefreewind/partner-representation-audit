#!/usr/bin/env python3
"""Freeze Phase 2S inputs and label-blind protein-pair inference requests."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/phase2S'
OUT.mkdir(parents=True, exist_ok=True)

def read(name):
    return pd.read_csv(ROOT/name, sep='\t', dtype=str, keep_default_na=False)

def write(frame, name):
    frame.to_csv(OUT/name, sep='\t', index=False, lineterminator='\n', na_rep='NA')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()

def main():
    sources = [
        'results/phase0R/gene_disease_variant_counts_exact_rcv.tsv',
        'results/phase0R/mechanism_gene_states_exact_strict.tsv',
        'results/phase1A/protein_mapping.tsv', 'results/phase1A/ppi_provenance.tsv',
        'results/phase1A/interface_mapping.tsv',
        'results/phase1A/partner_specific_matrix_structural.tsv',
        'results/phase1A/coverage_by_gene.tsv',
        'results/phase1A/mutpred_training_overlap.tsv',
        'results/phase1B0/structural_42_genes.tsv',
        'results/phase1B0/analysis_specification.md',
        'results/phase1B0/mutpred_training_overlap.tsv',
        'results/phase1B0R3/mutpred_ppi_v1_scores.tsv',
        'results/phase1B0R3/coverage_by_gene.tsv',
        'results/phase1B0R3/differential_score_coverage.tsv',
        'results/phase1B0R3/frozen_input_hash_validation.tsv',
        'reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md',
        'metadata/PHASE1B0_INPUT_FREEZE_MANIFEST.tsv',
        'data/raw/hgnc_non_alt_loci_2026-09-27.tsv',
    ]
    audit=read('metadata/PHASE1B0_INPUT_FREEZE_MANIFEST.tsv')
    for r in audit.itertuples(index=False):
        assert sha(ROOT/r.path)==r.SHA256, 'INPUT_FREEZE_VIOLATION: '+r.path
    original=ROOT/'reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md'
    assert '`NO_GO_EDGETIC_REPRESENTATION`' in original.read_text()
    c=read('results/phase1A/coverage_by_gene.tsv')
    u2=set(c.loc[~c.mechanism_stratum.isin(['KNOWN_INCOMPATIBLE','CONFIRMED_INCOMPATIBLE']),'gene'])
    u1=set(read('results/phase1B0/structural_42_genes.tsv').gene)
    assert len(c)==137 and len(u2)==132 and len(u1)==42 and u1<=u2
    universes=pd.DataFrame({'gene':sorted(u2)})
    universes['U1_STRUCTURAL_42']=universes.gene.isin(u1)
    universes['U2_STRICT_132']=True
    write(universes,'frozen_gene_universes.tsv')
    groups=read('results/phase0R/gene_disease_variant_counts_exact_rcv.tsv')
    groups=groups.loc[groups.strict_disease_ge10.eq('YES') & groups.gene.isin(u2)]
    write(groups,'frozen_disease_groups.tsv')
    strict=groups.groupby('gene').MONDO.agg(set).to_dict()
    variants=read('results/phase1A/protein_mapping.tsv')
    variants=variants.loc[variants.gene.isin(u2)].copy()
    variants['strict_MONDO_groups']=[ ';'.join(sorted(set(m.split(';'))&strict.get(g,set()))) for g,m in zip(variants.gene,variants.MONDOs)]
    variants=variants.loc[variants.strict_MONDO_groups.ne('')].copy()
    assert not variants[['gene','VariationID']].duplicated().any()
    variants['strict_group_n']=variants.strict_MONDO_groups.str.count(';')+1
    variants['multilabel']=variants.strict_group_n.gt(1)
    variants['frozen_primary_single_label']=variants.strict_group_n.eq(1)
    write(variants,'frozen_variants.tsv')
    write(variants.loc[variants.multilabel], 'multilabel_registry.tsv')
    pp=read('results/phase1A/ppi_provenance.tsv')
    pp=pp.loc[pp.gene.isin(u2) & pp.included_in_primary_cohort.eq('YES')].copy()
    assert not pp[['gene','partner']].duplicated().any()
    pp['evidence_stratum']=['EXPERIMENTAL' if e.lower()=='yes' else 'PIONEER_SUPPORTED' if p.lower()=='yes' else 'PREDICTION_ONLY_CREDIBLE' for e,p in zip(pp.experimental_support,pp.PIONEER_coverage)]
    hgnc=read('data/raw/hgnc_non_alt_loci_2026-09-27.tsv')
    acc={}
    for r in hgnc.itertuples(index=False):
        ids=[x for x in re.split(r'[|;,\s]+',r.uniprot_ids.strip('"')) if x]
        if len(ids)==1: acc[r.symbol]=ids[0]
    for gene,g in variants.groupby('gene'):
        ids=set(g.canonical_UniProt_accession)-{'','NA','nan'}
        if len(ids)==1: acc[gene]=next(iter(ids))
    pp['query_accession']=pp.gene.map(acc).fillna('')
    pp['partner_accession']=pp.partner.map(acc).fillna('')
    pp['accession_mapping_status']= ['CANONICAL_SINGLE_ACCESSION' if a and b and '-' not in a and '-' not in b else 'UNRESOLVED_CANONICAL_ACCESSION' for a,b in zip(pp.query_accession,pp.partner_accession)]
    write(pp,'frozen_partner_sets.tsv')
    requests=pp[['gene','partner','query_accession','partner_accession','accession_mapping_status','tier','evidence_stratum']].copy()
    write(requests,'protein_pair_requests.tsv')
    needed=set(pp.query_accession)|set(pp.partner_accession)
    needed.discard('')
    seqrows=[]
    seqdir=ROOT/'data/interim/phase1A_uniprot_sequences'
    cached_names={entry.name for entry in seqdir.iterdir()}
    for a in sorted(needed):
        path=seqdir/f'{a}.fasta'
        sequence=''.join(line.strip() for line in path.read_text().splitlines() if not line.startswith('>')) if path.name in cached_names else ''
        if sequence: sources.append(str(path.relative_to(ROOT)))
        seqrows.append({'UniProt_accession':a,'length':len(sequence) if sequence else 'NA','sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest() if sequence else 'NA','source':str(path.relative_to(ROOT)) if sequence else 'UNIPROT_CANONICAL_RETRIEVAL_PENDING','isoform_status':'CANONICAL','status':'AVAILABLE_FROZEN' if sequence else 'SEQUENCE_RETRIEVAL_PENDING'})
    write(pd.DataFrame(seqrows),'sequence_manifest.tsv')
    manifest=[]
    for n in sorted(set(sources)):
        path=ROOT/n
        manifest.append({'path':n,'size_bytes':path.stat().st_size,'sha256':sha(path),'status':'FROZEN'})
    mf=pd.DataFrame(manifest)
    write(mf,'input_freeze_manifest.tsv')
    mf.to_csv(ROOT/'metadata/PHASE2S_INPUT_FREEZE_MANIFEST.tsv',sep='\t',index=False)
    old=read('results/phase1B0R3/coverage_by_gene.tsv')
    assert len(old)==42 and old.eligible_for_rebuilt_primary_cohort.eq('YES').sum()==3
    write(old,'mutpred_frozen_coverage_by_gene.tsv')
    # Sparse long coordinates preserve per-gene frozen partners without a
    # cross-gene dense matrix containing meaningless gene/partner cells.
    a=read('results/phase1B0R3/mutpred_ppi_v1_scores.tsv')
    a=a.loc[a.gene.isin(u2)].copy()
    write(a[['gene','VariationID','partner','score','score_status','failure_reason']], 'representation_matrix_mutpred.tsv')
    print(json.dumps({'U1':len(u1),'U2':len(u2),'frozen_variants':len(variants),'single_label_primary':int(variants.frozen_primary_single_label.sum()),'frozen_gene_partner_pairs':len(pp),'available_sequences':sum(r['status']=='AVAILABLE_FROZEN' for r in seqrows),'required_sequences':len(seqrows),'original_verdict':'NO_GO_EDGETIC_REPRESENTATION'},ensure_ascii=False),flush=True)

if __name__=='__main__': main()
