"""Finalize only after all recovery and pilot evidence passes; no benchmark."""
from pathlib import Path
import json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';L=R/'logs/phase2SR'
def read(n):return pd.read_csv(P/n,sep='\t',keep_default_na=False)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
freeze=read('input_freeze_manifest.tsv')
for r in freeze.itertuples(index=False):assert digest(R/r.path)==r.SHA256,r.path
seq=read('canonical_sequence_recovery.tsv');assert len(seq)==7300 and seq.status.isin(['EXACT','OBSOLETE_REPLACED']).mean()>=.95
legacy=read('legacy_denominator_reproduction.tsv');assert (legacy.ge3_pct>=70).sum()==9
assert sorted(legacy.loc[legacy.balanced,'gene'])==['FGFR2','MEN1','TUBA1A']
assert 'PASS:' in (L/'legacy_independent_audit.log').read_text()
assert 'PILOT_VALIDATION_PASS' in (L/'pilot_final_audit.log').read_text()
for m in ['sppider','xpair']:
    assert read(m+'_example_validation.tsv').result.eq('PASS').all()
    assert 'environment_import_pass": true' in (L/('environment_'+m+'_validation.log')).read_text()
    mf=read(m+'_encoder_manifest.tsv');mf['offline_load']='PASS';mf['example']='PASS';mf['status']='OFFLINE_VALIDATED';mf.to_csv(P/(m+'_encoder_manifest.tsv'),sep='\t',index=False)
common=read('common_primary_membership.tsv');n=int(common.included_in_common_primary_denominator.sum());assert n==9978
pilot=[]
for m in ['sppider','xpair']:
    d=read('pilot_predictions_'+m+'.tsv');pilot.append({'method':m,'selected_pairs':len(d),'scored_pairs':int(d.status.eq('SCORED').sum()),'statuses':json.dumps(d.status.value_counts().to_dict())})
report=R/'reports/PHASE2SR_OFFLINE_ENCODER_AND_DENOMINATOR_RECOVERY.md'
old=report.read_text();legacy_section=old.split('## Exact legacy reproduction\n\n',1)[1].split('## Canonical sequence recovery',1)[0]
s='# Phase 2S-R offline encoder and denominator recovery\n\n## Gate\n\n`GO_RESUME_PHASE2S`\n\nBoth exact encoders loaded offline and both official example workflows passed. Canonical sequence recovery, exact legacy denominator reproduction and the disease-blind technical pilot passed. Recovery ends here. The full Phase 2S benchmark has not started.\n\n'
s+='## Preserved state\n\nAll '+str(len(freeze))+' frozen files were rehashed and remain unchanged, including original results and neural heads. U1=42, U2=132 and 15,942 frozen gene–partner pairs remain fixed.\n\n## Exact legacy reproduction\n\n'+legacy_section
s+='## Canonical sequences\n\n'+f'{len(seq):,} requested; {int(seq.status.isin(["EXACT","OBSOLETE_REPLACED"]).sum()):,} recovered. '+json.dumps(seq.status.value_counts().to_dict())+'. All 267 prior sequences retain identical hashes.\n\n'+read('sequence_recovery_coverage.tsv').to_markdown(index=False)+'\n\n'
s+='UniProtKB release 2026_03 (2026-09-02) and authoritative batch responses are preserved. Partners without frozen accessions remain unavailable; no isoform or accession identity was substituted.\n\n## Exact models and runtime\n\n'+pd.DataFrame(json.loads((P/'encoder_download_provenance.json').read_text())['encoders']).to_markdown(index=False)+'\n\n'
s+='The user authorized local download after switching networks on 2026-10-02. Official revision, byte sizes, SHA256 and tokenizer/config hashes are recorded in the manifests.\n\nIndependent Python 3.12.13 environments reside at `/Users/zy/.codex/tools/phase2s_sppider` and `/Users/zy/.codex/tools/phase2s_xpair`. Both imports and package consistency passed. Package pins are recorded in environment_sppider.txt and environment_xpair.txt. exFAT package copying failed; a workspace sparse-image mount was denied, so reusable user-level SSD environments were installed under the existing environment authorization. Task outputs remain in the workspace.\n\n'
s+='SPPIDER v2.0-manuscript uses unchanged ESM2, heads, MAX_TOKENS=1024 and STRIDE=512. No installed SPPIDER source patch was required. X-PAIR 0.1.0 uses the bundled checkpoint, Ankh T5EncoderModel/AutoTokenizer with dimension 1536 and the frozen 50–2000 length rule. Only its encoder path and local_files_only were patched; the exact diff is preserved. No unsupported protein was truncated.\n\n'
s+='## Example and pilot validation\n\nSPPIDER used its official TP53 three-partner FASTAs and validated all six receptor/peptide tracks. X-PAIR ran the documented official prediction workflow on those same three pairs; no executable example FASTA was bundled in its documentation. Both output lengths matched WT sequences, all probabilities were finite and within [0,1].\n\n'
s+='The 20 genes were fixed before inference: ten deterministic U1 genes and ten additional U2 genes. The first three lexicographic frozen partners per gene were selected for technical checks only. ABCA4 has no frozen partners and contributes no pilot pairs; it remains selected and unavailable. Total pilot pairs=57. No full-universe coverage or rescue conclusion is inferred from this subset.\n\n'+pd.DataFrame(pilot).to_markdown(index=False)+'\n\n'
s+='SPPIDER retains complete float32 receptor and peptide residue arrays in keyed npz caches. Its original window merge is used; short and long sequences, both endpoints and internal coverage were checked without padding or truncation. X-PAIR raw six-decimal tracks are archived before the unchanged official CLI rounds its final output to three decimals and removes its transient file. This wrapper preserves output only; inference and source remain unchanged. A second pass reuses the same embeddings and must reproduce every exported probability exactly. No full-precision tensor equivalence beyond the official six-decimal export is claimed.\n\n'
s+='Variant positions were extracted from cached WT tracks only when the frozen canonical mapping and reference residue matched. No mutant embeddings were generated. The ≥3 calculation counts distinct frozen partners with available scores and remains separate from pair-space availability. No coordinate-offset or sequence-mismatch failure remained. Evidence: sppider_window_validation.tsv, pilot_variant_scores_*.tsv, pilot_variant_partner_count_validation.tsv and inference_cache_manifest.tsv.\n\n'
s+='## Stop\n\nOriginal NO_GO_EDGETIC_REPRESENTATION was not altered.\n\nNo sequence-model rescue verdict was issued.\n\nDisease-label phenotype association testing was NOT started.\n\nWait for authorization before the full Phase 2S benchmark.\n'
report.write_text(s)
(P/'recovery_delivery_audit.json').write_text(json.dumps({'frozen_files_checked':len(freeze),'unchanged':True,'requested_accessions':len(seq),'common_primary':n,'pilot':pilot,'verdict':'GO_RESUME_PHASE2S','full_benchmark_started':False},indent=2)+'\n')
print('GO_RESUME_PHASE2S',flush=True)
