from pathlib import Path
import hashlib,json,subprocess
import pandas as pd
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR'
def wr(d,n):d.to_csv(P/n,sep='\t',index=False)
def rd(n):return pd.read_csv(P/n,sep='\t',keep_default_na=False)
def main(verdict):
    if (R/'external_models/esm2_t33_650M_UR50D/OFFICIAL_FILE_MANIFEST.json').exists():
        print('Verified encoder acquisition has started; old absence-only finalizer is superseded. Preserve manifests and report for updated ingestion.')
        return
    for method,repo,directory in [('sppider','facebook/esm2_t33_650M_UR50D','esm2_t33_650M_UR50D'),('xpair','ElnaggarLab/ankh-large','ankh-large')]:
        wr(pd.DataFrame([{'repo_id':repo,'local_directory':str(R/'external_models'/directory),'revision':'NA','weight_filename':'NA','byte_size':'NA','SHA256':'NA','config_SHA256':'NA','tokenizer_SHA256':'NA','status':'AWAITING_EXTERNAL_TRANSFER'}]),method+'_encoder_manifest.tsv')
        wr(pd.DataFrame([{'test':x,'result':'FAIL','execution_status':'NOT_RUN_ENCODER_ABSENT','reason':'Official encoder files must be transferred from another machine'} for x in ['offline_load','official_example','residue_length','finite_probabilities']+(['short_window','long_window','receptor_track','peptide_track'] if method=='sppider' else [])]),method+'_example_validation.tsv')
        env=R/('phase2s_'+method)
        proc=subprocess.run(['/Users/zy/.codex/tools/bin/uv','pip','freeze','--python',str(env/'bin/python')],capture_output=True,text=True)
        (P/('environment_'+method+'.txt')).write_text('ENVIRONMENT='+str(env)+'\nPYTHON=3.12.13\n'+proc.stdout+'\n'+proc.stderr)
    for n in ['pilot_predictions_sppider.tsv','pilot_predictions_xpair.tsv']:
        wr(pd.DataFrame(columns=['gene','partner','residue','probability','status']),n)
    wr(pd.DataFrame(columns=['method_version','encoder_model_SHA256','query_sequence_SHA256','partner_sequence_SHA256','cache_path','status']),'inference_cache_manifest.tsv')
    s=rd('canonical_sequence_recovery.tsv');valid=s.status.isin(['EXACT','OBSOLETE_REPLACED']);ok=set(s.loc[valid,'requested_accession'])
    partners=pd.read_csv(R/'results/phase2S/frozen_partner_sets.tsv',sep='\t',keep_default_na=False)
    u=pd.read_csv(R/'results/phase2S/frozen_gene_universes.tsv',sep='\t');u1=set(u.loc[u.U1_STRUCTURAL_42,'gene'])
    seqcov=[]
    for name,pp in [('U1',partners.loc[partners.gene.isin(u1)]),('U2',partners)]:
        required=set(pp.query_accession)|set(pp.partner_accession);required.discard('')
        complete=pp.query_accession.isin(ok)&pp.partner_accession.isin(ok)
        seqcov.append({'universe':name,'distinct_accessions':len(required),'recovered_accessions':len(required&ok),'accession_pct':100*len(required&ok)/len(required),'frozen_PPI_pairs':len(pp),'complete_pairs':int(complete.sum()),'pair_completeness_pct':100*complete.mean()})
    wr(pd.DataFrame(seqcov),'sequence_recovery_coverage.tsv')
    exact=rd('legacy_denominator_reproduction.tsv');names=exact.loc[exact.balanced.eq(True),'gene'].tolist();assert sorted(names)==['FGFR2','MEN1','TUBA1A']
    assert 'PASS:' in (R/'logs/phase2SR/legacy_independent_audit.log').read_text()
    common=rd('common_primary_membership.tsv');count=int(common.included_in_common_primary_denominator.eq(True).sum())
    freeze=rd('input_freeze_manifest.tsv');original='reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md'
    assert hashlib.sha256((R/original).read_bytes()).hexdigest()==freeze.loc[freeze.path.eq(original),'SHA256'].iloc[0]
    text='# Phase 2S-R offline encoder and denominator recovery\n\n## Gate\n\n`'+verdict+'`\n\nBoth official encoders await external transfer. No Hugging Face acquisition was attempted in this recovery phase on this computer. Offline loading and examples have not been executed; FAIL denotes a gate not passed, not an observed model inference defect. No pilot or full representation benchmark ran.\n\n'
    text+='## Preserved state\n\nU1=42, U2=132 and 15,942 frozen gene–partner pairs are unchanged. Current Phase 2S files, heads, checkpoint, original results and manifests were SHA256-frozen before recovery in `metadata/PHASE2SR_INPUT_FREEZE_MANIFEST.tsv`. Original outputs were not rewritten.\n\n'
    text+='## Exact legacy reproduction\n\nAn independent recomputation used original scores only and traced `scripts/phase1b0r3_finalize.py::make_membership`. It reproduced every gene primary denominator, ≥3 partner count, rounded coverage, balance class and final gate. A separate audit reproduced all 109 disease-group denominators, ≥3 counts, rounded percentages and multilabel exclusions. Results: **9 genes ≥70%; 3 balanced: FGFR2, MEN1, TUBA1A. LEGACY_DENOMINATOR_EXACT_PASS=YES.**\n\n'
    text+='The legacy rule first constructs disease membership from Phase 1A `primary_cohort_eligible=YES`. Its gene-level multilabel exclusion is then restricted to variants appearing in original score rows. Consequently 674 multilabel variants absent from score rows remain in legacy primary denominators (1,555 disease-membership rows). This exactly explains why an ordinary full-cohort single-label registry disagreed. We preserve this behavior explicitly rather than retrospectively fix the original result. Group balance uses original strict counts and excludes multilabel variants for the retained-group ≥10 gate, as in the original script.\n\n'
    text+=f'COMMON_PRIMARY_MEMBERSHIP contains **{count:,} included unique gene–variant units**. U1 retains exact immutable legacy membership. U2 outside U1 uses frozen mapping eligibility and one strict label, fixed before inference. Both models will evaluate the same frozen identities; no future method-specific membership or structure filter is permitted. This U2 extension is explicit and does not claim to have historical legacy scores. Pair-space denominators remain separate from variant ≥3 representation.\n\n'
    text+='## Canonical sequence recovery\n\n'
    text+=f'Requested {len(s):,}; recovered {int(valid.sum()):,} ({100*valid.mean():.4f}%). Status counts: '+json.dumps(s.status.value_counts().to_dict())+'.\n\n'
    text+=pd.DataFrame(seqcov).to_markdown(index=False)+'\n\n'
    text+='UniProtKB batch searches record complete authoritative JSON, including primary/secondary accession relationships and canonical sequences, under `external_models/canonical_sequences/`. Secondary-to-primary replacement is accepted only when explicitly present in UniProt metadata. Frozen cached sequences are SHA256 checked. NOT_FOUND records with access-pending reasons are technical retrieval failures, not established biological absence.\n\n'
    text+='## Encoder and environment contracts\n\nSPPIDER v2.0-manuscript requires facebook/esm2_t33_650M_UR50D (33 layers, hidden dimension 1280, maximum positions 1026). MAX_TOKENS=1024 and STRIDE=512 remain unchanged. The CLI already supports --esm-model-name; no code patch was needed yet. X-PAIR 0.1.0 requires official ElnaggarLab/ankh-large, T5EncoderModel plus AutoTokenizer, hidden dimension 1536; its unchanged source contract is in `xpair_encoder_contract.md`. Checkpoint length limits remain 50–2000; no truncation. No unverified path patch has been applied.\n\n'
    text+='UniProt release headers from an authoritative probe are preserved in `logs/phase2SR/uniprot_release_headers.txt` (2026_03, released 2026-09-02).\n\nSeparate Python 3.12.13 environments are `phase2s_sppider` and `phase2s_xpair`; actual installed packages are saved in environment text files. Install diagnostics are in logs/phase2SR. Dependencies require import and consistency checks before example execution. SPPIDER has no esm package requirement; X-PAIR uses its own declared dependencies.\n\n'
    text+='## External download and transfer instructions\n\nRun the following **on another machine with working official Hugging Face access**, after copying `scripts/phase2SR/download_encoders_external.py` there. Do not execute encoder download on this project Mac.\n\n```bash\npython3 -m venv encoder-transfer-venv\nsource encoder-transfer-venv/bin/activate\npython -m pip install huggingface_hub==1.0.9\npython download_encoders_external.py --output ./external_models\n```\n\n'
    text+='The script is available at ['+str(R/'scripts/phase2SR/download_encoders_external.py')+']('+str(R/'scripts/phase2SR/download_encoders_external.py')+'). It resolves an official commit once per model, records it in OFFICIAL_REVISION.json and pins all subsequent download/resume calls to that commit. It downloads one official weight format plus config/tokenizer support, then records every file size and actual SHA256 in OFFICIAL_FILE_MANIFEST.json, comparing official LFS hashes where provided. It checks encoder architecture in config.json. Reuse the same directories when resuming. Do not replace weights with third-party conversions or quantized files.\n\nTransfer the **entire two directories**, including both provenance JSON files, to:\n\n```text\n'+str(R/'external_models/esm2_t33_650M_UR50D')+'\n'+str(R/'external_models/ankh-large')+'\n```\n\n'
    text+='Environment installation was launched in two separate processes. Until `environment_sppider_validation.log` and `environment_xpair_validation.log` explicitly report environment_import_pass=true, environment completion is pending. This does not establish offline model loading or example PASS.\n\n'
    text+='After transfer, rehash locally, verify config architecture and tokenizer files, inspect the exact revision, apply only a verified local-path patch, and run examples with HF_HUB_OFFLINE=1 and TRANSFORMERS_OFFLINE=1. Validate both SPPIDER tracks and short/long windows, plus X-PAIR lengths and finite residue tracks. Pilot genes are deterministic first 10 U1 and first 10 U2 outside U1, stored before inference. Pilot outputs and cache manifest remain headers only because gates have not passed.\n\n'
    text+='GO_RESUME_PHASE2S requires all recovery gates and pilot checks. At that point stop and request authorization before full Phase 2S.\n\nOriginal NO_GO_EDGETIC_REPRESENTATION was not altered.\n\nNo sequence-model rescue verdict was issued.\n\nDisease-label phenotype association testing was NOT started.\n'
    (R/'reports/PHASE2SR_OFFLINE_ENCODER_AND_DENOMINATOR_RECOVERY.md').write_text(text)
    print(json.dumps({'requested':len(s),'recovered':int(valid.sum()),'pct':100*valid.mean(),'unresolved':int((~valid).sum()),'common_primary':count,'verdict':verdict}))
if __name__=='__main__':
    import sys;main(sys.argv[1])
