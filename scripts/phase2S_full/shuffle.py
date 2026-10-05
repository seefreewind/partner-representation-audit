from common import *
import numpy as np,subprocess,os,importlib.util
from scipy.stats import spearmanr,pearsonr
plan=rd('results/phase2S_full/shuffled_control_plan.tsv');C=O/'shuffle_cache';C.mkdir(exist_ok=True);rows=[]
spec=importlib.util.spec_from_file_location('sppider',R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/cli/sppider_seq.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
args=s.build_parser().parse_args(['--query','unused','--partners','unused','--output-dir',str(C),'--esm-model-name',str(R/'external_models/esm2_t33_650M_UR50D'),'--receptor-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_rec_run27_best.pt'),'--peptide-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_pep_run07_best.pt'),'--device','cpu'])
s.setup_runtime(args);s.load_ppi_models()
pairs=rd('results/phase2S_full/pair_status_sppider.tsv').set_index(['gene','partner'])
def compare(r,a,b,precision):
    assert len(a)==len(b) and np.isfinite(b).all()
    return dict(method=r.method,gene=r.gene,partner=r.partner,shuffled_partner_accession=r.shuffled_partner_accession,true_partner_length=r.true_partner_length,shuffled_length=r.shuffled_length,residue_count=len(a),Pearson=pearsonr(a,b).statistic if np.ptp(a)>0 and np.ptp(b)>0 else np.nan,Spearman=spearmanr(a,b).statistic if np.ptp(a)>0 and np.ptp(b)>0 else np.nan,mean_absolute_track_change=np.mean(np.abs(a-b)),RMS_track_change=np.sqrt(np.mean((a-b)**2)),status='EVALUATED',precision=precision,no_phenotype_labels=True)
for r in plan.loc[plan.method.eq('sppider')].itertuples(index=False):
    info=pairs.loc[(r.gene,r.partner)]
    if info.status!='SCORED':rows.append(dict(method='sppider',gene=r.gene,partner=r.partner,status='GENUINE_PAIR_UNAVAILABLE'));continue
    key=hashlib.sha256(f'{r.query_accession}|{r.shuffled_partner_accession}|{info.encoder_sha256}'.encode()).hexdigest();file=C/(key+'.npy')
    if file.exists():b=np.load(file)
    else:
        b=s.predict_query_given_context(seq(r.query_accession),seq(r.shuffled_partner_accession),s.receptor_site_model,cache_query=True,cache_context=True);np.save(file,b)
    a=np.load(R/info.score_array_path)['receptor'];rows.append(compare(r,a,b,'float32'));wr(pd.DataFrame(rows),'shuffled_partner_control.tsv');s.embedding_cache.clear()
d=plan.loc[plan.method.eq('xpair')];D=C/'xpair';D.mkdir(exist_ok=True)
import gc
s.embedding_cache.clear();s.esm_model=None;s.receptor_site_model=None;s.peptide_site_model=None;gc.collect()
pd.DataFrame(dict(id1=d.query_accession,id2=d.shuffled_partner_accession)).drop_duplicates().to_csv(D/'pairs.tsv',sep='\t',index=False)
ids=sorted(set(d.query_accession)|set(d.shuffled_partner_accession));(D/'sequences.fasta').write_text(''.join(f'>{a}\n{seq(a)}\n' for a in ids))
raw=D/'predictions_six_decimal.csv'
if not raw.exists():
    env=dict(os.environ,HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD='1',XPAIR_ANKH_LOCAL_DIR=str(R/'external_models/ankh-large'),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4')
    with (L/'xpair_shuffle.log').open('a') as f:subprocess.run(['/Users/zy/.codex/tools/phase2s_xpair/bin/python','scripts/phase2SR/xpair_preserve_output.py','--predict_pairs_path',str(D/'pairs.tsv'),'--predict_fasta_path',str(D/'sequences.fasta'),'--predict_output_path',str(D/'predictions.csv'),'--emb_dir_prediction',str(O/'embedding_cache/xpair'),'--num_workers','0','--batch_size_prediction','1','--interaction_threshold','0'],env=env,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
lookup={(r.id1,r.id2):np.fromstring(r.interface1_probs,sep=' ') for r in pd.read_csv(raw).itertuples(index=False)};pairs=rd('results/phase2S_full/pair_status_xpair.tsv').set_index(['gene','partner'])
for r in d.itertuples(index=False):
    info=pairs.loc[(r.gene,r.partner)]
    if info.status!='SCORED':rows.append(dict(method='xpair',gene=r.gene,partner=r.partner,status='GENUINE_PAIR_UNAVAILABLE'));continue
    a=np.load(R/info.score_array_path)['query'];b=lookup[(r.query_accession,r.shuffled_partner_accession)];rows.append(compare(r,a,b,'official_six_decimal_export'))
wr(pd.DataFrame(rows),'shuffled_partner_control.tsv');print('SHUFFLE_CONTROL_COMPLETE',flush=True)
