"""Serial, locked offline recovery pilot. Stops before full Phase 2S."""
from pathlib import Path
import os,subprocess,json,fcntl,datetime,traceback
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';L=R/'logs/phase2SR'
lock=(P/'pilot_runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
S='/Users/zy/.codex/tools/phase2s_sppider/bin/python';X='/Users/zy/.codex/tools/phase2s_xpair/bin/python'
env=dict(os.environ,HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD='1',XPAIR_ANKH_LOCAL_DIR=str(R/'external_models/ankh-large'),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MPLCONFIGDIR=str(P/'matplotlib_cache'))
def state(stage,status,**kwargs):
    data=dict(stage=stage,status=status,pid=os.getpid(),updated_at=datetime.datetime.now().astimezone().isoformat(),**kwargs)
    tmp=P/'pilot_runner_status.tmp';tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(P/'pilot_runner_status.json');print(json.dumps(data),flush=True)
def run(stage,cmd):
    state(stage,'RUNNING')
    with (L/(stage+'.log')).open('a') as stream:
        stream.write('\nRUN '+datetime.datetime.now().astimezone().isoformat()+'\n');stream.flush()
        subprocess.run(cmd,cwd=R,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
    state(stage,'PASS')
def xp(stage,out):
    (P/out).mkdir(exist_ok=True)
    run(stage,[X,'scripts/phase2SR/xpair_preserve_output.py','--predict_pairs_path',str(P/'pilot_inputs/xpair_pairs.tsv'),'--predict_fasta_path',str(P/'pilot_inputs/xpair_sequences.fasta'),'--predict_output_path',str(P/out/'predictions.csv'),'--emb_dir_prediction',str(P/'xpair_pilot/embeddings'),'--num_workers','0','--batch_size_prediction','1','--interaction_threshold','0'])
try:
    run('sppider_example_audit',[S,'scripts/phase2SR/audit_examples.py','sppider'])
    run('xpair_example_audit',[X,'scripts/phase2SR/audit_examples.py','xpair'])
    run('legacy_independent_audit',[S,'scripts/phase2SR/audit_legacy.py'])
    run('sppider_pilot',[S,'scripts/phase2SR/sppider_pilot.py'])
    xp('xpair_pilot','xpair_pilot')
    xp('xpair_pilot_repeat','xpair_pilot_repeat')
    run('pilot_final_audit',[S,'scripts/phase2SR/audit_pilot.py'])
    run('recovery_gate_finalize',[S,'scripts/phase2SR/finish_recovery_gate.py'])
    state('RECOVERY_COMPLETE','GO_RESUME_PHASE2S',full_benchmark_started=False)
except Exception as e:
    state('RECOVERY_STOPPED','NEEDS_DIAGNOSIS',error=str(e));traceback.print_exc();raise
