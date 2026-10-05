from common import *
import os,subprocess,time,fcntl,traceback
lock=(O/'full_pipeline.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
S='/Users/zy/.codex/tools/phase2s_sppider/bin/python';env=dict(os.environ,HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MPLCONFIGDIR=str(O/'matplotlib_cache'))
def alive(pid,script):
    r=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True);return r.returncode==0 and script in r.stdout
retry=[]
def run(name):
    checkpoint('pipeline',name,pid=os.getpid())
    with (L/(name+'.log')).open('a') as f:
        f.write('\nRUN '+now()+'\n');f.flush();subprocess.run([S,'-u',str(R/'scripts/phase2S_full'/(name+'.py'))],env=env,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
try:
    import argparse
    stages=['extract','coverage','structural_labels','validation','specificity','shuffle','failures_cost','figures','finalize']
    parser=argparse.ArgumentParser();parser.add_argument('--resume-from',choices=stages);args=parser.parse_args()
    if not args.resume_from:
        checkpoint('pipeline','SPPIDER_RUNNING',pid=os.getpid())
        while True:
            statefile=O/'sppider_state.json';state=json.loads(statefile.read_text()) if statefile.exists() else {}
            if state.get('stage')=='COMPLETE':break
            if alive(state.get('pid',22251),'sppider_full.py'):time.sleep(45);continue
            run('sppider_full');break
        run('xpair_full')
    for name in stages[stages.index(args.resume_from) if args.resume_from else 0:]:run(name)
    checkpoint('pipeline','COMPLETE',pid=os.getpid(),phase2P_started=False)
except Exception as e:checkpoint('pipeline','NEEDS_DIAGNOSIS',pid=os.getpid(),error=repr(e));traceback.print_exc();raise
