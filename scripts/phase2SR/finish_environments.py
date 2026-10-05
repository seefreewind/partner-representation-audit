"""Finish already launched installs; never load or download model weights."""
from pathlib import Path
import subprocess,time
R=Path(__file__).resolve().parents[2]
for method,pid in [('sppider',71689),('xpair',71688)]:
    while True:
        command=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True).stdout
        if 'uv pip install' not in command or ('phase2s_'+method) not in command:break
        time.sleep(30)
    with (R/'logs/phase2SR'/('environment_'+method+'_validation.log')).open('w') as log:
        a=subprocess.run(['/Users/zy/.codex/tools/bin/uv','pip','check','--python',str(R/('phase2s_'+method)/'bin/python')],stdout=log,stderr=log,cwd=R)
        b=subprocess.run([str(R/('phase2s_'+method)/'bin/python'),str(R/'scripts/phase2SR/verify_environment.py'),method],stdout=log,stderr=log,cwd=R)
        print('package_check_exit',a.returncode,'import_check_exit',b.returncode,file=log,flush=True)
subprocess.run(['python3',str(R/'scripts/phase2SR/finalize_recovery.py'),'TECHNICAL_HOLD_BOTH_ENCODERS'],cwd=R)
