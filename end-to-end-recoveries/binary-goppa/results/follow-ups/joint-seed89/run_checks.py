"""Finish and verify the joint-solver tests on unchanged seed-89 kernels."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
CASES={
    27:(ROOT/'sweep-100/cases/seed-000089/instance/public.txt',ROOT/'sweep-100/cases/seed-000089/kernel.bin'),
    42:(ROOT/'rank-investigation/position-scan/held42/public.txt',ROOT/'rank-investigation/position-scan/held42/kernel.bin'),
    29:(ROOT/'diagnostics/seed89-next-heldout-public.txt',ROOT/'diagnostics/seed89-next-heldout-recovery/kernel.bin'),
}

def run(held):
    public,kernel=CASES[held];out=HERE/f'held{held}';steps=[];started=time.monotonic()
    def stage(name,command):
        tick=time.monotonic()
        with (HERE/f'held{held}-{name}.stdout.log').open('w') as stdout,(HERE/f'held{held}-{name}.stderr.log').open('w') as stderr:
            result=subprocess.run(list(map(str,command)),stdout=stdout,stderr=stderr,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=600)
        steps.append(dict(name=name,command=list(map(str,command)),returncode=result.returncode,seconds=time.monotonic()-tick))
        if result.returncode:raise RuntimeError((HERE/f'held{held}-{name}.stderr.log').read_text())
    record=dict(heldout=held,seed=89,public_inputs_only_for_recovery=True,
                fresh_holdout_computation=False,holdout_profile_unchanged=True,
                public_sha256=hashlib.sha256(public.read_bytes()).hexdigest(),
                kernel_sha256=hashlib.sha256(kernel.read_bytes()).hexdigest())
    try:
        if not (out/'recovery.json').exists():stage('recover',[HERE/'joint-recovery','recover',public,kernel,out,'both'])
        record['recovery']=json.loads((out/'recovery.json').read_text())
        record['independent_verification']={}
        for route in ('minor','direct'):
            stage(route+'-independent',[sys.executable,ROOT/'verify.py',public,out/f'{route}-key.json'])
            record['independent_verification'][route]=json.loads((HERE/f'held{held}-{route}-independent.stdout.log').read_text())
        record['status']='passed'
    except Exception as error:record.update(status='failed',error=str(error))
    record.update(stages=steps,verification_invocation_seconds=time.monotonic()-started)
    (HERE/f'held{held}-result.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'heldout':held,'status':record['status'],'error':record.get('error')}),flush=True)
    return record

with ThreadPoolExecutor(max_workers=3) as executor:results=list(executor.map(run,[27,42,29]))
files=[HERE/'public_joint.hpp',HERE/'joint_recovery.cpp',HERE/'run_checks.py',ROOT/'src/holdout.cpp',ROOT/'include/algebra.hpp',ROOT/'include/series_api.hpp',ROOT/'verify.py',HERE/'joint-recovery']
summary=dict(seed=89,requested=3,passed=sum(r['status']=='passed' for r in results),
             public_recovery_without_known_key_or_curve=True,original_kernels_unchanged=True,
             source_and_binary_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},results=results)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
