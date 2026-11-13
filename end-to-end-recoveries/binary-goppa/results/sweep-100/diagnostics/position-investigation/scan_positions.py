"""Compare all non-information holdouts and a changed-basis control on seed 89."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import time

root=Path(__file__).resolve().parents[1]
base=root/'rank-investigation/position-scan'
base.mkdir(exist_ok=False)
words=(root/'sweep-100/cases/seed-000089/instance/public.txt').read_text().split()
header=list(map(int,words[1:9]));columns=list(map(int,words[9:]))
key=root/'diagnostics/seed89-next-heldout-recovery/recovery/direct-key.txt'
positions=[j for j,c in enumerate(columns) if c.bit_count()!=1]
jobs=[(f'held{j}',j,columns) for j in positions if j not in (27,29)]

# Pick another information set from retained columns and invert its binary matrix.
pivots={};selected=[]
for j in reversed(range(64)):
    if j==27:continue
    v=columns[j]
    while v:
        p=v.bit_length()-1
        if p not in pivots:
            pivots[p]=v;selected.append(j);break
        v^=pivots[p]
    if len(selected)==28:break
assert len(selected)==28
aug=[sum(((columns[selected[j]]>>a)&1)<<j for j in range(28))|(1<<(28+a)) for a in range(28)]
for a in range(28):
    p=next(i for i in range(a,28) if (aug[i]>>a)&1)
    aug[a],aug[p]=aug[p],aug[a]
    for i in range(28):
        if i!=a and (aug[i]>>a)&1:aug[i]^=aug[a]
inverse=[r>>28 for r in aug]
changed=[sum(((r&c).bit_count()&1)<<a for a,r in enumerate(inverse)) for c in columns]
assert [changed[j] for j in selected]==[1<<a for a in range(28)]
jobs.append(('held27-rebased',27,changed))
(base/'specification.json').write_text(json.dumps({'positions':positions,'known_positions':[27,29],'new_jobs':len(jobs),'workers':36,'rebase_information_set':selected,'rebase_matrix_rows':inverse,'seed':89},indent=2)+'\n')

def run(spec):
    name,position,cols=spec;out=base/name;out.mkdir()
    h=list(header);h[5]=position
    public=out/'public.txt';public.write_text('BGPK1\n'+' '.join(map(str,h))+'\n'+' '.join(map(str,cols))+'\n')
    started=time.monotonic();record={'name':name,'position':position,'stages':[]}
    commands=[('prepare',[root/'build/toy-goppa','prepare',public,out/'holdout.txt']),
              ('holdout',[root/'build/holdout',out/'holdout.txt',out/'kernel.bin','3']),
              ('analyze',[root/'rank-investigation/analyze',public,out/'kernel.bin',key,out/'analysis'])]
    try:
        for stage,command in commands:
            t=time.monotonic()
            with (out/f'{stage}.stdout.log').open('w') as stdout,(out/f'{stage}.stderr.log').open('w') as stderr:
                p=subprocess.run(list(map(str,command)),stdout=stdout,stderr=stderr,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=900)
            record['stages'].append({'stage':stage,'seconds':time.monotonic()-t,'returncode':p.returncode})
            if p.returncode:raise RuntimeError((out/f'{stage}.stderr.log').read_text())
        record.update(status='completed',analysis=json.loads((out/'analysis/report.json').read_text()),holdout=json.loads((out/'holdout.stdout.log').read_text()))
    except Exception as error:record.update(status='error',error=str(error))
    record['seconds']=time.monotonic()-started
    (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    return record

results=[];started=time.monotonic()
with ThreadPoolExecutor(max_workers=36) as executor:
    for future in as_completed([executor.submit(run,spec) for spec in jobs]):
        record=future.result();results.append(record)
        print(json.dumps({'name':record['name'],'status':record['status'],'analysis':record.get('analysis'),'completed':len(results),'total':len(jobs)}),flush=True)
        (base/'summary.json').write_text(json.dumps({'completed':len(results),'requested':len(jobs),'seconds':time.monotonic()-started,'results':sorted(results,key=lambda r:r['name'])},indent=2)+'\n')
