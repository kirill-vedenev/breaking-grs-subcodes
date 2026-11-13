"""Does removing the other exceptional point restore the missing polar equations?"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json,os,subprocess,time
root=Path(__file__).resolve().parents[1]
base=root/'rank-investigation/double-holdout';base.mkdir(exist_ok=False)
words=(root/'sweep-100/cases/seed-000089/instance/public.txt').read_text().split()
header=list(map(int,words[1:9]));cols=list(map(int,words[9:]));key=root/'diagnostics/seed89-next-heldout-recovery/recovery/direct-key.txt'
def run(other):
    out=base/f'27-{other}';out.mkdir();start=time.monotonic();steps=[]
    retained=[c for j,c in enumerate(cols) if j not in (27,other) and c.bit_count()!=1]
    operator=out/'holdout.txt';operator.write_text(f'28 5 {len(retained)}\n'+''.join(f'{c} 4\n' for c in retained))
    for j in (27,other):
        h=list(header);h[5]=j;(out/f'public-{j}.txt').write_text('BGPK1\n'+' '.join(map(str,h))+'\n'+' '.join(map(str,cols))+'\n')
    commands=[('holdout',[root/'build/holdout',operator,out/'kernel.bin',3]),('filter',[root/'rank-investigation/filter-curve',out/'public-27.txt',out/'kernel.bin',key,out/'curve-kernel.bin'])]
    commands.extend((f'analyze-{j}',[root/'rank-investigation/analyze',out/f'public-{j}.txt',out/'curve-kernel.bin',key,out/f'held{j}']) for j in (27,other))
    status='passed'
    for name,cmd in commands:
        t=time.monotonic()
        with (out/f'{name}.stdout.log').open('w') as stdout,(out/f'{name}.stderr.log').open('w') as stderr:
            result=subprocess.run(list(map(str,cmd)),stdout=stdout,stderr=stderr,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=900)
        steps.append({'stage':name,'seconds':time.monotonic()-t,'returncode':result.returncode})
        if result.returncode:status='error';break
    record={'omitted_positions':[27,other],'status':status,'seconds':time.monotonic()-start,'steps':steps}
    for j in (27,other):
        p=out/f'held{j}/report.json'
        if p.exists():record[f'held{j}']=json.loads(p.read_text())
    (out/'result.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[42,29]))
