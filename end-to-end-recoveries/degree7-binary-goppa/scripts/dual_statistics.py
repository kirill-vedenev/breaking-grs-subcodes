"""Public-only exhaustive dual weight counts; these never filter the cohort."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'scripts'))
from algebra import rref,kernel,restrict


def masks(rows): return [sum(int(x)<<j for j,x in enumerate(row)) for row in rows]

def sums(basis):
    out=np.zeros(1,dtype=np.uint64)
    for row in basis:out=np.concatenate([out,out^np.uint64(row)])
    return out


def main():
    start=time.monotonic()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='Fresh diagnostic directory')
    args=parser.parse_args();output=args.output.resolve()
    target=ROOT/'build/dual-profile'
    if not target.is_file():raise FileNotFoundError('Run make cpu first')
    output.mkdir(parents=True,exist_ok=False)
    reports=[]
    for trial in json.loads((ROOT/'TRIALS.json').read_text())['trials']:
        p=ROOT/trial['inputs'];public=json.loads((p/'public.json').read_text())
        original=json.loads((p/'original-public.json').read_text())
        P=masks(public['generator']);full=masks(original['parity_check'])
        dual=[restrict(row,public['original_indices']) for row in full]
        assert len(rref(full,64)[0])==len(rref(dual,54)[0])==24
        assert rref(dual,54)[0]==rref(kernel(P,54),54)[0]
        path=output/(trial['trial_id']+'-dual-input.txt')
        path.write_text(f"24 64 54 {trial['heldout_shortened']}\n"+''.join(f'{a} {b}\n' for a,b in zip(full,dual)))
        native=json.loads(subprocess.check_output([target,path],text=True))
        left,right=sums(dual[:12]),sums(dual[12:])
        distribution=np.zeros(55,dtype=np.int64)
        for first in range(0,len(left),128):
            words=(left[first:first+128,None]^right[None,:]).reshape(-1)
            distribution+=np.bincount(np.bitwise_count(words),minlength=55)
        assert distribution.tolist()==native['shortened_weight_distribution']
        assert int(distribution.sum())==2**24
        # Appending the public infinity column gives the other coset of the
        # extended dual: complement the finite vector and append a 1.
        assert all(((row.bit_count()&1)==((public['infinity_mask']>>j)&1)) for j,row in enumerate(P))
        low=int(distribution[1:13].sum());inf=int(distribution[43:].sum())
        report=dict(trial_id=trial['trial_id'],key_seed=trial['key_seed'],
                    public_input_sha256=hashlib.sha256((p/'public.json').read_bytes()).hexdigest(),
                    public_only=True,private_audit_read=False,key_filtered_or_replaced=False,
                    enumerated_dual_vectors=2**24,independent_enumerations_agree=True,
                    minimum_dual_weight=next(i for i in range(1,55) if distribution[i]),
                    finite_nonzero_weight_at_most_12=low,
                    additional_weight_at_most_12_involving_infinity=inf,
                    all_point_weight_at_most_12=low+inf,
                    actual_Weis_failure_demonstrated=False,
                    finite_weight_distribution=distribution.tolist(),native=native)
        reports.append(report)
        print(json.dumps({k:v for k,v in report.items() if k not in ['finite_weight_distribution','native']}),flush=True)
    (output/'dual-statistics.json').write_text(json.dumps(dict(status='passed',trials=reports,seconds=time.monotonic()-start),indent=2)+'\n')


if __name__=='__main__': main()
