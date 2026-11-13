#!/usr/bin/env python3
"""Verify saved recovered keys using public JSON only; imports no recovery code.

Lagrange dual multipliers give parity equations. Frobenius closure of their
row space gives the exact subfield-code dimension for alternant instances.
"""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, MatrixSpace, matrix, prod


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(case):
    public_path=case/'public.json'
    data=json.loads(public_path.read_text())
    p=data['characteristic']; modulus=data['modulus']; degree=len(modulus)-1
    if degree==1:E=GF(p)
    else:E=GF(p**degree,'a',modulus=PolynomialRing(GF(p),'z')(modulus))
    def dec(value):
        if degree==1:return E(value)
        return sum(E((value//(p**i))%p)*E.gen()**i for i in range(degree))
    def mat(rows):
        A=matrix(E,rows)
        return MatrixSpace(E,A.nrows(),A.ncols(),implementation='generic')(A.list()) if p!=2 and degree>1 else A
    Y=mat([[dec(a) for a in row] for row in data['matrix']])
    q,m,n,k,D=[data[s] for s in ('q','m','n','k','D')]
    assert Y.dimensions()==(k,n) and Y.rank()==k
    assert all(a**q==a for a in Y.list())
    result=dict(label=data['label'],seed=data['seed'],public_sha256=digest(public_path),public_only=True,independent_of_recovery=True,routes={})
    for route in ('minor','direct'):
        path=case/(route+'-key.json')
        if not path.exists():continue
        key=json.loads(path.read_text())
        support=[dec(a) for a in key['support']]; multiplier=[dec(a) for a in key['multiplier']]
        assert len(support)==len(multiplier)==n and len(set(support))==n and all(multiplier)
        derivative=[prod(support[j]-support[l] for l in range(n) if j!=l) for j in range(n)]
        def parity(weights):
            beta=[1/(weights[j]*derivative[j]) for j in range(n)]
            return mat([[beta[j]*support[j]**i for j in range(n)] for i in range(n-D-1)])
        H=parity(multiplier)
        assert H.rank()==n-D-1 and (Y*H.transpose()).is_zero()
        check=dict(status='passed',key_sha256=digest(path),distinct_finite_support=True,public_grs_containment=True)
        if m>1:
            closure=mat([[a**(q**i) for a in row] for i in range(m) for row in H])
            assert closure.rank()==n-k
            check.update(exact_public_subfield_code_equality=True,subfield_parity_rank=int(closure.rank()))
        changed=list(multiplier);changed[0]*=next(a for a in E if a not in (0,1))
        assert not (Y*parity(changed).transpose()).is_zero()
        check['tampered_multiplier_rejected']=True
        result['routes'][route]=check
    record=json.loads((case/'case.json').read_text())
    if record['status']=='passed':assert len(result['routes'])==(2 if p==2 else 1)
    result['status']='passed' if result['routes'] else 'no-key'
    (case/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('root',type=Path);args=ap.parse_args()
    results=[]
    for path in sorted(args.root.glob('cases/*/seed-*/case.json')):
        result=verify(path.parent);results.append(result)
        print(json.dumps(dict(label=result['label'],seed=result['seed'],status=result['status'])),flush=True)
    summary=dict(cases=len(results),cases_with_keys=sum(r['status']=='passed' for r in results),routes_verified=sum(len(r['routes']) for r in results),
                 verifier_sha256=digest(Path(__file__)),public_only=True,independent_of_recovery=True,records=results)
    (args.root/'independent-verification.json').write_text(json.dumps(summary,indent=2)+'\n')
