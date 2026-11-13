"""Independent degree-seven Taylor and series checks against Sage polynomials."""
from sage.all import GF, PolynomialRing, vector, prod
from pathlib import Path
from itertools import combinations
from math import comb
from random import Random
import argparse
import json
import os
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'recovery'))
from extend_jet import Cache


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError('test output exists')
    args.output.mkdir(parents=True)
    rng=Random(2026101611);k,d,count,depth=9,7,11,16
    E=GF(64,'a',modulus=PolynomialRing(GF(2),'z')([1,1,0,0,0,0,1]));R=PolynomialRing(E,'u')
    mons=sorted(sum(1<<i for i in c) for c in combinations(range(k),d))
    panel=np.array([rng.getrandbits(count) for _ in mons],dtype='<u8')
    panel.tofile(args.output/'panel.bin');np.array(mons,dtype='<u8').tofile(args.output/'monomials.bin')
    point=rng.getrandbits(k)
    subprocess.run([str(ROOT/'build/extract-jets'),str(args.output/'panel.bin'),
                    str(args.output/'monomials.bin'),str(k),str(d),str(point),str(count),
                    str(args.output/'taylor')],check=True)
    taylor_checks=0
    for u in range(d):
        derivatives=sorted(sum(1<<i for i in c) for c in combinations(range(k),u))
        got=np.fromfile(args.output/'taylor'/f'hasse-{u}.u64le',dtype='<u8')
        for derivative,value in zip(derivatives,got):
            expected=0
            for mon,word in zip(mons,panel):
                if derivative & mon == derivative and not ((mon^derivative)&~point):expected^=int(word)
            assert expected==int(value)
            taylor_checks+=count
    os.environ['TII_JET_SERIES_LIBRARY']=str(ROOT/'build/libseries.so')
    cache=Cache(k,d,E,depth,count,2,args.output/'panel.bin')
    arcs=[[E.zero()]*(depth+1) for _ in range(k)];checks=0
    def expected(r):
        values=[prod(R(arcs[a]) for a in range(k) if mon>>a&1)[r] for mon in mons]
        return vector(E,[sum((v for v,w in zip(values,panel) if int(w)>>j&1),E.zero()) for j in range(count)])
    try:
        for r in range(3):
            jet=vector(E,[E.from_integer(rng.randrange(64)) for _ in range(k)]);cache.set(r,jet)
            for a in range(k):arcs[a][r]=jet[a]
        for r in range(3):assert cache.evaluate(r)==expected(r);checks+=1
        for r in range(3,depth,2):
            assert cache.pair(r)==[expected(r),expected(r+1)];checks+=2
            odd=vector(E,[E.from_integer(rng.randrange(64)) for _ in range(k)])
            even=vector(E,[E.from_integer(rng.randrange(64)) for _ in range(k)])
            got=cache.patch(r,odd,even)
            for a in range(k):arcs[a][r],arcs[a][r+1]=odd[a],even[a]
            assert got==[expected(r),expected(r+1)];checks+=2
            assert cache.pair(r)==got;checks+=2
    finally:cache.close()
    report=dict(status='passed',degree=d,k=k,field_order=64,arbitrary_series=True,
                exact_taylor_bits_checked=taylor_checks,exact_series_field_elements_checked=checks*count,
                oracle='independent Sage polynomial products',public_end_to_end_recovery=False)
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':main()
