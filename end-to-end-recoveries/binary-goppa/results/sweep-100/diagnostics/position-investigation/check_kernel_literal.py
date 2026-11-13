"""Independent literal Hasse verification of every saved kernel polynomial."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
import multiprocessing
from pathlib import Path
import struct
import time

def check_point(item):
    position,point=item
    rows=[bytearray(STRIDE) for _ in DERIV]
    for c,(monomial,bits) in enumerate(MONS):
        missing=monomial&~point;count=missing.bit_count()
        if count>3:continue
        rest=[bit for bit in bits if not missing&bit]
        for extra_degree in range(4-count):
            for extra in combinations(rest,extra_degree):
                rows[INDEX[missing|sum(extra)]][c//8]|=1<<(c%8)
    checked=0
    for derivative,row in zip(DERIV,rows):
        literal=int.from_bytes(row,'little')
        for i,polynomial in enumerate(KERNEL):
            if (literal&polynomial).bit_count()&1:
                raise ValueError(f'position {position}, derivative {derivative}, polynomial {i}')
            checked+=1
    return dict(position=position,checked=checked,status='passed')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('public',type=Path);parser.add_argument('kernel',type=Path)
    parser.add_argument('--report',type=Path,required=True);parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args();started=time.monotonic()
    w=args.public.read_text().split();m,t,n,k,mod,held,d,s=map(int,w[1:9]);points=list(map(int,w[9:]));assert s==4
    blob=args.kernel.read_bytes();assert blob[:4]==b'WHK1'
    kk,dd,N,M=struct.unpack('<4I',blob[4:20]);assert (kk,dd)==(k,d)
    STRIDE=(N+7)//8;assert len(blob)==20+M*STRIDE
    KERNEL=[int.from_bytes(blob[20+i*STRIDE:20+(i+1)*STRIDE],'little') for i in range(M)]
    pivots={}
    for polynomial in KERNEL:
        while polynomial:
            p=polynomial.bit_length()-1
            if p not in pivots:pivots[p]=polynomial;break
            polynomial^=pivots[p]
    assert len(pivots)==M
    MONS=[(sum(1<<a for a in c),tuple(1<<a for a in c)) for c in combinations(range(k),d)]
    assert len(MONS)==N
    DERIV=[sum(1<<a for a in c) for r in range(s) for c in combinations(range(k),r)]
    INDEX={mask:i for i,mask in enumerate(DERIV)}
    nonunit=[(j,p) for j,p in enumerate(points) if j!=held and p.bit_count()!=1]
    units=[j for j,p in enumerate(points) if j!=held and p.bit_count()==1]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('fork')) as executor:
        results=list(executor.map(check_point,nonunit))
    record=dict(status='passed',heldout=held,kernel_polynomials=M,kernel_rank=len(pivots),
                retained_nonunit_points=len(nonunit),retained_unit_points=len(units),
                unit_conditions='Every squarefree quintic has multiplicity at least four at every unit vector.',
                full_hasse_scalar_conditions_checked=sum(r['checked'] for r in results),
                hasse_orders=[0,1,2,3],uses_m4ri=False,imports_native_recovery=False,
                method='literal Hasse monomial formula, Python integer dot products and XOR elimination',
                seconds=time.monotonic()-started,points=results)
    args.report.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='points'}),flush=True)
