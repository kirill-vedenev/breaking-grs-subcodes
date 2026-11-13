"""Literal Hasse derivatives and packed GF(2) dots; no native recovery imports."""
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('public',type=Path)
parser.add_argument('kernel',type=Path)
parser.add_argument('export',type=Path)
parser.add_argument('--report',type=Path)
args=parser.parse_args();started=time.monotonic()
w=args.public.read_text().split();assert w[0]=='BGPK1'
m,t,n,k,mod,held,d,s=map(int,w[1:9]);points=list(map(int,w[9:]));point=points[held]
blob=args.kernel.read_bytes();assert blob[:4]==b'WHK1'
kk,dd,N,M=struct.unpack('<4I',blob[4:20]);assert (kk,dd)==(k,d)
stride=(N+7)//8;assert len(blob)==20+M*stride
kernel=[int.from_bytes(blob[20+i*stride:20+(i+1)*stride],'little') for i in range(M)]
derivatives=[0]+[1<<a for a in range(k)]+[(1<<a)|(1<<b) for a,b in combinations(range(k),2)]
index={mask:i for i,mask in enumerate(derivatives)}
rows=[bytearray(stride) for _ in derivatives]
count=0
for c,coordinates in enumerate(combinations(range(k),d)):
    monomial=sum(1<<a for a in coordinates);missing=monomial&~point
    missing_count=missing.bit_count();count+=1
    if missing_count>2:continue
    rest=[1<<a for a in coordinates if not missing&(1<<a)]
    for extra_degree in range(3-missing_count):
        for extra in combinations(rest,extra_degree):
            deriv=missing|sum(extra)
            rows[index[deriv]][c//8]|=1<<(c%8)
assert count==N
taylor=[int.from_bytes(row,'little') for row in rows]
native=[]
for name in ('J.txt','Q.txt'):
    text=(args.export/name).read_text().split();assert text[0]=='GF64M1'
    nr,nc=map(int,text[1:3]);values=list(map(int,text[3:]));assert nr==M and len(values)==nr*nc
    native.append([values[r*nc:(r+1)*nc] for r in range(nr)])
checks=0
for i,polynomial in enumerate(kernel):
    assert (polynomial&taylor[0]).bit_count()%2==0
    actual=[(polynomial&row).bit_count()%2 for row in taylor[1:]]
    assert actual==native[0][i]+native[1][i],(i,'Taylor mismatch')
    checks+=len(actual)
result={'status':'passed','heldout':held,'polynomials':M,'monomials':N,
        'derivative_coefficients_checked':checks,'constant_coefficients_checked':M,
        'native_gradient_and_quadratic_coefficients_match':True,
        'method':'literal squarefree Hasse derivative formula and Python integer AND/popcount',
        'imports_native_recovery':False,'uses_m4ri':False,
        'kernel_sha256':hashlib.sha256(blob).hexdigest(),'seconds':time.monotonic()-started}
if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
