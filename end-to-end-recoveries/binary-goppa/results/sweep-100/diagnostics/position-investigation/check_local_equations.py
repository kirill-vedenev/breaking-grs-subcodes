"""Independent finite-field checks of the exported local equations."""
from pathlib import Path
import json

base=Path(__file__).resolve().parent
def mult(a,b):
    c=0
    while b:
        if b&1:c^=a
        b>>=1;a<<=1
        if a&64:a^=67
    return c
M=[[mult(a,b) for b in range(64)] for a in range(64)]
def power(a,n):
    c=1
    while n:
        if n&1:c=M[c][a]
        a=M[a][a];n>>=1
    return c
inv=[0]+[power(a,62) for a in range(1,64)]
def rref(A):
    A=[list(row) for row in A];n=0
    for col in range(len(A[0])):
        pivot=next((i for i in range(n,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[n],A[pivot]=A[pivot],A[n]
        A[n]=[M[inv[A[n][col]]][x] for x in A[n]]
        for i in range(len(A)):
            if i!=n and A[i][col]:A[i]=[a^M[A[i][col]][b] for a,b in zip(A[i],A[n])]
        n+=1
        if n==len(A):break
    return A[:n]
def read(path):
    w=path.read_text().split();assert w[0]=='GF64M1'
    nr,nc=map(int,w[1:3]);values=list(map(int,w[3:]));assert len(values)==nr*nc
    return [values[i*nc:(i+1)*nc] for i in range(nr)]

pairs=[(a,b) for a in range(6) for b in range(a,6)]
reports=[]
for held in (27,29):
    root=base/'rank-investigation'/f'held{held}'
    J=read(root/'J.txt');Q=read(root/'branch-quadrics.txt');S=read(root/'stationary-combinations.txt')
    assert len(J)==len(Q)==484 and all(set(row)<={0,1} for row in S)
    SQ=[]
    for combination in S:
        row=[0]*21
        for coefficient,q in zip(combination,Q):
            if coefficient:row=[a^b for a,b in zip(row,q)]
        SQ.append(row)
    R=rref(Q)
    monomials=[]
    for row in R:
        support=[(pairs[c],v) for c,v in enumerate(row) if v]
        assert len(support)==1 and support[0][1]==1
        monomials.append(support[0][0])
    expected=[(a,b) for a in range(6) for b in range(a+1,6) if held==29 or (b-a)%6 not in (1,5)]
    assert monomials==expected
    report={'heldout':held,'all_polynomials':len(Q),'zero_gradient_subspace_dimension':len(S),
            'all_polynomial_quadratic_rank':len(R),'stationary_quadratic_rank':len(rref(SQ)),
            'gradient_rank':len(rref(J)),'quadratic_monomial_basis':monomials}
    if held==27:
        cubic=[]
        for edge in range(6):
            C=read(root/f'edge-{edge}-cubic.txt')
            assert all(row[0]==row[3]==0 for row in C)
            assert rref(C)==[[0,1,0,0],[0,0,1,0]]
            cubic.append({'edge':[edge,(edge+1)%6],'homogeneous_equations':['a^2*b','a*b^2'],'rank':2})
        report['third_order_normalized_jet_constraints']=cubic
    reports.append(report)
(base/'independent-local-equations.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
