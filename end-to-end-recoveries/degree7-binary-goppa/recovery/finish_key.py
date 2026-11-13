"""Restore all original coordinates from either publicly recovered support.

The public infinity column fixes an affine chart. No generation secret,
other route's recovered support, or known Goppa polynomial is an input.
"""
from sage.all import GF, PolynomialRing, matrix, vector, prod, gcd
from pathlib import Path
import argparse
import hashlib
import json
import time


def enc(x):
    return int(x.to_integer())


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, default=int)+'\n')


def binary_expansion(A):
    m=A.base_ring().degree()
    return matrix(GF(2),A.nrows()*m,A.ncols(),
                  lambda r,j:(enc(A[r//m,j])>>(r%m))&1)


def multiplier_matrix(G, support, r):
    E=G.base_ring()
    V=matrix(E,r,len(support),lambda a,j:support[j]**a)
    return matrix(E,G.nrows()*r,len(support),
                  lambda a,j:G[a//r,j]*V[a%r,j])


def multipliers(G, support, r):
    A=multiplier_matrix(G,support,r)
    W=A.right_kernel_matrix()
    if W.nrows()!=1 or any(not x for x in W.row(0)):
        raise ValueError('public multiplier system does not have one full-support solution')
    nu=W.row(0)/W[0,0]
    assert not A*nu
    return nu, A.rank()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--public',type=Path,required=True)
    ap.add_argument('--support',type=Path,required=True)
    ap.add_argument('--original',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();begin=time.monotonic()
    if args.output.exists():
        raise FileExistsError('key output already exists')
    args.output.mkdir(parents=True)
    public=json.loads(args.public.read_text());candidate=json.loads(args.support.read_text())
    raw=args.original.read_bytes();original=json.loads(raw)
    digest=hashlib.sha256(raw).hexdigest()
    assert digest==public['public_sha256']
    assert original['parameters']==public['original_parameters']
    assert original['field_modulus']==public['field_modulus']
    assert candidate['original_indices']==public['original_indices']
    m,t,n,kfull=[original['parameters'][a] for a in ('m','t','n','k')]
    ns,k=[public['parameters'][a] for a in ('n','k')]
    D=ns-2*t-1;r=2*t
    E=GF(2**m,'a',modulus=PolynomialRing(GF(2),'x')(public['field_modulus']))
    dec=E.from_integer;R=PolynomialRing(E,'x');x=R.gen()
    P=matrix(GF(2),public['generator']);H=matrix(GF(2),original['parity_check'])
    G=matrix(GF(2),original['generator']);keep=public['original_indices'];removed=public['shortened_original_indices']
    assert P.dimensions()==(k,ns) and G.dimensions()==(kfull,n)
    assert G.rank()==kfull and H.rank()==n-kfull and not G*H.transpose()
    assert sorted(keep+removed)==list(range(n)) and len(set(keep+removed))==n
    assert H.matrix_from_columns(keep).right_kernel()==P.row_space()
    beta=[None if v is None else dec(v) for v in candidate['support']]
    assert len(beta)==ns and len(set(beta))==ns
    # First obtain any chart finite on the retained support.
    preliminary_pole=next((dec(i) for i in range(E.order()) if dec(i) not in beta),None)
    if None in beta:
        if preliminary_pole is None:raise ValueError('no affine chart for retained support')
        alpha=[E.zero() if z is None else 1/(z+preliminary_pole) for z in beta]
    else:
        alpha=list(beta)
    PE=P.change_ring(E)
    nu,rank=multipliers(PE,alpha,r)
    parity=matrix(E,r,ns,lambda a,j:nu[j]*alpha[j]**a)
    assert binary_expansion(parity).right_kernel()==P.row_space()
    Pi=prod(x-a for a in alpha);partials=[Pi//(x-a) for a in alpha]
    F=[sum((PE[a,j]*nu[j]*partials[j] for j in range(ns)),R.zero()) for a in range(k)]
    assert max(f.degree() for f in F)<=D and gcd(F).degree()==0
    for j in range(ns):
        scale=1/(nu[j]*Pi.derivative()(alpha[j]))
        assert scale*vector(E,[f(alpha[j]) for f in F])==PE.column(j)
    w=sum(PE.columns(),vector(E,k))
    assert w==vector(E,[(public['infinity_mask']>>a)&1 for a in range(k)])
    anchor=next(a for a,z in enumerate(w) if z)
    matches=[]
    for z in [dec(i) for i in range(E.order())]+[None]:
        value=vector(E,[f[D] if z is None else f(z) for f in F])
        if value[anchor] and value==value[anchor]/w[anchor]*w:
            matches.append(z)
    if len(matches)!=1 or matches[0] in alpha:
        raise ValueError('public infinity has no unique unused curve parameter')
    pole=matches[0]
    if pole is not None:alpha=[1/(a+pole) for a in alpha]
    nu,rank=multipliers(PE,alpha,r)
    assert rank==ns-1
    chart=dict(public_only=True,public_infinity_unique=True,
               public_infinity_parameter_before_final_chart=None if pole is None else enc(pole),
               projective_candidates_checked=int(E.order()+1),
               retained_support=[enc(a) for a in alpha],retained_multiplier_rank=int(rank))
    dump(args.output/'affine-chart.json',chart)
    powers=matrix(E,r,ns,lambda a,j:alpha[j]**a)
    support=[None]*n
    for j,a in zip(keep,alpha):support[j]=a
    candidates=[dec(i) for i in range(E.order()) if dec(i) not in set(alpha)]+[None]
    extension_records=[]
    for j in removed:
        positions=keep+[j]
        Gj=H.matrix_from_columns(positions).right_kernel_matrix()
        assert Gj.nrows()==k+1
        extra=next(row for row in Gj.rows() if row[-1]).change_ring(E)
        extra/=extra[-1]
        embedded=PE.augment(matrix(E,k,1))
        assert embedded.stack(matrix(E,[extra])).row_space()==Gj.change_ring(E).row_space()
        moments=powers*vector(E,[extra[a]*nu[a] for a in range(ns)])
        pivot=next((a for a,v in enumerate(moments) if v),None)
        if pivot is None:raise ValueError('zero extension moments')
        accepted=[]
        for z in candidates:
            target=vector(E,[int(a==r-1) if z is None else z**a for a in range(r)])
            scale=target[pivot]/moments[pivot]
            if scale and scale*moments==target:accepted.append(z)
        if len(accepted)!=1 or accepted[0] is None:
            raise ValueError(f'coordinate {j} has no unique finite extension')
        support[j]=accepted[0]
        extension_records.append(dict(original_position=j,accepted=enc(accepted[0]),
                                      projective_candidates=len(candidates),unique=True))
    assert all(a is not None for a in support) and len(set(support))==n
    full_nu,full_rank=multipliers(G.change_ring(E),support,r)
    full_parity=matrix(E,r,n,lambda a,j:full_nu[j]*support[j]**a)
    assert binary_expansion(full_parity).row_space()==H.row_space()
    square=R.lagrange_polynomial([(a,1/v) for a,v in zip(support,full_nu)])
    if square.degree()!=2*t or square.derivative():
        raise ValueError('reciprocal multipliers are not a degree-2t square')
    g=R([square[2*a]**(2**(m-1)) for a in range(t+1)]).monic()
    assert g.degree()==t and g.is_irreducible()
    assert square==square.leading_coefficient()*g*g
    gp=matrix(E,t,n,lambda a,j:support[j]**a/g(support[j]))
    assert binary_expansion(gp).row_space()==H.row_space()
    key=dict(schema='equivalent-binary-goppa-key-v1',parameters=original['parameters'],
             field_modulus=public['field_modulus'],public_sha256=digest,
             support=[enc(a) for a in support],g=[enc(a) for a in g.list()],
             dual_grs_multipliers=[enc(a) for a in full_nu],public_only_recovery=True)
    dump(args.output/'equivalent-key.json',key)
    report=dict(status='passed',public_only=True,secret_reference_read=False,
                full_support_count=n,shortened_support_count=ns,restored_positions=len(removed),
                extension_tests=extension_records,full_multiplier_rank=int(full_rank),
                exact_original_public_code_equality=True,goppa_degree=t,goppa_irreducible=True,
                original_public_sha256=digest,support_source_sha256=hashlib.sha256(args.support.read_bytes()).hexdigest(),
                seconds=time.monotonic()-begin)
    dump(args.output/'recovery-report.json',report)
    print(json.dumps(report,default=int),flush=True)


if __name__=='__main__':main()
