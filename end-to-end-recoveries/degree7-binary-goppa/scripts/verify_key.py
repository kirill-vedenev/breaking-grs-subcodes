#!/usr/bin/env python3
"""Independent standard-library certificate check; imports no recovery code."""
import argparse
import hashlib
import json
from pathlib import Path


def verify(public_path,key_path):
    raw=public_path.read_bytes();public=json.loads(raw);key=json.loads(key_path.read_text())
    p=public['parameters'];m,t,n,k=[p[a] for a in ('m','t','n','k')]
    q=1<<m;modulus=sum(int(a)<<i for i,a in enumerate(public['field_modulus']))
    if key['parameters']!=p or key['field_modulus']!=public['field_modulus']:
        raise ValueError('key parameter mismatch')
    if key['public_sha256']!=hashlib.sha256(raw).hexdigest():raise ValueError('public digest mismatch')
    def multiply(a,b):
        z=0
        for bit in range(m):
            if b>>bit&1:z^=a<<bit
        for bit in reversed(range(m,2*m-1)):
            if z>>bit&1:z^=modulus<<(bit-m)
        return z
    def power(a,e):
        z=1
        while e:
            if e&1:z=multiply(z,a)
            a=multiply(a,a);e>>=1
        return z
    def evaluate(g,a):
        z=0
        for c in reversed(g):z=multiply(z,a)^c
        return z
    def trim(a):
        a=list(a)
        while a and not a[-1]:a.pop()
        return a
    def remainder(a,b):
        a=trim(a);inv=power(b[-1],q-2)
        while len(a)>=len(b):
            shift=len(a)-len(b);scale=multiply(a[-1],inv)
            for i,c in enumerate(b):a[shift+i]^=multiply(c,scale)
            a=trim(a)
        return a
    def gcd(a,b):
        while b:a,b=b,remainder(a,b)
        return a
    def reduced(rows,width):
        basis={}
        for row in rows:
            while row:
                pivot=row.bit_length()-1
                if pivot in basis:row^=basis[pivot]
                else:
                    basis[pivot]=row;break
        for pivot in sorted(basis):
            for other in basis:
                if other>pivot and basis[other]>>pivot&1:basis[other]^=basis[pivot]
        return tuple(basis[p] for p in sorted(basis))
    support,g=key['support'],key['g']
    if len(support)!=n or len(set(support))!=n or len(g)!=t+1 or g[-1]!=1:
        raise ValueError('support or polynomial shape')
    if any(type(a)!=int or a<0 or a>=q for a in support+g):raise ValueError('field encoding')
    # Rabin irreducibility test, independent of the generation implementation.
    divisors=[p for p in range(2,t+1) if t%p==0 and all(p%d for d in range(2,p))]
    h=[0,1]
    for step in range(1,t+1):
        for _ in range(m):
            sq=[0]*(2*len(h)-1)
            for i,c in enumerate(h):sq[2*i]=multiply(c,c)
            h=remainder(sq,g)
        if step in [t//p for p in divisors]:
            difference=h+[0]*max(0,2-len(h));difference[1]^=1
            if len(gcd(g,trim(difference)))!=1:raise ValueError('reducible Goppa polynomial')
    if h!=[0,1]:raise ValueError('Goppa Frobenius test failed')
    rows=[0]*(m*t)
    for j,a in enumerate(support):
        denominator=evaluate(g,a)
        if not denominator:raise ValueError('Goppa support root')
        v=power(denominator,q-2)
        if multiply(v,denominator)!=1:raise ValueError('field inversion failure')
        for degree in range(t):
            for bit in range(m):rows[degree*m+bit]|=((v>>bit)&1)<<j
            v=multiply(v,a)
    expected=[sum(int(x)<<j for j,x in enumerate(row)) for row in public['parity_check']]
    actual=reduced(rows,n);target=reduced(expected,n)
    if actual!=target or len(target)!=n-k:raise ValueError('recovered Goppa code differs from public code')
    return dict(status='passed',standard_library_only=True,recovery_modules_imported=False,
                secret_reference_read=False,full_support_count=n,distinct_finite_support=True,
                goppa_degree=t,goppa_irreducible=True,binary_parity_rank=len(actual),
                exact_original_public_code_equality=True,public_sha256=hashlib.sha256(raw).hexdigest(),
                key_sha256=hashlib.sha256(key_path.read_bytes()).hexdigest())


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('public',type=Path);ap.add_argument('key',type=Path)
    ap.add_argument('--report',type=Path)
    args=ap.parse_args();report=verify(args.public,args.key)
    data=json.dumps(report,indent=2)+'\n'
    if args.report:args.report.write_text(data)
    print(data)


if __name__=='__main__':main()
