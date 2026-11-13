#!/usr/bin/env python3
"""Generate one fixed public Goppa instance; keep the secret outside the package."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
from random import Random
from algebra import Field, rref, kernel, restrict, columns, retained_orders, choose

ROOT = Path(__file__).resolve().parents[1]
SEED = 2026101601


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True)+'\n').encode()


def generate(seed=SEED):
    rng = Random(seed)
    m, t, n, k, degree, multiplicity = 6, 4, 64, 30, 7, 6
    field = Field(m, 67)
    draws = 0
    while True:
        g = [rng.randrange(field.q) for _ in range(t)]+[1]
        draws += 1
        if field.irreducible(g):
            break
    support = list(range(n))
    rng.shuffle(support)
    beta = [field.inv(field.evaluate(g, a)) for a in support]
    parity = field.binary_parity(support, beta, t)
    if len(parity) != m*t:
        raise ValueError('generated key has unexpected rank; seed is not replaced')
    original_generator, info = rref(kernel(parity, n), n)
    removed = info[:n-m*t-k]
    keep = [j for j in range(n) if j not in removed]
    hs = [restrict(row, keep) for row in parity]
    generator, information = rref(kernel(hs, len(keep)), len(keep))
    assert len(generator) == k and len(keep) == 54
    points = columns(generator, len(keep))
    assert all(points[information[a]] == 1 << a for a in range(k))
    held = next(j for j, point in enumerate(points)
                if j not in information and point.bit_count() >= degree)
    infinity = 0
    for point in points:
        infinity ^= point
    full = dict(schema='binary-goppa-public-v1',
                parameters=dict(m=m,t=t,n=n,k=n-m*t),
                field_modulus=[(67 >> i) & 1 for i in range(m+1)],
                parity_check=[[row >> j & 1 for j in range(n)] for row in parity],
                generator=[[row >> j & 1 for j in range(n)] for row in original_generator])
    full_bytes = encoded(full)
    digest = hashlib.sha256(full_bytes).hexdigest()
    constraint_points = [dict(original=keep[j],mask=point,multiplicity=multiplicity)
                         for j, point in enumerate(points) if j != held]
    constraint_points.append(dict(original=-1,mask=infinity,multiplicity=multiplicity))
    active = [p for p in constraint_points if p['mask'].bit_count() >= degree-multiplicity+1]
    levels = retained_orders(degree,multiplicity)
    assert levels == [3,5]
    N = comb(k,degree)
    J = len(active)*sum(comb(k,u) for u in levels)
    nnz = sum(choose(p['mask'].bit_count(), degree-u)*comb(k-degree+u,u)
              for p in active for u in levels)
    assert J > N
    public = dict(schema='degree7-k30-public-v1',challenge='degree7-k30',
                  parameters=dict(m=m,t=t,n=len(keep),k=k,degree=degree),
                  original_parameters=full['parameters'],public_sha256=digest,
                  field_modulus=full['field_modulus'],original_indices=keep,
                  shortened_original_indices=removed,information_set=information,
                  information_set_original=[keep[j] for j in information],
                  heldout=[held],heldout_original=keep[held],heldout_mask=points[held],
                  minimum_elements=256,default_elements=512,points=constraint_points,
                  infinity_mask=infinity,infinity_multiplicity=multiplicity,
                  multiplicities=[0 if j==held else multiplicity for j in range(len(keep))],
                  generator=[[row >> j & 1 for j in range(len(keep))] for row in generator])
    public['instance_id'] = hashlib.sha256(encoded(public)).hexdigest()
    op = ['TII_HOLDOUT_V1',f'instance_id {public["instance_id"]}',f'k {k}',f'degree {degree}',
          'minimum 256',f'heldout_original {keep[held]}',f'heldout_mask {points[held]}',
          f'points {len(constraint_points)}']
    op += [f'{p["original"]} {p["mask"]} {multiplicity}' for p in constraint_points]
    request = ['mceliecex-holdout-cuda-request-v1',f'public_sha256 {digest}',
               f'k {k}',f'degree {degree}','words 8','repetitions 1','seed 2026101602',
               'emit_output 0',f'expected_columns {N}',f'expected_rows {J}',
               f'expected_nonzeros {nnz}',f'point_count {len(active)}']
    request += [f'point {p["mask"]} {len(levels)} '+ ' '.join(map(str,levels)) for p in active]
    depth = 2*(len(keep)-2*t-1)
    report = dict(schema='degree7-k30-preparation-v1',seed=seed,
                  key_generator='Python random.Random; monic quartic rejection sampling for irreducibility',
                  irreducible_polynomial_draws=draws,geometric_seed_selection=False,
                  original_shape=[n,n-m*t],shortened_shape=[len(keep),k],
                  shortened_original_indices=removed,heldout_original=keep[held],
                  heldout_shortened=held,heldout_selection='first non-information column of weight at least d',
                  original_public_sha256=digest,instance_id=public['instance_id'],
                  N=N,J=J,nonzeros=nnz,active_points=len(active),retained_orders=levels,
                  block_bits=512,requested_samples=256,sequence_terms=2*((N+511)//512)+256,
                  fft_length=1 << (J-2).bit_length(),fft_batch_words=2,
                  csr_bytes=4*nnz+8*(J+1),three_blocks_and_range_bytes=(3*N+J)*64,
                  fft_explicit_buffer_bytes=128*((1 << (J-2).bit_length())+2)*8,
                  memory_excludes_cufft_plan_workspace=True,
                  forced_zeros=sum(public['multiplicities'])+multiplicity,
                  composition_degree_bound=degree*(len(keep)-2*t-1),
                  ordinary_vanishing_condition_holds=sum(public['multiplicities'])+multiplicity>degree*(len(keep)-2*t-1),
                  direct_jet_depth=depth,minor_contact_order=depth-len(keep)+1,
                  solver_seed=2026101602,compressor_seed=2026101603,
                  solver='existing homogeneous CUDA/CADO Wiedemann, not the proved sampler',
                  relation_to_weis='same m=6,t=4,n=64,k_short=30,d=7,s=6 family; not his exact key',
                  kernel_solved=False,recovery_tested=False)
    audit = dict(schema='generation-audit-only',excluded_from_server_package=True,
                 seed=seed,field_modulus=67,m=m,t=t,support=support,g=g)
    return {'original-public.json':full_bytes,'public.json':encoded(public),
            'operator.txt':('\n'.join(op)+'\n').encode(),
            'request.txt':('\n'.join(request)+'\n').encode(),
            'experiment.json':encoded(report)}, audit


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=ROOT/'inputs')
    ap.add_argument('--audit-output',type=Path,default=ROOT/'runs/generated-audit')
    ap.add_argument('--check',action='store_true')
    args=ap.parse_args();files,audit=generate()
    if args.check:
        for name,data in files.items():
            if (args.output/name).read_bytes()!=data:
                raise ValueError('fixed-input reproduction mismatch: '+name)
        print(json.dumps(dict(status='passed',fixed_inputs_reproduced=True,files=len(files))))
        return
    args.output.mkdir(parents=True,exist_ok=True)
    if any((args.output/name).exists() for name in files):
        raise FileExistsError('inputs already exist; use --check to verify them')
    args.audit_output.mkdir(parents=True,exist_ok=True)
    if (args.audit_output/'secret.json').exists():
        raise FileExistsError('audit secret already exists')
    for name,data in files.items():(args.output/name).write_bytes(data)
    (args.audit_output/'secret.json').write_bytes(encoded(audit))
    print(files['experiment.json'].decode())


if __name__=='__main__':main()
