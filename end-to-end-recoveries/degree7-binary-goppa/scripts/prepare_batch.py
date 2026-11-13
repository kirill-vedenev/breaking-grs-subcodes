"""Generate the four precommitted keys; never select or replace by outcomes."""
from pathlib import Path
import argparse
import shutil
import hashlib
import json
from algebra import Field, retained_orders
from constraint_profile import profile
from generate_instance import generate, encoded

ROOT = Path(__file__).resolve().parents[1]


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def certificate(public, audit):
    """Generation-secret diagnostic only; its result never selects a key."""
    F = Field(6, 67)
    tab = [[F.mul(a, b) for b in range(64)] for a in range(64)]
    inverse = [0] + [F.inv(a) for a in range(1, 64)]
    def trim(p):
        p = list(p)
        while p and not p[-1]: p.pop()
        return p
    def add(a, b):
        p = list(a) + [0] * max(0, len(b)-len(a))
        for i, c in enumerate(b): p[i] ^= c
        return trim(p)
    def mul(a, b):
        if not a or not b: return []
        p = [0] * (len(a)+len(b)-1)
        for i, c in enumerate(a):
            for j, e in enumerate(b): p[i+j] ^= tab[c][e]
        return trim(p)
    def divrem(a, b):
        a = trim(a); q = [0] * max(0, len(a)-len(b)+1)
        while len(a) >= len(b):
            off = len(a)-len(b); c = tab[a[-1]][inverse[b[-1]]]; q[off] = c
            for j, e in enumerate(b): a[off+j] ^= tab[c][e]
            a = trim(a)
        return trim(q), a
    def derivative(p): return trim([p[i] if i & 1 else 0 for i in range(1, len(p))])
    def rank(rows):
        rows = [list(r) for r in rows]; used = 0
        for j in range(len(rows[0])):
            at = next((i for i in range(used, len(rows)) if rows[i][j]), None)
            if at is None: continue
            rows[used], rows[at] = rows[at], rows[used]
            lookup = tab[inverse[rows[used][j]]]; rows[used] = [lookup[c] for c in rows[used]]
            for i in range(len(rows)):
                if i != used and rows[i][j]:
                    lookup = tab[rows[i][j]]
                    rows[i] = [a ^ lookup[b] for a, b in zip(rows[i], rows[used])]
            used += 1
        return used
    alpha = [audit['support'][i] for i in public['original_indices']]
    g2 = mul(audit['g'], audit['g']); Pi = [1]
    for a in alpha: Pi = mul(Pi, [a, 1])
    partials = [divrem(Pi, [a, 1])[0] for a in alpha]
    curve = []
    for row in public['generator']:
        numerator = []
        for bit, p in zip(row, partials):
            if bit: numerator = add(numerator, p)
        q, r = divrem(numerator, g2)
        assert not r and len(q) <= 46
        curve.append(q)
    assert max(map(len, curve)) == 46
    for p in curve:
        assert add(mul(Pi, derivative(p)), mul(derivative(Pi), p)) == mul(g2, mul(p, p))
    leading = [p[45] if len(p) > 45 else 0 for p in curve]
    infinity = [(public['infinity_mask'] >> i) & 1 for i in range(30)]
    assert rank([leading, infinity]) == 1
    factor = [1]
    for a, s in zip(alpha, public['multiplicities']):
        for _ in range(s): factor = mul(factor, [a, 1])
    degree = 7 * 45 - public['infinity_multiplicity'] - (len(factor)-1)
    assert degree == 3
    product = mul(Pi, factor)
    mods = [divrem(derivative([0]*j + product), g2)[1] for j in range(degree+1)]
    constraint_rank = rank([[p[i] if i < len(p) else 0 for p in mods] for i in range(8)])
    return dict(status='certified' if constraint_rank == 4 else 'not_certified',
                curve_vanishing_certified=constraint_rank == 4, forced_zeros=312,
                composition_degree=315, value_residual_degree=3, value_residual_coefficients=4,
                hidden_constraint_rank=constraint_rank, residual_nullity=4-constraint_rank,
                all_curve_differential_identities_verified=True,
                generation_secret_used_only_for_this_diagnostic=True,
                key_selected_or_replaced_by_this_result=False,
                method='g^2 divides (Pi*F(C))prime; exact GF(64) polynomial and linear arithmetic')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Fresh generation directory')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(ROOT/'PLAN.json', output/'PLAN.json')
    plan = json.loads((ROOT/'PLAN.json').read_text())
    audit_root = output/'generation-audits'
    audit_root.mkdir()
    trials = []
    for number, seed in enumerate(plan['key_seeds'], 1):
        label = f'trial-{number:02d}'; dest = output/'trials'/label/'inputs'
        dest.mkdir(parents=True, exist_ok=False)
        try:
            base, audit = generate(seed)
            public, config = json.loads(base['public.json']), json.loads(base['experiment.json'])
            original_id = public.pop('instance_id')
            eligible = [p for p in public['points'] if p['original'] >= 0 and p['mask'].bit_count() > 1]
            lowered = [p['original'] for p in eligible[:12]]
            assert len(lowered) == 12
            for point in public['points']:
                if point['original'] in lowered: point['multiplicity'] = 5
            for i, original in enumerate(public['original_indices']):
                if original in lowered: public['multiplicities'][i] = 5
            public.update(challenge='degree7-k30-312-four/'+label, trial_id=label,
                          constraint_revision='roots312-four-v1', baseline_instance_id=original_id,
                          lowered_original_positions=lowered, key_seed=seed)
            public['instance_id'] = digest(encoded(public))
            counts = profile(public)
            assert counts['forced_zeros'] == 312 and counts['N'] == 2035800
            assert counts['J'] == 2136372
            cert = certificate(public, audit)
            config.update(counts)
            config.update(schema='degree7-roots312-four-trial-v1', trial_id=label,
                          instance_id=public['instance_id'], constraint_revision='roots312-four-v1',
                          baseline_instance_id=original_id, lowered_original_positions=lowered,
                          lowered_position_rule=plan['profile_rule'], ordinary_vanishing_condition_holds=False,
                          curve_vanishing_certified=cert['curve_vanishing_certified'],
                          vanishing_criterion='Recorded diagnostic; an uncertified trial is retained, not replaced.',
                          vanishing_preparation_record='vanishing-record.json',
                          requested_samples=plan['sample_size'], block_bits=plan['block_bits'],
                          solver_seed=plan['solver_seed'], compressor_seed=plan['compressor_seed'],
                          expected_gradient_rank=23, expected_quadratic_rank=15, expected_continuation_rank=28,
                          multiplicity_histogram={'5':12,'6':12},
                          multiplicity_histogram_scope='24 explicit non-unit constraints including infinity',
                          retained_orders_by_multiplicity={'5':[3,4],'6':[3,5]},
                          geometric_seed_selection=False, no_failed_seed_replacement=True,
                          parameter_choice_uses_prior_recovery_diagnosis=True,
                          runtime_recovered_key_input=False)
            config.pop('retained_orders', None)
            config['csr_bytes'] = 4*counts['nonzeros'] + 8*(counts['J']+1)
            config['three_blocks_and_range_bytes'] = (3*counts['N']+counts['J'])*64
            config['fft_length'] = 1 << (counts['J']-2).bit_length()
            config['fft_explicit_buffer_bytes'] = 128*(config['fft_length']+2)*8
            op = ['TII_HOLDOUT_V1',f'instance_id {public["instance_id"]}','k 30','degree 7','minimum 256',
                  f'heldout_original {public["heldout_original"]}',f'heldout_mask {public["heldout_mask"]}',
                  f'points {len(public["points"])}']
            op += [f'{p["original"]} {p["mask"]} {p["multiplicity"]}' for p in public['points']]
            active = [p for p in public['points'] if p['mask'].bit_count() >= 8-p['multiplicity']]
            req = ['mceliecex-holdout-cuda-request-v1',f'public_sha256 {public["public_sha256"]}',
                   'k 30','degree 7','words 8','repetitions 1',f'seed {plan["solver_seed"]}','emit_output 0',
                   f'expected_columns {counts["N"]}',f'expected_rows {counts["J"]}',
                   f'expected_nonzeros {counts["nonzeros"]}',f'point_count {len(active)}']
            for p in active:
                levels = retained_orders(7,p['multiplicity'])
                req.append(f'point {p["mask"]} {len(levels)} '+ ' '.join(map(str,levels)))
            files = {'public.json':encoded(public),'original-public.json':base['original-public.json'],
                     'experiment.json':encoded(config),'operator.txt':('\n'.join(op)+'\n').encode(),
                     'request.txt':('\n'.join(req)+'\n').encode()}
            record = dict(schema='degree7-four-vanishing-record-v1', instance_id=public['instance_id'],
                          trial_id=label, **cert, input_sha256={name:digest(raw) for name,raw in files.items()},
                          recovered_key_or_curve_in_runtime_inputs=False, certificate_failure_filters_trial=False)
            files['vanishing-record.json'] = encoded(record)
            for name, raw in files.items(): (dest/name).write_bytes(raw)
            (audit_root/(label+'.json')).write_bytes(encoded(audit))
            item = dict(trial_id=label,key_seed=seed,status='prepared',
                        instance_id=public['instance_id'], original_public_sha256=public['public_sha256'],
                        inputs=str(dest.relative_to(output)),
                        input_sha256={name:digest(raw) for name,raw in files.items()},
                        heldout_original=public['heldout_original'],heldout_shortened=public['heldout'][0],
                        lowered_original_positions=lowered,curve_vanishing_certified=cert['curve_vanishing_certified'],
                        hidden_constraint_rank=cert['hidden_constraint_rank'],nonzeros=counts['nonzeros'])
        except Exception as error:
            item = dict(trial_id=label,key_seed=seed,status='generation_or_preparation_failed',error=str(error),replaced=False)
        trials.append(item)
        print(json.dumps(item),flush=True)
    manifest = dict(schema='degree7-four-public-input-manifest-v1',plan_sha256=digest((ROOT/'PLAN.json').read_bytes()),
                    trials=trials, no_outcome_based_replacement=True, private_audits_excluded_from_transport=True)
    (output/'TRIALS.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__': main()
