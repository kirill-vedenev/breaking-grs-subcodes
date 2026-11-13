"""Validate a precommitted public trial without generating or reading any key."""
from pathlib import Path
import argparse
import hashlib
import json
from constraint_profile import profile

ROOT = Path(__file__).resolve().parents[1]


def validate(inputs):
    if not __debug__: raise RuntimeError('Do not use Python -O')
    inputs = Path(inputs).resolve()
    registry = json.loads((ROOT/'TRIALS.json').read_text())
    assert hashlib.sha256((ROOT/'PLAN.json').read_bytes()).hexdigest() == registry['plan_sha256']
    config = json.loads((inputs/'experiment.json').read_text())
    matches = [t for t in registry['trials'] if t['trial_id'] == config['trial_id'] and t['status'] == 'prepared']
    if len(matches) != 1: raise ValueError('Trial not in the fixed cohort')
    trial = matches[0]
    for name, digest in trial['input_sha256'].items():
        if hashlib.sha256((inputs/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pinned trial input changed: '+name)
    public = json.loads((inputs/'public.json').read_text())
    assert public['instance_id'] == config['instance_id'] == trial['instance_id']
    counts = profile(public)
    for name in ['N','J','nonzeros','forced_zeros']:
        assert counts[name] == config[name]
    assert counts['forced_zeros'] == 312
    record = json.loads((inputs/'vanishing-record.json').read_text())
    assert record['curve_vanishing_certified'] == trial['curve_vanishing_certified']
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path)
    args=parser.parse_args()
    registry=json.loads((ROOT/'TRIALS.json').read_text())
    paths=[args.inputs] if args.inputs else [ROOT/t['inputs'] for t in registry['trials'] if t['status']=='prepared']
    for path in paths:
        record=validate(path)
        print(json.dumps(dict(status='passed',trial_id=record['trial_id'],curve_vanishing_certified=record['curve_vanishing_certified'])))


if __name__=='__main__': main()
