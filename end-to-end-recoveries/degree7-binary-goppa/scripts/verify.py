"""Check the distributed files, all four saved kernels, and all eight full keys."""
import argparse
import json
from pathlib import Path

from artifact import ROOT, TRIALS, sha, write
from materialize import verify_storage
from verify_key import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    checked = 0
    for line in (ROOT / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(None, 1)
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or sha(path) != digest:
            raise ValueError(f"File checksum mismatch: {name}")
        checked += 1
    results = []
    for trial in TRIALS:
        meta = verify_storage(trial)
        public = ROOT / "trials" / trial / "inputs/original-public.json"
        keys = {}
        for route in ("minor", "direct"):
            candidate = ROOT / "results" / trial / "recovery" / route / "full-key/equivalent-key.json"
            keys[route] = verify(public, candidate)
        results.append(dict(trial=trial, kernel_sha256=meta["sha256"],
                            polynomials=meta["polynomials"], keys=keys))
        print(f"{trial}: kernel bytes and both full keys verified", flush=True)
    report = dict(status="passed", checksum_entries=checked, kernels=4, full_keys=8,
                  kernel_membership_recomputed=False, recovery_rerun=False,
                  trials=results)
    if args.report:
        write(args.report, report)
    print(json.dumps({k: v for k, v in report.items() if k != "trials"}, indent=2))


if __name__ == "__main__":
    main()
