#!/usr/bin/env python3
"""Independently recheck all public keys from a completed sweep.

The verifier reads public inputs, recovered keys, and experiment records. It
does not read generated secret keys or import the native recovery program.
Raw kernels are hashed when present; their algebraic checks are recorded by
the separate native stages rather than repeated here.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from verify import verify


def check(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_batch(output):
    started = time.monotonic()
    manifest = json.loads((output / "manifest.json").read_text())
    summary = json.loads((output / "summary.json").read_text())
    seeds = manifest["seeds"]
    check(seeds and len(seeds) == len(set(seeds)), "invalid manifest seed list")
    expected_dirs = {f"seed-{seed:06d}" for seed in seeds}
    actual_dirs = {path.name for path in (output / "cases").iterdir() if path.is_dir()}
    check(actual_dirs == expected_dirs, "case directories differ from the requested seed set")
    passed, failed, verified_keys, hashed_kernels = 0, 0, 0, 0
    failures, verification, dimensions, public_hashes = [], [], [], []
    for seed in seeds:
        case = output / "cases" / f"seed-{seed:06d}"
        record = json.loads((case / "case.json").read_text())
        check(record["seed"] == seed, f"seed {seed}: case identity mismatch")
        check(record["status"] in ("passed", "failed"), f"seed {seed}: unfinished case")
        if "holdout" in record:
            dimensions.append(str(record["holdout"]["kernel_dimension"]))
        if "public_sha256" in record:
            check(sha256(case / "instance/public.txt") == record["public_sha256"],
                  f"seed {seed}: public hash mismatch")
            public_hashes.append(record["public_sha256"])
        kernel = case / "kernel.bin"
        if kernel.exists() and "kernel_sha256" in record:
            check(sha256(kernel) == record["kernel_sha256"], f"seed {seed}: kernel hash mismatch")
            hashed_kernels += 1
        if record["status"] == "failed":
            failed += 1
            failures.append(dict(seed=seed, stage=record["failed_stage"], error=record["error"]))
            continue
        check(all(stage["returncode"] == 0 and not stage["timed_out"] for stage in record["stages"]),
              f"seed {seed}: a passing record contains an unsuccessful stage")
        expected_stages = ["generate", "prepare", "holdout", "recover", "minor-verify",
                           "minor-independent", "minor-audit", "direct-verify",
                           "direct-independent", "direct-audit"]
        check([stage["name"] for stage in record["stages"]] == expected_stages,
              f"seed {seed}: missing, duplicated, or reordered stage")
        public = case / "instance/public.txt"
        for route in ("minor", "direct"):
            result = verify(public, case / "recovery" / f"{route}-key.json")
            original = record["independent_verification"][route]
            check(result["public_sha256"] == original["public_sha256"],
                  f"seed {seed}: independent public hash mismatch")
            check(result["recovered_key_sha256"] == original["recovered_key_sha256"],
                  f"seed {seed}: recovered key hash mismatch")
            verification.append(dict(seed=seed, route=route, **result))
            verified_keys += 1
        passed += 1
    check(summary["requested"] == len(seeds), "requested count mismatch")
    check(summary["completed"] == len(seeds), "incomplete sweep")
    check(summary["pending"] == 0 and summary["all_requested_completed"], "pending cases remain")
    check(summary["passed"] == passed and summary["failed"] == failed, "outcome counts mismatch")
    check(summary["failures"] == sorted(failures, key=lambda r: r["seed"]), "failure list mismatch")
    check(summary["kernel_dimensions"] == dict(Counter(dimensions)), "kernel histogram mismatch")
    expected_records = [dict(seed=seed, status=json.loads((output / "cases" / f"seed-{seed:06d}" /
                        "case.json").read_text())["status"], record=f"cases/seed-{seed:06d}/case.json")
                        for seed in sorted(seeds)]
    check(summary["records"] == expected_records, "summary record index mismatch")
    return dict(status="passed", requested=len(seeds), recovery_passed=passed, recovery_failed=failed,
                independent_keys_reverified=verified_keys, kernel_hashes_checked=hashed_kernels,
                unique_public_inputs=len(set(public_hashes)), duplicate_public_inputs=len(public_hashes)-len(set(public_hashes)),
                public_only=True, generated_secrets_read=False,
                independent_of_native_recovery=True, failures=failures,
                seconds=time.monotonic() - started, verification=verification)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = verify_batch(args.output)
    if args.report:
        args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "verification"}, indent=2))


if __name__ == "__main__":
    main()
