"""Reproduce selected cases from the fixed 820-trial cohort."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from common import ROOT, REPOSITORY, GENERIC, WILD, cohort, digest, write


def seed_list(text):
    values = []
    for part in text.split(","):
        if "-" in part:
            low, high = map(int, part.split("-", 1))
            if low > high:
                raise ValueError("Reversed seed range")
            values.extend(range(low, high + 1))
        else:
            values.append(int(part))
    if not values or len(set(values)) != len(values):
        raise ValueError("Seed selection must be nonempty and contain no duplicates")
    return sorted(values)


def invoke(command, log, env, timeout):
    with log.open("w") as stream:
        process = subprocess.Popen(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT,
                                   env=env, cwd=ROOT, start_new_session=True)
        try:
            return process.wait(timeout=timeout), False
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            return process.returncode, True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sets", help="Comma-separated labels; default: every set")
    parser.add_argument("--seeds", help="Subset of recorded seeds, e.g. 1 or 1-10")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--sage", default="sage")
    parser.add_argument("--wild-threads", type=int, default=1)
    parser.add_argument("--matrix-gib", type=float, default=16)
    parser.add_argument("--timeout", type=float, default=1800, help="Seconds per recovery or Schur-check process")
    parser.add_argument("--keep-kernels", action="store_true")
    parser.add_argument("--skip-schur", action="store_true", help="Omit the additional GRS shortening checks")
    parser.add_argument("--list", action="store_true", help="Print the selected cohort without running")
    args = parser.parse_args()
    if min(args.jobs, args.wild_threads, args.matrix_gib, args.timeout) <= 0:
        parser.error("Resource limits must be positive")
    data = cohort()
    known = {p["label"]: p for p in data["sets"]}
    labels = args.sets.split(",") if args.sets else list(known)
    if len(set(labels)) != len(labels) or any(label not in known for label in labels):
        parser.error("Unknown or duplicate parameter label; use --list")
    seeds = seed_list(args.seeds) if args.seeds else None
    selected = []
    for label in labels:
        p = known[label]
        chosen = p["seeds"] if seeds is None else seeds
        if not set(chosen).issubset(p["seeds"]):
            parser.error(f"Seed selection is outside the recorded cohort for {label}")
        selected.extend((p, seed) for seed in chosen)
        print(f"{label}: {len(chosen)} seeds, {chosen[0]}–{chosen[-1]}, [{p['n']},{p['k']}]", flush=True)
    print(f"Total: {len(selected)} trials", flush=True)
    if args.list:
        return
    if args.output is None:
        parser.error("--output is required for a run")
    sage = shutil.which(args.sage) if any(p["family"] != "wild-goppa" for p, _ in selected) else None
    if any(p["family"] != "wild-goppa" for p, _ in selected) and not sage:
        parser.error("SageMath was not found; pass --sage /path/to/sage")
    if any(p["family"] == "wild-goppa" for p, _ in selected):
        for name in ("wild-goppa", "holdout"):
            if not (WILD / "build" / name).is_file():
                parser.error("Build end-to-end-recoveries/wild-goppa first")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source_files = [p for p in (ROOT / "scripts").glob("*.py")]
    source_files += [GENERIC / "reference.py", WILD / "src/wild_goppa.cpp",
                     WILD / "src/holdout.cpp", WILD / "src/series.cpp", WILD / "verify.py"]
    source_files += list((WILD / "include").glob("*.hpp"))
    manifest = dict(cohort_sha256=digest(ROOT / "cohort.json"), requested=len(selected),
                    cases=[dict(label=p["label"], seed=seed) for p, seed in selected],
                    sources={str(p.relative_to(REPOSITORY)): digest(p) for p in source_files},
                    source_scope="Existing arithmetic implementations plus portable cohort runners",
                    skip_schur=args.skip_schur, kernel_retention=args.keep_kernels,
                    workers=args.jobs, timeout_seconds=args.timeout,
                    failure_policy="Keep failures and timeouts; no seed replacement or retry")
    write(out / "manifest.json", manifest)

    def run_case(p, seed):
        case = out / "cases" / p["label"] / f"seed-{seed:06d}"
        case.parent.mkdir(parents=True, exist_ok=True)
        log = out / "logs" / f"{p['label']}-{seed:06d}.log"
        log.parent.mkdir(exist_ok=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", DOT_SAGE=str(out / ".sage" / p["label"] / str(seed)))
        if p["family"] == "wild-goppa":
            command = [sys.executable, "-B", ROOT / "scripts/wild_worker.py",
                       "--threads", args.wild_threads, "--matrix-gib", args.matrix_gib]
            if args.keep_kernels:
                command.append("--keep-kernel")
        else:
            command = [sage, "-python", ROOT / "scripts/generic_worker.py"]
        command += ["--label", p["label"], "--seed", seed, "--case", case]
        started = time.monotonic()
        code, timeout = invoke(command, log, env, args.timeout)
        record = json.loads((case / "case.json").read_text()) if (case / "case.json").exists() else {}
        result = dict(label=p["label"], seed=seed, returncode=code, timed_out=timeout,
                      status="passed" if code == 0 and record.get("status") == "passed" else "failed")
        if result["status"] == "passed" and p["family"] == "grs" and not args.skip_schur:
            command = [sage, "-python", ROOT / "scripts/check_shortening.py", case,
                       "--orderings", p["shortening_orderings"], "--output", case / "schur-shortening.json"]
            code, timeout = invoke(command, log.with_suffix(".schur.log"), env, args.timeout)
            schur = json.loads((case / "schur-shortening.json").read_text()) if (case / "schur-shortening.json").exists() else {}
            result["schur"] = dict(returncode=code, timed_out=timeout,
                                   all_match_random=schur.get("all_match_random"),
                                   profiles_checked=schur.get("profiles_checked"))
            if code or timeout or not schur.get("all_match_random"):
                result.update(status="failed", failed_stage="schur-shortening")
        result["seconds"] = time.monotonic() - started
        write(case / "run-status.json", result)
        return result

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_case, p, seed): (p, seed) for p, seed in selected}
        for future in as_completed(futures):
            p, seed = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = dict(label=p["label"], seed=seed, status="failed", error=str(error))
            results.append(result)
            ordered = sorted(results, key=lambda x: (x["label"], x["seed"]))
            write(out / "summary.json", dict(requested=len(selected), completed=len(results),
                  passed=sum(r["status"] == "passed" for r in results),
                  failed=sum(r["status"] != "passed" for r in results),
                  pending=len(selected) - len(results), results=ordered))
            print(json.dumps(result), flush=True)
    if any(r["status"] != "passed" for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
