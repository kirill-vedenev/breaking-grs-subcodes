#!/usr/bin/env python3
"""Run reproducible public-input-only binary Goppa recovery across fresh seeds.

Requires the two compiled programs, Python 3, and GNU /usr/bin/time (Linux).
Every requested seed is recorded, including generation and recovery failures.
No failed seed is replaced. Resume skips only cases with a completed case.json.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parent
SOURCES = ["Makefile", "run.sh", "verify.py", "sweep.py", "src/holdout.cpp",
           "src/toy_goppa.cpp", "include/algebra.hpp", "include/series_api.hpp"]


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def capture(command):
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    return dict(command=command, returncode=result.returncode,
                stdout=result.stdout.strip(), stderr=result.stderr.strip())


def seeds_from(text):
    result = []
    for item in text.split(","):
        ends = item.split("-")
        if len(ends) == 1:
            result.append(int(item))
        elif len(ends) == 2:
            start, end = map(int, ends)
            if end < start:
                raise ValueError("descending seed range")
            result.extend(range(start, end + 1))
        else:
            raise ValueError("invalid seed range")
    if not result or min(result) < 0 or max(result) >= 2**64 or len(set(result)) != len(result):
        raise ValueError("seeds must be distinct unsigned 64-bit integers")
    return result


class StageFailure(RuntimeError):
    pass


def run_case(seed, output, timeout, matrix_gib):
    case = output / "cases" / f"seed-{seed:06d}"
    record = dict(seed=seed, status="running", started_at=utcnow(), stages=[])
    started = time.monotonic()
    stage_name = "initialize"
    case.mkdir(parents=True, exist_ok=False)
    save(case / "progress.json", record)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", LC_ALL="C")

    def stage(name, command, stdout_name=None):
        nonlocal stage_name
        stage_name = name
        record["current_stage"] = name
        save(case / "progress.json", record)
        stdout_path = case / (stdout_name or f"{name}.stdout.log")
        stderr_path = case / f"{name}.stderr.log"
        resource_path = case / f"{name}.resources.json"
        timed = ["/usr/bin/time", "-f",
                 '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M}',
                 "-o", str(resource_path), *map(str, command)]
        begin = time.monotonic()
        timed_out = False
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            process = subprocess.Popen(timed, stdout=stdout, stderr=stderr,
                                       cwd=ROOT, env=env, start_new_session=True)
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                code = process.wait()
        resources = None
        if resource_path.exists():
            for line in reversed(resource_path.read_text().splitlines()):
                try:
                    resources = json.loads(line)
                    break
                except json.JSONDecodeError:
                    continue
        entry = dict(name=name, command=list(map(str, command)), returncode=code,
                     timed_out=timed_out, seconds=time.monotonic() - begin,
                     resources=resources, stdout=stdout_path.name, stderr=stderr_path.name)
        record["stages"].append(entry)
        save(case / "progress.json", record)
        if code or timed_out:
            detail = stderr_path.read_text(errors="replace")[-4000:]
            raise StageFailure(f"{name}: {'timeout' if timed_out else f'exit {code}'}: {detail}")

    def read(relative):
        return json.loads((case / relative).read_text())

    try:
        program = ROOT / "build/toy-goppa"
        stage("generate", [program, "generate", seed, case / "instance"])
        public = case / "instance/public.txt"
        record["public_sha256"] = digest(public)
        record["generation"] = read("instance/generation.json")
        stage("prepare", [program, "prepare", public, case / "holdout.txt"])
        stage("holdout", [ROOT / "build/holdout", case / "holdout.txt",
                          case / "kernel.bin", matrix_gib], "holdout-report.json")
        record["holdout"] = read("holdout-report.json")
        record["kernel_sha256"] = digest(case / "kernel.bin")
        stage("recover", [program, "recover", public, case / "kernel.bin",
                          case / "recovery", "both"])
        record["recovery"] = read("recovery/recovery.json")
        record["independent_verification"] = {}
        record["secret_audits"] = {}
        for route in ("minor", "direct"):
            key = case / "recovery" / f"{route}-key"
            stage(f"{route}-verify", [program, "verify", public, key.with_suffix(".txt")],
                  f"{route}-verification.json")
            stage(f"{route}-independent", [sys.executable, ROOT / "verify.py", public,
                                          key.with_suffix(".json")], f"{route}-independent.json")
            record["independent_verification"][route] = read(f"{route}-independent.json")
            # Generated secret data are read only after public recovery and verification.
            stage(f"{route}-audit", [program, "audit", public, case / "instance/audit/secret.txt",
                                    key.with_suffix(".txt")], f"{route}-audit.json")
            record["secret_audits"][route] = read(f"{route}-audit.json")
        stage_name = "validate-reports"
        recovery = record["recovery"]
        checks = [record["holdout"]["full_hasse_verified"],
                  recovery["status"] == "passed", recovery["public_inputs_only"],
                  recovery["full_hasse_verified"], recovery["kernel_independence_verified"],
                  recovery["jet"]["full_substitution"], recovery["jet"]["normalized"],
                  recovery["geometry"]["coherent_jet"],
                  recovery["geometry"]["holdout_curve_identities_certified"],
                  not recovery["geometry"]["uses_original_secret"],
                  recovery["routes_mobius_equivalent"], recovery["tampered_multiplier_rejected"]]
        for route in ("minor", "direct"):
            verify = record["independent_verification"][route]
            audit = record["secret_audits"][route]
            checks.extend([recovery[route]["exact_public_goppa_equality"],
                           recovery[route]["exact_public_alternant_equality"],
                           verify["status"] == "passed", verify["exact_public_code_equality"],
                           verify["public_sha256"] == record["public_sha256"],
                           verify["recovered_key_sha256"] == digest(case / "recovery" / f"{route}-key.json"),
                           audit["status"] == "passed", audit["global_support_equivalence"]])
        if not all(checks):
            raise StageFailure("a required report check was false")
        record["status"] = "passed"
    except Exception as error:
        record.update(status="failed", failed_stage=stage_name,
                      error=f"{type(error).__name__}: {error}")
    record.update(finished_at=utcnow(), total_wall_seconds=time.monotonic() - started)
    record.pop("current_stage", None)
    save(case / "case.json", record)
    save(case / "progress.json", record)
    return record


def distribution(values):
    return dict(count=len(values), minimum=min(values), median=statistics.median(values),
                maximum=max(values), mean=statistics.mean(values)) if values else None


def summarize(output, manifest, records, started):
    ordered = sorted(records, key=lambda record: record["seed"])
    passed = [record for record in ordered if record["status"] == "passed"]
    failures = [dict(seed=r["seed"], stage=r["failed_stage"], error=r["error"])
                for r in ordered if r["status"] != "passed"]
    holds = [r["holdout"] for r in ordered if "holdout" in r]
    stage_names = sorted({s["name"] for r in ordered for s in r["stages"]})
    result = dict(schema_version=1, updated_at=utcnow(), parameters=manifest["parameters"],
                  requested=len(manifest["seeds"]), completed=len(ordered), passed=len(passed),
                  failed=len(failures), pending=len(manifest["seeds"]) - len(ordered),
                  all_requested_completed=len(ordered) == len(manifest["seeds"]),
                  failures=failures, failure_stages=dict(Counter(r["stage"] for r in failures)),
                  invocation_elapsed_seconds=time.monotonic() - started,
                  kernel_dimensions=dict(Counter(str(h["kernel_dimension"]) for h in holds)),
                  timing_scope="Wall timings under concurrent load, not isolated benchmarks",
                  case_wall_seconds=distribution([r["total_wall_seconds"] for r in ordered]),
                  stage_wall_seconds={name: distribution([s["seconds"] for r in ordered
                                     for s in r["stages"] if s["name"] == name]) for name in stage_names},
                  max_stage_rss_kib=max([s["resources"]["max_rss_kib"] for r in ordered
                                        for s in r["stages"] if s["resources"]] or [0]),
                  records=[dict(seed=r["seed"], status=r["status"],
                                record=f"cases/seed-{r['seed']:06d}/case.json") for r in ordered])
    save(output / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", default="1-100", help="comma-separated integers or inclusive ranges")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=3600, help="seconds per stage")
    parser.add_argument("--matrix-gib", type=float, default=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    seeds = seeds_from(args.seeds)
    if args.workers < 1 or args.timeout <= 0 or args.matrix_gib <= 0:
        parser.error("workers, timeout, and matrix-gib must be positive")
    for file in ["build/toy-goppa", "build/holdout"]:
        if not (ROOT / file).is_file():
            parser.error(f"build first: {file} is missing")
    if not Path("/usr/bin/time").exists() or platform.system() != "Linux":
        parser.error("stage resource measurements require GNU /usr/bin/time on Linux")
    output = args.output.resolve()
    fingerprint = {file: digest(ROOT / file) for file in [*SOURCES, "build/toy-goppa", "build/holdout"]}
    specification = dict(seeds=seeds, source_and_binary_sha256=fingerprint,
                         timeout_per_stage_seconds=args.timeout, matrix_allocation_cap_gib=args.matrix_gib,
                         parameters=dict(m=6, t=6, n=64, k=28, D=51, modulus=67, d=5, s=4))
    records = []
    if args.resume:
        manifest = json.loads((output / "manifest.json").read_text())
        if any(manifest.get(key) != value for key, value in specification.items()):
            parser.error("resume specification or source/binary fingerprint differs")
        for seed in seeds:
            case = output / "cases" / f"seed-{seed:06d}"
            if (case / "case.json").is_file():
                record = json.loads((case / "case.json").read_text())
                if record["seed"] != seed or record["status"] not in ("passed", "failed"):
                    parser.error(f"invalid completed case record: {case}")
                records.append(record)
            elif case.exists():
                parser.error(f"incomplete case needs inspection before resume: {case}")
    else:
        output.mkdir(parents=True, exist_ok=False)
        manifest = dict(specification, created_at=utcnow(), workers=args.workers,
                        host=platform.node(), platform=platform.platform(), cpu_count=os.cpu_count(),
                        python=sys.version, single_thread_per_process=True,
                        source_directory=str(ROOT), output_directory=str(output),
                        compiler=capture(["c++", "--version"]),
                        m4ri=capture(["pkg-config", "--modversion", "m4ri"]),
                        linked_libraries=capture(["ldd", str(ROOT / "build/holdout")]),
                        hardware=capture(["lscpu"]), memory=capture(["free", "-b"]),
                        failure_policy="retain all requested seeds; no substitutions or automatic retries",
                        recovery_inputs="public.txt and newly computed kernel.bin only",
                        audit_scope="generated secret read only in separate post-recovery audits")
        save(output / "manifest.json", manifest)
    started = time.monotonic()
    done = {r["seed"] for r in records}
    summarize(output, manifest, records, started)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        jobs = {executor.submit(run_case, seed, output, args.timeout, args.matrix_gib): seed
                for seed in seeds if seed not in done}
        for future in as_completed(jobs):
            record = future.result()
            records.append(record)
            summary = summarize(output, manifest, records, started)
            event = dict(at=utcnow(), seed=record["seed"], status=record["status"],
                         seconds=record["total_wall_seconds"], completed=summary["completed"],
                         passed=summary["passed"], failed=summary["failed"])
            if "error" in record:
                event["error"] = record["error"]
            with (output / "events.jsonl").open("a") as stream:
                stream.write(json.dumps(event) + "\n")
            print(json.dumps(event), flush=True)
    result = summarize(output, manifest, records, started)
    print(json.dumps({key: result[key] for key in ("requested", "completed", "passed", "failed")}), flush=True)
    return 0 if not result["failed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
