"""One fresh wild Goppa recovery using the existing native implementation."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

from common import WILD, parameters, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--matrix-gib", type=float, default=16)
    parser.add_argument("--keep-kernel", action="store_true")
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args()
    P = parameters(args.label)
    if P["family"] != "wild-goppa" or args.seed not in P["seeds"]:
        parser.error("Choose a wild Goppa case from cohort.json")
    if args.threads < 1 or args.matrix_gib <= 0:
        parser.error("Thread and matrix limits must be positive")
    program, holdout = WILD / "build/wild-goppa", WILD / "build/holdout"
    for path in (program, holdout):
        if not path.is_file():
            raise FileNotFoundError(f"Build end-to-end-recoveries/wild-goppa first: {path}")
    case = args.case.resolve()
    case.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE", OMP_MAX_ACTIVE_LEVELS="1")
    record = dict(label=args.label, seed=args.seed, parameters=P, status="running",
                  recovery_inputs="public.txt and a freshly computed kernel", stages=[])
    started, phase = time.monotonic(), "generate"

    def stage(name, command):
        nonlocal phase
        phase = name
        record["current_stage"] = name
        write(case / "progress.json", record)
        tick = time.monotonic()
        with (case / (name + ".log")).open("w") as log:
            if name == "holdout":
                with (case / "holdout.json").open("w") as output:
                    result = subprocess.run(list(map(str, command)), env=env, stdout=output, stderr=log)
            else:
                result = subprocess.run(list(map(str, command)), env=env, stdout=log, stderr=subprocess.STDOUT)
        record["stages"].append(dict(name=name, seconds=time.monotonic() - tick,
                                     returncode=result.returncode))
        if result.returncode:
            raise RuntimeError(f"{name} exited with {result.returncode}; see {name}.log")

    try:
        stage("generate", [program, "generate", P["t"], P["n"], P["d"], args.seed, case / "instance"])
        if args.generate_only:
            record["status"] = "generated"
        else:
            public, kernel = case / "instance/public.txt", case / "kernel.bin"
            stage("holdout", [holdout, public, kernel, 256, args.matrix_gib])
            record["holdout"] = json.loads((case / "holdout.json").read_text())
            stage("recover", [program, "recover", public, kernel, case / "recovery", "both"])
            record["recovery"] = json.loads((case / "recovery/recovery.json").read_text())
            spec = importlib.util.spec_from_file_location("wild_public_verifier", WILD / "verify.py")
            verifier = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(verifier)
            checks = {}
            for route in ("minor", "direct"):
                stage(route + "-verify", [program, "verify", public, case / "recovery" / (route + "-key.txt")])
                checks[route] = verifier.verify(public, case / "recovery" / (route + "-key.json"))
                stage(route + "-audit", [program, "audit", public, case / "instance/audit/secret.txt",
                                        case / "recovery" / (route + "-key.txt")])
            record.update(status="passed", independent_verification=checks)
    except Exception as error:
        record.update(status="failed", failed_stage=phase, error=f"{type(error).__name__}: {error}",
                      traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        if not args.keep_kernel and (case / "kernel.bin").exists():
            (case / "kernel.bin").unlink()
    record.pop("current_stage", None)
    record["total_wall_seconds"] = time.monotonic() - started
    write(case / "case.json", record)
    write(case / "progress.json", record)
    print(json.dumps(dict(label=args.label, seed=args.seed, status=record["status"])), flush=True)
    return 0 if record["status"] in ("passed", "generated") else 1


if __name__ == "__main__":
    raise SystemExit(main())
