"""Replay both recovery routes from a saved or newly computed public kernel."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import time

from artifact import ROOT, TRIALS, environment, execute, sha, write
from materialize import materialize


def main():
    if not __debug__:
        raise RuntimeError("Do not run recovery with Python -O")
    from sage.all import version
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trial", choices=TRIALS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, help="Use a newly computed THK1 kernel")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    for binary in ("build/extract-jets", "build/libseries.so", "holdout/build/holdout-cpu"):
        if not (ROOT / binary).is_file():
            raise FileNotFoundError(f"Run make cpu first: {binary}")
    started = time.monotonic()
    out = materialize(args.trial, args.output, panel=True, kernel=args.kernel)
    inputs = out / "inputs"
    public = json.loads((inputs / "public.json").read_text())
    config = json.loads((inputs / "experiment.json").read_text())
    p, held = public["parameters"], public["heldout"][0]
    env = environment(out, args.threads)
    env["TII_JET_SERIES_LIBRARY"] = str(ROOT / "build/libseries.so")
    logs, common = out / "logs", out / "common"
    common.mkdir()
    shutil.copyfile(inputs / "public.json", common / "public.json")
    panel = out / "polynomials/kernel-panel.u64le"
    index = out / "polynomials/index.json"

    def call(script, *argv, name=None):
        return execute([sys.executable, "-B", ROOT / script, *argv],
                       logs / ((name or Path(script).stem) + ".log"), env=env)

    try:
        execute([ROOT / "holdout/build/holdout-cpu", "verify", inputs / "operator.txt",
                 out / "kernel.thk1", out / "kernel-verification.json",
                 "--threads", args.threads], logs / "kernel-verification.log", env=env)
        verified = json.loads((out / "kernel-verification.json").read_text())
        if not (verified["all_original_orders_verified"] and
                verified["independent_rank"] == verified["verified_elements"] == 256 and
                verified["heldout_nonzero_values"] == 0 and
                verified["instance_id"] == public["instance_id"]):
            raise ValueError("Kernel verification failed")
        call("recovery/verify_input_encoding.py", out / "kernel.thk1", panel,
             common / "encoding-verification.json")
        execute([ROOT / "build/extract-jets", panel,
                 out / "polynomials/monomials-colex.u64", p["k"], p["degree"],
                 public["heldout_mask"], 256, common], logs / "extract-jets.log", env=env)
        call("recovery/check_criterion.py", inputs / "public.json", index, common)
        call("recovery/extend_jet.py", inputs / "public.json", panel, common,
             "--depth", config["direct_jet_depth"], "--threads", args.threads)
        routes = {}
        for route in ("minor", "direct"):
            dest = out / route
            dest.mkdir()
            for name in ("public.json", f"deep-jets-{held}.json"):
                shutil.copyfile(common / name, dest / name)
            if route == "minor":
                call("recovery/bootstrap.py", dest, "--position", held,
                     "--depth", config["minor_contact_order"], name="minor-bootstrap")
                boot = json.loads((dest / "bootstrap.json").read_text())
                if boot["minor"]["rank"] != 41 or not boot["recovered_subfield_code_equals_public_code"]:
                    raise ValueError("Minor finishing failed")
                call("recovery/prepare_retained.py", dest, name="minor-prepare")
                support = dest / "retained-support.json"
            else:
                call("recovery/extract_support.py", inputs / "public.json", dest,
                     name="direct-extract")
                support = dest / "direct-support.json"
                if json.loads(support.read_text())["report"]["status"] != "verified":
                    raise ValueError("Direct finishing failed")
            call("recovery/finish_key.py", "--public", inputs / "public.json",
                 "--support", support, "--original", inputs / "original-public.json",
                 "--output", dest / "full-key", name=route + "-finish-key")
            call("scripts/verify_key.py", inputs / "original-public.json",
                 dest / "full-key/equivalent-key.json", "--report",
                 dest / "independent-verification.json", name=route + "-verify-key")
            routes[route] = json.loads((dest / "independent-verification.json").read_text())
        report = dict(status="passed", trial=args.trial, public_only=True,
                      generation_secret_read=False, saved_keys_or_jets_used=False,
                      kernel_sha256=sha(out / "kernel.thk1"), sage=str(version()),
                      routes=routes, seconds=time.monotonic() - started)
        write(out / "summary.json", report)
        print(f"{args.trial}: both full recovery routes passed", flush=True)
    except Exception as error:
        write(out / "summary.json", dict(status="failed", trial=args.trial,
              error=str(error), seconds=time.monotonic() - started))
        raise


if __name__ == "__main__":
    main()
