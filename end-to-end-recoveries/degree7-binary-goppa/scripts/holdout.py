"""Compute a fresh degree-seven kernel with the archived CUDA/CADO algorithms."""
import argparse
import json
from pathlib import Path
import shutil
import time

from artifact import ROOT, TRIALS, environment, execute, sha, write
from validate_inputs import validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trial", choices=TRIALS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lingen", type=Path, required=True,
                        help="Pinned, patched CADO lingen_b64 executable")
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--ram-gb", type=int, default=64)
    args = parser.parse_args()
    if args.threads < 1 or args.ram_gb < 1:
        parser.error("Thread and memory limits must be positive")
    worker = ROOT / "holdout/build/square-worker"
    verifier = ROOT / "holdout/build/holdout-cpu"
    lingen = args.lingen.resolve()
    for path in (worker, verifier, lingen):
        if not path.is_file():
            raise FileNotFoundError(path)
    inputs = ROOT / "trials" / args.trial / "inputs"
    validate(inputs)
    config = json.loads((inputs / "experiment.json").read_text())
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = environment(out, args.threads)
    started = time.monotonic()
    state, parts = "-", []
    try:
        for first in range(0, config["sequence_terms"], 1024):
            last = min(first + 1024, config["sequence_terms"])
            dest = out / "segments" / f"{first:07d}-{last:07d}"
            dest.parent.mkdir(exist_ok=True)
            execute([worker, "segment", inputs / "request.txt", state, last - first,
                     dest, config["fft_batch_words"], config["solver_seed"],
                     config["compressor_seed"]], out / "logs" / f"sequence-{first:07d}.log", env)
            result = json.loads((dest / "result.json").read_text())
            if (result["start_step"], result["end_step"], result["columns"], result["rows"]) != (
                    first, last, config["N"], config["J"]):
                raise ValueError("GPU sequence metadata mismatch")
            part = dest / "sequence.raw.bin"
            if part.stat().st_size != (last - first) * config["block_bits"] ** 2 // 8:
                raise ValueError("GPU sequence length mismatch")
            parts.append(part)
            state = dest / "state.bin"
        sequence = out / "sequence.raw.bin"
        with sequence.open("wb") as stream:
            for part in parts:
                with part.open("rb") as source:
                    shutil.copyfileobj(source, stream, 4 << 20)
        generator = out / "square.gen"
        execute([lingen, "prime=2", "m=512", "n=512", f"afile={sequence}",
                 f"input_length={config['sequence_terms']}", f"ffile={generator}",
                 "rhs=none", f"thr={args.threads}x1", f"max_ram={args.ram_gb}",
                 "tuning_thresholds=recursive:100000,notiming:1"],
                out / "logs/generator.log", env, cwd=out)
        if not generator.is_file() or generator.stat().st_size == 0:
            raise ValueError("CADO returned no generator")
        dest = out / "reconstruction"
        execute([worker, "recover", inputs / "request.txt", generator, dest,
                 config["fft_batch_words"], config["solver_seed"],
                 config["compressor_seed"], 1, config["requested_samples"]],
                out / "logs/reconstruction.log", env)
        kernel = dest / "kernel.thk1"
        execute([verifier, "verify", inputs / "operator.txt", kernel,
                 out / "verification.json", "--threads", args.threads],
                out / "logs/verification.log", env)
        report = json.loads((out / "verification.json").read_text())
        if not (report["instance_id"] == config["instance_id"] and
                report["all_original_orders_verified"] and
                report["independent_rank"] == report["verified_elements"] == 256 and
                report["heldout_nonzero_values"] == 0):
            raise ValueError("Fresh kernel verification failed")
        write(out / "summary.json", dict(status="kernel_verified", trial=args.trial,
              kernel="reconstruction/kernel.thk1", kernel_sha256=sha(kernel),
              solver_seed=config["solver_seed"], compressor_seed=config["compressor_seed"],
              public_only=True, seconds=time.monotonic() - started))
        print(f"{args.trial}: fresh kernel verified; run scripts/recover.py next")
    except Exception as error:
        write(out / "summary.json", dict(status="failed", trial=args.trial,
              error=str(error), seconds=time.monotonic() - started))
        raise


if __name__ == "__main__":
    main()
