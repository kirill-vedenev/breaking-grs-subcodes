"""Separate public-only holdout diagnostic for failed baseline seed 89."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

root = Path(__file__).resolve().parents[1]
diagnostics = root / "diagnostics"
case = root / "sweep-100/cases/seed-000089"
original = case / "instance/public.txt"
words = original.read_text().split()
header = list(map(int, words[1:9]))
columns = list(map(int, words[9:]))
old_held = header[5]
new_held = next(i for i in range(old_held + 1, len(columns)) if columns[i].bit_count() != 1)
header[5] = new_held
public = diagnostics / "seed89-next-heldout-public.txt"
public.write_text("BGPK1\n" + " ".join(map(str, header)) + "\n" + " ".join(map(str, columns)) + "\n")
output = diagnostics / "seed89-next-heldout-recovery"
record = dict(seed=89, separate_from_baseline_sweep=True, automatic_retry_in_baseline=False,
              original_heldout=old_held, alternate_heldout=new_held,
              selection="next non-unit public column; selected without generated secret data",
              public_generator_unchanged=True,
              original_public_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
              alternate_public_sha256=hashlib.sha256(public.read_bytes()).hexdigest(), stages=[])
started = time.monotonic()


def stage(name, command):
    begin = time.monotonic()
    with (diagnostics / f"seed89-alternate-{name}.stdout.log").open("w") as stdout, \
         (diagnostics / f"seed89-alternate-{name}.stderr.log").open("w") as stderr:
        result = subprocess.run(list(map(str, command)), cwd=root, stdout=stdout, stderr=stderr)
    record["stages"].append(dict(name=name, command=list(map(str, command)),
                                 returncode=result.returncode, seconds=time.monotonic()-begin))
    return result.returncode == 0


ok = stage("recovery", [root / "run.sh", public, output, "both"])
if ok:
    for route in ("minor", "direct"):
        ok = stage(route + "-independent", [sys.executable, root / "verify.py", public,
                   output / "recovery" / f"{route}-key.json"]) and ok
        ok = stage(route + "-audit", [root / "build/toy-goppa", "audit", public,
                   case / "instance/audit/secret.txt", output / "recovery" / f"{route}-key.txt"]) and ok
record.update(status="passed" if ok else "failed", total_seconds=time.monotonic()-started)
if (output / "recovery/recovery.json").exists():
    record["recovery"] = json.loads((output / "recovery/recovery.json").read_text())
if (output / "holdout-report.json").exists():
    record["holdout"] = json.loads((output / "holdout-report.json").read_text())
(diagnostics / "seed89-alternate-result.json").write_text(json.dumps(record, indent=2)+"\n")
print(json.dumps(record), flush=True)
sys.exit(0 if ok else 1)
