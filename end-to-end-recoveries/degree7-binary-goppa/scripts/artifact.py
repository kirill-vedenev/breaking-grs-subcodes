"""Small file and process helpers for the portable artifact commands."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
TRIALS = tuple(f"trial-{i:02d}" for i in range(1, 5))


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            value.update(block)
    return value.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def execute(command, log, env=None, cwd=None):
    command = list(map(str, command))
    log = Path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    print(f"Running {log.stem}", flush=True)
    with log.open("w") as stream:
        process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                 env=env, cwd=ROOT if cwd is None else cwd)
    record = dict(command=command, returncode=process.returncode,
                  seconds=time.monotonic() - start)
    write(log.with_suffix(".execution.json"), record)
    if process.returncode:
        raise RuntimeError(f"Command failed ({process.returncode}); see {log}")
    return record


def environment(output, threads=None):
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE="1", DOT_SAGE=str(Path(output) / ".sage"))
    if threads is not None:
        env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
                   OMP_MAX_ACTIVE_LEVELS="1")
    return env
