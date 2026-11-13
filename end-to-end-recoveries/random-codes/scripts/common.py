"""Shared paths, cohort validation, and result writing."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent.parent
GENERIC = REPOSITORY / "toy-notebooks"
WILD = ROOT.parent / "wild-goppa"


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=int) + "\n")
    temporary.replace(path)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def cohort():
    data = json.loads((ROOT / "cohort.json").read_text())
    labels = [p["label"] for p in data["sets"]]
    if len(set(labels)) != 9:
        raise ValueError("Expected nine distinct parameter sets")
    for p in data["sets"]:
        expected = list(range(1, 21 if p["label"] == "wild-b" else 101))
        if p["seeds"] != expected:
            raise ValueError("The recorded cohort's seed range changed")
    if sum(len(p["seeds"]) for p in data["sets"]) != 820:
        raise ValueError("Expected 820 prescribed trials")
    return data


def parameters(label):
    return next(p for p in cohort()["sets"] if p["label"] == label)
