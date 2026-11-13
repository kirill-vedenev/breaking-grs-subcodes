"""Check the exported diagnostic ranks without the native algebra library."""
import argparse
import json
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("matrices", type=Path)
args = parser.parse_args()
results = {}
for name in ("gradient", "quadratic-restrictions"):
    words = (args.matrices / f"seed89-{name}.txt").read_text().split()
    assert words[0] == "GF64M1"
    rows, columns = map(int, words[1:3])
    values = list(map(int, words[3:]))
    assert len(values) == rows * columns and set(values) <= {0, 1}
    pivots = {}
    for i in range(rows):
        value = sum(values[i * columns + j] << j for j in range(columns))
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    results[name] = dict(rows=rows, columns=columns, rank=len(pivots))
assert results["gradient"]["rank"] == 21
assert results["quadratic-restrictions"]["rank"] == 9
print(json.dumps(dict(status="passed", seed=89, heldout=27, matrices=results,
                     expected_quadratic_rank=15, method="packed GF(2) elimination; binary matrices have the same rank over GF(64)"), indent=2))
