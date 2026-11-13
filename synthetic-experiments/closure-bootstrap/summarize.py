"""Summarizes the JSON-lines records of closure_chain.py.

Usage: python3 summarize.py RESULTS_DIR

Reads keys.jsonl (or keys.jsonl.gz) and writes summary.md and summary.json to
the same directory. Uses only the Python standard library.
"""

import collections
import gzip
import json
import os
import statistics
import sys


def load(path):
    if os.path.exists(path):
        f = open(path)
    elif os.path.exists(path + ".gz"):
        f = gzip.open(path + ".gz", "rt")
    else:
        return []
    with f:
        return [json.loads(line) for line in f if line.strip()]


def closures_label(r):
    if r["family"] != "goppa":
        return "2D+1-n"
    if r["mode"] == "proved":
        return "n-4t-1"
    return {0: "n-5t-1", 2: "n-5t+1"}.get(r.get("extra_closures", 0),
                                         "n-5t-1+%d" % r.get("extra_closures", 0))


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "results"
    rows = load(os.path.join(d, "keys.jsonl"))
    groups = collections.OrderedDict()
    order = {"goppa": 0, "alternant": 1}
    for r in sorted(rows, key=lambda r: (order[r.get("family", "goppa")], r["m"], r["n"], r["t"],
                                         r["mode"], -r.get("extra_closures", 0), r["seed"])):
        r.setdefault("family", "goppa")
        key = (r["family"], r["m"], r["n"], r["t"], closures_label(r))
        groups.setdefault(key, []).append(r)
    out = []
    lines = ["| Family | `(m,n,t)` | `k` | Closures | Positions | Keys | Recovered | Binary quadrics | "
             "Other outcomes | Audits passed | Median seconds |",
             "|---|---|---:|---|---:|---:|---:|---|---|---:|---:|"]
    for key, rs in groups.items():
        seeds = [r["seed"] for r in rs]
        assert len(seeds) == len(set(seeds)), key
        ok = [r for r in rs if "error" not in r]
        status = collections.Counter(r.get("status", "error") for r in rs)
        recovered = status.pop("recovered", 0)
        quadrics = collections.Counter(r.get("binary_quadrics") for r in ok)
        audits = sum(1 for r in ok if r.get("audit_support_moebius")
                     and len(r.get("audit_tangent_conjugates", {})) == 1)
        row = dict(family=key[0], m=key[1], n=key[2], t=key[3], closures=key[4],
                   k=sorted(set(r["k"] for r in ok)),
                   positions=sorted(set(r["closure_positions"] for r in ok)),
                   keys=len(rs), recovered=recovered,
                   binary_quadrics={str(a): b for a, b in sorted(quadrics.items())},
                   predicted_binary_quadrics=sorted(set(r["predicted_binary_quadrics"] for r in ok
                                                        if r.get("predicted_binary_quadrics") is not None)),
                   other_outcomes=dict(status),
                   first_position_passed=sum(1 for r in ok if "jet_position" in r
                                             and len(r.get("separation_attempts", [])) == 1),
                   inconsistent_continuation_orders=sum(r.get("continuation_inconsistent_orders") or 0
                                                        for r in ok),
                   gradient_rank_failures=sum(len(r.get("gradient_rank_failures", [])) for r in ok),
                   audits_passed=audits,
                   median_seconds=round(statistics.median(r["total_seconds"] for r in ok), 1))
        out.append(row)
        lines.append("| %s | (%d,%d,%d) | %s | `%s` | %s | %d | %d | %s | %s | %d | %s |" % (
            "Binary Goppa" if key[0] == "goppa" else "Binary alternant", key[1], key[2], key[3],
            ",".join(map(str, row["k"])), key[4], ",".join(map(str, row["positions"])),
            row["keys"], recovered,
            ", ".join("%s (%d keys)" % (a, b) for a, b in row["binary_quadrics"].items()),
            ", ".join("%s: %d" % (a, b) for a, b in row["other_outcomes"].items()),
            audits, row["median_seconds"]))
    total = sum(r["keys"] for r in out)
    text = ["# Summary", "", "%d keys in %d parameter sets and closure counts." % (total, len(out)), "",
            "\"Recovered\" means that the support was recovered and the public code embeds in",
            "`GRS_{D+1}` on it. \"Audits passed\" counts keys whose recovered tangent planes all lie",
            "on one Frobenius conjugate and whose recovered support is a Möbius image of a",
            "conjugate of the secret support. Times are for one process and include the audits.", ""]
    text += lines
    with open(os.path.join(d, "summary.md"), "w") as f:
        f.write("\n".join(text) + "\n")
    with open(os.path.join(d, "summary.json"), "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
