"""One fresh GRS/alternant recovery using the paper's generators and panel rules."""
import argparse
import json
from pathlib import Path
import random
import sys
import time
import traceback

from common import GENERIC, parameters, write


def main():
    if not __debug__:
        raise RuntimeError("Do not run the experiments with Python -O")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()
    P = parameters(args.label)
    if P["family"] not in ("grs", "alternant") or args.seed not in P["seeds"]:
        parser.error("Choose a GRS/alternant case from cohort.json")
    case = args.case.resolve()
    case.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    record = dict(label=args.label, seed=args.seed, parameters=P, status="running",
                  recovery_inputs="public matrix and public parameters",
                  supplied_secret_tangents=False, supplied_secret_jets=False, stages=[])
    phase = "imports"

    def begin(name):
        nonlocal phase
        phase = name
        record["current_stage"] = name
        write(case / "progress.json", record)
        return time.monotonic()

    def end(tick):
        record["stages"].append(dict(name=phase, seconds=time.monotonic() - tick))

    try:
        tick = begin("imports")
        from sage.all import GF
        from sage.version import version as sage_version
        sys.path.insert(0, str(GENERIC))
        import reference as r
        from verify_generic import verify
        record["sage_version"] = sage_version
        end(tick)
        tick = begin("generate")
        q, m, n, k, D, d, s, delta = [P[key] for key in ("q", "m", "n", "k", "D", "d", "s", "delta")]
        if P["family"] == "grs":
            Y, audit_secret = r.generate_grs(q, n, k, D, args.seed)
        else:
            Y, audit_secret = r.generate_alternant(q, m, n, P["r"], args.seed)
        E = GF(q ** m, "a")
        characteristic = int(E.characteristic())
        assert Y.dimensions() == (k, n)
        position = next(j for j in range(n) if j not in Y.pivots())

        def enc(value):
            value = E(value)
            return int(value) if E.degree() == 1 else sum(
                int(c) * characteristic ** i for i, c in enumerate(value.polynomial().list()))

        public = dict(label=args.label, seed=args.seed, q=q, m=m, n=n, k=k, D=D, d=d, s=s,
                      position=position, characteristic=characteristic,
                      modulus=[int(c) for c in E.modulus().list()],
                      matrix=[[enc(a) for a in row] for row in Y])
        write(case / "public.json", public)
        end(tick)
        tick = begin("holdout")
        assert s * (n - 1) > d * D
        K_all, monomials, holdout = r.holdout_kernel(Y, position, d, s)
        assert K_all.nrows() > 0 and holdout["heldout_values_zero"]
        if P["family"] == "grs" and K_all.nrows() > P["panel_size"]:
            rng, elements = random.Random(0), list(E)
            while True:
                mixing = r.mat(E, P["panel_size"], K_all.nrows(),
                               [rng.choice(elements) for _ in range(P["panel_size"] * K_all.nrows())])
                if mixing.rank() == P["panel_size"]:
                    break
            K = mixing * K_all
        elif P["panel_size"] is None:
            K = K_all
        else:
            K = K_all[:min(P["panel_size"], K_all.nrows())]
        assert K.rank() == K.nrows()
        holdout.update(polynomials_used=int(K.nrows()), panel_selection=P["panel_selection"])
        record["holdout_report"] = holdout
        end(tick)
        tick = begin("branches")
        J, Q, pairs = r.local_maps(K, monomials, Y.column(position))
        branches, stationary, report = r.separate_branches(J, Q, pairs, Y.column(position), E, m)
        record["branch_report"] = report
        end(tick)
        tick = begin("jet")
        depth = delta * D - 1
        jets, report = r.continue_jet(K, monomials, J, Q, pairs, stationary,
                                     Y.column(position), branches[0], depth)
        if P["family"] == "grs":
            assert report["first_full_span_order"] is not None and report["first_full_span_order"] < depth
            assert (K_all * r.series_products(jets, monomials, depth)).is_zero()
            report["full_kernel_substitution"] = True
        record["jet_report"] = report
        end(tick)
        keys, boot = {}, None
        # The historical alternant notebooks ran minor before direct; the larger
        # GRS runner ran direct first. Preserve those orders and panel choices.
        route_order = P["routes"] if P["family"] == "alternant" else list(reversed(P["routes"]))
        for route in route_order:
            tick = begin(route)
            if route == "minor":
                boot = r.bootstrap(Y, E, jets, D, form_degree=2)
                key = r.minor_finish(Y, E, boot, D, expected_dimension=D)
                record["bootstrap"] = boot["report"]
            else:
                key = r.direct_finish(Y, E, jets, position, D, form_degree=delta)
            record[route] = key["report"]
            record[route + "_public_check"] = r.audit_key(Y, key, D)
            assert r.reject_tampered_key(Y, key, D)
            keys[route] = key
            write(case / (route + "-key.json"), dict(support=[enc(a) for a in key["support"]],
                  multiplier=[enc(a) for a in key["multiplier"]], report=key["report"]))
            end(tick)
        tick = begin("independent-public-verification")
        write(case / "case.json", record)
        independent = verify(case)
        assert independent["status"] == "passed" and len(independent["routes"]) == len(P["routes"])
        record["independent_verification"] = independent
        end(tick)
        tick = begin("geometry-audit")
        record["geometry_audit"] = r.audit_geometry(
            Y, audit_secret, K, monomials, position, jets, q, bootstrap_data=boot, branches=branches)
        for route, key in keys.items():
            record[route + "_secret_check"] = r.equivalent_support(audit_secret["support"], key["support"], E, q)
        if len(keys) == 2:
            record["route_agreement"] = r.equivalent_support(keys["minor"]["support"], keys["direct"]["support"], E, q)
        end(tick)
        record["status"] = "passed"
    except Exception as error:
        record.update(status="failed", failed_stage=phase, error=f"{type(error).__name__}: {error}",
                      traceback=traceback.format_exc())
        traceback.print_exc()
    record.pop("current_stage", None)
    record["total_wall_seconds"] = time.monotonic() - started
    write(case / "case.json", record)
    write(case / "progress.json", record)
    print(json.dumps(dict(label=args.label, seed=args.seed, status=record["status"])), flush=True)
    return 0 if record["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
