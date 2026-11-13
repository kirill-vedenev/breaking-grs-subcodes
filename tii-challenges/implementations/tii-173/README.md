# TII-173 implementation

This implementation recovered all 96 support positions and a degree-eight
Goppa polynomial. The hold-out sample contains 256 degree-five
forms in 35 variables. Finishing uses an order-296 normalized jet, direct
extraction on the 91 retained columns, and restoration of five removed columns.
The field is GF(128), defined by `z^7+z+1` (`0x83`).

Run the commands from this directory and choose fresh output directories.
Recovery was not rerun when preparing this distribution.

## Verify the saved key

Verify the saved key with Python 3, without rerunning recovery:

```sh
python3 recovery/src/independent_verify_key.py \
  ../../results/tii-173/recovered-sk.json \
  --public ../../results/tii-173/pk_McEliece_173.txt \
  --report tii173-key-verification.json
```

## Recover from the saved hold-out sample

Requirements: SageMath 10.7, NumPy, M4RI, GNU Make, and a C++17 compiler.
The native Makefile uses Sage's include/library paths; `SAGE` and
`SAGE_LOCAL` may be set explicitly for a nonstandard installation.

```sh
OUT="$PWD/runs/recovery"
python3 ../../kernels/materialize.py ../../kernels/tii-173 --output "$OUT"
make -C recovery/src libseries.so libincremental.so
export DOT_SAGE="$OUT/.sage"

sage -python recovery/src/extend_normalized.py "$OUT" --position 0 --depth 296
sage -python recovery/src/direct_jet_recovery.py "$OUT" --position 0 \
  --output "$OUT/direct-recovery.json"
sage -python recovery/src/deshorten.py "$OUT" \
  --public-key ../../results/tii-173/pk_McEliece_173.txt
python3 recovery/src/independent_verify_key.py "$OUT/equivalent-key.json" \
  --public ../../results/tii-173/pk_McEliece_173.txt \
  --report "$OUT/independent-full-key-verification.json"
sage -python recovery/src/verify_chain.py "$OUT"
```

The chain audit reconstructs the curve, checks every hold-out composition
through degree 370, verifies the seven branches and normalized continuation,
and checks coefficient-level jet coherence. The new recovered key is written
to `runs/recovery/equivalent-key.json`.

## Compute a new hold-out sample

Requirements: CUDA development toolkit, C++17 compiler, and two visible GPUs.
The saved run used two H200 NVL GPUs, CUDA 13.3, and the exact command parameters
below. The matrix occupies about 13.93 GiB in total; checkpointing requires
about 60 GB of free disk. The two-GPU elimination runs on one host.

```sh
make -C holdout cpu cuda CUDA_ARCH=90
KERNEL_OUT="$PWD/runs/holdout"
holdout/build/holdout-cuda solve ../../kernels/tii-173/operator.txt \
  "$KERNEL_OUT" --count 256 --seed 17320260930 --devices 0,1 \
  --panel-bits 8 --max-gib 64 --checkpoint-every 2048
holdout/build/holdout-cpu verify ../../kernels/tii-173/operator.txt \
  "$KERNEL_OUT/kernel-colex.bin" "$KERNEL_OUT/verification.json" --threads 16
holdout/build/holdout-cpu export "$KERNEL_OUT/kernel-colex.bin" "$KERNEL_OUT/kernel.bin"
cp ../../kernels/tii-173/public.json "$KERNEL_OUT/public.json"
```

Require the verification report to accept 256 independent forms and every
original Hasse order before using this new sample. Set `OUT` to this new
directory and start the previous section at the native-library build.
`kernel-colex.bin` is THK1; export explicitly permutes it to WHK1.

The prepared public input is included. To regenerate the shortening and
constraints from the original public key without solving for a kernel:

```sh
python3 holdout/scripts/prepare.py --output runs/prepared
```

This source is the original preparation implementation with relocated paths.
The [kernel format](../../kernels/README.md),
[saved results](../../results/tii-173/README.md), and
[source provenance](../PROVENANCE.json) specify the exact stored inputs.

## Recorded minor finishing

The minor finishing run for TII-173 used the same public matrix and coherent
jet as the direct run. At bootstrap depth 58, it recovered 489 quadrics and a full
GRS [91,66] minor code, then recovered all 91 retained support
coordinates. Exact shortened public-code equality passed. Its internal
finishing time was 2.412 s, excluding Sage startup and the shared jet.

The recovered supports agree up to a Möbius transformation, with Frobenius
power 0. Full-key extension was not rerun for the minor route. The evidence
is in [minor-finishing results](../../results/tii-173/evidence/minor-finishing).
The original bootstrap algorithm is in `recovery/minor/bootstrap.py`; the
bundled engine contains the unchanged mathematical helpers used by that run.

To compare the saved supports with SageMath:

```sh
sage -python recovery/minor/compare_finishings.py \
  --public ../../kernels/tii-173/public.json \
  --minor ../../results/tii-173/evidence/minor-finishing/bootstrap-public-data.json \
  --direct-key ../../results/tii-173/recovered-sk.json \
  --report tii173-finishing-comparison.json
```

Use a fresh report path. The expected entry is `tii-173` with
`moebius_equivalent: true` and `frobenius_powers_matching: [0]`.
The saved full key is used only by this comparison, after recovery.

To repeat minor finishing from the saved public jet, use a fresh directory:

```sh
MINOR_OUT="$PWD/runs/minor-finishing"
mkdir -p runs
mkdir "$MINOR_OUT"
cp ../../kernels/tii-173/public.json "$MINOR_OUT/public.json"
cp ../../results/tii-173/evidence/deep-jets-0.json "$MINOR_OUT/deep-jets-0.json"
export DOT_SAGE="$MINOR_OUT/.sage"
sage -python recovery/minor/bootstrap.py "$MINOR_OUT" --position 0 --depth 58
```

This reads only `public.json` and `deep-jets-0.json` and writes
`bootstrap.json` and `bootstrap-public-data.json`. The results distributed
here come from the original run; this finishing was not repeated during packaging.

## Optional audits of saved results

The saved-key check above validates the recovered Goppa key. These optional
audits check the hold-out sample and intermediate algebraic results against
the paper's claims. They read saved artifacts and write to a fresh temporary
directory without rerunning recovery.

Run from this `tii-173/` implementation directory. Requirements are Python
3.10+, a C++17 compiler, and SageMath 10.7 with NumPy for the deeper audits, plus M4RI for the native substitution helper.
Run Python/Sage without `-O`, which disables the assertions used by these audits.

```sh
TII_REPO="$(cd ../.. && pwd)"
TII_CHECKS="$(mktemp -d "${TMPDIR:-/tmp}/tii-verification.XXXXXX")"
export PYTHONDONTWRITEBYTECODE=1
unset PYTHONOPTIMIZE
export DOT_SAGE="$TII_CHECKS/.sage"
mkdir "$TII_CHECKS/build"

(cd "$TII_REPO/results" && shasum -a 256 -c SHA256SUMS)
(cd "$TII_REPO/kernels" && shasum -a 256 -c SHA256SUMS)
(cd "$TII_REPO/implementations" && shasum -a 256 -c SHA256SUMS)
```


### Saved kernel

```sh
c++ -O3 -std=c++17 -pthread \
  "$TII_REPO/implementations/tii-173/holdout/src/cpu.cpp" \
  -o "$TII_CHECKS/build/holdout-cpu"
python3 "$TII_REPO/kernels/materialize.py" "$TII_REPO/kernels/tii-173" \
  --output "$TII_CHECKS/tii173"
"$TII_CHECKS/build/holdout-cpu" verify "$TII_REPO/kernels/tii-173/operator.txt" \
  "$TII_CHECKS/tii173/kernel.thk1" "$TII_CHECKS/tii173-hasse.json" --threads 4
```

Expect `all_original_orders_verified: true` and `independent_rank: 256`.
This audit uses a different thread count from the original timed run.

Check that the archived THK1 and WHK1 encodings contain exactly the same forms:

```sh
"$TII_CHECKS/build/holdout-cpu" export "$TII_CHECKS/tii173/kernel.thk1" \
  "$TII_CHECKS/tii173-export.whk1"
cmp "$TII_CHECKS/tii173-export.whk1" "$TII_CHECKS/tii173/kernel.bin"
```

`cmp` should exit successfully without output.

### Saved jet and curve

The audit expects the key to be named `equivalent-key.json` beside
`public.json`, `kernel.bin`, and the saved jet/direct-extraction outputs.
Assemble those names in the already restored temporary directory. The native
library is built beside temporary copies of the original source because
this implementation loads the library from its source directory.

```sh
cp "$TII_REPO/results/tii-173/recovered-sk.json" \
  "$TII_CHECKS/tii173/equivalent-key.json"
cp "$TII_REPO/results/tii-173/evidence/deep-jets-0.json" \
  "$TII_REPO/results/tii-173/evidence/direct-recovery.json" \
  "$TII_CHECKS/tii173/"

mkdir "$TII_CHECKS/tii173-code"
cp "$TII_REPO/implementations/tii-173/recovery/src/"*.py \
  "$TII_REPO/implementations/tii-173/recovery/src/native.cpp" \
  "$TII_REPO/implementations/tii-173/recovery/src/Makefile" \
  "$TII_CHECKS/tii173-code/"
make -C "$TII_CHECKS/tii173-code" libseries.so
sage -python "$TII_CHECKS/tii173-code/verify_chain.py" "$TII_CHECKS/tii173"
```

Sage's M4RI headers and library must be available to the supplied Makefile;
`SAGE` and `SAGE_LOCAL` can be passed explicitly to `make` for a nonstandard
installation. This target builds only the substitution library.

The output `tii173/chain-verification.json` should certify all 256 hold-out
polynomials as exact curve identities through composition degree 370,
seven Frobenius branches, all continuation ranks equal to 33, and all
10,395 coordinate-jet coefficients coherent. It also checks support-only
extraction from the saved prefix through order 148. It does not compute a
new hold-out sample, extend the jet, or restore a new full key.
