# End-to-end sweeps: code and seeds

Code and seeds for the paper's 600 GRS/alternant and 220 wild Goppa end-to-end
trials. [cohort.json](cohort.json) lists all 820 cases with their seeds and
settings. Historical kernels, recovered keys, and run
outputs are omitted.

The runners use `../../toy-notebooks` for GRS and alternant arithmetic and
`../wild-goppa` for wild Goppa recovery. Keep these directories in the
repository layout. [PROVENANCE.json](PROVENANCE.json) identifies those sources and the
new run scripts.

## Fixed cohort

Here q is the code's base field size, m the extension degree, and D the
ambient polynomial degree bound. Both routes are run in characteristic two;
only direct recovery is run in odd characteristic.

| Label | Family | q | m | [n,k] | D | (d,s) | Seeds |
|---|---|---:|---:|---|---:|---|---|
| grs-q256-k16 | GRS subcode | 256 | 1 | [48,16] | 27 | (3,2) | 1–100 |
| grs-q256-k20 | GRS subcode | 256 | 1 | [72,20] | 42 | (3,2) | 1–100 |
| grs-q257-k16 | GRS subcode | 257 | 1 | [48,16] | 27 | (3,2) | 1–100 |
| grs-q257-k20 | GRS subcode | 257 | 1 | [72,20] | 42 | (3,2) | 1–100 |
| 02-alternant-char2 | Alternant | 8 | 2 | [60,20] | 39 | (3,2) | 1–100 |
| 04-alternant-odd | Alternant | 11 | 2 | [45,15] | 29 | (3,2) | 1–100 |
| wild-a | Wild Goppa, t=4 | 4 | 3 | [63,27] | 46 | (4,3) | 1–100 |
| wild-b | Wild Goppa, t=3 | 4 | 3 | [55,28] | 42 | (5,4) | 1–20 |
| wild-c | Wild Goppa, t=4 | 4 | 3 | [64,28] | 47 | (4,3) | 1–100 |

The 820 trials exclude the seed-zero pilots, the earlier small GRS examples
over GF(32) and GF(19), and the k=12 controls. The binary Goppa 100-seed sweep
is already in the separate binary Goppa collection and brings the paper's
random-code count to 920.

## Requirements

The runner needs Python 3.10 or later and a POSIX system. GRS and alternant
cases use SageMath 10.7 and the unchanged [reference implementation](../../toy-notebooks/reference.py). Use the same
Sage version to reproduce generation and finite-field ordering.

Wild Goppa cases use the C++17/M4RI programs already in
[wild-goppa](../wild-goppa/README.md). Build them there
before running those cases:

```sh
make -C ../wild-goppa
```

If M4RI is outside the system paths, use that collection's `M4RI_PREFIX`
option. The large `wild-b` cases need an OpenMP-enabled M4RI build. On Linux
with GCC and M4RI development headers installed:

```sh
(cd ../wild-goppa && ./build-openmp.sh)
make -B -C ../wild-goppa \
  M4RI_OMP="$(cd ../wild-goppa/deps/m4ri-omp && pwd)"
```

The helper pins the M4RI source revision. Each `wild-b` packed matrix uses
about 3.93 GiB before solver workspace; set the number of concurrent jobs to
fit available RAM. The two smaller wild sets use about 0.10 and 0.13 GiB.

## Run

Run commands from this directory. To list the cases without running them:

```sh
python3 -B scripts/run.py --list
```

Run one GRS case, including its Schur-square checks:

```sh
python3 -B scripts/run.py --sets grs-q256-k16 --seeds 1 --output runs/one-grs
```

Run all six GRS/alternant sets:

```sh
python3 -B scripts/run.py \
  --sets grs-q256-k16,grs-q256-k20,grs-q257-k16,grs-q257-k20,02-alternant-char2,04-alternant-odd \
  --jobs 4 --output runs/generic-600
```

Run the 220 wild Goppa cases after building the native programs:

```sh
python3 -B scripts/run.py --sets wild-a,wild-b,wild-c \
  --jobs 2 --wild-threads 8 --output runs/wild-220
```

Omit `--sets` to run all 820 cases. `--sage /path/to/sage` selects the Sage
launcher. `--seeds 1-5` selects a subset of the recorded seeds; the runner
rejects seeds outside the chosen cohort. `--timeout` sets the time limit in
seconds for each recovery or Schur-check process (default 1800).

Use a fresh output directory for each invocation. The runner records every
case, including failures and timeouts, without replacing seeds or retrying.
`summary.json` reports completed, passing, failed, and pending cases.
`manifest.json` records the exact cohort selection and source hashes.
Stage logs, public matrices, and newly recovered keys stay under the chosen
output directory, which should be outside the distributed source files.

The generic worker keeps its hold-out polynomials in memory. The wild worker
uses a temporary kernel file and deletes it at the end of the case; pass
`--keep-kernels` to retain those newly computed files. No stored experimental
kernels or keys are inputs to either worker.

## Algorithm settings

The listed seed is passed to the reference generator for GRS and alternant
codes, or to the native generator for wild Goppa codes. Each run holds out the
first position outside the public information set. The original generation
checks remain in place; recovery outcomes do not filter the keys.

For the larger GRS codes, the sample is 48 independent linear combinations
of the full hold-out kernel, using Python `random.Random(0)` for the mixing
matrix and retrying it until it has full row rank. If the kernel has at most
48 rows, all are used. This public sampling rule does not inspect the secret
curve. The characteristic-two alternant set uses the first 48 kernel-basis
rows; the odd alternant set uses the full basis. Each wild run exports 256
independent polynomials using the unchanged native sampler.

The generic worker recovers a jet of order `2D−1`; both finishings use
quadratic forms. The GRS worker also checks that the recovered jet satisfies
the whole hold-out kernel, beyond the selected sample. Independent public
verification checks containment in the recovered GRS code and exact
subfield-code equality for alternant cases. Separate post-recovery audits
check the curve geometry and support equivalence. The generated secret is
used only for generation and those audits, never as an input to recovery.

The wild worker runs generation, hold-out computation, both finishings,
native public verification, independent Python public verification, and
separate support-equivalence audits. It recovers a radical whose third and
fourth powers both reproduce the public GF(4) code.

## Schur-square and shortening checks

The four GRS sets also run the original public-data checker on each generated
code and its dual. For each side it uses eight shuffled information-set
orderings, with Python random seed `9417 + instance_seed`, and checks every
nonzero shortening dimension. It compares square ranks with
`min(length, dimension*(dimension+1)/2)`. If the sampled products fall short of this bound, the
checker computes the rank using every product. The information sets themselves
are sampled; the check does not cover all coordinate subsets.

The result is `schur-shortening.json` in each GRS case directory. To rerun
that check alone:

```sh
sage -python scripts/check_shortening.py \
  runs/one-grs/cases/grs-q256-k16/seed-000001 \
  --orderings 8 --output runs/one-grs-shortening.json
```

`--skip-schur` skips this additional check when only recovery is needed.
The independent generic key checker can also be run on an existing batch:

```sh
sage -python scripts/verify_generic.py runs/generic-600
```

## Sources

`check_shortening.py` and the independent generic verifier are the original
sources recovered from the experiment records. The run scripts were rewritten
for this repository using the recorded settings.

[PROVENANCE.json](PROVENANCE.json) records the source hashes and these changes.
`SHA256SUMS` covers this collection's code, seed manifest, and
instructions. Check it with `shasum -a 256 -c SHA256SUMS`.
