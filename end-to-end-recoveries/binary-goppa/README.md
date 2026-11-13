# Toy binary Goppa key recovery

C++ implementation for unshortened binary Goppa parameters
`m=6, t=6, n=64, k=28`, over GF(64) with modulus `x^6+x+1`.
The ambient degree bound is `D=n-2t-1=51`.

The program computes a public hold-out kernel, separates the six Frobenius
branches, continues a binary-normalized jet, and recovers an equivalent full key.

| Finishing method | Jet prefix | Output |
|---|---:|---|
| Quadratic bootstrap, aligned tangents, minor code, Sidelnikov–Shestakov | 39 | Full GRS `[64,45]` minor code and complete key |
| Direct extraction through a square-root quotient | 101 | Complete key from quadratic sections |

Each method returns 64 distinct finite support elements, a compatible GRS
multiplier, and a monic irreducible degree-six Goppa polynomial. The recovered
alternant and Goppa codes are checked for exact equality with the public binary
code. No shortening or support extension is needed.

## Build

Build with a C++17 compiler and M4RI. Extension-field arithmetic uses C++
lookup tables. The recovery executable runs without Sage or Python.

If M4RI is available through `pkg-config`:

```sh
make
make check
```

Otherwise supply the installation prefix containing `include/m4ri` and the library:

```sh
make M4RI_PREFIX=/path/to/m4ri
make M4RI_PREFIX=/path/to/m4ri check
```

Compiler and dependency details and source hashes are recorded in
[results/build-info.json](results/build-info.json).

The packed hold-out matrix occupies about 1.32 GiB, with additional storage needed
by the solver. `TOY_MATRIX_GIB` sets the matrix-allocation cap, which defaults to
3 GiB; it does not bound total process memory.

## Run

Use a fresh output directory:

```sh
./run.sh instances/a/public.txt runs/demo both
```

The final argument is `minor`, `direct`, or `both`. This command computes a new
hold-out kernel and verifies the recovered key. Minor-only mode computes a
normalized jet through order 40 and uses its order-39 prefix. Direct mode computes
through order 102 and uses only the order-101 prefix for its sections.

To run recovery from a saved kernel:

```sh
./build/toy-goppa recover instances/a/public.txt runs/demo/kernel.bin runs/replay both
```

Kernel independence and every retained Hasse condition are checked before the
polynomials are used.

## Fresh-key end-to-end sweeps

On Linux with GNU `/usr/bin/time`, run a consecutive seed batch with:

```sh
python3 sweep.py --seeds 1-100 --workers 24 --output runs/sweep-100
```

Each seed generates a fresh public code, computes and checks a fresh hold-out
kernel, recovers both keys, checks their public geometric certificates, and runs
the independent Python verifier on each key. Separate secret-key audits run
only after public recovery and verification. The recovery algorithms receive
only the public matrix and the newly computed kernel.

The seed list is fixed before execution. Failed generation, rank checks,
verification, and timeouts remain in the denominator. Failed seeds are neither
replaced nor retried automatically. `summary.json` counts all completed and
pending cases; `cases/seed-*/case.json` retains every stage, command, exit status,
wall time, peak resident memory, and verification report. Raw kernels, logs,
public inputs, generated audit keys, and recovered keys remain beside the case
records. The manifest records source and binary hashes, hardware, and library
version. Sweep timings are measured with concurrent jobs and are not isolated
runtimes.

`--workers` limits simultaneous cases and every subprocess uses one OpenMP
thread. Each hold-out matrix alone is about 1.32 GiB; leave additional memory for
the solver and other stages. `--timeout` defaults to 3,600 seconds per stage,
and `--matrix-gib` defaults to a 3 GiB matrix allocation cap. These sweeps vary
the seed while keeping the parameters above fixed.

To resume the same specification, add `--resume`. Completed cases, including
failures, are kept. The runner stops if it finds an interrupted case or changed
source or binary hashes, leaving the existing records intact.

After the batch completes, reverify every recovered key and check that every
requested seed has a saved outcome:

```sh
python3 verify_sweep.py runs/sweep-100 --report runs/sweep-100/validation.json
```

This uses the independent public verifier and checks input, key, and available
kernel hashes. Its `status` describes the integrity of the evidence;
`recovery_passed` and `recovery_failed` give the actual experimental outcomes.

## Input instances

[instances/a](instances/a) and [instances/b](instances/b) contain fixed test
inputs with the same parameters. Each has a public generator in `public.txt` and
an audit key in `audit/secret.txt`. Recovery reads only the public file. The
included seed and generator settings for instance B are in its `generation.json`.

Generate another key with an integer seed:

```sh
./build/toy-goppa generate 1 runs/new-key
./run.sh runs/new-key/public.txt runs/new-key-recovery both
```

Generation uses `std::mt19937_64`, an explicit support shuffle, and a degree-six
irreducibility test. It checks the public dimension and the binary square-free
Goppa identity. Failed hypotheses stop the run.

## Verify a saved key

The C++ verifier uses the public input alone:

```sh
./build/toy-goppa verify instances/a/public.txt results/a/recovery/direct-key.txt
```

The independent verifier uses Python's standard library and imports none of the
recovery implementation:

```sh
python3 verify.py instances/a/public.txt results/a/recovery/direct-key.json
python3 verify.py instances/a/public.txt results/a/recovery/minor-key.json
```

An optional audit compares the recovered key with the generated secret:

```sh
./build/toy-goppa audit instances/a/public.txt instances/a/audit/secret.txt results/a/recovery/direct-key.txt
```

It checks that one global Frobenius power and one Möbius transformation relate
the supports. Secret data enter only generation and this separate audit.

## Exact checks

1. Hold-out: 98,280 squarefree degree-five monomials and multiplicity four
   at retained points. Orders 1 and 3 suffice for construction. All orders 0
   through 3 are verified, with unit-column conditions automatic.
2. Branches: gradient rank 21, nullity 7, and quadratic restriction rank 15.
   All six branches have continuation rank 26.
3. Jets: even coefficients are squares of their half-order coefficients.
   Odd coefficients satisfy the regular and next-even-order equations. A fresh
   full substitution verifies the completed series against every kernel polynomial.
4. Minor finishing: 406 quadratic unknowns, constraint rank 97, and 309
   quadrics; tangent rank 26 at all 64 positions; 378 minor rows spanning the
   full GRS `[64,45]` code.
5. Direct finishing: quadratic sections with contact orders 102 and 100;
   square roots of their public evaluation ratios recover the support.
6. Certificates: exact public-code equality, a reconstructed curve of degree
   at most 51, coefficient-level jet coherence, all branch planes, and global
   tangent alignment. The two recovered supports agree up to a Möbius map.

The initial ordinary root count is 252, below the degree bound 255. Verified
coherence through order 102 supplies 103 additional roots at the held-out point,
so `252+103=355>255`. Together with the Hasse checks and full substitution, this
certifies the returned curve identities using the recovered public representation.

## Results

| Instance | Kernel dimension | Hold-out and verification | Recovery and verification | Outcome |
|---|---:|---:|---:|---|
| A | 552 | 186.6 s | 13.9 s | Both methods pass |
| B | 414 | 194.8 s | 13.1 s | Both methods pass |

Recovery timings include the kernel recheck, branch recovery, continuation,
both finishing methods, and public geometric certificates. Optional secret audits
are excluded. The two individual toy runs provide no statistical success estimate.

### 100-seed sweep

The consecutive seeds 1–100 were run from generation through a fresh hold-out
kernel, both finishing routes, independent public verification, and separate
secret audits, at the fixed parameters above. No seeds were discarded or
replaced.

| Check | Outcome |
|---|---:|
| Generated distinct public inputs | 100 / 100 |
| Fresh kernels passing every retained Hasse condition | 100 / 100 |
| Full recovery through both routes on the initial hold-out | 99 / 100 |
| Independently verified keys from those successful cases | 198 / 198 |
| Failed seed recovered with a different public hold-out | 1 / 1 |

Seed 89 fails the quadratic branch-separation criterion at zero-based position
27: its gradient rank is 21, but the restricted quadratic rank is 9 instead of
15. An independent binary row reduction confirms both ranks. A separate
follow-up changes only the held-out column to position 29, the next non-unit
public column. A newly computed kernel then has quadratic rank 15, and both
routes recover independently verified keys. The original sweep succeeded on
99/100 keys. The separate follow-up recovered the remaining key
with one additional hold-out computation.

The initial kernels have dimensions 372–706. Every successful initial recovery
has six branches, continuation rank 26, bootstrap constraint rank 97, minor
rank 45, and full jet coefficient span first reached at order 33. All returned
keys pass exact alternant and Goppa code equality, coefficient-level coherence,
global tangent alignment, and Frobenius/Möbius audits.

The initial sweep took 1,707.36 seconds (28.46 minutes) with 24 simultaneous
single-threaded cases on an AMD EPYC-Genoa VM with 96 logical CPUs and 188 GiB
RAM. Median hold-out time was 338.59 seconds; median successful recovery time was
11.35 seconds. Maximum recorded resident memory for one stage was 1,417,884 KiB
(1.35 GiB). These wall times were measured under concurrent load. The separate
seed-89 follow-up took 317.99 seconds and is excluded from the initial sweep's
time and success rate.

The [summary](results/sweep-100/summary.json),
[per-seed table](results/sweep-100/cases.csv), and
[independent batch audit](results/sweep-100/raw-validation.json) record all requested
seeds, including the failure. The evidence directory contains
all 100 public inputs, every recovered JSON key, detailed case reports, and the
[seed-89 diagnostics](results/sweep-100/diagnostics/seed89-alternate-result.json).
The full kernels and intermediate matrices are retained in the experiment
archives; all 100 initial kernel hashes were checked after transfer.
File paths in the saved records are relative to the directory in which the
experiment was run, and time stamps were removed from them.

Recheck the compact evidence with:

```sh
python3 verify_sweep.py results/sweep-100
python3 results/sweep-100/diagnostics/check_ranks.py results/sweep-100/diagnostics
python3 verify.py results/sweep-100/diagnostics/seed89-next-heldout-public.txt results/sweep-100/diagnostics/seed89-minor-key.json
python3 verify.py results/sweep-100/diagnostics/seed89-next-heldout-public.txt results/sweep-100/diagnostics/seed89-direct-key.json
```

A subsequent [local investigation](results/sweep-100/diagnostics/position-investigation/conclusions.json)
first checked all 36 hold-outs outside the fixed information set, then changed
the information set to check the remaining 28 positions. Positions 27 and 42
have quadratic rank 9; the other 62 positions have rank 15. A new systematic basis at
position 27 reproduces rank 9. In coordinates along the six actual Frobenius
tangents, the missing terms are exactly the six neighboring products
`z_s z_(s+1)`, with indices modulo six. All nine non-neighboring products
remain. Thus the quadratic zero set in the projective tangent quotient is a
hexagon of six lines, containing 384 GF(64) directions rather than the six
true directions. The six lost terms form one Frobenius orbit and cannot be treated as
independent rank defects.

All 484 kernel polynomials are used at position 27. Restricting to its
463-dimensional zero-gradient subspace gives the same quadratic span.
Independent Python recomputations match every gradient and quadratic
coefficient at positions 27, 29, 42, and the rebased position 27. A separate
Python calculation verifies all 62,390,020 retained Hasse conditions and
independence of the 484 original kernel vectors; PLUQ elimination independently
confirms the original hold-out matrix's nullity 484.

On every false line `v = a tau_s + b tau_(s+1)`, the third-order equations
with the binary normalization `v_2 = v_1^2`, after eliminating the unknown
third coefficient using zero-gradient polynomials, span `a^2 b = a b^2 = 0`.
They remove all 378 spurious directions. This diagnostic uses the verified
recovered tangent basis; it does not implement a public cubic solver.
The usual continuation matrix also drops to rank 24 at both exceptional
positions, so bypassing the quadratic rank check alone would not repair the
existing recovery method.

The rank loss belongs to the local evaluation of the degree-five hold-out
interpolation space, which depends on the omitted position. The recovered
curve has a unique unramified preimage at both positions 27 and 29, and its
measured Frobenius jet-span dimensions through orders 1, 3, and 5 are 7, 13,
and 19 at all 64 positions. These checks describe the local rank defect. A
general closed-form criterion in the Goppa polynomial and support value for
predicting exceptional positions is still missing.

The tangents do not descend to a proper subfield: after
quotienting by the radial direction and normalizing projectively, every tangent
requires GF(64), including positions 27 and 42. None is defined over GF(2),
GF(4), or GF(8). The
[subfield report](results/sweep-100/diagnostics/position-investigation/tangent-subfields.json)
includes explicit Frobenius tests at all 64 positions.

The defect depends on the interpolation restrictions. Using the verified curve,
the full 3,497-dimensional binary space of global cubic identities gives
quadratic rank 15 at every position. Separately, omitting both 27 and 42 restores
rank 15 at 27, but so does omitting 27 and ordinary position 29. In the latter
tests, only genuine curve identities are retained: eight series coefficients
at position 27, together with the 62 retained multiplicity-four points, give
`248+8=256>255`. This diagnostic filters with the previously recovered key,
so it cannot serve as a public recovery procedure. These controls exclude a
proper-subfield explanation and a necessary pairwise interaction between the
two exceptional positions. A closed-form criterion predicting their exact
locations from the Goppa polynomial remains unresolved.

### Joint tangent and third-coefficient solving on the original kernels

The [joint-solver experiment](results/follow-ups/joint-seed89/summary.json) recovers
seed 89 at both exceptional positions, 27 and 42, using their original uniform
multiplicity-four kernels. Position 29 is a normal control. Recovery uses the
existing kernels and constraint profiles, without a known key or curve. The original quadratic-only executable is retained alongside
this experimental implementation.

Write `h_r` for the degree-r Taylor part of a kernel polynomial at the held-out
point, and `B_h` for the polar form of `h_2`. Under binary normalization,
`v_2 = v_1^2` componentwise. The order-three equation is

```text
J_h v_3 + C_h(v_1) = 0,
C_h(v) = h_3(v) + B_h(v, v^2).
```

Taking the zero-gradient linear combinations of the same kernel polynomials
eliminates `v_3` and gives public cubic equations `C_h(v_1)=0`. Combined with
the original quadrics, their degree-three ideal has rank 50 among the 56 cubic
monomials in the six tangent-quotient coordinates. Its six-dimensional quotient
allows an eigenvalue computation to recover all six branches using only public
data. The cubic equations are constructed from the actual order-three series
coefficients; all returned directions are checked by fresh substitution.

Orders three and four alone still give rank 24 in `v_3` at the exceptional
positions. Apart from the two usual scalar/parameter freedoms, two unwanted
degrees of freedom remain. For zero-gradient polynomials the order-five
equation contains no `v_5`; its part depending on `v_3` is
`D C_h(v_1)[v_3]`. These rows raise the continuation rank to 26. More generally,
at every odd order `r >= 3`, the zero-gradient equations at order `r+2` supply
the same additional linear matrix in `v_r`. The next even coefficient is
already determined by binary normalization. The same augmented linear system
therefore determines every subsequent odd coefficient.

| Heldout | Original quadratic rank | Continuation through next order | Continuation with extra stationary order | Full recovery |
|---:|---:|---:|---:|---|
| 27 | 9 | 24 | 26 | Both routes pass |
| 42 | 9 | 24 | 26 | Both routes pass |
| 29 | 15 | 26 | 26 | Both routes pass |

All three runs reach jet order 102 and pass fresh full substitution against
every original kernel polynomial, binary normalization, public code equality,
jet coherence, and global geometry checks. The order-103 check is a stationary
lookahead check; the saved complete jet is through order 102. All six recovered
keys pass the independent Python verifier, and the kernel hashes match the
previous experiments. The recovery stages take about 14 seconds per case on
the measured server, excluding the already completed kernel acquisition.

[Validation](results/follow-ups/joint-seed89/validation.json) and the experimental
[public joint solver](results/follow-ups/joint-seed89/public_joint.hpp) are saved
with the results. The projected order-five conditions resolve the failure that
remains after joint solving through order four.

[results/](results) contains complete recovered keys, intermediate matrices,
and verification records. Minor-only and direct-only modes on instance A return
the same keys as combined mode. Mutated kernel and multiplier inputs are rejected.

See [FORMATS.md](FORMATS.md) for all encodings.
