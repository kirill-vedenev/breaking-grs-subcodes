# Experimental evidence

The toy and TII recovery runs start from public matrices and parameters;
generated secret data, where available, are used only for separate audits of
those runs. Four collections test individual steps using secret-derived data:
tangent planes for minor-code fullness, curve polynomials for the
direct-extraction existence condition, one local jet for synthetic direct
recovery, and Frobenius closures of the tangent planes for the closure
bootstrap.

## GRS subcodes and alternant codes

The [notebook collection](toy-notebooks/README.md) contains four fixed examples.
Each computes its hold-out kernel, recovers a coherent jet, and extracts a support
and multiplier.

| Code | Base field | Ambient field | Jet order | Finishing | Result |
|---|---|---|---:|---|---|
| GRS subcode `[24,6]` | GF(32) | GF(32) | 13 | Minor and direct | Pass |
| Alternant `[60,20]` | GF(8) | GF(64) | 77 | Minor and direct | Pass |
| GRS subcode `[17,4]` | GF(19) | GF(19) | 14 | Direct | Pass |
| Alternant `[45,15]` | GF(11) | GF(121) | 57 | Direct | Pass |

All examples satisfy the ordinary interpolation bound `s(n-1)>dD`. Audits check
the curve identities, recovered branches, coefficient-level jet coherence, and
global Frobenius/Möbius equivalence of the supports. Both characteristic-two
examples verify exact minor-GRS equality and agreement of the finishing routes.
Both alternant examples verify exact equality of the recovered subfield code
with the public code. Corrupted multipliers are rejected.

The [notebook result records](toy-notebooks/results/README.md) contain the
parameters, ranks, timings, and individual checks.

## TII challenge recoveries

The [experiment specification](tii-challenges/EXPERIMENTS.md) records the shortening,
systematic coordinate sets, run settings, and full pipeline for each code.
Its [hardware and timing tables](tii-challenges/EXPERIMENTS.md#5-measured-hardware-and-runtimes)
give the measured durations and distinguish complete stages from component
timers and sums across separate runs.

The [TII results](tii-challenges/results/README.md) contain full equivalent binary Goppa
keys recovered from the original public matrices and single-position hold-out
samples. No published secret challenge answers were inputs.

| Instance | Original code | Goppa degree | Saved hold-out forms | Jet order | Full support |
|---|---|---:|---:|---:|---:|
| [TII-173](tii-challenges/results/tii-173/README.md) | `[96,40]` | 8 | 256 quintics in 35 variables | 296 | 96 |
| [TII-249](tii-challenges/results/tii-249/README.md) | `[235,107]` | 16 | 280 quintics in 50 variables | 290 | 235 |
| [TII-253](tii-challenges/results/tii-253/README.md) | `[214,142]` | 9 | 280 sextics in 46 variables | 82 | 214 |

TII-173 and TII-249 were also finished through quadratic bootstrapping and
minor codes on the same saved jets. Their minor ranks are 66 and 129, and the
retained supports match the direct results up to a Möbius transformation.
These additional runs did not repeat full-key extension. TII-253 has full
minor and direct recoveries. The direct run extends the shared
public jet to order 198 and produces an independently verified 214-position
degree-nine key. The [coverage table](tii-challenges/EXPERIMENTS.md#finishing-coverage) and saved `evidence/minor-finishing/`
outputs distinguish these scopes.

Each key retains its exact extension-field polynomial and original column
order. Saved independent verifiers establish equality of the reconstructed
binary Goppa parity-check row space with the original public check space,
distinct support, support avoidance, and irreducibility of the recovered
polynomial. TII-249 and TII-253 also retain the exhaustive coordinate-extension
candidate records. The larger-instance hold-out reports certify sampled rank
and every original Hasse order; later geometry certificates are stored with
the recovery evidence.

[Kernel storage](tii-challenges/kernels/README.md) preserves the exact accepted samples.
[Source provenance and instructions](tii-challenges/implementations/README.md) identify the implementations used. The original results
are preserved alongside the additional TII-253 direct run, which reused the
hold-out sample.

## Binary Goppa codes

The [C++ collection](end-to-end-recoveries/binary-goppa/README.md) contains two fixed
unshortened instances and a 100-seed end-to-end sweep with
`(m,t,n,k)=(6,6,64,28)` and ambient degree bound `D=51`.

| Instance | Hold-out dimension | Branches | Continuation rank | Minor code | Direct extraction |
|---|---:|---:|---:|---|---|
| [A](end-to-end-recoveries/binary-goppa/results/a/recovery/recovery.json) | 552 | 6 | 26 | Full GRS `[64,45]` | Pass |
| [B](end-to-end-recoveries/binary-goppa/results/b/recovery/recovery.json) | 414 | 6 | 26 | Full GRS `[64,45]` | Pass |

Both finishing routes recover all 64 support positions, a compatible multiplier,
and a monic irreducible degree-six Goppa polynomial. Independent public verifiers
check exact code equality. Separate audits check equivalence with the generated
secret key. Minor-only and direct-only runs on instance A return the same keys
as combined mode.

Minor finishing uses a jet prefix through order 39; binary direct finishing uses
order 101 and takes a square root of the quotient.

The initial ordinary root count, 252, is below the degree bound 255. After
recovery, exact coefficient-level coherence through order 102 and verified Hasse
conditions certify the hold-out identities with `252+103=355>255` roots counted
with multiplicity. This certificate uses the reconstructed curve and public data.

The [validation record](end-to-end-recoveries/binary-goppa/results/validation.json) includes
arithmetic checks, independent key verification, and rejection of corrupted
kernel and multiplier inputs. The fixed examples do not assert that the rank
hypotheses hold for every key.

The [100-seed sweep](end-to-end-recoveries/binary-goppa/results/sweep-100/summary.json)
generates seeds 1–100 and computes a fresh kernel for each.
All 100 kernels pass the full Hasse checks; 99 cases pass both recovery routes
and all public geometric certificates. Independent local verification rechecks
all 198 returned keys and all 100 kernel hashes. No failed seed is omitted or
replaced. Kernel dimensions range from 372 to 706.

Seed 89 fails at quadratic branch separation: at held-out position 27 its
restricted quadratic rank is 9 instead of 15, despite gradient rank 21. The
exported matrices confirm these ranks under independent binary row reduction.
A [separate follow-up](end-to-end-recoveries/binary-goppa/results/sweep-100/diagnostics/seed89-alternate-result.json)
uses position 29, the next non-unit public column, with the same public
generator. Its new kernel has quadratic rank 15 and both routes recover keys
that pass independent verification and secret-support audits. The baseline
remains 99/100, with one additional successful hold-out computation recorded
separately. All 100 sampled keys are therefore eventually recovered.

The baseline takes 28.46 minutes with 24 single-threaded workers on a 96-vCPU
AMD EPYC-Genoa VM. Median hold-out and successful recovery times under that
concurrent load are 338.59 and 11.35 seconds. The
[case table](end-to-end-recoveries/binary-goppa/results/sweep-100/cases.csv) and
[batch audit](end-to-end-recoveries/binary-goppa/results/sweep-100/raw-validation.json)
retain every outcome.

The [subsequent position analysis](end-to-end-recoveries/binary-goppa/results/sweep-100/diagnostics/position-investigation/conclusions.json)
finds the same rank-nine defect at positions 27 and 42 among all 64 positions
(using a second information set for the original 28 systematic positions),
and reproduces it at position 27 after a systematic basis change.
Independent literal Hasse computations, native-coefficient comparisons, and
alternative PLUQ elimination support the rank result. The six missing mixed
quadratic terms form the neighboring Frobenius orbit. Third-order normalized
jet constraints remove the resulting six lines of spurious tangent directions;
the ordinary continuation matrix nevertheless has rank 24 instead of 26 at
these positions.
All projective tangent directions require GF(64); no proper-subfield descent
occurs. Global cubic identities of the recovered curve give rank 15 everywhere,
and removing one additional ordinary-position constraint also restores the
missing equations among genuine curve identities. The exact locations of the
two exceptions have not been characterized by a closed-form condition on the
Goppa polynomial.

A [public joint-solver follow-up](end-to-end-recoveries/binary-goppa/results/follow-ups/joint-seed89/summary.json)
also recovers seed 89 at positions 27 and 42 using the original failing kernels,
with position 29 as a control. Third-order compatibility cubics separate all six
branches. Orders three and four still leave continuation rank 24 at the two
exceptions; stationary order-five conditions raise it to 26. The same extra
lookahead restores the linear recurrence at every odd order. All three runs
reach a coherent jet through order 102 and pass both finishing routes; all six
keys pass independent public verification. No known curve or key is supplied
to this recovery, and no hold-out kernel is changed or recomputed.

## Wild Goppa recovery over GF(4)

The [wild Goppa collection](end-to-end-recoveries/wild-goppa/README.md) recovers complete
keys for GF(4) codes `[63,27]`, `[55,28]`, and `[64,28]`, with support in GF(64).
Every run computes one hold-out kernel and continues one coherent branch jet.
Stationary quadratic rank three separates the three Frobenius branches; the
fixed continuation matrix has rank `k-2`.

Both finishing methods pass in all three cases. Bootstrapping yields full minor
GRS codes of dimensions 46, 42, and 47 respectively. Direct extraction uses
ordinary evaluation ratios, with jet orders 91, 83, and 93. The recovered radical
defines both a cube and a fourth-power Goppa description exactly equal to the
public GF(4) code. Independent verifiers check all six recovered keys.

The [result table](end-to-end-recoveries/wild-goppa/results/README.md) gives kernel dimensions,
exported panel sizes, and measured runtimes. The large hold-out uses OpenMP-enabled
M4RI; the public recovery and coefficient-level coherence checks use no generated
secret data.

## Degree-seven toy recoveries

The [degree-seven collection](end-to-end-recoveries/degree7-binary-goppa/README.md) contains four
fresh binary Goppa codes with `(m,t)=(6,4)`, shortened from `[64,40]` to
`[54,30]`. The four seeds and the mixed multiplicity rule were fixed before
execution; the earlier pilot is excluded and no failed key was replaced.
The saved 256-polynomial sample for each trial is included, together with the
public inputs, native and Sage sources, both full recovered keys, and the
original stage logs and certificates.

All four samples satisfy every original Hasse condition and have gradient
rank 23 and restricted quadratic rank 15. Both routes recover all 64 support
positions and an irreducible degree-four Goppa polynomial, independently
verified against the original public code. Standard continuation reaches order
90 from the public first jet. The minor code has dimension 41 in every trial.
The [summary](end-to-end-recoveries/degree7-binary-goppa/results/summary.json) gives each outcome
and the measured kernel-stage sums, about 43 minutes per trial on one H200.
Those stages overlapped between trials, so the sums are not isolated wall-time
benchmarks.

The experiment used the practical CUDA/CADO solver; the proved uniform sampler
was not used. The instructions cover verification of saved files, CPU recovery
from the saved polynomials, and a fresh GPU kernel computation. The four
successful runs do not establish a general success probability.

## Reproduction code and seeds for the 820 additional trials

The [end-to-end sweep collection](end-to-end-recoveries/random-codes/README.md) provides code
and fixed seed lists for the paper's 600 GRS/alternant and 220 wild Goppa
trials. Every set uses seeds 1–100 except the larger `[55,28]` wild Goppa set,
which uses 1–20. The cohort excludes seed-zero pilots, earlier smaller GRS
examples, and the k=12 controls from the larger-GRS investigation.

The runners reuse the existing arithmetic implementations and include the
public Schur-square and ordinary-shortening checks for the four GRS sets. The
collection contains code and seeds; run outputs and hold-out kernels are
omitted. Each new run records failures and timeouts without replacing seeds.
The fixed examples above have their own saved results.

## Minor-code fullness

The [minor-code fullness collection](synthetic-experiments/minor-code-fullness/README.md)
tests 54 parameter sets in characteristic two, using tangent planes computed
from generated secret keys. It checks the minor code and subsequent support
recovery; it does not compute the tangent planes from a public hold-out kernel.
The 44 smaller sets have 100 codes each, and the 10 shortened Classic McEliece
sets have 20 codes each, for 4,600 trials.

| Family | Sets | Codes | Full minor GRS code and successful recovery |
|---|---:|---:|---:|
| GRS subcodes | 7 | 700 | 700 / 700 |
| Alternant | 11 | 1,100 | 1,099 / 1,100 |
| Wild Goppa | 15 | 1,500 | 1,500 / 1,500 |
| Binary Goppa, smaller sets | 11 | 1,100 | 1,100 / 1,100 |
| Classic McEliece, shortened | 10 | 200 | 200 / 200 |
| Total | 54 | 4,600 | 4,599 / 4,600 |

All 4,599 trials with two-dimensional tangent planes at every position yield
a full minor GRS code. Sidelnikov–Shestakov recovers an equivalent support,
and independent checks verify containment of the public code in the recovered
ambient GRS code. For shortened codes, these checks cover the retained
coordinates.

The remaining trial, `alternant-a3` trial 85, has a one-dimensional tangent
space at position 40 (zero-based). That column of the minor matrix is zero,
and the complete minor-code rank is 50 rather than 51. This tangent failure
is retained in the results; no trial raised an error or timed out. The
[per-set summary](synthetic-experiments/minor-code-fullness/results/summary.md) gives all
ranks, outcomes, and timings.

## Synthetic direct recovery

The [synthetic collection](synthetic-experiments/direct-recovery/README.md) runs the direct
recovery from one local jet computed with the secret key, at a random position,
with a random reparametrization and on a random Frobenius branch. The recovery
uses only the public generator matrix and the jet, and every recovered key is
verified against the public code and audited against the secret support.

| Family | Fields | Sets | Keys recovered |
|---|---|---:|---:|
| Binary Goppa, square-root variant | GF(2^8) to GF(2^13), up to shortened `mceliece6688128` | 4 | 160 / 160 |
| Wild Goppa | q = 3, 4, 5 | 6 | 360 / 360 |
| Generic alternant, quadratic and cubic forms | q = 2, 3, 4, 5 | 7 | 540 / 540 |

The binary Goppa control with adjacent orders fails on all 20 keys because no
second form exists, as the paper predicts. A SageMath cross-check recomputed
36 recovered keys independently and rejected corrupted ones. The
[summary](synthetic-experiments/direct-recovery/results/summary.md) contains all tables.

## Closure bootstrap

The [closure collection](synthetic-experiments/closure-bootstrap/README.md) starts from the
Frobenius closures of the tangent planes at a random subset of the positions,
computed with the secret key. It solves for the binary quadrics, separates the
branches at one position, continues a local jet there, and then runs
bootstrapping and the minor-code recovery. The recovery uses only the public
generator matrix and the closures. A key counts as recovered when the public
code embeds in the GRS code on the recovered support, and every recovered key
is audited against the secret support.

| Family | Positions | Sets | Keys recovered |
|---|---|---:|---:|
| Binary Goppa, `m = 6, 7, 8` | `n − 5t + 1` | 5 | 500 / 500 |
| | `n − 5t − 1` | 5 | 496 / 500 |
| | `n − 4t − 1` | 5 | 500 / 500 |
| Shortened `mceliece348864` and `mceliece6688128` | `n − 5t + 1` | 2 | 40 / 40 |
| | `n − 5t − 1` | 2 | 39 / 40 |
| Binary alternant | `2D + 1 − n` | 3 | 300 / 300 |

Each of the five keys lost at `n − 5t − 1` positions had one binary quadric
that does not vanish on the curve. Three further binary Goppa sets with fewer
binary quadrics than the gradient rank `k − m − 1` requires gave 0 / 300, as
expected. The [summary](synthetic-experiments/closure-bootstrap/results/summary.md)
contains all rows.

## Existence condition for binary direct extraction

The [existence-condition collection](synthetic-experiments/goppa-product-rank/README.md)
generates binary Goppa keys and computes `W`, the span of the products
`f_a f_b` of the secret curve polynomials. Equality `dim W = 2D−t+1` implies
that the quadratic forms of the binary direct extraction exist at every
position. The rank computation stops once this upper bound is reached.

| Keys | Count | `dim W = 2D−t+1` | Positions with both forms |
|---|---:|---:|---:|
| Classic McEliece, shortened with the GIAJS and Weis dimensions | 10,000 | 10,000 | 80,000 / 80,000 |
| One-point shortenings of 200 of these keys, 20 positions each | 4,000 | 4,000, at the bound of the shortened code | — |
| Classic McEliece `(m,t)`, `k` near the counting threshold | 400 | `min(C(k+1,2), 2D−t+1)` for all | — |
| Moderate shortened codes, `m` from 9 to 13 | 33,000 | 32,993 | all on keys at the bound |
| Unshortened codes with the TII challenge parameters | 3,000 | 3,000 | 24,000 / 24,000 |
| Small codes above the counting threshold | 98,008 | 94,065 | — |

In the small codes, every deficit at positive slack equals the number of linear
conditions from double and base points of the curve, including the point at
infinity. These points are repeated or zero columns of the public generator
matrix extended by its parity column. The remaining 26 deficient keys have zero
slack and one extra dependency each. The seven deficient moderate keys have
`t ≤ 20` and `k` at most three above the counting threshold: five have a double
point and two have zero slack. A SageMath cross-check recomputed 67 keys from
their coefficients and agreed on every rank and membership. For a key at the
bound, the code shortened at a position reaches its own bound if and only if
the quadratic forms vanishing on the curve have gradient rank `k−2` at that
position, so the 4,000 shortened codes give this rank at 4,000 positions. The
[summary](synthetic-experiments/goppa-product-rank/results/summary.md) contains all
tables.
