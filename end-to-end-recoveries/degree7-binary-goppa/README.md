# Degree-seven toy key-recovery experiments

This directory contains the four degree-seven trials in *Extending
Distinguishing to Key Recovery for Subfield Subcodes of GRS codes*: their
public inputs, 256 hold-out polynomials per trial,
source code, full keys recovered by both routes, and verification and timing
records. The earlier pilot is excluded.

## Experiment and results

Each instance is a binary Goppa code over GF(64), with a random monic
irreducible quartic and a shuffled full-field support. It is shortened from
[64,40] to [54,30], so D = 45. The hold-out forms are homogeneous multilinear
polynomials of degree seven in 30 variables. The first twelve retained finite
non-unit constraint columns receive multiplicity five; the other retained
points, including infinity, receive six. The held-out position receives zero.
The unit columns satisfy the required conditions automatically. This gives
312 imposed zeros, 2,035,800 unknowns, and 2,136,372 reduced constraint rows.

[PLAN.json](PLAN.json) records the four seeds, multiplicity rule, and solver
seeds; [TRIALS.json](TRIALS.json) records the input hashes and coordinate choices.
All four keys were retained regardless of their ranks, curve vanishing, dual
weights, or recovery outcomes.

| Trial | Key seed | Gradient / quadratic rank | Minor / direct full key | Kernel stages, minutes | Finite dual words of weight 1–12 |
|---|---:|---|---|---:|---:|
| trial-01 | 2026101701 | 23 / 15 | Pass / pass | 43.56 | 373 |
| trial-02 | 2026101702 | 23 / 15 | Pass / pass | 43.03 | 360 |
| trial-03 | 2026101703 | 23 / 15 | Pass / pass | 43.20 | 370 |
| trial-04 | 2026101704 | 23 / 15 | Pass / pass | 43.41 | 350 |

All four samples have 256 independent polynomials satisfying every original
Hasse condition. Each gives six Frobenius conjugate branches and continuation
rank 28. Continuation starts from the recovered first jet and reaches order 90.
Minor finishing uses its order-37 prefix, obtains tangent rank 28 at all 54
retained positions, and constructs a [54,41] GRS minor code. Both routes then
recover all 64 original coordinates and an irreducible degree-four Goppa
polynomial. The independent verifier checks equality with the original public
code. Both routes share the same polynomial sample and continued jet.

[results/summary.json](results/summary.json) contains the numerical results. [results/dual-statistics.json](results/dual-statistics.json)
contains the two independent exhaustive dual-word enumerations. The four
trials establish recovery on these instances only. They do not establish a
general success rate or compare attacks on the same keys.

The kernel computation used the practical CUDA/CADO block Wiedemann solver
with standard projection, a 512-bit block, 8,210 sequence terms, solver seed
2026101602, and compressor seed 2026101603. The proved uniform sampler from
the companion note was not used.

## Verify the saved artifact

Run commands from this directory. Python 3.10 or later is sufficient:

```sh
python3 -B scripts/verify.py
```

This checks every file in SHA256SUMS, all input identities, concatenated kernel
hashes and headers, and all eight full keys against their original public
codes. It does not repeat the Hasse tests or recovery. An optional report can
be saved with `--report runs/verification.json`. The checksum-only equivalent
is `shasum -a 256 -c SHA256SUMS`.

To verify one recovered key independently:

```sh
python3 -B scripts/verify_key.py trials/trial-01/inputs/original-public.json \
  results/trial-01/recovery/direct/full-key/equivalent-key.json
```

This verifier uses Python's standard library and imports no recovery code or
generation secrets.

## Recover again from the saved polynomials

Requirements: a C++17 compiler, GNU Make, and SageMath 10.7 with its NumPy.
No GPU, CADO, or M4RI development installation is needed for this replay.
Use fresh output directories; the scripts refuse to overwrite an existing run.

```sh
make -j 4 cpu
sage -python scripts/recover.py trial-01 --output runs/replay-01 --threads 4
```

Repeat for trial-02, trial-03, and trial-04, with a fresh output directory for
each. The command reconstructs the coefficient matrix and panel, checks every
original Hasse order and sample independence with the native verifier, extracts
the Hasse tensors, checks branch separation, continues the jet, runs both
finishing routes, restores the shortened coordinates, and independently checks
both full keys. The archived keys, jets, and local-rank reports are not inputs.

The new keys are written to `runs/replay-01/minor/full-key/equivalent-key.json`
and `runs/replay-01/direct/full-key/equivalent-key.json`. The run's summary and
stage logs are saved beside them. An equivalent recovered key need not have
the same bytes as the archived key. Allow about 1 GiB of working disk per
replay. The recorded local finishing took about 40 seconds per trial; fresh
replays also include kernel reconstruction and verification.

The native series and Taylor implementations have an independent small
Sage-polynomial check:

```sh
sage -python tests/test_native.py --output runs/native-check
```

## Recompute a hold-out sample

This part requires a CUDA development toolkit with cuFFT, an appropriate
NVIDIA GPU, and the pinned, patched CADO-NFS `lingen_b64`. The original run
used one H200 NVL and CUDA 12.8. Build the worker for your device:

```sh
make cpu
make gpu CUDA_ARCH=90
```

CADO must be at revision `2d98fe176c11342cfa5d6c0d795a992cbe433c11`.
Its upstream URL, revision, and patch hash are recorded in
[holdout/patches/cado.json](holdout/patches/cado.json). On Linux, use a GNU
C/C++ toolchain with OpenMP, CMake, and GMP development files; the original
build also had FLINT, MPFR, and HWLOC development packages. For an unpatched
checkout at this revision:

```sh
CADO_SRC=/path/to/cado-nfs
patch -d "$CADO_SRC" -p1 < holdout/patches/cado-coefficient-parallel.patch
cmake -S "$CADO_SRC" -B "$CADO_SRC/build-degree7" \
  -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF
cmake --build "$CADO_SRC/build-degree7" --target lingen_b64 -j 4
LINGEN="$CADO_SRC/build-degree7/linalg/bwc/lingen_b64"
python3 -B scripts/holdout.py trial-01 --output runs/new-kernel-01 \
  --lingen "$LINGEN" --threads 32 --ram-gb 64
```

The driver uses the original 1,024-term segmentation, seeds, block size,
CADO parameters, and reconstruction orientation. It saves sequence, generator,
reconstruction, and verification logs. It processes one trial at a time and
retains the replay files; allow at least 8 GiB of disk per trial. A fresh
kernel can be passed through the public recovery pipeline:

```sh
sage -python scripts/recover.py trial-01 --output runs/new-recovery-01 \
  --kernel runs/new-kernel-01/reconstruction/kernel.thk1 --threads 4
```

The arithmetic sources are unchanged. The GPU run script was rewritten for
this directory and has not been used to repeat the kernel computation. The
original script is in the source archive.

## Regenerate the public inputs and dual-weight counts

The generator reproduces all four trials from the fixed seeds. Generation
secrets go into a separate `generation-audits` directory under the requested
output and are never inputs to kernel computation or recovery:

```sh
python3 -B scripts/prepare_batch.py --output runs/regenerated
python3 -B scripts/validate_inputs.py --inputs runs/regenerated/trials/trial-01/inputs
```

The input records include the curve-vanishing checks made during generation.
These checks did not filter the four keys or supply recovery with a jet or
curve. The flags `prepared`, `kernel_solved: false`, and
`recovery_tested: false` in these unchanged records describe their creation time;
completed outcomes are in `results/`.

To repeat the exhaustive dual-weight counts using both C++ and NumPy:

```sh
make cpu
sage -python scripts/dual_statistics.py --output runs/dual-weights
```

This enumerates all 2^24 dual words per trial and checks that both independent
implementations give the same distribution. It requires NumPy 2 or later.

## Checks of this distribution

[validation/PACKAGING-CHECKS.json](validation/PACKAGING-CHECKS.json) records the checks performed on
this distribution. All four CPU recoveries were rerun from the packaged kernel
chunks, starting without saved jets or keys; both routes passed in every case.
All public input files and the trial manifest regenerated byte for byte. The
independent native Taylor/series test and both exhaustive dual enumerations
passed. New verification records are in [validation/](validation); the
original experiment results in `results/` are unchanged apart from the file
paths described under [File formats and provenance](#file-formats-and-provenance).
GPU kernel generation was not rerun.

## File formats and provenance

| Path | Contents |
|---|---|
| `trials/trial-*/inputs/` | Original and shortened public matrices, exact operator/request files, seeds, preparation diagnostics |
| `trials/trial-*/kernel/` | Saved coefficient chunks, their storage manifest, original polynomial index |
| `results/trial-*/holdout/` | Original stage logs, measurements, kernel and branch certificates |
| `results/trial-*/recovery/` | Original jets, both recovered keys, intermediate outputs, independent checks, and stage logs |
| `results/batch/` | Original hardware description and cohort execution state |
| `holdout/`, `recovery/` | Numerical sources used in the completed experiment |
| `scripts/` | Generation, verification, reconstruction, and portable run commands |
| `provenance/` | Original aggregate reports, source manifest, and complete original source archive |

Each file stores 256 independent forms sampled from the hold-out kernel. They
do not span the whole kernel. Its THK1 header is 20 bytes: four ASCII magic
bytes followed by four unsigned little-endian 32-bit values
`(k,d,N,count)`. The payload has 256 rows of `ceil(N/8)`
bytes, least-significant bit first. Monomials are in colex order: the index of
`{a1 < ... < ad}` is `sum(binomial(ai,i), i=1..d)`. Variables correspond to the
exact generator rows in `inputs/public.json`. The complete file has
65,145,620 bytes and is split into consecutive raw chunks of at most 32 MiB.
Each `storage.json` records offsets, sizes, and hashes.

To reconstruct a kernel without running recovery:

```sh
python3 -B scripts/materialize.py trial-01 --output runs/input-01
```

Adding `--panel` under `sage -python` also recreates the 8-word-wide transposed
panel and colex monomial masks, checking their original hashes. Higher Hasse
tensors, the panel, monomial masks, GPU sequence states, and CADO generators
are reproducible working files and are omitted from the saved distribution.
All coordinate indices are zero based; `original_indices` maps the shortened
column order back to the original public matrix. Field elements use the
polynomial basis modulo `x^6+x+1` (integer modulus 67).

In 81 original records under `results/` and `provenance/`, absolute file paths
were made relative to the original experiment directory. The SHA-256 hashes
and byte counts that 13 further records store of these files were updated to
the packaged files. Date and time fields were removed from seven records
and from `PLAN.json`, with stored hashes updated in the same way. No other content of the records
was changed.

[PROVENANCE.json](PROVENANCE.json) records each copied file's original path and hash, along
with the packaging changes. Concatenating the kernel chunks restores the
original bytes. Recovery and hold-out arithmetic sources are unchanged.
Paths and omitted working files mentioned in old records refer to the original
experiment directory and are needed only to trace the original runs. [provenance/original-source.tar.gz](provenance/original-source.tar.gz) is the unchanged source
archive, including patched CADO sources and their licenses. Its SHA256SUMS
checks the extracted archive; this directory has its own manifest.

The original kernel-stage durations sum GPU sequence, CPU generator, GPU
reconstruction, and verification times. One CPU stage overlapped another
trial's GPU stage, so their sums are not isolated wall-clock benchmarks. The
four-trial queue took 130.21 minutes. Peak child-process RSS was about 4.05
GiB, and peak sampled device-wide memory was about 14.62 GiB. These
measurements do not determine the minimum hardware requirements. The host had
256 GiB of RAM, and CADO used 32 threads with a 64 GiB setting.

The repository's MIT license applies to the author's code and documentation.
The Saarinen worker retains its own MIT notice under `holdout/third_party/`.
The CADO patch and archived CADO sources retain their upstream licenses,
including the notices for bundled dependencies.
