# Key recovery from first-order Frobenius closures

These experiments recover a support from the Frobenius closures `E_i` of
the tangent planes at a set `S` of positions. They use the bootstrapping
variant in Section 5.4 of the paper ("Other variants of bootstrapping") and
the recovery methods of Sections 4 to 6.

We computed the closures from the secret key, so these synthetic tests begin
after closure recovery and do not compute a hold-out kernel. The paper
reports their results in the "Closures" rows of its synthetic-experiments
table.

## Codes

- Binary Goppa codes `Γ(α, G)` over `F_2` with `K = F_{2^m}`, `G` a
  uniformly random monic irreducible polynomial of degree `t` and `α` a random
  subset of `K` of size `n`. Here `D = n − 2t − 1`.
- Shortened Classic McEliece sizes: binary Goppa codes with
  `(m, n, t) = (12, 922, 64)` and `(13, 1888, 128)`, the lengths
  `n = k + mt` of `mceliece348864` and `mceliece6688128` shortened to
  `k = 154` and `k = 224`. A shortened Goppa code is the Goppa code on the
  retained support, so the codes are generated directly at this length. They
  are not produced by the Classic McEliece key generator.
- Binary alternant codes `GRS_{n−t}(α, λ) ∩ F_2^n` of degree `t` with
  uniformly random nonzero multipliers. Here `D = n − t − 1`.

The public matrix `Y` is the reduced generator matrix, with columns `y_j`.

## The closures

With `ν` the multiplier of the dual of the ambient GRS code, the vector

```
ŷ_j = Σ_{l ≠ j} ν_l y_l / (α_j − α_l)
```

spans the tangent plane at `y_j` together with `y_j`. The closure `E_j` is
the span over `F_2` of `y_j` and the binary coordinates of `ŷ_j`; it has
dimension `m + 1`. The recovery receives a reduced basis of `E_i` for `i` in
a random set `S`, and nothing else derived from the secret.

Number of positions:

| Label | `|S|` | Meaning |
|---|---|---|
| `n-4t-1` | `2D + 1 − n` | `n + |S| > 2D`, so step A below gives exactly the binary quadrics of the curve |
| `n-5t-1` | `2D + 1 − t − n` | the smallest number allowed by the `t` hidden constraints of binary Goppa codes |
| `n-5t+1` | `2D + 3 − t − n` | the same with two positions in reserve |
| `2D+1-n` | `2D + 1 − n` | alternant codes |

## Recovery

Inputs: `Y`, `m`, `D`, the Goppa degree `t` for a binary Goppa code, and the
closures.

| Step | Computation | Rank condition checked |
|---|---|---|
| A | Binary quadrics `h` with `h(y_j) = 0` for all `j` and `∇h(y_i)·v = 0` for `v ∈ E_i`, `i ∈ S` | |
| B | At one random position: the kernel of the gradients, the restriction of the quadrics to it, and one tangent vector `v_1` of one Frobenius conjugate branch | nullity `m + 1`, quadratic rank `C(m,2)` |
| C | Local jet of order `2D − n + 1` by continuation on the binary quadrics. Binary Goppa: the normalized expansion `v_{2j} = v_j^{∘2}`. Alternant: the look-ahead recursion, with the quadrics singular at the point | rank `k − 2`; every order consistent |
| D | Quadrics over `K` with `h(y_j) = 0` for all `j` and `T^{r+1} ∣ h(V(T))` | rank `2D + 1 − t` (binary Goppa), `2D + 1` (alternant) |
| E | Tangent plane at every public column, as the kernel of the gradients | rank `k − 2` at every column |
| F | Minor code, Sidelnikov–Shestakov, multiplier of the public code on the recovered support `α′` | minor code of dimension `D − t` (binary Goppa), `D` (alternant) |

For alternant codes step B uses the quadrics singular at the point. For
binary Goppa codes it uses all binary quadrics, since `h(v_1) = 0` for every
binary quadric `h` of the curve.

Step D uses all monomials. For the tangent-plane computation in step E, we
chose `k + 40` random elements of that solution space. Step F uses a random
subset of the minor-matrix rows and `src/minor.py` from the minor-code
fullness collection.

A key counts as recovered when Sidelnikov–Shestakov returns a support
`α′`, the minor code is a GRS code on it, and the public code embeds in
`GRS_{D+1}(α′, λ′)` for a recovered multiplier `λ′`; the embedding is checked
on all rows. The coordinates removed by shortening are not restored.

## Audits

After recovery, two audits use the secret key to check the following
properties. Their results are not used by the recovery algorithm.

- Every recovered tangent plane contains the same Frobenius conjugate of the
  true tangent vector.
- The recovered support is a Möbius image of a Frobenius conjugate of the
  secret support (all cross-ratios agree).

Every recovered key passed both audits.

## Results

Binary Goppa codes, 100 keys per parameter set and number of positions:

| `(m,n,t)` | `k` | Binary quadrics | `n-5t+1` | `n-5t-1` | `n-4t-1` |
|---|---:|---:|---:|---:|---:|
| (6,64,6) | 28 | 149 | 100/100 | 99/100 | 100/100 |
| (7,112,10) | 42 | 370 | 100/100 | 97/100 | 100/100 |
| (7,122,11) | 45 | 457 | 100/100 | 100/100 | 100/100 |
| (8,192,16) | 64 | 1,007 | 100/100 | 100/100 | 100/100 |
| (8,256,24) | 64 | 751 | 100/100 | 100/100 | 100/100 |

Shortened Classic McEliece sizes, 20 keys per row and number of positions:

| Parameter set | `(m,n,t)` | `k` | Positions | Binary quadrics | `n-5t+1` | `n-5t-1` |
|---|---|---:|---:|---:|---:|---:|
| `mceliece348864` | (12,922,64) | 154 | 603 / 601 | 3,812 | 20/20 | 19/20 |
| `mceliece6688128` | (13,1888,128) | 224 | 1,249 / 1,247 | 7,113 | 20/20 | 20/20 |

Binary alternant codes, 100 keys per parameter set:

| `(m,n,t)` | `k` | Positions | Binary quadrics | Recovered |
|---|---:|---:|---:|---:|
| (6,64,5) | 34 | 53 | 213 | 100/100 |
| (7,120,9) | 57 | 101 | 826 | 100/100 |
| (8,240,20) | 80 | 199 | 1,408 | 100/100 |

In every recovered key the first position tried in step B passed both rank
conditions, every order of the continuation was consistent, the gradient rank
was `k − 2` at every public column, and the number of binary quadrics of a
binary Goppa code was `C(k+1,2) − (n+1) − m(n−5t−2)`.

The `n-5t-1` experiment used exactly `n − 5t − 1` positions and failed on
five of the 540 binary Goppa keys. Each of these keys had one additional
binary quadric compared with the computation using more positions. This
quadric did not vanish on the curve, and step B or step C failed. With two
positions in reserve all 540 keys were recovered. Every choice of the number
of positions used the same seeds and therefore the same keys.

Step B needs gradient rank `k − m − 1`. Three further binary Goppa sets
have fewer binary quadrics than this required rank. No key was recovered
in these sets:

| `(m,n,t)` | `k` | Binary quadrics (most keys) | `k − m − 1` | Recovered |
|---|---:|---:|---:|---:|
| (6,56,6) | 20 | 9 | 13 | 0/100 |
| (7,99,10) | 29 | 6 | 21 | 0/100 |
| (8,241,24) | 49 | 31 | 40 | 0/100 |

The [summary](results/summary.md) has one row per parameter set and number of
positions.

## Not covered

- Closures obtained from hold-out kernels: every closure here comes from the
  secret key.
- Direct recovery from the continued jet; only the minor-code method is run.
- Wild Goppa codes and odd characteristic.

## Build and run

The scripts require SageMath 10.7. They run under `sage -python` and use
Sage's NumPy; the summary requires Python 3. No compilation is needed.

One key:

```sh
sage -python src/closure_chain.py --m 6 --n 64 --t 6 --extra 2 --seeds 1 --output runs/one.jsonl --verbose
```

Options: `--family goppa|alternant`, `--mode saving|proved` (for binary Goppa
codes, `n − 5t − 1` or `n − 4t − 1` positions), `--extra E` (`E` further
positions), `--seeds a-b`. A seed determines the key, the positions and all
random choices of the recovery.

All parameter sets:

```sh
./run.sh                           # everything; JOBS processes in parallel
SKIP_CM=1 KEYS_TOY=10 ./run.sh     # ten keys per small set, no Classic McEliece sizes
python3 summarize.py results       # rebuild results/summary.md from the records
```

`run.sh` overwrites `results/keys.jsonl.gz`; set `OUT` to write elsewhere.

Time on one core: 0.6 to 12 s per binary Goppa key with `m ≤ 8`, 1 to 40 s per
alternant key, about 7 minutes at `(12,922,64)` and about one hour at
`(13,1888,128)`.

A rerun on another platform reproduced every recorded rank and outcome. The
selected Frobenius conjugate can differ, because it depends on the order in
which the roots of a polynomial are returned.

## Files

- `src/closure_chain.py`: key generation, closures, recovery and audits.
- `src/minor.py`: copy of `../minor-code-fullness/src/minor.py`
  (minor code checks, Sidelnikov–Shestakov, GRS multipliers).
- `run.sh`: all parameter sets; `summarize.py`: the summary tables.
- `results/keys.jsonl.gz`: one JSON record per key (2,180 keys).
- `results/summary.md`, `results/summary.json`: per-set tables.
- `results/build-info.json`: software versions and source hashes.
