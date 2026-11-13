# Toy end-to-end key recovery

Four SageMath Jupyter notebooks start with a freshly generated public matrix,
compute a hold-out kernel, recover a coherent local jet, and return an equivalent
support and multiplier.

| Notebook | Public code | Ambient field | `D` | Hold-out `(d,s)` | Finishing |
|---|---|---|---:|---|---|
| [01-grs-char2.ipynb](01-grs-char2.ipynb) | `[24,6]` over GF(32) | GF(32) | 7 | `(3,1)` | Minor and direct |
| [02-alternant-char2.ipynb](02-alternant-char2.ipynb) | `[60,20]` over GF(8) | GF(64) | 39 | `(3,2)` | Minor and direct |
| [03-grs-odd.ipynb](03-grs-odd.ipynb) | `[17,4]` over GF(19) | GF(19) | 5 | `(3,1)` | Direct only |
| [04-alternant-odd.ipynb](04-alternant-odd.ipynb) | `[45,15]` over GF(11) | GF(121) | 29 | `(3,2)` | Direct only |

All four use generation seed 0. As fixed examples, they do not estimate
success probability. In each case `s(n-1)>dD`, so curve vanishing follows from
the ordinary interpolation bound. The alternant examples recover the Frobenius
branches from the public equations, without supplied tangents or jets.

## Read and run

The notebooks include their executed outputs. Open one with a SageMath kernel
and run all cells, or execute the complete collection from this directory:

```sh
sage -python execute.py
```

To execute one notebook:

```sh
sage -python execute.py 02-alternant-char2.ipynb
```

The notebooks were tested with SageMath 10.7. The runner uses `nbformat`,
`nbclient`, and `jupyter_client`, and discovers a registered SageMath kernel automatically;
`--kernel NAME` selects one explicitly. Every notebook gets a fresh kernel.
The runner saves outputs only after successful execution and writes a summary
with content hashes under [results/](results).
[results/PROSE-REVISION.json](results/PROSE-REVISION.json) records later Markdown edits; the
saved execution hashes identify the notebook versions that were run.

The implementation is in [reference.py](reference.py). Run the notebooks from
this folder. They require neither a C++ compiler nor precomputed kernels.

## What is checked

1. Public hold-out. Impose the retained Hasse conditions, verify the kernel
   residual, and check zero evaluation at the held-out point. Retained unit
   columns are handled exactly by the permitted monomial exponents.
2. Local branch recovery. Measure gradient rank and restricted stationary
   quadratic rank, and recover the projective branches by the quadratic quotient
   eigenvalue method. Select one branch using only these public equations.
3. Coherent continuation. Check that the fixed regular/next-order matrix has
   exactly the point/tangent plane as its kernel. Continue past saturation of the
   coefficient span and independently substitute the final series into every
   selected polynomial.
4. Characteristic-two minor finishing. Recover the full extension-field
   quadratic ideal from public values and a jet prefix, check tangent rank at
   every position, form all square-root minors, and verify exact equality with
   the recovered GRS code.
5. Direct finishing. Solve for the high-contact sections using exactly the
   sufficient jet prefix of order `delta*D-1`. Extract support ratios, handle
   infinity, and solve for a public-code multiplier. The odd GRS example uses
   cubic sections; the other three use quadrics.
6. Independent audits. After recovery, use the generated secret to check the
   polynomial identities, all branch planes, a coherent scalar/parameter series,
   and global Frobenius/Möbius equivalence of the recovered support. For the minor
   route, compare the recovered ideal with the full secret-curve ideal and check
   all tangent planes against one branch. The two characteristic-two finishing
   routes must agree up to a common Möbius transformation.

Both alternant cases also verify exact equality between the recovered
subfield code and the public code. A deliberately corrupted multiplier must fail
the public verification in every notebook.

The generated secret is kept in a separate variable for the final audits; the
recovery functions do not receive it. Each notebook shows the inputs and ranks
at every stage. Recovered supports and multipliers remain in `direct` and,
where applicable, `minor`; the result records are written to `results/`.

## Scope and interpretation

These examples cover generic GRS subcodes and generic alternant codes.
Characteristic two alone does not imply the special binary-Goppa relation
`v_(2r)=v_r^2`: these notebooks use the general stationary continuation, and their
direct extraction uses adjacent orders and ordinary ratios. Binary Goppa recovery
is implemented in the [C++ collection](../end-to-end-recoveries/binary-goppa/README.md).

In the `[24,6]` GRS example, the public points already determine the quadratic
ideal, so the bootstrap uses order zero. The `[60,20]` alternant example requires
order 19 to select one branch from the recovered jet. Its direct finishing
uses order 77. The odd examples use orders 14 and 57 respectively.

The notebooks recover every coordinate of the displayed codes without
preliminary shortening. Restoration of removed coordinates and performance
at cryptographic sizes are outside these tests.
