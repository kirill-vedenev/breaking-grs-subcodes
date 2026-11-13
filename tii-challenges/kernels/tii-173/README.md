# TII-173 hold-out sample

The sample contains 256 independent degree-5 polynomials in
35 variables, represented by their binary coefficients on 324,632
squarefree monomials. The held-out position is original column 0.
The WHK1 format uses lexicographic tuples.
The sample spans only part of the hold-out kernel.

- `public.json`: exact variable basis, shortened public matrix, field modulus,
  original coordinate map, held-out point, and multiplicities.
- `operator.txt`: the complete public Hasse constraints.
- `verification.json`: saved independent membership and rank check.
- `storage.json`: dimensions, ordered chunks, original SHA-256, and byte sizes.
- `kernel.thk1`: the original accepted colex encoding for the standalone
  CPU Hasse verifier. It contains the same 256 forms as `kernel.whk1`.

The original packed file has 10,388,244 bytes and SHA-256
`8eb72336771ffc13fc91f9846c0a93e57029b69fc56a9973eb26c0a0aad18f18`. It is stored in 1 file, listed
in `storage.json`; concatenation restores the original coefficient bytes.
The [root storage specification](../README.md) defines the header, bit order,
monomial order, panel layout, and coordinate conventions.

From `kernels/`, restore the input into a fresh directory:

```sh
python3 materialize.py tii-173 --output tii173-input
```

The output contains `kernel.bin` (WHK1, for recovery), `kernel.thk1` (THK1,
for Hasse verification), and `public.json`. Both archived encodings are
authenticated by `storage.json` and `SHA256SUMS`. The THK1 SHA-256 is
`a1572050cf858742af5c496f0416b7bb56b8f44c52b817848904a30fa144d6e6`.
The [optional implementation audits](../../implementations/tii-173/README.md#optional-audits-of-saved-results)
give CPU-only membership/rank checks and an exact THK1-to-WHK1 comparison.

The [implementation](../../implementations/tii-173/README.md)
contains direct verification and recovery commands. The
[final key and later audits](../../results/tii-173/README.md) are stored
separately. The hold-out report's historical `curve_identity_certified` flag
describes that stage only.
