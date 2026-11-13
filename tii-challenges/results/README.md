# TII challenge results

These equivalent binary Goppa decoding keys were recovered from public data.
Each instance includes the original public input, the recovered key, the
extension-field definition, and verification reports. The recovered key may
differ from the challenge author's original key.

Each instance README gives a short Python command to verify its saved key.
Optional intermediate audits are documented in the corresponding
[implementation README](../implementations/README.md).

| Instance | Original code | Field | Field-defining polynomial | Goppa degree | Result |
|---|---|---|---|---:|---|
| [TII-173](tii-173/README.md) | `[96,40]` | GF(128) | `z^7 + z + 1` (`0x83`) | 8 | Full equivalent key |
| [TII-249](tii-249/README.md) | `[235,107]` | GF(256) | `z^8 + z^4 + z^3 + z^2 + 1` (`0x11d`) | 16 | Full equivalent key |
| [TII-253](tii-253/README.md) | `[214,142]` | GF(256) | `z^8 + z^4 + z^3 + z^2 + 1` (`0x11d`) | 9 | Full equivalent key |

Both direct and minor finishings were tested on the retained TII-173 and
TII-249 codes; their archived minor results are in `evidence/minor-finishing/`.
TII-253 has a full key from each route. The direct-route key
(`direct-recovered-sk.json`) was recovered independently.
For TII-173 and TII-249, full-key extension was run only after direct finishing. See
[finishing coverage](../EXPERIMENTS.md#finishing-coverage).

## Files and conventions

- `pk_McEliece_N.txt`: original binary parity-check rows, followed by the
  field-modulus coefficient list. The public code is the binary kernel of
  this matrix, not its row space.
- `recovered-sk.json`: byte-for-byte copy of the completed equivalent key.
  `support[j]` is the field element at original public column `j`;
  `g` lists the monic Goppa polynomial coefficients in increasing degree.
- `field.json`: explicit field polynomial, its coefficient list and integer
  encoding, original parameters, and coordinate convention.
- `verification.json`: the saved independent full-key verification report.
- `timings.json`: recorded hardware, elapsed stage times, timing boundaries,
  and the hashes of their original measurement records.
- `evidence/`: saved jets, intermediate algebraic outputs, and verification
  certificates. These are outputs, not inputs to fresh recovery.
- `source.json`: original public-input provenance. Historical scope labels
  describe the input-acquisition stage, not the final recovery status.

An integer `a` encodes `sum_i ((a >> i) & 1) z^i` modulo the field polynomial.
Thus both `field_modulus` and `g` are low-degree-first coefficient lists,
but entries of `field_modulus` are binary and entries of `g` are encoded
extension-field elements. Support and multiplier arrays use the original
public column order, so no permutation needs to be inferred.

`dual_grs_multipliers` stores `nu`, with check rows
`nu_j * alpha_j^r`, `r=0..2t-1`, normalized by `nu[0]=1`.
`ambient_grs_multipliers` stores `lambda`, with generator rows
`lambda_j * alpha_j^r`, `r=0..n-2t-1`, where
`lambda_j = 1/(nu_j * product_{i!=j}(alpha_j + alpha_i))`.
When present, `goppa_parity_multipliers` stores `1/g(alpha_j)` for the
degree-`t` Goppa check; this is a different convention from `nu`.
The pair `(support,g)` and its field definition suffice to describe the key.

The standard-library verifiers in the corresponding
[implementation folders](../implementations/README.md) check exact
equality of the reconstructed binary parity-check row space with the original
public matrix, distinct support, support avoidance, and irreducibility of `g`.
Each instance README gives its verification command.

The original hold-out samples and their storage format are in
[kernels](../kernels/README.md). File checksums and
`PROVENANCE.json` in the implementation root record the original source
paths and hashes. The original runs are preserved, including the additional
TII-253 direct finishing, which reused the hold-out sample.
Some historical reports contain local paths from those runs; the reproduction
commands do not depend on those paths. Dates were removed from the saved
records and from the recorded source names.
From this directory, `shasum -a 256 -c SHA256SUMS` checks every distributed file.
