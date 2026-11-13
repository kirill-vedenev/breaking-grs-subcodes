# End-to-end recovery on generated codes

These experiments start recovery from generated public matrices. They compute
hold-out polynomials, recover local jets, and finish by the minor-code or
direct method. Generation secrets are kept for separate verification and
diagnostics, as described in each collection.

| Collection | Codes and runs | Material included |
|---|---|---|
| [Binary Goppa](binary-goppa/README.md) | Two fixed examples and a 100-seed sweep at `[64,28]` | C++ code, public inputs, recovered keys, and follow-up records |
| [Wild Goppa](wild-goppa/README.md) | Three fixed examples over GF(4) | C++ code, public inputs, recovered keys, and verification records |
| [Degree-seven binary Goppa](degree7-binary-goppa/README.md) | Four `[64,40]` codes shortened to `[54,30]` | Code, saved hold-out samples, both recovered keys per trial, and reproduction commands |
| [Random-code sweeps](random-codes/README.md) | 600 GRS/alternant and 220 wild Goppa trials | Code and fixed seed lists; run outputs are generated when the experiments are repeated |

The SageMath [toy notebooks](../toy-notebooks/README.md) give four smaller
worked examples. The tests in [synthetic-experiments](../synthetic-experiments/README.md)
start from secret-derived jets, tangent planes, or related data and test the
later recovery steps in isolation.
