# Synthetic experiments

These experiments use generated secret keys to supply tangent planes, local
jets, or Frobenius closures, or to compute the product space of the hidden
curve. They test individual recovery steps and rank conditions without
computing hold-out kernels.

| Collection | Input computed from the secret key | What is tested |
|---|---|---|
| [Minor-code fullness](minor-code-fullness/README.md) | Aligned tangent planes at every position | Minor-code dimension, Sidelnikov–Shestakov recovery, and containment of the public code in the recovered GRS code |
| [Direct recovery](direct-recovery/README.md) | One local jet, with random reparametrization and scaling | Direct support and multiplier recovery, followed by public-code verification |
| [Closure bootstrap](closure-bootstrap/README.md) | Frobenius closures of tangent planes at selected positions | Recovery of a local jet, bootstrapping, aligned tangents, and minor-code recovery |
| [Goppa product rank](goppa-product-rank/README.md) | Products of the secret curve polynomials | The sufficient condition for binary direct recovery and its relation to minor-code fullness and gradient rank |

The recovery algorithms receive the stated inputs, public code, and public parameters.
The product-rank experiments compute ranks directly from the secret curve.
For shortened codes, the minor-code and closure experiments verify the
retained support; they do not restore the removed coordinates.

Each collection documents its parameters, commands, saved results, and
checks. Run its commands from that collection's directory.
