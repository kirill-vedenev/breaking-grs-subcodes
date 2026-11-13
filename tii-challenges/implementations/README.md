# TII key-recovery implementations

This directory contains the sources and stage commands used to recover
TII-173, TII-249, and TII-253. The original results are included without a new
run. TII-253 also has an additional direct finishing run, using its
saved public jet prefix and hold-out sample.

Each instance README starts with a Python command that verifies the saved
key and ends with optional audits of the kernel, jet, and curve geometry.

| Instance | Hold-out computation | Finishing |
|---|---|---|
| [TII-173](tii-173/README.md) | Exact binary panel elimination on two GPUs | Direct and minor finishing on the same order-296 jet; full key via direct finishing and deshortening |
| [TII-249](tii-249/README.md) | GPU sequence → CPU generator → GPU reconstruction | Direct and minor finishing on the same order-290 jet; full key via direct finishing and IKZ extension |
| [TII-253](tii-253/README.md) | GPU sequence → CPU generator → GPU reconstruction | Minor finishing from an order-82 jet; direct finishing after continuation to 198; both produce full keys |

TII-173 and TII-249 also contain the recorded `recovery/minor/` implementation.
Those runs reached equivalent retained supports; full-key extension was not
repeated separately. TII-253 has separately completed direct and minor finishings. See the
[finishing coverage](../EXPERIMENTS.md#finishing-coverage) for exact scope.

Each instance includes all its computational sources. `holdout/` contains
native kernel computation and public-input preparation; `recovery/` contains
jet recovery, finishing, and independent verifiers. For the larger hold-outs,
GPU sequence generation, CPU polynomial generation, and GPU reconstruction
are separate commands. Run or schedule them as needed; server launchers,
transfers, supervisors, download scripts, and multi-server orchestration are
not included.

The commands reference the sibling [kernels](../kernels/README.md)
and [results](../results/README.md) folders. Use fresh output folders.
Recovery starts from public inputs and either the saved hold-out sample or
a newly computed sample. Only the verification commands read saved keys and
evidence.

## Dependencies and source fidelity

The completed CPU recoveries used SageMath 10.7 with NumPy and a C++17 compiler.
TII-173's native formal-substitution helpers also link M4RI from Sage. The
independent key verifiers use only Python's standard library. Hold-out GPU code
requires a CUDA development toolkit; the original runs used H200 hardware.
TII-249 and TII-253 additionally use cuFFT and the pinned, patched CADO-NFS
`lingen_b64` described in their hold-out instructions.

[PROVENANCE.json](PROVENANCE.json) maps every copied source and archived
artifact to its original path and SHA-256. Computational source is unchanged
except for documented path relocation, a fresh output argument for TII-253
preparation, and removal of an obsolete module search path. TII-253's public
preparation helper contains six unchanged functions extracted from the original
module; its unused command-line entry point is omitted. Build files and
documentation are adapted for this layout. The kernel materializer is new
file-conversion tooling that checks the reconstructed bytes against their
original hashes.

MIT notices for the reused Saarinen code are retained under each hold-out's
`third_party/saarinen/`. The CADO patch, pinned revision, and upstream license
are retained under `holdout/patches/`. Build the binaries locally; precompiled
machine-specific binaries are not included. This layout has passed checks of syntax, archive integrity, and the saved
keys, but has not been used for
a new GPU-to-key run. The additional TII-253 direct finishing in
`recovery/direct/` uses the sources from the completed public-input run.

`PACKAGING-CHECKS.json` records the integrity and saved-key checks performed
for this layout, including the saved-artifact audits in the per-instance
READMEs. These checks did not rerun recovery. `SHA256SUMS` covers the
distributed files under this root, excluding Finder metadata, caches, and
generated build/run files. Run `shasum -a 256 -c SHA256SUMS` from this directory to check them. The results
and kernel collections each have their own checksum file.
