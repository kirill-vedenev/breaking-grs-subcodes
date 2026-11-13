# Supplementary material

This repository contains the programs, experimental artifacts, and a companion note for the paper
[*Extending Distinguishing to Key Recovery for Subfield Subcodes of GRS codes*](https://eprint.iacr.org/2026/1747),
Cryptology ePrint Archive, Paper 2026/1747.

## Contents

| Folder | Contents |
|---|---|
| [toy-notebooks](toy-notebooks/README.md) | Four SageMath notebooks illustrating key recovery for GRS subcodes and alternant codes in characteristic two and odd characteristic. Includes executed outputs, a reference implementation, and verification records. |
| [tii-challenges](tii-challenges/README.md) | Implementations and artifacts for TII-173, TII-249, and TII-253: public keys, saved hold-out polynomial samples, recovered equivalent keys, verification reports, and timing records. |
| [end-to-end-recoveries](end-to-end-recoveries/README.md) | Recovery experiments on generated codes, starting from public matrices. Contains binary and wild Goppa examples with saved recovery artifacts, and code and seeds for GRS/alternant and wild Goppa sweeps. |
| [synthetic-experiments](synthetic-experiments/README.md) | Programs and saved results testing minor-code fullness, direct recovery, closure bootstrapping, and Goppa product ranks. These experiments use auxiliary data computed from generated secret keys. |
| [wiedemann-cost-estimate](wiedemann-cost-estimate/README.md) | Companion note on uniform kernel sampling with block Wiedemann, including Boolean gate, storage, and memory-access cost estimates for Classic McEliece. PDF and LaTeX source. |

[EVIDENCE.md](EVIDENCE.md) summarizes the experiments, outcomes, and verification
checks. Dependencies, build steps, and running instructions are documented in
the README files of the corresponding folders.

## Citation

If you use this code or these results, please cite the paper.

```bibtex
@misc{cryptoeprint:2026/1747,
  author       = {Kirill Vedenev},
  title        = {Extending Distinguishing to Key Recovery for Subfield
                  Subcodes of {GRS} codes},
  howpublished = {Cryptology {ePrint} Archive, Paper 2026/1747},
  year         = {2026},
  url          = {https://eprint.iacr.org/2026/1747}
}
```

## License

The code, documentation and results in this repository are released under the
[MIT License](LICENSE), with these exceptions:

- `tii-challenges/implementations/*/holdout/third_party/saarinen/` contains code by
  Markku-Juhani O. Saarinen from
  [tii254-artifact](https://github.com/mjosaarinen/tii254-artifact), under his
  MIT License in the `LICENSE.txt` beside it.
- `tii-challenges/implementations/*/holdout/patches/cado-coefficient-parallel.patch`
  modifies [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs) and is
  distributed under the GNU Lesser General Public License, version 2.1, as
  CADO-NFS itself; see `CADO-COPYING` beside it.
- `end-to-end-recoveries/degree7-binary-goppa/holdout/` retains the same Saarinen and CADO
  notices beside its worker and patch. Its `provenance/original-source.tar.gz`
  also contains the original CADO source tree and bundled dependencies, under
  the upstream licenses retained in that archive.

Some programs link against M4RI and M4RIE, which are licensed under the GNU
General Public License, version 2 or later, and some scripts run under
SageMath. These libraries are not included. Compiled programs that contain
them must be distributed under the terms of the GPL.
