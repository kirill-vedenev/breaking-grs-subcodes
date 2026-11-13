# Sampling cost estimates

[*Sampling Binary Kernels with Block Wiedemann*](main.pdf) is a companion note to
[*Extending Distinguishing to Key Recovery for Subfield Subcodes of GRS codes*](https://eprint.iacr.org/2026/1747).

The note gives a sampler with exactly uniform accepted outputs and a proved
rejection bound for every fixed binary matrix. It counts the Boolean gates,
storage, and memory access needed to sample polynomials from the mixed hold-out
systems at the Classic McEliece parameters. The conditions and costs of subsequent
key recovery are treated in the main paper.

- [main.pdf](main.pdf): compiled note.
- [main.tex](main.tex): LaTeX source.
- [refs.bib](refs.bib): bibliography.

To rebuild with LuaLaTeX, Biber, and latexmk, run from this folder:

```sh
latexmk -lualatex -outdir=build main.tex
```

The rebuilt PDF is written to `build/main.pdf`.
