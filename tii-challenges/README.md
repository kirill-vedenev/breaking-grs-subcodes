# TII challenge recoveries

Sources and saved results for TII-173, TII-249, and TII-253:

- [Experiments](EXPERIMENTS.md): shortening choices, coordinate conventions,
  parameters, timings, and recovery steps.
- [Results](results/README.md): original public keys, recovered equivalent
  keys, and verification reports.
- [Hold-out kernels](kernels/README.md): saved polynomial samples and their
  storage format.
- [Implementations](implementations/README.md): code, build requirements,
  saved-key checks, and reproduction commands.

Start with an implementation README to verify a saved key with Python 3.
Recovering a key again from its saved kernel requires SageMath and a C++
compiler. Recomputing the large hold-out kernels also requires the GPU tools
listed in the corresponding instructions.
