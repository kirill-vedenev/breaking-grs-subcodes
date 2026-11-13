# Recorded recoveries

| Instance | Public code | Kernel dimension | Panel | Jet order | Minor rank | Kernel time | Recovery time |
|---|---|---:|---:|---:|---:|---:|---:|
| A | [63,27] | 5867 | 256 | 91 | 46 | 8.77 s | 3.29 s |
| B | [55,28] | 22223 | 256 | 83 | 42 | 279.83 s | 38.67 s |
| C | [64,28] | 7844 | 256 | 93 | 47 | 13.06 s | 2.38 s |

Both finishing methods passed on all three instances. Each recovered key
passed independent public verification and comparison with the generated key.
Each successful run uses one held-out position and one computed kernel. Recovery times include complete
panel verification and geometric certificates.

The recorded large-instance kernel uses an OpenMP-enabled M4RI build with a
16-thread limit. The two small-instance timings use the system M4RI build.
The driver defaults to the available CPU count; its thread limit is configurable.

Minor-only and direct-only checks on instance A reproduce the combined-mode keys.
These fixed examples provide no statistical success estimate.

With all 96 available threads, the OpenMP configuration check recomputed
instance A in 2.12 seconds and produced a byte-identical kernel.
The driver defaults to the available CPU count and rejects a large allocation
when linked to single-threaded M4RI, unless explicitly overridden.
