# Summary

2180 keys in 25 parameter sets and closure counts.

"Recovered" means that the support was recovered and the public code embeds in
`GRS_{D+1}` on it. "Audits passed" counts keys whose recovered tangent planes all lie
on one Frobenius conjugate and whose recovered support is a Möbius image of a
conjugate of the secret support. Times are for one process and include the audits.

| Family | `(m,n,t)` | `k` | Closures | Positions | Keys | Recovered | Binary quadrics | Other outcomes | Audits passed | Median seconds |
|---|---|---:|---|---:|---:|---:|---|---|---:|---:|
| Binary Goppa | (6,56,6) | 20 | `n-5t-1` | 25 | 100 | 0 | 9 (96 keys), 10 (4 keys) | branch_separation_failed: 100 | 0 | 0.1 |
| Binary Goppa | (6,64,6) | 28 | `n-4t-1` | 39 | 100 | 100 | 149 (100 keys) |  | 100 | 0.6 |
| Binary Goppa | (6,64,6) | 28 | `n-5t+1` | 35 | 100 | 100 | 149 (100 keys) |  | 100 | 3.4 |
| Binary Goppa | (6,64,6) | 28 | `n-5t-1` | 33 | 100 | 99 | 149 (99 keys), 150 (1 keys) | branch_separation_failed: 1 | 99 | 0.6 |
| Binary Goppa | (7,99,10) | 29 | `n-5t-1` | 48 | 100 | 0 | 6 (98 keys), 7 (2 keys) | branch_separation_failed: 100 | 0 | 0.4 |
| Binary Goppa | (7,112,10) | 42 | `n-4t-1` | 71 | 100 | 100 | 370 (100 keys) |  | 100 | 2.5 |
| Binary Goppa | (7,112,10) | 42 | `n-5t+1` | 63 | 100 | 100 | 370 (100 keys) |  | 100 | 15.6 |
| Binary Goppa | (7,112,10) | 42 | `n-5t-1` | 61 | 100 | 97 | 370 (97 keys), 371 (3 keys) | continuation_inconsistent: 1, branch_separation_failed: 2 | 97 | 2.5 |
| Binary Goppa | (7,122,11) | 45 | `n-4t-1` | 77 | 100 | 100 | 457 (100 keys) |  | 100 | 3.4 |
| Binary Goppa | (7,122,11) | 45 | `n-5t+1` | 68 | 100 | 100 | 457 (100 keys) |  | 100 | 20.0 |
| Binary Goppa | (7,122,11) | 45 | `n-5t-1` | 66 | 100 | 100 | 457 (100 keys) |  | 100 | 3.3 |
| Binary Goppa | (8,192,16) | 64 | `n-4t-1` | 127 | 100 | 100 | 1007 (100 keys) |  | 100 | 11.6 |
| Binary Goppa | (8,192,16) | 64 | `n-5t+1` | 113 | 100 | 100 | 1007 (100 keys) |  | 100 | 62.7 |
| Binary Goppa | (8,192,16) | 64 | `n-5t-1` | 111 | 100 | 100 | 1007 (100 keys) |  | 100 | 11.4 |
| Binary Goppa | (8,241,24) | 49 | `n-5t-1` | 120 | 100 | 0 | 31 (99 keys), 32 (1 keys) | branch_separation_failed: 100 | 0 | 2.2 |
| Binary Goppa | (8,256,24) | 64 | `n-4t-1` | 159 | 100 | 100 | 751 (100 keys) |  | 100 | 10.9 |
| Binary Goppa | (8,256,24) | 64 | `n-5t+1` | 137 | 100 | 100 | 751 (100 keys) |  | 100 | 59.4 |
| Binary Goppa | (8,256,24) | 64 | `n-5t-1` | 135 | 100 | 100 | 751 (100 keys) |  | 100 | 10.7 |
| Binary Goppa | (12,922,64) | 154 | `n-5t+1` | 603 | 20 | 20 | 3812 (20 keys) |  | 20 | 884.6 |
| Binary Goppa | (12,922,64) | 154 | `n-5t-1` | 601 | 20 | 19 | 3812 (19 keys), 3813 (1 keys) | branch_separation_failed: 1 | 19 | 1175.1 |
| Binary Goppa | (13,1888,128) | 224 | `n-5t+1` | 1249 | 20 | 20 | 7113 (20 keys) |  | 20 | 5149.5 |
| Binary Goppa | (13,1888,128) | 224 | `n-5t-1` | 1247 | 20 | 20 | 7113 (20 keys) |  | 20 | 7287.4 |
| Binary alternant | (6,64,5) | 34 | `2D+1-n` | 53 | 100 | 100 | 213 (100 keys) |  | 100 | 2.1 |
| Binary alternant | (7,120,9) | 57 | `2D+1-n` | 101 | 100 | 100 | 826 (100 keys) |  | 100 | 16.0 |
| Binary alternant | (8,240,20) | 80 | `2D+1-n` | 199 | 100 | 100 | 1408 (100 keys) |  | 100 | 54.4 |
