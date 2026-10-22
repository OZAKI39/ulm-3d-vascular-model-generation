# Constrained optimization results

| Method | Observations | Complete | COM-feasible | Raw hits | First hit | Charged wall (s) | Ask/tell (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sobol | 32 | 32 | 12 | 1 | 26 | 98.200 | 0.000 |
| Ax_BoTorch | 32 | 32 | 23 | 9 | 10 | 106.082 | 10.303 |

Ax reached a discovery-feasible T1 point at observation 10, Sobol at 26. Ax produced 9 raw target hits, Sobol 1, with about 10.3 s optimizer overhead. These are discoveries on S2_r0; independent target confirmation is separately shown below.

| Target | Frozen source | Discovery aspect change | COM (nm) | Independent confirmations | Final test raw hits | Confirmed |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | extension_028 | 2.9974% | 2.1325 | 2/3 | 2/3 | False |
| T2 | ax_14 | 0.4096% | 1.2628 | 2/3 | 3/3 | False |
| T3 | baseline_OFF | 0.0000% | 0.0000 | 0/3 | 0/3 | False |
| T4 | recovery_05 | 3.0208% | 0.3866 | 0/3 | 0/3 | False |

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
