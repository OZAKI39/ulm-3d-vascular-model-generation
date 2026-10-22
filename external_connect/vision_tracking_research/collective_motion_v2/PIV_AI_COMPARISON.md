# PIV versus AI: retain PIV

Frozen final test: 1,536 triples = 1,320 VALID_GLOBAL_MOTION + 54 NON_AFFINE_MOTION + 162 NOT_OBSERVABLE. Single source recording. MAE is over all finite estimates in class0, including rejected predictions; unavailable estimates are counted separately and as failure in overall global recovery. PIV has 1,114 finite class0 estimates, whereas every AI has 1,320. Units: px/s, px/s, 1/s, 1/s, 1/s, 1/s.

| Method | MAE denominator | Six-parameter MAE (vx, vy, div, vort, s1, s2) | Accepted / 1536 | Failed accepted |
|---|---:|---|---:|---:|
| PIV single pass | 1114 | 1.894, 1.698, 0.020, 0.039, 0.009, 0.010 | 1111 | 85 |
| Direct6 gray | 1320 | 11.420, 7.994, 0.057, 0.327, 0.045, 0.030 | 0 | 0 |
| Direct6 + heat | 1320 | 11.394, 8.011, 0.061, 0.332, 0.040, 0.031 | 0 | 0 |
| Direct6 heat only | 1320 | 13.416, 9.526, 0.059, 0.300, 0.044, 0.031 | 0 | 0 |
| Direct6 + heat robust | 1320 | 11.554, 7.937, 0.057, 0.288, 0.043, 0.026 | 0 | 0 |
| GeoMotion | 1320 | 10.833, 7.597, 0.057, 0.282, 0.040, 0.023 | 0 | 0 |

PIV accepted 1,111: 1,110 class0 (84 over tolerance) plus one wrongly accepted non-affine case. Thus joint accepted risk is 85/1111 = 7.65%; class0 accepted error risk is 84/1110 = 7.57%. This is the unmodified PIV quality gate, not a 5%-risk guarantee. AI had no feasible validation threshold with at least 20 accepted samples and <=5% joint error; all formal AI fields are rejected. **0/0 is undefined risk, not zero-error success.**

| Method | Error risk at 25% diagnostic coverage (384/1536) |
|---|---:|
| PIV single pass | 0.52% |
| Direct6 gray | 80.73% |
| Direct6 + heat | 86.98% |
| Direct6 heat only | 83.07% |
| Direct6 + heat robust | 79.43% |
| GeoMotion | 86.46% |

These curves retrospectively compare fixed rankings at matched coverage; they do not authorize deployment thresholds selected from test. PIV ranks within its existing accepted quality set; AI ranks algebra-valid predicted-global candidates before the selected reject-all threshold. Test rankings demonstrate that merely relaxing AI rejection would admit many errors.

Two-pass PIV and CPD were rerun only on the first 112 preregistered test samples (four material frames ×28 cases), not all1536: multipass86 accepted/7 failed; CPD41 accepted/10 failed. These subset results are never substituted for full-test denominators. F10U historical216-window results remain unchanged and separate.

| Method | Pure model p50/p95 ms | Full pipeline p50/p95/max ms | Consumed /120 | Fresh six-parameter Hz |
|---|---:|---:|---:|---:|
| PIV single pass | N/A | 34.35, 44.29, 47.32 | 114 | 25.24 |
| Direct6 gray | 2.45, 3.25 | 19.72, 25.16, 32.23 | 120 | 0.00 |
| Direct6 + heat | 3.24, 4.39 | 16.08, 21.13, 22.98 | 120 | 0.00 |
| GeoMotion | 2.94, 3.42 | 18.92, 23.56, 26.06 | 120 | 0.00 |

Identical raw ROI, actual PTS, hardware and software environment; each pipeline has its own preprocessing. AI internally resizes to128² and maps the full affine generator back to raw pixels; PIV uses raw resolution. The timing comparison includes these costs. It is a short same-recording causal replay, not a repeated-session benchmark. True camera end-to-end latency is NOT_MEASURED.
