# OpenPIV validation

Installed OpenPIV 0.25.3 (GPLv3), NumPy/SciPy backend. Rust API absent; no Rust timing claim. Single pass: 48 px windows, 24 overlap. Official two-pass: 64→32 px, 32→16 overlap; peak-to-peak threshold 1.2, original texture SD >=2, local median outlier rejection. Invalid vectors may be replaced inside the official deformation predictor but are excluded from published final vectors and fitting. Coarse flags are used internally and final invalid-vector flags are retained; final multipass numerical SNR is not exposed by the official wrapper and is not invented. Current adapter also exposes coarse-pass SNR.

Development seed 3101 compared raw, local Gaussian-background subtraction, causal median fixed-background suppression (up to 5 prior frames + current), candidate enhancement and official multipass. Selected `multipass` using the preregistered error/availability rule. The first engineering masked-array failure was fixed and rerun on development only, preserved as `benchmark_development_api_failure.json`. Test thresholds unchanged.

| Method | Endpoint tests | All-six output fraction | Accepted failures | Failure / accepted truth outputs | Non-affine correct rejection |
| --- | --- | --- | --- | --- | --- |
| PIV_multipass | 216 | 0.6991 | 21 | 0.1429 | 0.6667 |
| CPD | 216 | 0.5417 | 12 | 0.1081 | 0.5 |

| Parameter | MAE | 95th abs error | Sign fraction (nonzero truth) | Nominal 95% interval coverage |
| --- | --- | --- | --- | --- |
| translation_vx | 1.668 | 11.51 | 1 | 0.6259 |
| translation_vy | 0.9699 | 7.528 | 1 | 0.7279 |
| divergence | 0.02002 | 0.06067 | 0.9455 | 0.4218 |
| vorticity | 0.06342 | 0.1723 | 0.974 | 0.3741 |
| shear1 | 0.007145 | 0.03034 | 1 | 0.7007 |
| shear2 | 0.005327 | 0.03294 | 1 | 0.7687 |

These overall errors include difficult cases and only accepted truth-bearing outputs. 4 accepted local-opposite cases had no single affine truth and are counted separately; false accepted estimates were not discarded. A small residual does not prove material motion: stationary background may produce a coherent but wrong fit. Fully raw images failed the chosen correlation-quality gate; preprocessing is necessary for this recording and these images. Candidate and fixed-background versions sometimes suppress slow-moving objects too.

Clean-case mean absolute errors, six parameters in order:

| Case | MAE vx, vy, delta, omega, shear1, shear2 |
| --- | --- |
| translation | 0.2255, 0.2428, 0.002623, 0.00186, 0.001458, 0.001278 |
| rotation_cw | 0.2227, 0.1985, 0.03065, 0.03457, 0.00395, 0.002589 |
| rotation_ccw | 0.1331, 0.2006, 0.03193, 0.03344, 0.002827, 0.003095 |
| contraction | 0.1467, 0.2212, 0.01726, 0.003252, 0.002226, 0.002718 |
| expansion | 0.08291, 0.1427, 0.03314, 0.004194, 0.001768, 0.002238 |
| stretch | 0.1216, 0.1379, 0.01442, 0.003583, 0.0101, 0.001199 |
| shear | 0.2141, 0.06504, 0.003774, 0.01653, 0.001668, 0.007907 |

Fixed PDMS-like retained background, dropout, sparse and locally opposing motions reveal failures. Sparse cases are not universally mathematically unobservable: enough accidental window texture can still pass support tests, which is a limitation, not a reason to mark them successful. The supplementary local-defocus/distractor checks are post-benchmark diagnostics and did not change gates. Nominal interval coverage fails calibration; do not feed these values into a controller as certified measurements.
