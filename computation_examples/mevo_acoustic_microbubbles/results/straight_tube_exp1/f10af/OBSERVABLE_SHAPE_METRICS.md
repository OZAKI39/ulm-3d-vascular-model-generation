# Observable shape metrics

The exact F10Y distribution() implementation is vendored with provenance. Its numerical coordinates are micrometres, NOT calibrated camera pixels. 3D truth and XY projection are distinguished. Added: support-corrected structure factor, orientation, spacing, covered periods, occupied bands, second nonparallel peak, local components, psi4/psi6, radial/angle ring checks, histograms and rotation decomposition. A fixed 1.0 x 0.8 mm observation window is used; the structure-factor detector is not a classifier for the compact cloud. Planar rigid motion has affine vorticity = 2 Omega_xy; the requested Omega_z uses the full 3D denominator, so finite Z spread slightly changes that equality. 512 independent random null clouds had zero registered stripe/lattice positives; this finite test is not universal specificity.

## Evidence

- `METRIC_NULL_CALIBRATION.json`
- `TEST_RESULTS.json`
- `scripts/straight_tube_exp1/f10af/vendor/f10y_aggregation.py`
- `series/`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
