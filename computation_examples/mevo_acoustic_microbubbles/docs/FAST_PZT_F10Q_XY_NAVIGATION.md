# F10Q independent XY navigation

A failure of the old strong control did not block B. Eight newly optimized controls were registered before simulation and independently completed 0.5 ms. The proposal optimizer reused the original cubic water-field cache, applied four independent complex source weights, checked actual phasor pressure, and used only a frozen-state COM predictor to propose controls. All success claims use actual 32-particle integration.

Exactly two controls were selected from measured short-window displacement and a documented pressure-risk cost: X_P1200 and Y_P1800. Each was freshly run to 2 ms at dt=1 µs, repeated at 0.5 µs, and repeated with 1.5 µm cache spacing instead of 3 µm. OFF was independently repeated in all three settings.

| Control | Final COM (µm) | Assigned error (µm) | Basic | COM | Shape | Graph | ID | Target |
|---|---|---:|---|---|---|---|---|---|
| OFF | (0, 0, 0) | 0.2 | PASS | PASS | PASS | PASS | PASS | False |
| X | (0.200896324, -1.09787383e-06, 1.25731453e-16) | 0.000896325125 | PASS | PASS | PASS | PASS | PASS | True |
| Y | (-2.96732256e-08, 0.051776285, 1.21430114e-14) | 0.148223715 | PASS | PASS | FAIL | PASS | FAIL | False |

Original targets remain X=(0.2,0,0) µm and Y=(0,0.2,0) µm, tolerance 0.075 µm. X reaches its original target in this same-model numerical prototype. Y does not. These are no-flow, homogeneous-water numerical controls, not calibrated experimental results.

## Sensitivity and shape

| Control / comparison | COM (nm) | ID max (nm) | Rg (nm) | Relative Q | Centered matched RMS/R |
|---|---:|---:|---:|---:|---:|
| NAV_X_DT05 | 0.0002255681 | 0.0008928674 | 8.339666e-07 | 3.903443e-09 | 3.514164e-07 |
| NAV_X_FINE | 0.005704564 | 0.005752851 | 2.817215e-07 | 1.276924e-09 | 1.708313e-08 |
| NAV_Y_DT05 | 0.0001001698 | 395.6084 | 96.41726 | 0.04042844 | 0.1377272 |
| NAV_Y_FINE | 0.003639943 | 764.6402 | 306.7297 | 0.1216458 | 0.3149122 |

X has a resolved but tiny Rg change of -0.29529204 nm; its X/Y ratio ends at 1.00000008337, so this is nearly pure translation, not substantial shape control. Its leading in-plane eigenvalues remain degenerate and a unique principal direction is NOT_IDENTIFIABLE.
Y coarse Rg changes by -127.979772 nm and X/Y ends at 0.992642044654, but shape differences of 96.4/306.7 nm and ID differences of 395.6/764.6 nm exceed original limits. This is an actual GROUP_POSITION pass with GROUP_SHAPE and PARTICLE_ID failures, not a synthetic metadata fixture. All initial 32-member graphs survive 2 ms; no physical fusion is implied.

X pressure is 690.39–1113.33 Pa and Y pressure 1742.13–1801.09 Pa across actual sampled steps. Initial proposal caps are not hard runtime amplitude clipping; actual pressure is re-evaluated at every step against the unchanged APECSS range. Source amplitude is not bubble pressure. All source controls and per-step pressure extrema are archived.

## Exact selected controls

Order: Y−, Y+, X−, X+. Frequency 2.5 MHz; all sources enabled; duration 2 ms. Full precision is in NAVIGATION_SELECTION.json and each INPUT.json.

| Control | Source amplitudes (Pa) | Phases (rad) |
|---|---|---|
| X | 638.918244486, 698.829667169, 1499.78201485, 1203.15906757 | -0.630736763718, -0.616128204014, -1.24064667724, 2.17225806568 |
| Y | 1499.99999916, 1499.9999997, 477.335749145, 477.336068802 | 2.54407616313, -0.598297263499, -0.544445470398, -0.54444548106 |

No two-segment/5/10/30 ms experiment was needed: X already reaches the short target, while Y shape/ID convergence should be repaired before claiming dependable shape control. Group COM admission would permit considering a longer experiment, but it is not a requirement or evidence of morphology validity. No budget was spent merely to fill the allowance.

CONTACT_MODEL_DEPENDENT is established for historical shape contrasts; cross-model dependence of these NEW controls is NOT_VERIFIED because they were not rerun with M0/B2. EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED, TRAINING_READY=false.
