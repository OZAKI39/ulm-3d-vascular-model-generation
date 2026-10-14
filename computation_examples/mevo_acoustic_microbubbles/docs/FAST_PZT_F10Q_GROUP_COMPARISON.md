# F10Q three-model common-window comparison

All controls and initial coordinates match bitwise. Only saved times 0 and 0.1 ms are exactly shared within the 0–0.174 ms legal intersection. No interpolation is used. Different contact and mobility approximations are compared, not reinterpreted as integration error.

| Pair | COM difference (nm) | Maximum ID difference (nm) | Rg difference (nm) | Relative Q difference |
|---|---:|---:|---:|---:|
| M0_VS_B2 | 1.24266805e-05 | 15.9548594 | 13.2558255 | 0.0045088271 |
| M0_VS_HARD | 5.80284105e-06 | 3.76645935 | 2.50914486 | 0.000928184264 |
| B2_VS_HARD | 1.82288379e-05 | 18.4450986 | 15.7649703 | 0.00545288315 |

The largest cross-model COM contrast is 1.822883789e−5 nm, while Rg differs by up to 15.76497033 nm. This supports limited translation consistency for this initial state, this control and these two saved states only. It does not establish general contact-model independence or any 2/30 ms agreement. Shape is model-dependent.

Within-model dt and cache differences at exactly the same saved times are separately stored in THREE_MODEL_COMPARISON.json. The original F10N full-window dt COM error is not used as the cross-model difference. Zero-flow OFF is exactly stationary in all inherited records and the new independent three-way OFF repeats; directional acoustic displacement equals displacement relative to that null.

Model-specific contact quality remains separate: M0 uses its original PBD budget, B2 its original soft-overlap budget, and hard spheres their own complementarity tests. There is no hard-sphere reaction complementarity metric for M0 or B2. All three remain exploratory and uncalibrated.
