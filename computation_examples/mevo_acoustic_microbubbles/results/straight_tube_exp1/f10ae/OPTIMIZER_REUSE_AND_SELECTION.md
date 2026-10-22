# Open-source optimizer reuse

Official [Ax Client](https://ax.dev/docs/tutorials/getting_started/) and its [outcome constraints](https://ax.dev/docs/recipes/outcome-constraints/) are used with BoTorch, not a home-written Gaussian process or acquisition function. The comparator is [SciPy Sobol](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html). Installed versions and licenses are in `evidence/DEPENDENCIES.json`.

Two continuous parameters: X pair common amplitude 200–1100 Pa-equivalent and relative phase 150–240 degrees. Both methods receive the same 8 starting observations and 32 total observations each, representing 56 distinct new full simulations. Failed solves are marked failed, with no endpoint loss. Ax models COM_nm ≤20 as an outcome constraint.

Nevergrad is installed but NOT_RUN: two mature methods sufficed for this bounded comparison. Winner expansion beyond the frozen initial comparison is NOT_RUN; budget is reserved for recovery and independent validation. One optimizer seed does not establish general superiority.

| Method | Observations | Complete | COM-feasible | Raw hits | First hit | Charged wall (s) | Ask/tell (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sobol | 32 | 32 | 12 | 1 | 26 | 98.200 | 0.000 |
| Ax_BoTorch | 32 | 32 | 23 | 9 | 10 | 106.082 | 10.303 |

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
