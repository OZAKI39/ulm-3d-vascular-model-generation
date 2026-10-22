# Common curved geometry validation

See [中文总报告](COMMON_CURVED_GEOMETRY_MAPPING_REVIEW_ZH.md) and [figures](OPEN_RESULTS.html).

Volume patch: validated in the native element harness. Overall: COMMON_GEOMETRY_FIX_NOT_EVALUABLE; Phase F weak-boundary mapping is blocked. No CFD, linear solve, time step or particle run is permitted by this study.

Reproduction uses the retained production source/objects identified in `SOURCE_IDENTITY.json`, `COMMON_GEOMETRY_PATCH_IDENTITY.json` and `BOUNDARY_BUILD_IDENTITY.json`. With the existing server Python environment, run `supervise.py preflight`, `baseline`, `build`, `validate`, `boundary`, `finalize` sequentially. Each stage writes only this directory; the build has no production solver main and the capture-only LinearAlgebra.solve throws. Do not rerun stages to overwrite frozen evidence; use a separate review directory and recheck identities. `closeout.py` runs only after the supervisor has exited and TASK is absent.

`COMMON_GEOMETRY_PATCH.diff` applies to exact retained source and adds `common_geometry.h`; the patch is opt-in with `SVMP_VALIDATION_COMMON_GEOMETRY`. It is not installed or production-certified. Tracing is a separate build transformation recorded in `build.py`; original and instrumented local R/K are checked. Source copies, objects and executables stay under ignored `runtime/` on the server; tracked identities and scripts make them reconstructible.

`ATTEMPT_01_*` retains the first Hessian gate failure and exact patch identity. Final results use the corrected common-inverse Hessian implementation with the original tolerance. Historical A-D1 checks were not reexecuted. Scientific boundary FAILs are preserved separately from engineering checks that correctly detect them.
