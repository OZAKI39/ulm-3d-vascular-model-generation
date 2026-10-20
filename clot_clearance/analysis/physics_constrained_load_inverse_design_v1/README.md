# PHYSICS_CONSTRAINED_CLOT_LOAD_INVERSE_DESIGN_V1

Active isolated run: `20261016T144922Z`. The accepted native library and previous
experimental directories are read-only. This prototype optimizes a **prescribed
unsteady plane Stokes flow**, used as a low-fidelity surrogate for streaming.
It does not calculate acoustic streaming from microbubble dynamics.

Responsibilities:

| Module | Executed responsibility |
|---|---|
| `traction_library.py` | Solve diffusion modes for moving-plane Couette startup; independent finite-difference stress check; real V3 face template |
| `traction_to_v3.py` | Original surface support, affine scalar weights, transpose power interpolation; callable frozen SI node-force provider and high-cost admission gate |
| `inverse_optimizer.py` | Optuna SQLite ask/subprocess/tell, persisted sampler, constrained BoTorch qLogEI proposals |
| `mechanical_evaluator.py` | Existing precracked specimen and unchanged native variational LPS; traction callback, event forecast, original topology and energy ledger |
| `mevo_feasibility.py` | Read actual historical k-Wave, current F10S phasor/trajectory/APECSS data; bounded feasibility classification |
| `run_inverse_design.py` | Freeze and staged resumable execution |
| `analyze.py` | Read actual outputs; report, scientific figures and final statuses |

The full V3 native run is deliberately gated. Its frozen load provider is
callable, but **a new full V3 dynamics execution is not validated in this round**.
Small-specimen observations must never be passed off as V3 chunk predictions.

Reproduction commands and exact installed dependencies are in the isolated
verification directory. Running with `--limit 4 --next` resumes the same study;
it does not erase failures. A failed initialization was retained, then its
identical candidate was retried through the native admitted `v3forces` port.

All results remain exploratory: material Gc, interface parameters, finite-size
convergence, and acoustic-to-streaming conversion are uncalibrated.
