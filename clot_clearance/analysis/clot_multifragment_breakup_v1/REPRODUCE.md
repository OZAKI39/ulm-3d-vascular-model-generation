# CLOT_MULTIFRAGMENT_BREAKUP_V1

This is an isolated extension of the latest accepted single-contact crushing implementation. It is not a new material calibration or validated bubble model. Do not rerun into an existing evidence directory.

## Environment and immutable inputs

Project root `/workspace/y_clot/project`; `source /workspace/y_clot/env.sh`. Use one rank and one thread (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`). The task is `clot_multifragment_breakup_v1`, run `20261017T110849Z`.

The baseline manifest checks all consumed historical states against the prior sealed delivery. The specimen is byte-copied from `verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z/provenance`. Its material, 1792 particles, 336 anchors, 95666 original bonds and 2736 original cuts remain unchanged. Historical `analysis/clot_multifragment_crushing_mechanism_v1/indenter.py`, `analysis/v3_softening_pilot_v1/native.py`, `analysis/softening_feasibility_v1/model.py` and `analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py` are required, with exact hashes in the manifests.

The isolated library is `verification/clot_multifragment_breakup_v1/20261017T110849Z/build/liblammps.so`. `build.py` compiles only the changed native force coordinator and reuses hash-verified previous objects. `BUILD_AUDIT.json` records compiler commands, object hashes and the new library hash. Prior libraries are untouched. To build on another machine, restore the recorded previous source/build dependency chain first; the small delivered source directory alone is not a standalone LAMMPS distribution.

## Original execution order

1. `prepare_baseline.py`: verify/read old evidence and freeze baseline. Refuses to overwrite.
2. `offline_crack_bridge_audit.py`: bounded historical diagnostics only. The corrected heterogeneous onset field is required; the initial erroneous threshold comparison is retained as a bookkeeping correction.
3. `multifragment_objective.py`, then `distributed_multisite_loading.py`: freeze independent objectives, external sites and both schedules before either trajectory. Refuse to overwrite.
4. `build.py`, `check_multisite.py`: isolated coordinator and direct interface checks. The documented failed stationary-center-crossing fixture is retained; successful checks use a short, smooth contact-entry interval and do not alter formal kinematics.
5. `manage.py simultaneous`, then after completion `manage.py staggered`: isolated foreground Supervisor jobs. Existing case histories/configurations are protected; do not invoke these against a completed run.

For an explicitly requested reproduction, create a new run identifier and output directory, preserve the same frozen specimen, objective, sites and schedule, and point the isolated library path there. Do not overwrite this run or silently retune parameters. The published two formal cases are the entire authorized computation for this task.

For direct invocation in the correctly isolated environment:

```bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_multifragment_breakup_v1/20261017T110849Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
python analysis/clot_multifragment_breakup_v1/run_multifragment_cases.py --name simultaneous
```

This command intentionally refuses to overwrite the already completed simultaneous history.

## Read-only reduction and validation

After both cases complete, these commands reconstruct analysis from actual saved states; they do not perform a new dynamics run:

```bash
python analysis/clot_multifragment_breakup_v1/fragment_topology_analysis.py
python analysis/clot_multifragment_breakup_v1/validate_outputs.py
python analysis/clot_multifragment_breakup_v1/figures.py
python analysis/clot_multifragment_breakup_v1/report.py
```

Re-running these writes derived artifacts and invalidates an existing delivery manifest; use a copied output directory if preserving the sealed delivery. Analysis reuses the prior reviewed failed-edge lineage/local-affine functions by AST extraction, without importing old run constants. Pinned cuts are diagnostics only and never fed back into dynamics. `validate_outputs.py` replays original-ID irreversible damage, recomputes topology whenever a full failure occurs, compares every accepted component membership, and checks source immutability, actual contact, key-state shape and energy bookkeeping. It does not repeat unrelated material/MPI/restart suites.

In this execution, the first completed trajectory was validated while the second frozen case ran: `validate_outputs.py --case simultaneous`, then `--case staggered`. Each proof records the validator hash and every consumed raw-evidence hash. The second invocation combines the two unchanged proofs; the combined status is PASS only when both cases are validated. A PARTIAL_ONE_CASE receipt is not a complete delivery. Postprocessing scripts can be finished after a trajectory starts; the formal runtime sources and frozen input protocols remain byte-identical, and any changed nonruntime files are disclosed in the proof.

The driver uses 0.25 µs base dt, a 0.03125 µs predictive damage subdivision floor, Δd≤0.1 and 1% normalized energy residual admission. A 1e−14 s end-time roundoff guard avoids a physically negligible extra terminal step; actual times are retained. Work uses accepted piecewise-linear particle/actuator paths and contact-entry/exit quadrature. Recorded 4-vs-8 point errors and ledger residuals bound the observed numerical error, not arbitrary long sphere-center-crossing paths.

## Delivery

The WSL-only `update_wsl_document.py` requires validated results, preserves all prior `clot_rupture.md` bytes, and writes before/after hashes plus a unified diff. Run it once after syncing the derived results, then sync its receipt and document backups to the isolated server evidence directory.

`seal.py` requires validation, a successful WSL update and both named Supervisor jobs EXITED. It records dependencies, environment, source diff and hashes. `local_readback.py` uses only the Python standard library to validate the selectively synchronized local delivery, HTML links, figures, 18 answers, 12 final status fields and the live WSL document hash.

Sync source, necessary FP64 states, compressed damage/component histories, scalar CSV, logs, figures and manifests. Exclude reproducible build `.o` and `liblammps.so*`; retain their hashes and exact build audit on both machines. No GitHub push or unrelated large dataset transfer is part of this task.
