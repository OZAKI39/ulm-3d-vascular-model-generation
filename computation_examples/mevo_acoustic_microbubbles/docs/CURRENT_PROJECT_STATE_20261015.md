# Current acoustic–microbubble project — F10S, 2026-10-15

Complete scientific delivery: `31a44d3f342e5d0bca55f7f9b83f4d32b8189696`.
Source worktree: `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10s_20261015/computation_examples/mevo_acoustic_microbubbles`.
Curated GitHub branch: `sync/mevo-acoustic-f10s-20261015`, based on the published F10Q synchronization `467b17411deeb13a77370565b7382f071b41a774`.

This supersedes the F10Q handoff for current development. Historical reports and failures remain unchanged. This synchronization packages data only; no scientific computation or hardware communication is run.

## Model and active code

Four independently controlled 7 mm PZT discs on X±/Y± faces of a 32.5×9×9 mm block; central 1.6 mm horizontal channel, two 1.5 mm top openings, both X ends closed. Accessible cavities are water-filled; background velocity is strictly zero. Acoustic model: `HOMOGENEOUS_WATER_CONTROL_PROXY`, 2.5 MHz. The 32 fixed-ID bubbles have radius 1.78 µm; the connected initial geometry is a prescribed demonstration, not acoustically formed in these runs.

The original APECSS response, F2 secondary-force hypothesis, pair lubrication, Sequential QP hard-sphere engine and GPU wall/force queries are retained. Source amplitude is Pa-equivalent, not calibrated PZT voltage or measured face pressure.

- [F10S active code map](F10S_ACTIVE_CODE_MAP.md), [complete Chinese report](FAST_PZT_F10S_DEVELOPMENT.md).
- [X pair 3 ms](F10S_X_PAIR_3MS.md), [independent X sources](F10S_INDEPENDENT_X_SOURCES.md), [phase-force stability](F10S_PHASE_TRAP_FEASIBILITY.md), [Y/PDMS scope](F10S_Y_PDMS_SCOPE.md).
- [F10R report](FAST_PZT_F10R_DEVELOPMENT.md) and [scheduled-control protocol](F10R_SIMULATION_CONTROL_PROTOCOL.md).
- Current orchestration: `scripts/straight_tube_exp1/f10s/`; preregistration: `configs/straight_tube_exp1/f10s/protocol.json`.
- Control interface: `src/mevo_mb/control/scheduled_pzt.py`; frozen physics: `src/mevo_mb/fast_sim/water_sequential.py` and `sequential_hard_sphere.py`.
- Protected original geometry, response and coarse/fine acoustic caches: `results/straight_tube_exp1/f10n/`.

## Actual outcomes and limits

F10S freshly integrates the full-precision F10R X pair to 3 ms: ΔX=0.197999447 µm, error=0.002000553 µm to the fixed (0.2,0,0) µm target. Original 0.075 µm and new 0.01 µm endpoint criteria both pass. Full-window dt/cache COM, shape, topology and fixed-ID checks pass; no tolerances were loosened.

Independent X− and X+ single sources (400 Pa-equivalent) move the identical initial group by ±0.015835804 µm in 2 ms; dt/cache checks pass. Rg decreases only ~0.29518 nm and aspect ratio barely changes. This demonstrates finite proxy-model translation, not arbitrary shape control or physical fusion. Eleven F10S trajectories complete their registered windows; longest physical time is 3 ms.

Six phase candidates have local restoring single-bubble force equilibria in the proxy. Pressure peaks and force zeros can differ by ~18.94 nm. Frozen internal-shape COM diagnostics predict X relaxation times ~0.368–2.473 s; these are not full particle trajectories. High-peak secondary forces increase compression risk. Long moving-trap trajectories were not run.

Y conclusions are `HOMOGENEOUS_WATER_PROXY_SPECIFIC`: the near-field zero is highly geometry/wavelength-sensitive. Conditional PDMS grids require ~60–110 GiB and lack validated 2.5 MHz material/coupling inputs; `PDMS_Y_FIELD_NOT_ADMITTED`, `PDMS_FIELD_NOT_RUN`. Historical F10Q Y shape/ID failures and strong-control rejection remain valid; F10R weak Y misses its target. No independent arbitrary Z navigation or experimental calibration is established.

Retain **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false**. No AFG1062, amplifier or Arduino communication has been developed in F10R/S.

## Evidence, environment and continuation

`results/straight_tube_exp1/f10r/` and `f10s/` retain original controls, all returned server files, trajectories, pressure/force/contact evidence, numerical failures/NOT_RUN records, source execution archives, final PNG/PDFs and visual reviews. [Current Chinese progress](MICROBUBBLE_AGGREGATION_PROGRESS_20261015.md) is a byte-identical copy of the live WSL document; historical before/after/diff audits remain available.

Existing server environment: `source /workspace/mevo_acoustic_microbubbles/activate.sh`. Dependency declarations and deployment evidence are inherited. Recorded server state is historical; a stale RUNNING field is not an instruction to restart anything. Future scientific runs require a new output directory and Supervisor label; never overwrite archived results. Prior development commit IDs are provenance identities, not necessarily GitHub-reachable commits; included source archives preserve the executed code.

Recommended next research: bounded X-distance extension and push/OFF or feedback stopping, with unchanged numerical budgets. Moving traps, complex shape control, Y-field improvements and hardware work require separate scientific development.

Packaging verification: `python scripts/verify_sync_20261015.py` (Python, NumPy, h5py; no GPU or solver call). See [scope, exclusions and hashes](../sync/20261015/README.md). Original full-delivery validators can reference deliberately omitted previews; use this package verifier for the curated Git tree.
