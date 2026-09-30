# Existing clot_clearance: failure and transport extension

This extends the existing prototype in this directory. The original NOSB-PD implementation, mild-case configuration, results, documentation and visualization are preserved. Start with [OPEN_FRAGMENTATION.html](OPEN_FRAGMENTATION.html) for the new movies, figures and numerical report; [OPEN_RESULTS.html](OPEN_RESULTS.html) remains the original mild-damage entry point.

The new **qualitative forced-failure verification; UNCALIBRATED** case demonstrates prescribed streaming-like loading → PD deformation and empirical cyclic damage → bond failure → detachment → fragment transport. It does not establish correct ultrasound frequency response, microbubble shape oscillations, clinical ablation rates, experimental fatigue lifetime or quantitative fragment sizes. No new BraVa CFD or population calculation was launched.

## Two separate demos

| Case | Configuration | Completed output | Meaning |
|---|---|---|---|
| `straight_pipe_damage_demo` | `configs/straight_pipe_damage_demo.json` (byte-identical copy of the original) | `results/straight_pipe_damage_demo` → `straight_pipe_demo` | Original 8,000-cycle mild damage, zero broken bonds and zero detachment |
| `streaming_fragmentation_demo` | `configs/streaming_fragmentation_demo.json` | `results/streaming_fragmentation_demo` → `../runs/fragmentation_pilot_003` | 25,000-cycle uncalibrated failure/transport verification |

The mild case was also rerun in `verification/fragmentation_extension/mild_reproduction`; every saved numerical array is bitwise equal to the original. Its internal legacy configuration name is intentionally unchanged. The completed fragmentation result is pilot 003, not the earlier displacement-guard failure in pilot 002.

## Build, test and reproduce

Use the existing shared environment. All commands below run from `/home/lzy/projects/clot_clearance`. They write only into new output/build directories; runners refuse an existing output directory.

```bash
cd /home/lzy/projects/clot_clearance
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD"
export NUMBA_CACHE_DIR="$PWD/build_fragmentation_replay/numba_cache"
export MPLCONFIGDIR="$PWD/build_fragmentation_replay/matplotlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
PD_PYTHON=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python
cmake -S . -B build_fragmentation_replay -DPython3_EXECUTABLE="$PD_PYTHON"
cmake --build build_fragmentation_replay
ctest --test-dir build_fragmentation_replay --output-on-failure
"$PD_PYTHON" -B -m pd_clot.runner --config configs/straight_pipe_damage_demo.json --output runs/mild_replay
"$PD_PYTHON" -B -m pd_clot.fragmentation --config configs/streaming_fragmentation_demo.json --output runs/fragmentation_replay
"$PD_PYTHON" -B scripts/audit_fragmentation.py --run runs/fragmentation_replay
"$PD_PYTHON" -B scripts/prove_fragment_transport.py --run runs/fragmentation_replay
"$PD_PYTHON" -B scripts/render_fragmentation.py --run runs/fragmentation_replay --output visualization/fragmentation_replay
```

The delivered CTest run passed all **26 pytest tests**, including the original 18. The selected simulation took approximately 192 seconds on the current CPU after compilation. Compilation and rendering add time. Existing build/test logs are in `logs/fragmentation/`; reproducibility is tied to `results/streaming_fragmentation_demo/IDENTITY.json` and the matching `provenance/fragmentation_extension/source_v1/` source copy.

## Extension boundaries and outputs

New modules `fragmentation.py`, `fragment_topology.py`, `fragment_fluid.py`, `fragment_mechanics.py` use the original geometry, analytical traction provider and full-rank NOSB-PD kernel. New tests and scripts are separate files. No artificial cleavage utility is imported by the new runner, and no particle is removed after failure or clearance.

All 360 particles remain simulated. Detached fragments retain their surviving PD bonds. At late time, however, 172 particles have no active bond, and 38 other particles have only rank-one or rank-two support. Those particles use an explicitly documented intrinsic-support approximation; the 172 isolated particles have zero internal bond force and still undergo fluid-driven motion. This is a significant limitation, not evidence of quantitatively valid three-dimensional clot fragments. See [FRAGMENTATION_MODEL.md](FRAGMENTATION_MODEL.md).

The result directory includes `states.npz`, `history.csv`, `history.json`, `components.json`, `lineage.json`, `failure_ledger.jsonl`, `clearance_ledger.json`, `particles.pvd` and `bonds.pvd`. Open the PVD files in ParaView for all 26 saved states. Every particle VTP contains damage, stable fragment ID, anchor reachability, fluid-exposure mask, velocity, displacement and deformation rank. Bond VTP stores damage, integrity and active status. The failure ledger records each bond's measured cycle amplitude and before/after damage.

The primary movies use actual positions, displacement factor 1, equal spatial units and a fixed camera. They show every particle, an active-bond spanning forest, a bounded sample of recently broken bonds, fixed anchors, the prescribed bubble marker and the clearance plane. The bubble marker is enlarged only for visibility. Both movies and GIFs are in `visualization/streaming_fragmentation_demo/`. Figures are supplied as PNG and PDF.

The complete Chinese numerical report and parameter changes are in [FRAGMENTATION_REPORT_ZH.md](FRAGMENTATION_REPORT_ZH.md) and [CONFIG_DIFF.json](provenance/fragmentation_extension/CONFIG_DIFF.json). Final delivery checks are recorded separately in [FINAL_QC.json](provenance/fragmentation_extension/FINAL_QC.json).

The next scientific step remains a resolved or experimentally validated microbubble-generated flow field and validated fluid-to-solid forcing. The current prescribed field, fatigue law, surface quadrature and relaxation time are verification surrogates.
