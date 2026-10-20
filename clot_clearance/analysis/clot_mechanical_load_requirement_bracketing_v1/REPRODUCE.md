# Reproduce the frozen mechanical campaign

Server root: `/workspace/y_clot/project`. Activate `/workspace/y_clot/env.sh`. Prepend `verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/build` to `LD_LIBRARY_PATH`. Use one MPI rank and one BLAS/OpenMP thread per case. Required native library SHA: `da705cdd2358eab259d9f99ab402e743f6a64c829959eb094c148f8492552ae7`.

The accepted run is `verification/clot_mechanical_load_requirement_bracketing_v1/20261016T202533Z`. Do not overwrite it. `DELIVERY_MANIFEST.json` and the local readback receipt identify the delivered files. `LOAD_CASES.json` indexes the actual requests, including failed attempts and the successful continuation.

## Physical protocol and numerical versions

Read `BASELINE_AND_SPECIMEN_FREEZE.json`. All amplitudes share the original 1792 material / 336 anchor specimen, material, original precut, load patch and 2 ms half-cosine rise + 0.5 ms hold. These are nominal reference-area dead tractions, not follower pressure, current-surface pure shear or validated bubble forces. Loaded original material nodes remain loaded after detachment.

- Original 1/10 Pa cases and the first failed 100 Pa attempt use `native_case.py`.
- `STARTUP_RECOVERY_FREEZE.json` records the resolved 0–8 μs startup. Its original driver is retained as `provenance/native_case_v2.before_prediction_fix.py`.
- `PREDICTOR_IDENTITY_FIX.json` records the algebraically equivalent half-velocity ordering of the independent VV predictor. The native integrator and the 2e-18 m check are unchanged. The intermediate driver is retained as `provenance/native_case_v2.before_timer_field_fix.py`.
- `REPLAY_TIMER_FIELD_FIX.json` separates native field8[7], a wall-clock profiling timer, from physical replay identity. Every remaining scalar, kinematic, force, history, graph and counter field is checked bitwise. The final driver is `native_case_v2.py`.

Each v2 case's `REQUEST.json` records its **executed driver SHA**, so choose the matching archived source for exact reproduction. The first replay preflight failed by comparing the profiling timer; its log and `FAILED_PREFLIGHT.json` remain. The successful 1000 Pa continuation replays the original 9802 accepted `run 1` calls, verifies identity, then advances to 2.5 ms. No native history or dissipation is reset.

## Replay one case without modifying accepted evidence

Create a fresh output root, for example `/workspace/y_clot/project/verification/clot_mechanical_replay/<new_run_id>`; create its `provenance` and `cases` directories. Copy the three accepted `provenance/reference.npz`, `reference.native`, and `body.data` files there byte-for-byte. Per-case relative symlinks expect this layout.

Copy the selected JSON from the accepted run's `requests/` directory to the fresh output root. Change **only its `case` path** to `<new_output_root>/cases/<case_name>`. Preserve `params`, including `replay_source` for an exact continuation reproduction. That source is read-only. These replay calls deliberately reuse the original frozen inputs and their SHA, while all outputs go to the fresh case directory.

Set `PYTHONPATH` to include `analysis/clot_mechanical_load_requirement_bracketing_v1` in addition to the existing environment. Run the source whose SHA matches the original case's `REQUEST.json`:

```bash
python /path/to/matching/native_case_source.py --request /path/to/new/request.json
```

A completed case returns its existing summary without rerunning; an interrupted directory is refused. Use a new directory for another attempt. For a genuinely new physical protocol, use a new task/run and a new freeze; do not silently modify this campaign's common constants or original material assets.

## Campaign order and postprocessing

`campaign.py` contains the initial shear decades. `campaign_recovered.py` contains the recovered 100/1000 Pa segment. `supplement.py` performs two observed-onset geometric midpoints, then independent local inward/outward comparisons at 100 and 1000 Pa. Their requests are registered serially before at most two native cases run concurrently. `chunk_midpoint.py` adds exactly one 316.227766 Pa shear case to classify the coarse no-detachment / detachment interval; it does not assume fragment size monotonicity.

Long runs use the archived supervisor configurations and wrappers. No application port or unrelated service is changed. Case resources measure native-process wall/CPU/RSS. One failed replay preflight preserves 249 s of native wall time from the LAMMPS log (1 s resolution), but its process CPU/RSS were not captured before the assertion; they are explicitly unknown, not imputed.

For the delivered run, `analyze.py`, `report.py`, and `validate_outputs.py` reduce accepted data and generate six main figure groups plus a full-domain projection. `OPEN_RESULTS.html` in the results folder links figures, CSV tables, review and status. `seal.py` runs only after native jobs have exited, the WSL document has been updated and delivery validation passes.

The full specimen reference is stored once. Accepted key states, scalar histories and compressed damage changes are retained. The successful continuation links unchanged prefix key states and uses atomic replacement for newly generated states, protecting the source. Failed candidates are never used as accepted visualization frames. Topology uses original live bonds and anchor reachability; contact proximity never creates a connected fragment. No additional native library, solver package, CFD run, material calibration or full V3 calculation was introduced.
