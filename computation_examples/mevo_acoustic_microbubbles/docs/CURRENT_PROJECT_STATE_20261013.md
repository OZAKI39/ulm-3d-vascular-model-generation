# Current acoustic–microbubble handoff — 2026-10-13

The active scientific delivery is F10K, completed on 2026-10-12 at
`850fd2df3225351f4fe9564f39cbe72a38662b30`, branch
`dev/fast-pzt-f10k-science-scope-scaling-20261012` in the local complete archive.
The GitHub branch for this curated synchronization is
`sync/mevo-acoustic-f10k-20261013`, based on the previously published
`c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0` F5 synchronization.

## What is current

The complete F6–F10K source/configuration/report changes are copied without
changing scientific code or results. Older Stage 4, Stage 5 and F1–F5 remain
historical evidence; their stage names, old hosts or archived RUNNING strings
are not instructions to resume computations.

Read [F10K development](FAST_PZT_F10K_DEVELOPMENT.md),
[scientific scope](FAST_PZT_F10K_SCIENTIFIC_SCOPE.md),
[bounded scaling](FAST_PZT_F10K_SCALING.md) and
[FINAL_GATE](../results/fast_pzt/f10k/analysis/FINAL_GATE.json).

- F5's 417.5 ms state is geometrically admissible within roundoff, but near-wall
  B2 is not admitted. It was not moved elsewhere to force a splitting result.
- A real interior F10D state initialized two new B2 counterfactual windows,
  dt10/5 µs and duration3 ms. Both show one31+1→32 graph merge at0.290 ms,
  persisting2.710 ms. This is not physical bubble fusion. Persistent splitting
  remains untested; instantaneous edge timing remains FAIL.
- Actual F10A saved trajectories cross the20 µm downstream plane at≈2.024 ms.
  O−OFF transverse increments are≈20.32/20.53 nm, below the original50 nm gate.
  The proposed fixed target is preregistered but not run; navigation is not established.
- Converged fixed-state N32/128/512 comparisons are retained. GPU advantages
  at128/512 apply to those regular synthetic networks, not arbitrary dense groups.
  N1024 was NOT_RUN because the predefined CPU-array memory budget was exceeded.
- Physical calibration and AI-scientific-data admission remain absent:
  **EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false**.

## Environments and paths

Local complete source archive:
`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10k_20261012/computation_examples/mevo_acoustic_microbubbles`.
This synchronization's separate worktree:
`/home/lzy/projects/computation_examples/_archives/mevo_acoustic_sync_20261013`.

The last verified F10K deployment is `/workspace/mevo_acoustic_microbubbles/f10k`
on the existing `mevo-new` server, using its parent's isolated `.venv`.
[F10K ENVIRONMENT.json](../results/fast_pzt/f10k/ENVIRONMENT.json) records the
actual runtime: Python3.12.3, Warp1.18.0+cu12, NumPy2.2.6, SciPy1.15.3 and
threadpoolctl3.6.0. The latter was added during F10K; other environment/driver
versions were not upgraded. Original compute/render locks and migration evidence
remain under `deployment/newserver_20261009/`. This sync did not alter environments.

The canonical original geometry/acoustic/CFD inputs are still under
`results/fast_pzt/f5/inputs/`; the33 MB acoustic basis is inherited from the
published parent and not uploaded as another copy. Server status files here are
dated records, not a live status service. Original F10K jobs had all exited at
the archived final verification; no job was started for this synchronization.

## Git snapshot versus full scientific archive

The copied compact results, reports, selected original trajectories and plots are
byte-identical to the complete source delivery. Every included and excluded delta
file has a source path, size, original Git blob ID and SHA-256 in
`sync/20261013/`. The five source-only F10K execution tarballs preserve exact
runtime source bytes without making the large local result history reachable
from the new branch. They include the two failed startup attempts' source commit.

Large repeated batch trajectories, full per-step solver matrices/observer traces,
raw radial waveforms and large derived arrays remain in the complete local archive.
They are not deleted. Existing full-delivery `ARTIFACT_MANIFEST.json`,
`VERIFICATION.json` and scientific gate files are preserved unchanged and may
reference omitted files. Their original PASS statements describe the original
complete delivery, not an assertion that this curated clone contains every file.

Use `python scripts/verify_sync_20261013.py` for this curated tree. Full-original
stage verification and some historical plot regeneration require restoring the
listed omitted data by hash. The F10K figure inputs and selected original HDF5
trajectories are included; no claim is made that all F6–F10J full campaigns can be
rerun from this subset. Do not modify historical scientific validators to conceal
these deliberate omissions.

The live Chinese progress document is copied byte-for-byte to
[MICROBUBBLE_AGGREGATION_PROGRESS_20261013.md](MICROBUBBLE_AGGREGATION_PROGRESS_20261013.md).
Its absolute historical paths remain provenance, not portable clone paths.
Next scientific priorities remain fixed-target control authority, a legal interior
persistent-split example, and strict edge-time sensitivity before production scaling.
