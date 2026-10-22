# F7 execution and readout

The F7 branch adds files only; F1–F6 defaults/source/results are unchanged.
Use an isolated worktree from F6 a8324747, not an older remote F5 branch.

Remote compute uses the existing `.venv`, Warp1.18+cu12, RTX4090, one supervisor
science job at a time. `remote.py` accepts only committed source and unique labels.
Completed supervisor jobs are EXITED (supervisorctl may return3); inspect the JSON
status and stderr rather than treating EXITED itself as failure.

Original executed sequence (with the per-stage commits in the report):

```bash
python3 -m scripts.fast_pzt_f7.remote fixed --label fixed_v1 --timeout 180
python3 -m scripts.fast_pzt_f7.remote replay --label pilot_v1 --timeout 180
python3 -m scripts.fast_pzt_f7.remote replay --label replay_v1 --timeout 300 --full
python3 -m scripts.fast_pzt_f7.remote control --label predictor_v1 --phase check
python3 -m scripts.fast_pzt_f7.remote control --label controller_pilot_v1 --phase pilot
python3 -m scripts.fast_pzt_f7.remote control --label controller_v1 --phase formal
# First feedback failure retained; a new commit/protocol then enabled SLSQP + pressure constraint.
python3 -m scripts.fast_pzt_f7.remote control --label controller_pilot_v2 --phase pilot
python3 -m scripts.fast_pzt_f7.remote control --label controller_v2 --phase formal
python3 -m scripts.fast_pzt_f7.remote control --label controller_half_v1 --phase half --dt 0.000005
python3 -m scripts.fast_pzt_f7.remote reachability --label reachability_v1
python3 -m scripts.fast_pzt_f7.remote control --label slow_v1 --phase slow --dt 0.0000025 --timeout 360
```

Do not rerun these names or claim the final protocol reproduces the first failed
policy. To reproduce history, check out its exact execution commit and run in a
new isolated result root; overwrite protection is intentional. The pipeline uses
named prior-admission gates; it is a bounded registered study, not a random-search
campaign scheduler. Deployment `--deploy-only` starts no computation.

Closed-data archival and local plots:

```bash
PYTHONPATH=src:. /path/to/python -m scripts.fast_pzt_f7.archive
PYTHONPATH=src:. /path/to/python -m scripts.fast_pzt_f7.analyze
PYTHONPATH=src:. /path/to/python -m scripts.fast_pzt_f7.figures
```

Archive verifies server file bytes and SHA256 and recorded Python source blobs
against each actual Git execution commit. Inspect all four PNGs and PDFs and
record review separately; generation alone is not visual acceptance.

After reviewing the four actual outputs, run `python -m scripts.fast_pzt_f7.finalize`
with the same `PYTHONPATH=src:.`. It verifies the existing review hashes, scientific
archive, original inputs, unchanged inherited sources and applied narrative update;
then closes delivery gates and refreshes the final artifact manifest. Re-running
analysis/figures resets delivery status and requires a new matching review before
finalization. The narrative update currently targets its explicitly recorded WSL path.

`SUMMARY.json` and HDF5 `quality.summary_json` retain navigation and all controls.
`quality/raw_F1_diagonal` is cumulative per saved frame: columns3/4/6 are accepted
overlap, summed contact correction, maximum contact-pass correction. The inherited
`cumulative_free_Stokes_displacement_m` key denotes actual free RPY path for M2;
read the explicit `free_path_definition`/mobility metadata, not that old name.
Instantaneous predictor velocities are pre-PBD; compare them with next-segment
realized centroid velocities in `analysis/plot_data.npz`.

Controller v2 enforces100Pa only at observed positions and the selected event;
this is not a continuous pressure lower-bound guarantee. SLSQP is limited to20
iterations; the legacy `optimizer_maxfun` field is used by v1 L-BFGS-B only.
All accepted and invalid-candidate counts/timings are preserved, with no hidden
Stokes fallback. V2 changes prediction, optimizer and excitation policy together;
do not label the old/new closed-loop result a single-factor mobility ablation.

For AI use: keep all paired no-flow/affine trajectories in the same group, including
F6 and failed variants. The CFD family is separate. Slow flow is an independent
physical validation *condition*, not an independent training split/initial family.
The inherited sample reader returns full future control history and end-of-run
quality for audit; never feed those future values as causal input features.
Connectivity is geometric, not coalescence; retain EXPLORATORY labels.

Slow-flow recovery records are distinct:

- `slow_v1`: observer-count wrapper exception; original attempt has no saved HDF5.
- `failure_record_v1`: injected exception check of matched-prefix preservation.
- `slow_v2`: unchanged science/control replay with true0–699.5ms prefix and700ms failure coordinates.
- `feasibility_v1`: exact failed-state check of one analytic source-amplitude projection, only within the original cap and only followed by a full GPU domain/force check.
- `slow_v3`: `control --phase slow --dt 0.0000025 --repair`, admitted only after that check. No added optimizer iterations or relaxed constraints.

The optional projection changes actual recorded source amplitudes, not particle
forces after evaluation. All previous failures stay distinct and retain their
execution identities. For incomplete trajectories, saved time and total integrated
physical steps can differ by the final failed observer interval; the actual failed
coordinates/time are in `FAILURE_PROPOSAL.json`, never presented as a successful
saved-frame continuation.
