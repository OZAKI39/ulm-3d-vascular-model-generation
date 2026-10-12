# F9 reproducibility

Run from the project root. Scientific GPU jobs use the existing remote Warp environment. Local analysis needs NumPy, h5py and Matplotlib; no acoustic/CFD recomputation.

```bash
export PYTHONPATH=src:.
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f9.analyze
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f9.scope_review
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f9.interface_check
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f9.figures
```

Open all regenerated PNGs/PDFs and renew VISUAL_REVIEW.json; old review hashes must not be reused. Analyze resets delivery-only gates to PENDING. Finalize runs only after archive, visual and progress-document checks.

GPU entry points: checks, wall_admission, campaign. Remote deploys committed source to `.execution/<SHA>` with unique Supervisor labels and refuses concurrent MEVO science. Use recorded IDENTITY/REGISTERED_PROTOCOL/REPRODUCTION files to reproduce actual runs. Preserve literal `1e-5`/`5e-6`: `10*1e-6` differs at the floating-point bit level. Y2_10 was conditionally admitted; Y2_5 remains unadmitted. No extra scientific run is needed for this delivery.

Archive pulls closed files and verifies server bytes and Git source identities. Server_record is read-only and saves environment and six completed process timings. Update_progress makes the requested one-time minimal WSL insertion with backup/diff/hashes. Finalize checks inherited source, history, closed archives, noninterference, figure/data provenance, AI restrictions and document update without running science.

Data: `motion_ledger/displacement_m` is `[accepted_step,9,3]`, terms named in attributes. `step_meta` holds start time, dt, accepted flag, maximum particle residual, PBD/wall mean path lengths and wall reconstruction residual. Sparse `wall_contact_events` columns are in attributes. All physical quantities are SI. `per_ID_cumulative_m` is `[32,9,3]`. Rejected intervals remain separate, never valid examples. Only fixed-step M0 single-group mode is admitted, with no general adaptive rollback/restart claim.

`src/mevo_mb/fast_sim/f9_samples.py:read_audit_sample` gives observed causal prefixes. Uncomputed historical ledgers are None. Final quality is explicitly audit-only, never causal input. All original-Y paired controls/timesteps/replays stay in one family; training_ready remains false.
