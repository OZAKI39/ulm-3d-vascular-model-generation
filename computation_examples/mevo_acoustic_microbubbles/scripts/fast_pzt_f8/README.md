# F8 reproduction and audit

Read `docs/FAST_PZT_F8_DEVELOPMENT.md` and `configs/fast_pzt/f8/protocol.json`. F8 adds files and uses immutable F7 inputs. Never overwrite a saved label.

The completed server jobs used the following commands from their recorded commits:

```bash
python3 -m scripts.fast_pzt_f8.remote ablation --label a_static_v1 --phase static
python3 -m scripts.fast_pzt_f8.remote ablation --label a_pilot_v1 --phase pilot
python3 -m scripts.fast_pzt_f8.remote ablation --label a_formal_v1 --phase formal
python3 -m scripts.fast_pzt_f8.remote y_navigation --label y_static_v1 --phase static
python3 -m scripts.fast_pzt_f8.remote y_navigation --label y_pilot_v1 --phase pilot
python3 -m scripts.fast_pzt_f8.remote y_navigation --label y_formal_v1 --phase formal
python3 -m scripts.fast_pzt_f8.remote contact_check --label contact_v1
python3 -m scripts.fast_pzt_f8.remote y_navigation --label y_half_v1 --phase half
```

These are provenance commands, not instructions to repeat completed work. For reproduction use the exact execution commits in the report and a new isolated results root; static/pilot admission and one major GPU job at a time still apply. The launcher rejects dirty executable source, existing output labels and concurrent mevo GPU jobs. Supervisor owns the processes.

WSL postprocessing, no new integration:

```bash
export PYTHONPATH=src:.
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f8.archive
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f8.analyze
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f8.figures
# Inspect PNG/PDF outputs and deliberately record their reviewed hashes.
# Finalization requires review and progress-document evidence.
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f8.finalize
```

Analysis resets delivery gates pending review; plotting does not manufacture approval. Finalization checks the actual reviewed hashes. All seven complete formal trajectories, five pilots, diagnostic states and logs are archived; build caches are excluded. No contact-candidate navigation or formal RPY Y run was admitted.

`src/mevo_mb/fast_sim/f8_samples.py:read_navigation_sample` adds causal task features to the inherited F6 interface. End-run quality remains audit-only; no independent AI validation set exists. Scientific classifications are separate from numerical quality.
