# Code and report synchronization

- Date: 2026-09-17
- Destination: `OZAKI39/ulm-3d-vascular-model-generation`
- Branch: `sync/ulm-microbubble-traj-gen-2d-20260917`
- Base: `main` at `8a264014a6181d23d30f3ad4447931ce182e1ac4`
- Source: local `hatimb-particle_flow_simulator/ulm_microbubble_traj_gen_2D`,
  companion vascular-generator code, and the two related root-level reports.
- Requested scope: code and reports; animations and simulation datasets excluded.

## Inventory

261 original files, 12,780,546 bytes (approximately 12.2 MiB):

| Type | Files |
| --- | ---: |
| Python source and existing tests | 184 |
| Markdown reports and documentation | 49 |
| YAML configuration and environment | 5 |
| Requirements text | 2 |
| Static PDF report figures | 10 |
| Static PNG report figures | 11 |

219 files come from the microbubble package, 40 from its companion generator,
and 2 are related reports (`report_sp.md`, `transport.md`). The largest original
file is 1,145,767 bytes. All original-file hashes are recorded in
[SYNC_MANIFEST.csv](SYNC_MANIFEST.csv). New synchronization documentation,
ignore rules, and the empty-data-directory README are additional to this count.

No GIF/video, NPZ/MAT trajectory or field arrays, VTK exports, vessel inputs,
interactive HTML outputs, raw reference-paper PDFs, Python/Numba caches, or
temporary diagnostics were added. Existing data already on the destination's
`main` branch is unaffected. Newly synchronized files use ordinary Git storage.

## Validation

- All 261 copied original files matched the source byte counts and SHA-256 hashes.
- All 261 staged Git blobs matched the same manifest; no line-ending changes
  were introduced. The staged change contains no excluded data or cache files.
- All 184 Python files parsed successfully with Python 3.14.
- Generation and visualization entry points both completed `--help` successfully
  from the isolated destination checkout.
- Seven existing test modules were run from that checkout: **82 passed,
  1 failed, 39 subtests passed** in 12.51 seconds. Command:

```bash
python -m pytest -q \
  ulm_microbubble_traj_gen_2D/test_files/test_generation_skip_render.py \
  ulm_microbubble_traj_gen_2D/test_files/test_field_reuse.py \
  ulm_microbubble_traj_gen_2D/test_files/test_particle_dynamics_config.py \
  ulm_microbubble_traj_gen_2D/test_files/test_molecular_config.py \
  ulm_microbubble_traj_gen_2D/test_files/test_progress.py \
  ulm_microbubble_traj_gen_2D/test_files/test_particle_inlet_flux_v14.py \
  ulm_microbubble_traj_gen_2D/test_files/test_molecular_binding.py
```

The failing case is
`GenerationSkipRenderCliTests.test_candidate_preparation_waits_for_selection_and_continues_generation`.
It reads the default configuration and resolves a real vascular model before
reaching its mocked workflow. It failed because this code-only checkout has no
matching `.swc`/`.vessels.npz` pair in `vessel_swc_models/`. Those inputs are
intentionally excluded by the requested scope. The test and numerical code are
preserved unchanged; this is not a claim that the full suite passes.

Full CFD and trajectory experiments were not rerun. DOLFINx is unavailable in
the validation interpreter, and experiment inputs are not included. Historical
path limitations and input requirements are documented in [README.md](README.md).
