# Current acoustic–microbubble project · 2026-10-10

This is the current entry point for the existing project. The active model is F5,
using the original Y vessel, stored acoustic basis/response table, and prescribed
steady CFD. This synchronization preserves the delivered numerical sources and
actual results; it introduces no new physical model or scientific trajectory.

## Identity and scientific results

| Item | Identity |
|---|---|
| Sync branch | `sync/mevo-acoustic-f5-migration-20261010` |
| F5 delivery source | `66d4cebec18555a26bc9a8035cc3082756d0e76f` |
| Actual F5 science execution | `ee6b35dce83bd3583457a7e85c163b6886e59ce5` |
| Original F4 comparison baseline | `e247fa0c61f668a9b0b2e26c909dd687eedde389` |
| New-server migration verification | 2026-10-09, `COMPLETE_VERIFIED` |

The formal result remains [F5 science_v1](../results/fast_pzt/f5/science_v1):
32 fixed-ID bubbles, radius 1.78 µm, fixed step 10 µs. Flow is ON throughout the
main 0–600 ms trajectory; ultrasound is ON for 0–300 ms and OFF for 300–600 ms.
Branch A is the exact main-run suffix. Branch B starts at the identical switch
state with ultrasound and flow both OFF. Full connectivity before release and
subsequent internal rearrangement are numerical observations, not claims of
physical bubble coalescence or fragmentation.

Read the [F5 development report](FAST_PZT_F5_FLOW_RELEASE_DEVELOPMENT.md),
[registered protocol](../results/fast_pzt/f5/ACOUSTIC_RELEASE_PROTOCOL.json),
[final gate](../results/fast_pzt/f5/FINAL_GATE.json), and
[source/input hashes](../results/fast_pzt/f5/INPUT_MANIFEST.json) for exact scope.
The frozen CFD uses water parameters, one-way prescribed flow, approximate local
mobility, PBD contact, and geometric wall exclusion. GPU agreement and rendering
tests do not establish local physical accuracy or convergence.

## Current server and environments

The working deployment is `/workspace/mevo_acoustic_microbubbles` on
`root@95.3.33.46:15283`; the user's existing local SSH alias is `mevo-new`.
The project shortcut `project` resolves to
`repo/computation_examples/mevo_acoustic_microbubbles` under that root.

```bash
ssh mevo-new
source /workspace/mevo_acoustic_microbubbles/activate.sh
# For visualization instead:
source /workspace/mevo_acoustic_microbubbles/activate-render.sh
```

| Environment | Verified stack |
|---|---|
| Compute `.venv` | Python 3.12.3; NumPy 2.2.6; SciPy 1.15.3; h5py 3.16.0; Warp 1.18.0+cu12; PyVista 0.46.5; VTK 9.5.2 |
| Render `.venv-render` | Python 3.13.11; NumPy 2.5.3; SciPy 1.18.1; PyVista 0.49.0; VTK 9.7.0 |

Both environments are independent of global `/venv/main`. The RTX 4090 host's
550.120 driver requires the official same-version CUDA12 Warp wheel instead of
the default CUDA13 wheel. The CUDA12 download URL and hash are pinned in the
[compute lock](../deployment/newserver_20261009/inventory/requirements-newserver.lock.txt).
The [render lock](../deployment/newserver_20261009/inventory/requirements-render.lock.txt)
preserves the newer PyVista API required by the renderer. No host driver change
was made. The archived deployment tools target the recorded server layout;
they are not generic installers and should not be run from an arbitrary clone.

The project-level relative F1/F4 basis/response links resolve to the exact F5
copies already committed in this repository. They reproduce the verified
new-server input layout without copying the 33 MB acoustic basis again.
The legacy `scripts/fast_pzt_f5/remote.py` still targets the old `vast4090` host;
do not use it for this new server. Use the current root and endpoint explicitly.

## Preserved migration verification

The [validation report](../deployment/newserver_20261009/validation/VALIDATION.json)
records 12 original bounded F5 GPU check groups passing on the new host, including
bitwise agreement with archived ON/OFF trajectory prefixes. Maximum sampled
GPU/CPU velocity difference was 1.734935221290501e-18 m/s. Complete dataset reads
covered 54 HDF5 file paths, and 672 saved F5 flow samples were checked against
the CPU sampler. Offscreen rendering used the RTX 4090, the delivered 1,201-frame
video fully decoded, and the H.264 encoder smoke check passed.

The first attempt's rendering failure and logs are retained under
`deployment/newserver_20261009/attempts/initial_gpu_env_render` alongside the
successful separate-render-environment results. These are migration/software
checks dated 2026-10-09, not new physical experiments run for this sync.

## Visualization sources and evidence

The previously local-only vessel renderer and mechanism/attraction renderers and
validators are now tracked under `scripts/stage4`. Their original source bytes
are retained. Selected original HDF5 trajectories, controls, identity/validation
JSON, logs, paper PNG/PDF exports and compact MP4s accompany them:

- [Vessel without PZT](../results/stage4/publication_figures/vessel_only_reference_v1)
- [Matched attraction and PZT translation illustrations](../results/stage4/visualization/aggregation_translation_PZT_v1_20261008)
- [Continuous initial → attraction ON → attraction OFF illustration](../results/stage4/visualization/aggregation_sequence_v2_20261008)

These earlier F3/F4 illustrations use their recorded controls; in the sequential
illustration, only secondary attraction is disabled while primary forcing remains.
They must not be relabeled as the F5 ultrasound-OFF/flow-ON release experiment.
Original manifests may list deliberately omitted previews, extracted frames or
duplicate export formats. Use the new [export inventory](../sync/20261010/README.md)
to determine what is present. The historical full-delivery validators expect those
complete original export folders; they are not drop-in validators for this curated copy.

## Storage and omissions

The new server's `/workspace` is container storage without a persistent host
volume. Stop/start and recycle/destroy have different retention behavior; do not
treat this deployment as permanent archival storage.

On 2026-10-10, the user explicitly requested deletion of the local complete
old-server backup and redundant local wheel/reference copies (about 62.1 GB).
See the [deletion receipt](../deployment/newserver_20261009/LOCAL_BACKUP_DELETION.json).
The active local project/worktrees, new-server files/environments, and small
migration records remain. Historical 2026-10-09 reports describing the then-existing
complete backup are retained as dated evidence; that backup is no longer available.

This branch is a curated continuation snapshot, not a complete historical archive.
Native CFD/APECSS/k-Wave builds, old raw response-node histories, large contact/event
streams, environment binaries, wheels, caches and frame caches are not added.
Current F5 uses the included stored fields/tables and does not need a fresh CFD
or acoustic solve. Arial is an external installed font, not redistributed here;
the recorded server retains its reference copy, and existing renderers provide
a fallback where implemented. No paused historical job is resumed by this sync.
