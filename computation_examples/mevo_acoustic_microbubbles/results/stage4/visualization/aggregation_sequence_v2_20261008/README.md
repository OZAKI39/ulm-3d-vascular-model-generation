# Initial state → Attraction ON → Attraction OFF

This revision uses **one uninterrupted 600 ms simulation**, not side-by-side simulations from separate initial states. The dispersed 32-bubble cloud evolves with the Secondary Bjerknes attraction enabled from 0 to 300 ms. At 300 ms only this force is disabled, and the same particles continue until 600 ms. Positions are not reset, replaced, interpolated or played backwards. PZT drive, Primary force, mobility, pair contact and geometric WALL constraints are retained.

- [Video](Video_attraction_initial_on_off_60fps.mp4): 1920 × 1080, 60 fps, 1201 frames, 20.017 s. Every frame is an actual sample, spaced by 0.5 ms of physical time. The switch occurs at video frame 600 (10.000 s).
- [Figure](Fig_attraction_initial_on_off.png): black background, 4800 × 2100 pixels, 300 dpi. [PDF](Fig_attraction_initial_on_off.pdf) and [SVG](Fig_attraction_initial_on_off.svg) are also supplied. The snapshots are at 0, 299.5 and 600 ms.

The video has a pure black RGB (0, 0, 0) background and retains the reference PZT materials, reflection, transparency and layout. The requested explanatory footers and the “Local views…” and “EXPLORATORY” overlays have been removed. Remaining labels identify the sequence, source layout, time, switch time, coordinates and bubble counts. All four PZT sources stay ON; the changing label concerns bubble–bubble attraction. The single moving bubble panel uses a fixed camera and actual relative sizes, without rigid cluster proxies.

## What the calculation actually shows

The dispersed cloud aggregates into a 32-bubble cluster before the switch. The cluster is still connected after attraction is switched off. Its radius of gyration is 4.9004847 µm at the switch and 4.9005592 µm at 600 ms. Switching off attraction is not the same experiment as starting dispersed with attraction disabled, and it does not introduce a force that separates an existing cluster. Therefore no artificial breakup, reverse playback or substitution of the earlier dispersed OFF trajectory is shown.

The original request's “initial → on → off” sequence is followed literally. Both the near-contact model and this switching experiment remain exploratory. Removing that label from the image does not change the scientific status. The earlier matched ON/OFF-from-dispersed comparison and its evidence are preserved in `aggregation_translation_PZT_v1_20261008`.

## Actual computation and identity

The new run is [data/run/trajectory.h5](data/run/trajectory.h5), with configuration, source identities, switch audit and numerical-quality summary alongside it. It uses the frozen F4 engine `e247fa0c61f668a9b0b2e26c909dd687eedde389`, original F4 WALL mesh, the same parent-wall initial cloud, true 1.78 µm bubble radius, 10 µs integration step and existing four-element PZT drive. Source centres are (0, ±4, 0) mm and (0, 0, ±4) mm, with 6 mm disc diameter. Display backing and rims follow the v6 reference and are illustrative housing geometry.

The source is `scripts/stage4/visual_sequential_ablation.py`. Only its isolated subclass changes the Secondary force schedule. At the switch it discards previously captured CUDA graphs so an ON kernel sequence cannot silently continue into the OFF segment. It verifies that switching leaves the positions bitwise unchanged and that the retained Primary force is nonzero while the Secondary addition is exactly zero. The run's 301 common 1 ms samples through 300 ms match the archived ON trajectory bitwise. It completed in 37.11 s of simulation-loop wall time and passed the inherited numerical-quality checks.

The server run is `/workspace/mevo_visual_sequence_20261008_v2/run`, managed by supervisor program `mevo_visual_sequence_20261008_v2` with finite timeout and automatic restart disabled. Existing engine files, scientific archives and the first visualization delivery are preserved.

## Reproduction

From `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles`, using the existing shared Python environment:

```bash
python -m scripts.stage4.visual_attraction_sequence
python -m scripts.stage4.validate_attraction_sequence
```

The rendering command reads existing data and does not start a scientific run. `--preview` produces selected preview frames; `--video-only` re-renders only the video. Offline rendering used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LP_NUM_THREADS=4`.

`FRAME_MAP.npz` maps each video frame to its actual source frame, time and Secondary-force state. `VIDEO_IDENTITY.json`, `FIGURE_IDENTITY.json` and `RENDER_SOURCE.json` retain display and source identities. `VALIDATION.json` records independent metric recomputation, switch continuity, preserved archive hashes and complete 60 fps video decoding. `ARTIFACT_MANIFEST.json` hashes the delivery files. Frozen script copies are in `source/`.
