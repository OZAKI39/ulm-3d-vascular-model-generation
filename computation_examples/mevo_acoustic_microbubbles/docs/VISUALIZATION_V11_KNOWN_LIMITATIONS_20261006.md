# v11 morphology movies: read-only interpretation notes

Recorded during the 2026-10-06 discussion. This note changes neither simulation nor rendering. It supplements the historical delivery's per-frame geometry/encoding checks, which do not establish temporal continuity of the display mapping.

Source: `results/stage4/lubrication_upgrade/visualization/B2_B3_visible_morphology_v11_20261006/RENDER_MANIFEST.json`. Both full 6 s inputs have `secondary_enabled: false`; they are not Stage 4M A2 results. Their 301 native saved states occur at 20 ms intervals. Each is held for four video frames (the endpoint for one), giving 1,201 frames, 60 fps and 0.30× playback, or 15 distinct states per video second.

Display rule in `scripts/stage4l/magnified_geometry.py`: `y_i = x_anchor + 24 * (x_i - x_anchor) + shift_group`. Radius is also scaled by 24. Each connected component uses its own anchor and a rigid translation to fit inside the lumen and separate inflated groups. There is no member repacking or group rotation, but component merges/splits can change this mapping discontinuously. Absolute displayed locations and displacements are not scientific measurements.

For B2, the largest native group at t=6 s contains 63 equal-radius bubbles (R=1.78 µm). Principal-axis spans of the **centers** are about 48.15 × 13.88 × 0.48 µm. This is a near-wall sheet, not a roughly spherical 3D cluster. Its primary forces point toward the wall; its displayed group is rigidly shifted about 0.188 mm. The ordered appearance occurs in the recorded model coordinates, whereas apparent detachment from the wall is a display limitation. Physical realism of this regular packing is not established by this observation.

Examples of display discontinuities, checked against native arrays, display arrays, group transforms and projection records:

| B2 particle ID | Native interval (s) | Video time (s) | Native step (µm) | Display step (µm) | Screen jump at 4K (px) | Membership change |
|---|---|---:|---:|---:|---:|---|
| 54 | 2.44–2.46 | 8.20 | 10.406 | 1,215.909 | 121.39 | 3 → 12 bubbles; anchor 2 → 82 |
| 2 | 3.04–3.06 | 10.20 | 9.729 | 915.505 | 91.42 | 4 → 52 bubbles; anchor 9 → 25 |
| 189 | 4.96–4.98 | 16.60 | 4.055 | 980.577 | 101.48 | singleton → 24 bubbles; anchor 189 → 76 |

The B2 active mask has zero inactive-to-active transitions and 90 removals. Thus the cited flashes are display jumps, not bubble respawning. The prior `DISPLAY_CHECKS.json` already records a maximum B2 display step of 1.216 mm but does not reject it as a temporal-continuity failure. Preserving a group's geometry within each frame does not guarantee truthful displayed motion across frames.

Future visualization work must explicitly address the conflict between independently magnified components, changing membership, exact within-group similarity and temporal continuity. Raising the encoded frame rate alone cannot repair an anchor jump. No scientific coordinate smoothing, synthetic packing or physical-accuracy claim is justified by these notes. No fix or new video was generated as part of the sync task.
