# Regularization development state

User prompt: `../../provenance/regularization_20260930/USER_PROMPT.txt`.

All work extends the existing `clot_clearance` project. Existing source/config/run/movie files must remain unchanged. Snapshot and 627 protected-file hashes saved in `provenance/regularization_20260930/`. No Git repository exists here. Baseline: 26 existing tests passed (26 s).

User authorized a demonstration fracture energy. Baseline target: **0.01 J/m²**, uncalibrated. Keep existing material constants, clot physical bounds, mean pipe flux, 0.01 m/s streaming velocity scale and 25,000 represented cycles; do not raise loading to force attractive fragmentation. New manufactured field center: **[-0.825, 0, -0.725] mm**, 0.2 mm upstream of the clot left face at x=-0.625 mm, centered vertically. Final movie: existing black/Arial regular/damage-color style, no arrows, no N header, no large blue bubble. Include a small upstream location marker if needed for the requested source position.

Phase order: 0 preservation/baseline; 1 energy accounting; 2 numerical coupon; 3 energy cyclic damage; 4 fragment diagnostics; 5 unified manufactured field; 6 run; 7 cycle study; 8 spacing study; 9 fragment drag; 10 importer; 11 reports and visuals. Run phase-specific checks before progressing. Do not proceed after unexplained energy growth/invalid J/NaN/graph inconsistency.

Current state: **PAUSED AT USER REQUEST**, to allow computer hibernation. Read `PAUSE_RESUME_ZH.md` and `PAUSED_STATE.json` before resuming. No jobs from this task remain active. Phases 0–5 passed their checks (Phase 2 coarse coupon only); Phase 6 runner was written but remains untested and unrun. No new regularized animation exists yet. Wait for the user's explicit resume instruction.

Explicit scientific limitation to preserve: replacement of independent 100 Pa synthetic traction by stress derived from the existing 0.01 m/s velocity scale changes the actual loading substantially. Absence of failure at the specified exposure is a result, not grounds to weaken D_break or tune the video. First-detachment acceptance criteria cannot be claimed if no detachment occurs.
