# Current workflow

- Start with README.md, CURRENT_WORKFLOW.json and CURRENT_RESULTS.md.
- User explicitly requested local cleanup on 2026-09-27. Obsolete worktrees, historical reports, logs and server setup copies were deleted. Do not restore old snapshots or require preservation of all Git-tracked history.
- Latest completed run: `particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/`, ROI-only-balanced-pressure-v1 flow, 1500 trajectories, nominal dt=0.0005 s. Use its isolated adapter and verification script; preserve the original size distribution and report actual outlet coverage.
- Preserve geometry/network inputs and the retained Network-H0/P9A4/CORE500 reference and its dt=0.001 s contract. The latest user authorization changes only the separately versioned new campaign. Do not use the geometry/RBC frozen_reference as production velocity input.
- Some source files have old stage names but are live dependencies or part of the mandatory 114-source identity contract. Preserve their bytes unless a separately authorized scientific change requires updating provenance.
- Keep trajectory completion hashes and source/event identity checks. Distinguish contact-supported stationary states from residence-time censoring.
- RBC displays are separate retained models; do not claim H0-coupled RBC production from their existence.
- For the latest batch use its `scripts/verify_collected.py`. `scripts/check_current.py` and `scripts/verify_current_data.py` still check the retained H0 reference. Full re-integration parity tests require --full-simulation.
- Completed development deliverables need readable figures (300 dpi PNG and PDF) and an HTML entry. A cleanup operation does not imply unfinished animations or production delivery gates have completed.
- Current server paths are in /home/lzy/projects/CURRENT_SERVER_PATHS.md. This cleanup was local only; do not delete server content without authorization.
