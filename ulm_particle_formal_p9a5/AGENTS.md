# Current workflow

- Start with README.md, CURRENT_WORKFLOW.json and CURRENT_RESULTS.md.
- User explicitly requested local cleanup on 2026-09-27. Obsolete worktrees, historical reports, logs and server setup copies were deleted. Do not restore old snapshots or require preservation of all Git-tracked history.
- Preserve current geometry/network inputs, NEW Network-H0 flow, P9A4 accepted ordering, P9A5 CORE500 and nominal dt=0.001 s. Do not substitute the historical default frozen_reference for NEW.
- Some source files have old stage names but are live dependencies or part of the mandatory 114-source identity contract. Preserve their bytes unless a separately authorized scientific change requires updating provenance.
- Keep trajectory completion hashes and source/event identity checks. Distinguish contact-supported stationary states from residence-time censoring.
- RBC displays are separate retained models; do not claim H0-coupled RBC production from their existence.
- Use scripts/check_current.py for current regressions and scripts/verify_current_data.py for existing formal trajectory hashes. Full re-integration parity tests require --full-simulation.
- Completed development deliverables need readable figures (300 dpi PNG and PDF) and an HTML entry. A cleanup operation does not imply unfinished animations or production delivery gates have completed.
- Current server paths are in /home/lzy/projects/CURRENT_SERVER_PATHS.md. This cleanup was local only; do not delete server content without authorization.
