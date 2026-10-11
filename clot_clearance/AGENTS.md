# Current clot research snapshot (2026-10-11)

Read `github_sync_20261011/README.md` for current scientific scope and export limitations. The new clot-resolved CFD and conservative nodal map are offline; historical fragmentation trajectories retain their original unobstructed background flow. Surface correction / energy failure is not production-certified.

Use new output directories for experiments or replay. Preserve accepted scientific data and existing SHA manifests. `pd_clot/portable_paths.py` resolves archived asset paths; `Y_CLOT_DEPENDENCY_ROOT` and `Y_CLOT_FONT_DIR` can override dependency/font roots.

The full native environment is isolated at `/workspace/y_clot` on the current server; activate `/workspace/y_clot/env.sh`. Clone copies can verify all 501 accepted frames with `github_sync_20261010/replay_compact.py`. Full native force audits and binaries were intentionally excluded from GitHub. Native source patch and version records are in `github_sync_20261010/native_patches/`.

After each development round, update the WSL file `/home/lzy/projects/clot_rupture.md` only with conclusions supported by implementation and numerical evidence (including reliable mismatches and explicit limitations), preserving its existing style and sections. Deliver the corresponding review figures and a static results index. Save a text backup, recheck the document hash before writing, and retain the diff.
