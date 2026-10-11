# Current server workspace

This project was migrated on 2026-10-09. Read `/workspace/y_clot/README_MIGRATION_ZH.txt` for current paths and verified scope. Activate `source /workspace/y_clot/env.sh` before Python, MPI, LAMMPS, CFD, or rendering work. All dependencies are scoped to `/workspace/y_clot`; do not rely on the old server paths in historical handoff documents.

Accepted scientific data are in `verification/y_clot_point_cloud_and_early_flow_coupling_v1` and the preserved preceding baselines. Use a new output directory for new experiments. Archived evidence keeps its original path strings; `pd_clot/portable_paths.py` relocates them at runtime. The current CFD source field is an unobstructed native Y flow, not clot-resolved two-way FSI.

Native LAMMPS source: `tools/lammps-contact-onset-v1`; shared build: `tools/lammps-native-y-wall-shared-v1-build`. Latest whole-Y renderer: `analysis/y_clot_whole_vessel_4k_v1/render.py`. New Git commits are migration source snapshots; original provenance lives in `/workspace/y_clot/migration/original` and SHA manifests.
