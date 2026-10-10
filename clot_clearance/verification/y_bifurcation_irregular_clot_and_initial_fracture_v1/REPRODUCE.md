# Saved-data reproduction
Run from the clot_clearance repository with /home/lzy/projects/computation_examples/_shared/python_environment/bin/python and OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.

Offline entrypoints (no solver invocation):
- python -m pd_clot.lammps_fragmentation.y_fragment_visualization figures
- python -m pd_clot.lammps_fragmentation.y_fragment_visualization movie
- python -m pd_clot.lammps_fragmentation.y_fragment_visualization paraview
- python -m pd_clot.lammps_fragmentation.y_fragment_delivery eventaudit
- python -m pd_clot.lammps_fragmentation.y_fragment_delivery aggregate

Actual executed mechanics sequence was morphology freeze, geometry build (one generic pre-dynamic wall-chord repair retained), cases freeze, one run0, C0 zero-load gate, C1, C2, C3. Each completed mechanics directory refuses recreation. No further dynamics is authorized by this delivery. Sources and exact input decks/requests are in source_snapshot and cases; source_before and CODE_CHANGES.diff document changes. No extra dependencies were installed.

Native trajectories have 234 frames (step 0 to 4660, every20). Every distinct native break step is also audited for transient >=200/300um components, without interpolating positions. All positions are SI, display only multiplies by1000 for mm. Reference normal projection is frozen; Gaussian envelopes use current xyz each step. MP4 is slow playback of 2.33ms; frame map records exact physical time.
