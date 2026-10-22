# Fast PZT aggregation F1

Base: reviewed R3 archive `5c0a0f10b19b376d9321510149cb6db50bcf4331`.
New branch: `dev/fast-pzt-aggregation-f1-20261007`. R3 worktree HEAD matches;
its untracked `SOURCE_COMMIT_R3.txt` is left untouched. Sparse checkout excludes
historical results from this new worktree. No earlier solver or gate is modified.

## Implementation

Reuse the existing complex acoustic basis and PZT amplitude/phase conventions,
APECSS adaptive response data, cycle-average harmonic force definitions,
Stokes drag, Stage 4L normal-resistance coefficient, and connected-component
cluster definition. Inspect actual table validity before selecting inputs.

Implement one explicitly selected NVIDIA Warp GPU backend in `fast_sim`.
Positions, eight complex volume harmonics, forces, velocities, contact workspaces,
and quality counters stay on the GPU. Full-pair forces preserve harmonic order
and attraction/repulsion signs. Recombine real basis grids only at control events;
query them by bounded trilinear interpolation. At fixed radius/frequency,
reduce adaptive cells to their exact piecewise-linear pressure sections without
resampling, extrapolation or new APECSS solves. Zero drive has exact zero response.

Use overdamped motion, optional local normal-mobility damping, and a small fixed
number of mobility-weighted PBD contact passes. This is a new exploratory motion
approximation, not the R3 coupled solve. Keep a classical far-field option and a
versioned continuous near-field regularization; preserve signed correlations.
Do not add adhesion, cluster locking, velocity clipping, walls or CFD.

Fixed macro steps of 1/2/5/10/20/50 microseconds are candidates, evaluated only
on one 32-particle formation scene. Report real internal steps and contact
corrections. Use bounded substeps only if evidence requires them. Align control
events exactly; never cross a control switch with old parameters. CUDA graphs
batch repeated kernels while retaining actual per-step quality counters.

## Ordered evidence and limits

1. GPU end-to-end smoke and tiny NumPy formula comparisons; two-bubble OFF,
   attraction, signed response test, and smooth near-field transition.
2. Fixed-seed dispersed 32-particle formation with Secondary ON/OFF, zero flow.
3. Limited step study, GPU kernel breakdown, one 200-particle/5-ms timing.
4. One real-basis formation and PZT control-switch trajectory plus fixed-control
   continuation if needed to distinguish control response from natural evolution.
5. One actual long simple run if predicted to fit a short budget, then eight
   independent batched samples and the corresponding serial throughput.

Stop each exploratory run after at most 180 wall seconds; compilation is separate.
At most two evidence-based numerical revisions of a failed step/contact candidate
before documenting the limit. No high-fidelity reference campaign. Save failed
attempts rather than relabeling them. Freeze numerical quality and statistical
step-selection criteria before running the step study.

Save real sampled positions/IDs/radii/control events in HDF5, cluster membership,
shape tensors/principal axes and observed lineage events. Four English figures
and one continuous 3D MP4 use actual positions and radii, with no interpolation
or stitched trajectories. Final gate separates implementation, aggregation,
control response, numerical quality, throughput and EXPLORATORY physics scope.

References: [Warp installation](https://nvidia.github.io/warp/latest/user_guide/installation.html),
[stable releases](https://github.com/NVIDIA/warp/releases),
[XPBD](https://mmacklin.com/xpbd.pdf),
[Doinikov & Bouakaz 2015](https://doi.org/10.1103/PhysRevE.92.043001),
[Fan et al. 2020](https://doi.org/10.1016/j.ultsonch.2019.104790).
The near-field regularization is not a calibrated implementation of either paper.

Prototype correction: finite local block-Jacobi iterations can violate the
net Stokes force balance of a lubrication component. F1 v2 adds GPU union-find
and a common velocity residual correction per current component, preserving
all within-component relative velocities. This is a three-component force
balance, not locking or a QP. Keep v1 records as superseded evidence; verify
the invariant and only the 1/20 microsecond representative endpoints again.
The original numerical-quality thresholds remain unchanged.
