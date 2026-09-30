# Model, topology and transport assumptions

This document applies only to `streaming_fragmentation_demo`. Existing `MODEL_ASSUMPTIONS.md` describes the preserved original prototype. The reference paper supplied by the user is an aspiration/CFD–PD study, not a calibration of this ultrasound fatigue or hydrodynamic relaxation model; see the existing [paper reading](references/PAPER_READING.md).

## Damage and failure

For each original undirected PD bond, the new state contains cumulative raw damage `D`, irreversible `active`, and effective integrity `g`. With `Q` the half peak-to-peak bond-stretch amplitude measured over the representative loading cycle, the existing empirical form is

```
DeltaD = DeltaN * C_damage * max(Q / Q_ref - 1, 0)**m_damage
D_new = clip(D_old + DeltaD, 0, 1)         # only for still-active bonds
active_new = active_old and D_new < D_break
g = (1 - D_new) if active_new else 0
```

The default threshold is 1; the selected uncalibrated case deliberately uses 0.25. Thus a broken bond may have raw `D < 1` while its effective integrity is exactly zero. Broken-bond damage freezes and bonds never heal. The damage plotted on particles is the original weighted neighborhood integrity loss, not simply a raw-bond-D average. Both raw bond damage and effective integrity are exported so this threshold distinction is inspectable.

The deformation-dependent amplitude, damage increment and threshold crossing are recorded separately for all 11,054 broken bonds. No cleavage, imposed gap, particle-position overwrite or manually assigned fragment IDs are used. Positions advance through velocities; the only positional constraint is the unchanged fixed base.

## NOSB-PD continuation after support loss

Damage weights are applied to the shape tensor, correspondence force state and stabilization pair interactions. The full-rank force calculation is the existing NOSB-PD kernel. Fixed-base reactions constrain the wall-connected material; surviving interactions inside a detached component continue to transmit force.

A heavily broken particle's shape tensor may have rank below three. Inverting it as though it were valid three-dimensional support is not justified. The new module explicitly diagonalizes the surviving weighted reference neighborhood, records its numerical rank, and uses only supported eigenvectors `U`. For rank `r = 1 or 2`, define `B = F_supported U`, `C = B.T B`, `J_r = sqrt(det(C))` and

```
W_r = G/2 * (trace(C) - r) - G*log(J_r) + lambda/2*log(J_r)**2
P_r = G*B + (lambda*log(J_r) - G)*B*inverse(C)
```

This is an **uncalibrated intrinsic correspondence continuation**, with the same nonaffine residual stabilization and pairwise action/reaction assembly. It is not a claim that a line or a sheet retains a physically valid three-dimensional clot stress. Rank-three points retain the original three-dimensional neo-Hookean expression and determinant guard. Rank-zero points have no active internal bonds or stress, but retain mass, position, velocity and hydrodynamic integration. There are no ghost bonds, hidden particle deletion or frozen unsupported fragments.

The final global ranks `(0, 1, 2, 3)` contain `(172, 18, 20, 150)` particles. Tests cover full-rank agreement with the old implementation, intrinsic energy/force gradients, rotation objectivity, a live rank-one pair and a disconnected particle. These tests establish implementation consistency, not experimental validity. The old runner's strict support checks remain untouched. In the new explicit rank policy, the legacy `minimum_neighbors` setting does not replace the rank assessment; the configured eigenvalue tolerance and retained-support condition guard determine admissible supported directions.

## Exposure and applied fluid forces

After each macro damage update the exposure detector compares weighted remaining support `sum(w_ij * g_ij * V_j)` with the complete intact cubic-lattice bulk support for the same spacing and horizon. A ratio below 0.8 flags surface exposure; anchors are excluded. This criterion responds to both gradual weakening and broken neighbors. It is an approximate surface criterion, not a resolved fluid interface.

An approximate outward normal is the negative normalized first directional moment of the remaining current bond neighborhood. Attached exposed mobile particles receive the existing analytical traction evaluated at their **current** positions and normals, multiplied by `h**2 * area_factor`. Normals below the configured moment tolerance do not receive surface traction. This surface area rule is a verification quadrature, not an exact reconstructed fracture area.

Every detached particle instead receives relaxation forcing, including the interior particles of a surviving fragment:

```
u_f(x,t) = u_Poiseuille(x) + u_vortex(x,t)
a_h = (u_f - v) / tau_h
f_h = rho_s * V * a_h
```

The localized divergence-free x–z vortex is proportional to `exp(-|x-bubble|**2 / (2*L**2)) * (-dz, 0, dx)/L`, with configured speed and the same prescribed slow temporal envelope as the loading proxy. It is not a resolved acoustic or bubble-shape solution. The velocity and surface-traction amplitudes are independently prescribed; they do not form a solved viscous stress/velocity boundary-value problem. Attached particles use traction and detached particles use relaxation, avoiding duplicate applied surface traction plus drag on the free material. Anchors receive neither relaxation updates nor displacement.

Two exact exponential relaxation half-steps, each with locally frozen fluid velocity, surround the mechanical step. This split approximation is replaceable through `FragmentFluid.velocity`, `surface_force` and `relaxation_half_step`. Detached numerical damping is zero; fluid relaxation still damps velocity relative to the prescribed flow. Attached numerical damping remains 3,500/s. The background velocity is analytically continued along the pipe axis; the 6 mm length is a reference for the original pressure field, not a particle-removal boundary.

There is no two-way flow feedback, resolved fluid obstruction, fragment–fragment contact, chemical dissolution or bubble dynamics. A radial safety guard aborts a particle outside the configured pipe radius; this is not a wall-contact model. The selected run stayed within the radial guard.

## Connectivity, identity and clearance

After every macro damage update, active edges define the graph. Components containing any base anchor are attached; all others are detached. Stable IDs are assigned algorithmically. On splitting, the anchor-bearing child keeps the parent's ID; if none carries an anchor, the largest child keeps it, with particle-index tie breaking. Other children receive monotonically increasing IDs. The full lineage is saved. No merges are possible without healing.

Every component has membership, volume-weighted COM, COM velocity, bounding box and volume. For travel, compare its current COM with the COM of the same material particles at their individual first-detachment positions. The reported distance is net displacement since detachment, not cumulative path length; singleton components are included. This definition avoids ID changes turning a split into apparent transport.

For each current detached component, its current COM and the previous saved positions of the **same current membership** are compared with `x_clearance`. An upstream-to-downstream crossing credits only members never credited before. Splitting, recrossing and persistent downstream residence cannot double count material. The ledger records every credited particle ID and volume. Crossing is detected at macro resolution (2 ms mechanical increments here); a crossing and recrossing between saved macro states would not be resolved. The plane at x = 1 mm is a diagnostic plane, not deletion or an outlet condition. All cleared material remains in the simulation.

## Time interpretation and evidence limits

The mechanical integration step is 2 microseconds. A 500 Hz representative loading proxy supplies one simulated cycle for each 1,000-cycle empirical damage jump; the first macro step includes one warm-up cycle. At N = 25,000, mechanical proxy time is 0.052 s. N divided by the metadata carrier frequency 1 MHz is 0.025 s; this is only cycle bookkeeping and is not the mechanical integration clock. Neither clock is a calibrated ultrasound treatment duration. Slow-cycle geometry can change during a damage update, so late, rapidly changing motion is not established steady acoustic fatigue.

This prototype proves only the numerical sequence from prescribed loading to damage, actual bond failure, detachment and subsequent fluid-driven motion. It does **not** prove correct ultrasound frequency response, microbubble shape oscillations, clinical ablation rate, calibrated fatigue lifetime or quantitative fragment-size distribution. No mesh/time-step convergence or experimental validation of this new failure case is claimed. The exact event cycles and fractions depend on these uncalibrated choices, including the failure threshold, geometry, support continuation, macro increment and relaxation time.
