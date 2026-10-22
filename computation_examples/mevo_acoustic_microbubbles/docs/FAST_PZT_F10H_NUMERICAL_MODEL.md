# F10H numerical modes and geometry consistency

F10G baseline: `8df3f063af6e0099ac630cd350bdd8a6a22ed4bd`. All additions are
opt-in files. The original F2/M0 default remains two lubrication iterations,
12 simultaneous Jacobi PBD passes, 10 us. The 17.8 nm geometric residual and
178 nm individual contact-iteration correction budgets are unchanged.

## Original GPU order and explicit four-pass mode

The original engine already has a bounded `mobility_iterations` parameter.
F3 allocation requires the unchanged F2 configuration; the diagnostic allocates
under that configuration, then selects 2 or 4 existing GPU kernel invocations.
No lubrication kernel, force kernel, CFD field or geometry file is rewritten.
`M0_LUB4_PBD12_DT10` differs from the two-pass control only in the iteration count.
`M0_LUB4_PBD12_DT5` keeps four lubrication and twelve PBD passes per real 5 us step.

Actual order: consistent acoustic primary response -> secondary pair force ->
CFD lookup and combined free Stokes velocity -> relaxed block Jacobi lubrication
-> ordered component-mean balance -> position prediction -> initial wall pass
-> 12 simultaneous PBD passes (wall pass after every third) -> single-pass budget
check -> original terminal clipping/final geometry/CFD-domain checks -> contact
residual check and acceptance -> accepted clock. This is not sequential PBD.
All original pressure, finite-value, wall and domain checks remain active.
The F10H wrapper checks quality after every physical step and stops at the first
rejection; rejected coordinates never replace the accepted state or advance time.
Wall corrections in these interior experiments are audited separately.

Every science attempt stores full-ID positions, raw/free/lubricated velocities,
forces and step times. The scalar ledger retains all PBD-round residuals and
single-round corrections; the last 32 attempts also retain complete round vectors.
Only valid accepted states enter morphology. Graph event records are filtered by
the exact accepted clock, not by a rejected proposal or an invented later frame.
The final accepted state between 50 us output frames remains in step evidence.

The resistance A, coefficient, cutoff and normalization are exactly the F10G
operator in `FAST_PZT_F10G_CONTACT_MODEL.md`. Its saved 2/4/8/direct static results
are reused. Offline audits substitute each actually executed velocity into that
same resistance; residuals and discrete-reference errors are not physical errors.
Four iterations passing contact is not an assertion of algebraic convergence.

## Geometry-consistent B2 diagnostic

Use the fixed F10G candidate `kc = 0.04641319518378658 N/m` (2kc), not a new fit.
During one physical step, A=A(x_old) and external acoustic/CFD inputs are frozen
at their original evaluation time. Let zeta=6*pi*mu*R and

    U_base = A^-1 U_free
    x_pred = x_old + dt U_base
    gap_ij(y) = |y_j-y_i| - 2R
    delta_ij(y) = max(-gap_ij(y), 0)
    n_ij(y) = (y_j-y_i)/|y_j-y_i|
    F_i(y) -= kc delta_ij(y) n_ij(y), F_j(y) += ...
    y = x_pred + dt A^-1 F(y)/zeta.

Contact normals, true center distances and active contacts are re-evaluated at
every nonlinear iterate and at the final returned coordinates. Separation means
exactly zero geometric spring force. Reaction magnitudes cannot be negative;
forces are central at the final coordinates, pair-balanced, and torque-balanced.
The force is a geometric overlap penalty, not a predictive nonpenetration
constraint. No force from the old linearized active set is carried across steps.
No extra PBD, RPY, adhesion, shell deformation or empirical force scaling is added.

In scaled displacement q=(y-x_pred)/R, with alpha=dt*kc/zeta, solve the stationary
equation of

    E(q) = 0.5 q^T A q + 0.5 alpha sum max(2-|Y_j-Y_i|,0)^2,
    Y = (x_old-origin + dt U_base)/R + q.

The gradient is A q - alpha B(Y)^T delta_scaled(Y). The bounded iteration uses
the positive Gauss-Newton metric A+alpha B_active^T B_active and Armijo decrease.
The full geometric energy is not claimed globally convex: its true Hessian also
has tangential curvature terms. Therefore a converged local stationary solution
is not a proof of uniqueness or the global energy minimum. The pre-registered
limits are 24 outer iterations and 20 line-search trials; failure is retained.
The displacement-unit residual and final position/force equation must be below
1e-13 m. All 496 pairs of the 32-bubble state are evaluated.

Independent references use a separately written pair-force loop and the finite-
difference Powell hybrid solver through [SciPy root](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html).
They do not reuse candidate gradient, B assembly, active-set iteration or metric.
Reference residual, library success flag, function-evaluation count and position
error are recorded separately; no bare solver-success claim replaces a residual.
Collinear two-bubble cases also check the analytical overdamped backward-Euler
formula. ID exchange and nonnegative/contact-only reactions are independently
checked. These references validate the chosen discrete law, not SonoVue physics.
SciPy's current online manual is newer than the installed server 1.15.3; environment
versions are archived. No library upgrade or new literature formula was needed.

The GPU adapter reuses original force/CFD/domain kernels and solves only contact
on the CPU. Any nonzero wall projection rejects this interior-only prototype.
No CUDA Graph or batch-production admission exists for this hybrid method.
Inherited PBD configuration fields are allocation compatibility; actual contact
identity is `B2_GEOM_K2_DT10/DT5`, with no executed PBD passes. Explicit HDF5
numerics and `NUMERICAL_IDENTITY.json` distinguish the executed method.

## Preregistered morphology, costs and stopping

`configs/fast_pzt/f10h/protocol.json` was committed before formal experiments.
Same-model comparisons use identical initial coordinates and control history;
external forces are evaluated by the same law at each evolving state. They are
not artificially forced to have equal force vectors along different trajectories.
Compare only coincident accepted physical times, with no temporal interpolation.
The fixed engineering limits are centroid 10 nm, fixed ID 178 nm, Kabsch 89 nm,
Rg 1%, tensor/aspect 2%. Largest-member identity and edge on/off histories are
also checked. Event times allow 10 us difference; two-ms persistence is evaluated
from true edge times, and right-censored events are explicitly not evaluated.
These are numerical comparison thresholds, not experimental accuracy claims.
PBD correction path and nonlinear iteration path are algorithmic workloads, not
the same physical observable. Between-model comparisons do not diagnose either
model's own temporal convergence.

Only two original short windows, at most one O/1102 10 ms LUB4 attempt and two
conditional DT5 windows are allowed. B2 uses small systems, saved32, then at most
O/1102 5 ms and a conditional same-model DT5 comparison; no B2 10 ms or B64.
All rejected attempts remain archived. New schemes and dt variants retain their
original split_group; no new independent seed or OFF/S/O triplet is invented.

Matched short cost runs use the same 0.2 ms input, checks and output requirements.
CPU wall time, integration/quality, diagnostic serialization, readback, morphology,
graph capture, HDF5 and complete process time are separated. CUDA events spanning
host-controlled loops include CPU gaps; those values are not pure kernel times.
The legacy recorder's completed-short-window samples/hour is not production yield.
No mode becomes the default or `training_ready` merely by crossing an old failure.
All results remain EXPLORATORY / PHYSICS_UNCALIBRATED.
