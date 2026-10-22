# F10P hard-sphere model and numerical audit

## Mathematical scope

Fixed-radius frictionless nonadhesive spheres, homogeneous-water control proxy, zero background flow. This is a change in contact model, not experimental calibration. It contains no B2 spring stiffness, penalty spring or after-the-fact PBD.

At the current geometry, reuse exactly the F10G dimensionless pair resistance A and Stokes factor ζ=6πμR. Let B contain all pair distance gradients (-n,+n), g the current gaps, s=dt v/R and b=dt F/(ζR). Solve

min_s (1/2)sᵀAs − bᵀs, subject to B s + g/R ≥ 0.

The independently checked KKT conditions are As−b−Bᵀλ=0, λ≥0, w=Bs+g/R≥0 and λ⊙w=0. Physical pair reactions are ζRλ/dt. Normals and resistance are frozen within a physical step and updated on the next step. A=I+Σ c BᵀB with c≥0 is SPD; each actual matrix must additionally pass symmetry and Cholesky checks. The *local velocity QP* is strictly convex. No claim is made about unique nonlinear multisphere final configurations, trajectories, contact force multipliers under redundant constraints, or experimental bubble deformation.

The plane constraint is conservative: |r+Δr| ≥ n·(r+Δr)=2R+g+n·Δr. Nevertheless the true updated Euclidean distance is independently measured for every pair. Local-linear complementarity and true-geometry reaction-gap products are both saved; the latter has a separately preregistered finite-step tolerance. No true-geometry gap is clipped or overwritten.

OSQP lower-bound dual signs are reversed to obtain repulsive λ. Solved inaccurate, time limit and iteration limit fail. A truly unconstrained resistance solution is used directly only when *all* plane constraints are already feasible, with exact zero reaction. This avoids numerical motion in OFF. Each QP setup starts fresh; caching refers to the acoustic basis, not reusing an old contact solution.

## Actual bounded solver revision

The first OSQP attempt was rejected at step4 (3 µs legal prefix): maximum3000 iterations, force residual2.55e−15 N. This failure remains archived. Revision2 solves the *same* local QP with at most96 equality working-set updates using SciPy Cholesky and rank-revealing least squares. It adds violated constraints, removes negative multipliers and rejects repeated working sets, inconsistent equalities or exhausted iterations. It is not Euclidean position projection and adds no physical stiffness. A frozen failed-input recheck and all small references passed before a new150-step and0.25 ms pilot. Original admission tolerances were not changed. Successful direct solves are labeled active_set_solved; the OSQP failure is never relabeled solved.

## Preregistered independent checks

- Original overlap limit 17.8 nm and contact correction limit 178 nm remain; hard-sphere geometry has the stricter 1 nm diagnostic.
- Mobility-weighted equation residual max_i |R A⁻¹(As−b−Bᵀλ)| ≤1e−13 m.
- Force residual ≤1e−16 N + 1e−8 × force scale, where scale=max(max particle external force, max absolute pair reaction,1e−18 N).
- Linear primal violation ≤1e−12 m; negative reaction ≤1e−16 N.
- max |reaction × linear gap| / force scale ≤1e−12 m; true-geometry counterpart ≤1e−11 m. Raw N·m values are also saved.
- The 178 nm contact correction check applies conservatively to the entire difference from the original resistance-only predictor. No intermediate OSQP iterate is treated as a physical substep or a PBD iteration.
- Original R3 GPU wall query checks predictor and accepted proposal. Any wall correction or scope error rejects this interior prototype, which is stricter than the original 178 nm wall allowance.
- Pressure/APECSS, frequency, source control, radius, cache interior and strict zero flow retain the original fail-closed interfaces.

Thresholds, target tolerances and calculation limits were committed in 36343716a before any F10P numerical result. Changing contact model does not authorize changing budgets.

## Foundation references and history

Two-sphere compression has zero normal closing velocity and reaction equal to the imposed compression. Separation/release/no-contact must have zero reaction. Four-sphere simultaneous contact and contact-network translation are compared with independent dual active-set enumeration of at most 64 subsets. Both use the same original resistance and force. B2 results use the unchanged geometric implicit spring and original stiffness for comparison only. Non-SPD, nonsymmetric and forced iteration-limit rejection fixtures are not scientific trajectories.

Historical F10O strong proposals at 0.175 ms and M0 trajectories remain untouched. Cross-model differences in contact and mobility approximations are not time-step error. The original F2 near-field force, lubricating-gap floor, acoustic response interpolation and real wall hydrodynamics remain uncalibrated.

## Primary implementation references

- [OSQP mathematical solver](https://osqp.org/docs/solver/index.html): convex QP interface and residuals.
- [OSQP status definitions](https://osqp.org/docs/interfaces/status_values.html): inaccurate and unsuccessful exit states.
- [OSQP settings](https://osqp.org/docs/interfaces/solver_settings.html): explicit tolerances, iteration/time limits and polishing.

These references support the numerical library usage, not validation of physical microbubble contact. OSQP 1.0.5 was installed only in the existing project Python environments; no CUDA driver or system service was changed. Final measured evidence is in the control-results report and analysis/FINAL_GATE.json.
