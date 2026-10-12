# F10G: overdamped contact model, reference and limitations

This is an opt-in CPU diagnostic model, not a production replacement. The archived
M0 default remains 10 µs, two relaxed lubrication iterations and twelve simultaneous
Jacobi PBD passes. Residual overlap is limited to 17.8 nm and every individual
contact position update to 178 nm. Acceptance, solver convergence, agreement between
models, temporal convergence, cost and experimental calibration are separate questions.

## Existing resistance, not an invented hydrodynamic model

For each distinct pair with n_ij = (x_j-x_i)/d, e=(d-2R)/R < 0.1, the actual
F4 active-mask kernel uses C_ij=c_ij n_ij n_ij^T, where

    phi(e) = 0.5/e - 0.45 log(e) - (3/56) e log(e)
    c_ij = 0.5 [phi(max(e,0.001)) - phi(0.1)]
    A_ii = I + sum_j C_ij
    A_ij = -C_ij
    A U = U_free = U_CFD + F_acoustic / zeta
    zeta = 6 pi mu R                       [N s/m]

A is dimensionless. For the saved equal-radius, all-active states, C_ij=C_ji is
positive semidefinite. Therefore v^T A v = sum_i |v_i|² + sum_(i<j)
c_ij [(v_i-v_j).n_ij]² > 0 for nonzero v. Symmetry and eigenvalues are also
checked numerically; the 96D matrix is not assumed to describe a full Stokes flow.
Its condition numbers are approximately 1514 and 1204. There are 42/35 floored pairs.

The original iteration is U_new = 0.5 [U_old + A_ii^-1
(U_free + sum_j C_ij U_old_j)], initialized at U_free. This is a block Jacobi
iteration, not a sequential Gauss–Seidel sweep. CPU reconstruction reproduces
actual saved GPU iterations to <1e-12 m/s; Cholesky and sparse direct solutions
independently agree. Residuals A U - U_free are recorded in m/s, and zeta times
that residual in N, with normalization by maximum free speed and separate block
preconditioning. None is an experimental velocity error.

Original prediction adds the common shift mean(U_free-U_local) independently
within each current lubrication component. This enforces the mean but does not
change relative velocities within a component. Exact A satisfies A times each
component-uniform translation = that translation; thus its exact solution already
has the free mean. The finite-iteration map followed by mean correction is not
claimed to be the symmetric inverse A^-1. Raw and balanced residuals are both
retained. Eight iterations are still far from the reference in some slow collective
modes; no additional Jacobi-count search was performed.

## Overdamped spring derivation and units

Let B contain a row (-n at i, +n at j) for each pair; B dx is the linearized
change in separation. For a repulsive scalar pair force f>=0, F_contact=B^T f.
It separates the centers; it is equal and opposite per pair. With fixed mobility M,

    U = U_base + M B^T f
    x_next = x_pred + dt M B^T f
    delta_linear = b - dt B M B^T f
    f = kc max(delta_linear, 0)

Here b=2R-d_initial-dt B U_base [m], kc [N/m], M [m/(N s)], and f [N].
For two bubbles, M_rel=B M B^T and

    delta_next = max(b,0) / (1 + dt kc M_rel).

The denominator is dimensionless. Noncontact produces zero reaction in this
linearized model; normal force cannot attract, and tangential velocity is unchanged
in the two-bubble normal-only reference. The implicit formula is checked against
an independent dual-force active-set enumeration, including ID exchange and the
hard-contact limit. The explicit comparator uses the spring at the old overlap;
negative signed predicted overlap denotes overshoot/separation, not negative force.

For B1, M=I/zeta and U_base=U_free. For B2, M=A^-1/zeta and
U_base=A^-1 U_free. Thus B2 satisfies A U = U_free + B^T f/zeta:
external forcing and contact forcing each pass through the same resistance once.
It does not damp the already-damped velocity again, add RPY, or stack PBD on top
of the spring. B2 is internally consistent with this approximate discrete resistance,
not a calibrated dense-bubble hydrodynamic theory.

The bounded CPU prototype solves a convex, piecewise-quadratic displacement problem:
q=dx_contact/R, b_scaled=b/R, alpha=dt kc/zeta,

    min_q  0.5 q^T A q + 0.5 alpha ||max(b_scaled - B q,0)||².

A=I for B1. A semismooth Newton step uses A+alpha B_active^T B_active,
an Armijo search and a maximum of 24 iterations. It checks the displacement-unit
mobility residual at 1e-13 m; 100 times tighter fixed-state tests are separate from
time-step convergence. The independent two/four-bubble reference enumerates at
most 64 active sets in force space, solves (I+dt kc B M B^T) f_scaled=b_scaled
on each subset and tests the unilateral conditions on all pairs. It is not the same
primal Newton algorithm relabeled as a reference. No global QP production backend
was added. All 496 pairs are represented in the 32-bubble diagnostic to avoid
missing contacts activated by a neighboring correction.

## This is not inertial XPBD or a measured shell-contact law

[Macklin et al., XPBD (2016)](https://mmacklin.com/xpbd.pdf), section 4 and equations
4–18, start from inertial Newton equations and scale compliance with dt². Their
position multipliers and inverse mass cannot be substituted for the force and
mobility above. This derivation uses dt M; it is called an overdamped spring
prototype, not validated XPBD microbubble contact.

[Marmottant et al. (2005)](https://liphy-annuaire.univ-grenoble-alpes.fr/pages_personnelles/philippe_marmottant/Publications_files/Marmottant2005_JASACoatedBubble.pdf)
models coating tension during radial oscillation through buckling, elasticity and
rupture. Its area compression modulus chi is not a local two-shell contact stiffness,
even though both may be expressed in N/m. No kc=chi identification was made.
The F10G candidate was registered before soft-contact outputs as the maximum
positive 0.5(F_i-F_j).n over saved lubrication-neighbor pairs in both original
states, divided by 8.9 nm. This gives kc=0.02320659759189329 N/m from
2.065387185678503e-10 N. Only 0.5/1/2 kc are scientific sensitivity cases;
10000 kc is solely a two-bubble hard-limit engineering test. This selection is a
numerical force/overlap scale, not endpoint fitting or experimental calibration.

## Frozen normals: a measured limitation

The implemented force uses the initial normal and linearized overlap for each
physical step. It does not solve the fully nonlinear equation with normals and
geometric overlap at x_next. Actual spherical overlap is always measured separately
and used for the unchanged 17.8 nm quality budget. By convexity of Euclidean distance,
the linearized overlap can exceed the actual geometric overlap. Consequently, a
small repulsion can remain after geometric separation. The short-window audit found
1 pair-step occurrence in O and 16 in preparation, at most 2.428/7.047 pN;
maximum linearized/geometric overlap differences were 0.137/0.190 nm. These are
stored explicitly, not hidden as floating-point noise. Contact normals also are not
exactly central at the final coordinates; the actual residual contact torque is audited.

Thus SOFT_CONTACT_*_PASS means the stated frozen-normal discrete problem passed
its engineering and original geometric-budget checks. GEOMETRIC_SPRING_LAW_EXACT
is false. Neither nonlinear shell contact nor physical shape is validated.
Normal refresh, true geometric complementarity and same-model time convergence
are appropriate next investigations before production integration.

## GPU integration scope and retained acceptance

The two conditional short trajectories use original GPU field/response/primary/
secondary/CFD kernels and CPU resistance/contact solution, with one attempted step
at a time. No approximate new force table, frozen acoustic force or cached flow is
introduced. Original field boundaries, pressure domain, CFD lookup, final geometry,
active-state commit and accepted clock remain checked. Any nonzero wall projection
rejects this interior-only prototype with diagnostic bit 1048576; all studied wall
corrections are zero. Soft-solver nonconvergence uses the same diagnostic bit.
Geometric contact budget failure uses 2048 before acceptance. Failure leaves position,
active state and accepted clock unchanged; the raw invalid proposal is retained.
There is no CUDA Graph for the CPU-contact adapter and no batch admission.

The inherited recorder allocates buffers under F2Physics and retains compatibility
fields such as numerical_mode=FULL_PARTICLE_BATCH and actual_graph_chunk=1.
These are not the executed contact-law identity. The explicit original_physics_mode
F10G_B2_K2_DT10, OBSERVER.json, SOFT_STEPS.json.gz and this model document define
the executed diagnostic. Graph replay is disabled; the graph preflight is a no-op.
Legacy valid_samples_per_gpu_hour means completed short diagnostic windows and
must not be used as production yield or full 10ms sample throughput.

## Literature access and appropriate use

- XPBD original full PDF read, especially the inertial derivation; conceptual contrast only.
- Marmottant original author-hosted full PDF read; radial shell response scope only.
- [Chèvremont, Bodiguel & Chareyre (2020)](https://www.sciencedirect.com/science/article/pii/S0032591020305064): publisher abstract/highlights available; full text fetch unavailable. Used only for its stated lubricated-contact motivation, not unpublished formulas or bubble stiffness. Its formulas are not required by this independently derived candidate.
- [Cundall & Strack (1979)](https://doi.org/10.1680/geot.1979.29.1.47): DOI access unavailable; historical DEM reference only, no claimed reading or formula attribution.
- [Wajnryb et al. (2013)](https://arxiv.org/abs/1307.7312): primary abstract read; collective RPY background, no RPY introduced in F10G.
- [SciPy linear algebra](https://docs.scipy.org/doc/scipy/reference/linalg.html): official solver documentation read; actual installed remote SciPy is 1.15.3, not the version on the current web page.
- [Warp 1.18 FAQ](https://nvidia.github.io/warp/v1.18/user_guide/faq.html): official versioned documentation read; actual environment is Warp 1.18.0+cu12.

No inaccessible paper is necessary to supply a missing implementation formula,
so no user PDF upload is needed. Source identities and derived numerical evidence,
not literature alone, support the engineering conclusions.

## Input pressure and deformation scope

At the two saved pre-failure states, the uncoupled local pressure amplitudes are
1355.863–1415.967 Pa (O/1102) and 1027.172–1069.624 Pa (preparation1104).
The archived pair-feedback effective maxima are about 2267.835 and 1709.586 Pa.
Actual input identities and original quality are in configs/fast_pzt/f10g/model_scope.json.
This is the original saved input scope, not a new pressure scan. The target 8.9 nm
is a total center-distance deficit (0.5% of R), not a measured deformation of each
shell. Fixed R spheres are retained; no shell bending, adhesion, coalescence,
friction, local contact-area or nonlinear deformation constitutive law is modeled.
