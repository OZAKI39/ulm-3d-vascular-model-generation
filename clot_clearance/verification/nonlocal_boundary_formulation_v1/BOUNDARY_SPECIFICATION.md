# NONLOCAL_BOUNDARY_FORMULATION_V1 — frozen boundary specification

Model: `PFPD_AT1_SPECTRAL_CANDIDATE_V1`. This namespace evaluates **intact mechanics only**. The historical BAQP evaluator, family construction, energy and mechanical solver sources are unchanged; their input hashes are checked in `HISTORICAL_INTEGRITY.json`. Phase, fracture, activation, stability and CFD are outside this run.

## Geometry and quadrature

The reference domain is a Cartesian block \(\Omega=\prod_{a=1}^3[l_a,u_a]\). `BlockBoundary` exposes its six faces, \(\Gamma_D\), \(\Gamma_N\), outward normals, tangent bases, nearest-face projection, distance and interaction-region masks. Face/edge/corner quadrature keeps each incident face separately; it does not invent a single corner normal. The current benchmark has \(\Gamma_D=\{z^-,z^+\}\) and the remaining four faces in \(\Gamma_N\). These are geometry parameters, not hard-coded particle IDs.

Main grid: \(19^3=6859\) particles, 690,020 undirected bonds, \(L=1.125\,\mathrm{mm}\), \(h=62.5\,\mu\mathrm m\), \(\delta=251.25\,\mu\mathrm m\), \(\delta/h=4.02\). Nodal volume is the tensor product of Cartesian trapezoidal widths; face, edge and corner nodes receive half, quarter and eighth cell volumes. Thus \(\sum_i V_i=L^3\). Surface quadrature independently uses trapezoidal areas.

The direct interaction region is \(I_\delta=\{X\in\Omega:d(X,\partial\Omega)<\delta\}\). The conversion's dependency support starts from these nodes and expands through two actual bond-graph hops. This follows the BAQP coefficient dependency envelope; it is not a physical collar prescription. In this N19 block that support covers all physical nodes. The geometric core \(d\ge\delta\) has 729 nodes; it is **not a fully closed three-hop bulk region**. No extrapolation of N19 core measurements to infinite bulk is made.

## Material and local surrogate

Every material call explicitly uses SVK, \(\mu=1000\,\mathrm{Pa}\), \(\lambda=0\), \(K=2\mu/3\). No implicit `Material()` defaults are used.

\[
E=\tfrac12(F^TF-I),\quad S=\lambda\operatorname{tr}(E)I+2\mu E,
\quad P=FS,\quad W=\tfrac\lambda2(\operatorname{tr}E)^2+\mu E:E.
\]

The analytic surrogate is \(u_L(X)=(F-I)X+c\), with normal \(F_N=\mathrm{diag}(1,1,1.003)\) or shear \(F_S=I+0.003e_x\otimes e_z\). These solve homogeneous local equilibrium with the corresponding boundary data. The nominal traction on **every** face is \(t=PN\), split into \(t_n=(t\cdot N)N\) and \(t_t=t-t_n\). In shear, lateral x-face tractions are nonzero and are included. A top-face-only shear would be a different problem.

Future traction inputs may use `pressure_shear_traction(p,N,tau)`, representing \(-pN+\tau\) with \(\tau\cdot N=0\). This helper supplies no fluid solver or adhesion law.

## Four treatments

| Treatment | Unknowns and loading | Meaning |
|---|---|---|
| B0 | Top/bottom physical displacement grips; classical \(PN\) nodal traction on the other faces | Legacy surface-style baseline; grip reactions are included in external work |
| B1 | Same grips; frozen \(f_i^{BC}=-r_i^{BAQP}(u_L)\) on the dependency support | Surrogate-to-volume-force conversion; main candidate |
| B2 | Prescribed \(u_L\) on all physical nodes with \(d<\delta\); remaining physical nodes free; no applied force | **LOCALITY-STRONG REFERENCE**, an explicitly in-domain constrained-support adaptation |
| B3 | None in V1 | Scabbia ghost control deferred until BAQP ghost coefficients and physical energy/flux accounting can be derived |

B1 uses the existing finite-domain BAQP internal **nodal force**, measured in N. The corresponding force density is \(b_i^{BC}=f_i^{BC}/V_i\), measured in N/m³. With internal force \(r=-\partial\mathcal E_h/\partial u\), the frozen-load potential is \(\Pi_h(u)=\mathcal E_h(u)-\sum_i f_i^{BC}\cdot u_i\). The load is not updated at mechanical iterations. B1 starts from a smooth nonzero perturbation of the surrogate, with amplitude \(0.02\times0.003L\), vanishing at all faces.

The construction identity \(r(u_L)+f^{BC}=0\) only checks implementation. Work is independently compared to analytic continuum face integrals, and ordering, force/moment, field recovery and resolution changes are reported. Affine energy consistency explains why affine virtual work is a necessary but restricted test; this run is not a theorem for nonaffine boundary conditions.

B2 constrains physical material, not exterior material. **No exterior/fictitious particles are added**, no exterior energy is integrated, and no extra physical DOFs are introduced. It is not proved equivalent to D'Elia's exterior-domain Dirichlet strategy. It tests what the one-horizon physical constraint can enforce under unchanged BAQP mechanics. Its finite-width free-core residual is allowed to reveal that limitation.

## Force, reaction and work accounting

For a mechanical solution, a constrained node has reaction \(R_i=-(r_i+f_i^{applied})\). The total external force used for balance/work is applied force plus reaction, so at constrained nodes it equals \(-r_i\). B1's converted load is also checked independently of reactions. B2 reactions include **all** constrained physical nodes, not just the outermost planes. A surface reaction and a distributed volume force need not have the same spatial allocation.

\[
F_h=\sum_i f_i,\quad F_*=\sum_{\text{faces }a}|\Gamma_a|PN_a=0,
\]
\[
M_h=\sum_i(X_i-X_c)\times f_i,\quad
M_*=\sum_a(X_a-X_c)\times(|\Gamma_a|PN_a),
\]
\[
\delta W_h=\sum_i f_i\cdot[A(X_i-X_c)+c],\quad
\delta W_*=|\Omega|P:A+c\cdot F_*.
\]

All nine Cartesian \(A_{ij}\) basis variations and three rigid translations are checked. About **reference** coordinates, the shear target moment is \(M_y=|\Omega|(P_{xz}-P_{zx})\); it is small but nonzero. It must not be confused with angular momentum expressed using current coordinates and symmetric Cauchy stress.

Predeclared gates are \(10^{-6}\): force normalized by \(\|P\|_F L^2\), reference moment by \(\|P\|_F L^3\), nonzero virtual work by its own absolute continuum target. Zero-work tests use the stable stress-volume variation scale. Global force is zero by symmetry, so its gate alone is weak. Independent face targets, the first-moment tensor, nonzero work variations and half-space load diagnostics are retained. Half-space forces have no automatic one-to-one identification with classical face tractions.

Normal/tangential load channels are evaluated through the same force assembly at frozen geometry with constant stress channels \(P_n=\mathrm{diag}(P)\), \(P_t=P-P_n\). This is a linear-in-stress diagnostic for Cartesian faces, not a constitutive superposition assumption or a curved-surface algorithm.

## Field and support diagnostics

Reported deformation errors use \(\|\bar F_i-F_*\|_F/\|F_*-I\|_F\), which is more stringent than normalization by \(\|F_*\|\). A separate predeclared known-field gate is \(10^{-5}\); the solver force gate is \(10^{-8}\). Core/boundary RMS, maximum, total energy, positive spectral energy maximum/P95, nodal positive-energy profiles, first interior displacement jumps, reactions and free-node residuals are saved. \(\Psi^+\) is an intact Green-strain spectral diagnostic only.

Support quantiles use cumulative \(\sum_i\|f_i^{BC}\|\) ordered by nearest-boundary distance. Both \(d/h\) and \(d/\delta\) are stored; d50/d90/d95/d99 are grid-resolved quantiles without interpolation. They describe this conversion's force distribution and are not automatically physical boundary-layer thicknesses.

## Existing solver and resources

The existing `regularization.equilibrium._scipy` function code and SciPy TRF/LSMR iteration are reused. The new adapter provides an analytic Jacobian-vector product of the unchanged intact SVK BAQP residual. Directional finite differences and reciprocal bilinear forms test that derivative. The one reporting adaptation skips dense condition-number SVD for a `LinearOperator`; JSON records null, not an invented condition estimate. There is no new Newton, ADR or PETSc implementation.

The legacy solver's `surface_traction_mean_Pa` and `surface_traction_max_Pa` fields are unpopulated placeholders in this nodal-force adapter. Their zeros do **not** mean zero prescribed classical shear. Authoritative tractions are the six analytic face targets and the saved external nodal forces. These legacy surface-only diagnostic fields are not used by any acceptance gate.

The normal cases used direct NumPy indexed aggregation. Subsequent cases use equivalent CSR incidence matrices only for bond-to-node sums inside the tangent action. N9 comparison against the archived original adapter gives a relative action difference of 3.29e-16; all focused tests pass. This changes arithmetic aggregation, not topology, constitutive response, force evaluation, solver iterations or tolerances. The original adapter is archived in `references/MECHANICS_ADAPTER_INITIAL.py.txt`; numerical/performance evidence is in `TANGENT_SPARSE_ASSEMBLY_CHECK.json`.

Serial scientific execution uses one CPU worker and one BLAS thread; no GPU or installations. Per-case cap: 80 residual evaluations, 120 LSMR iterations per linear solve, 1800 s wall time. A timed-out or force-unbalanced result is retained as incomplete. A first short attempt was interrupted to repair the condition-number reporting adapter before accepting any result; its log/resources remain separate.

The pre-run estimates, gates and two-resolution comparisons are frozen in `PREDECLARED_CONFIG.json`. After normal/shear gates, the fixed-L refinement compares N19/N21 at fixed delta and fixed delta/h. Only affine work/support evaluations are planned at N21, with no convergence or asymptotic-compatibility claim. Monitoring records minimum available system memory, peak process RSS and wall time; the required headroom is 20%.

Reproduction commands (from `clot_clearance`, with the project's existing Python environment):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B -m scripts.nonlocal_boundary_audit --stage references
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B -m scripts.nonlocal_boundary_audit --stage solves
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B -m scripts.nonlocal_boundary_audit --stage refinement
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python -B -m pytest tests/test_nonlocal_boundary_formulation_v1.py -q -p no:cacheprovider
```
