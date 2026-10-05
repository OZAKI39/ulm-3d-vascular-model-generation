# Stage 4L: scoped normal lubrication translation

This opt-in backend changes overdamped translation only. The frozen k-Wave basis/effective fields, APECSS response table, radii, control amplitudes/phases, CFD flow, classical secondary Bjerknes expression, and branch-classification planes are inputs and are unchanged. Existing results retain the label `LEGACY_STOKES_TRANSLATION_BASELINE`; no old run is relabelled or resealed.

Production description: **Stokes self-drag + normal particle-particle and particle-wall lubrication resistance**. Interface scope: **tangentially immobile-interface approximation for near-contact hydrodynamic resistance**. This is not validated full near-contact bubble physics.

## Resistance problem and conventions

For stationary walls, solve

\[
A U=b,\quad A=D+L_{pp}+L_{pw},\quad D_i=6\pi\mu R_i I,
\quad b_i=F_{acoustic,i}+D_i u_{flow,i}.
\]

For equal spheres, \(n=(x_j-x_i)/d\), \(h=d-2R\), and \(q=(U_i-U_j)\cdot n\). The force on i is \(-a_{pp}qn\), with the opposite force on j. It adds \(+a_{pp}nn^T\) to the two diagonal blocks and its negative to the two off-diagonal blocks. Common translation creates no pair drag.

\[
\phi_{pp}(\epsilon)=\frac{1}{2\epsilon}-\frac9{20}\ln\epsilon-\frac3{56}\epsilon\ln\epsilon,
\quad a_{pp}=3\pi\mu R[\phi_{pp}(\epsilon_{eff})-\phi_{pp}(0.1)].
\]

The leading coefficient multiplying **relative** normal velocity is \(3\pi\mu R^2/(2h)\). In symmetric squeeze with velocities +V and -V the relative speed is 2V: the individual force is therefore 2aV. The explicit factor of two is covered by the analytic tests. These project conventions must not be mixed with source conventions based on individual approach speed.

The wall gap is the center-to-wall distance minus R, with outward unit normal \(n_w\). Its force is \(-a_{pw}(U_i\cdot n_w)n_w\): tangential motion is unaffected.

\[
\phi_{pw}(\epsilon)=\frac1\epsilon-\frac15\ln\epsilon-\frac1{21}\epsilon\ln\epsilon,
\quad a_{pw}=6\pi\mu R[\phi_{pw}(\epsilon_{eff})-\phi_{pw}(0.1)].
\]

The leading wall coefficient is \(6\pi\mu R^2/h\). Both corrections are zero for \(h\ge0.1R\), and use \(\epsilon_{eff}=\max(h/R,0.001)\) inside activation. At R=1.78 µm these distances are 178 nm and 1.78 nm. The floor saturates a coefficient; it is **not** an equilibrium spacing, shell repulsion, adhesion, or a rigid-cluster constraint. Negative residual gaps also use the saturated coefficient until nonpenetration handles them.

The backend currently requires monodisperse radii (relative floating tolerance 8 machine epsilons). Other radii raise `UNEQUAL_RADIUS_LUBRICATION_NOT_IMPLEMENTED`. No silently substituted effective radius is used.

## Geometry and linear solve

`wall_gap_and_normal(x, r)` returns gap, outward normal and validity. The production Y geometry reuses the exact `vtkImplicitPolyDataDistance` level/gradient and its existing coordinate cache. No Warp or approximate surface is used. The cylinder implementation is for software fixtures. Curved-wall lubrication uses the local-plane approximation.

Assembly uses a SciPy CSC matrix and cKDTree pair candidates. Diagonal self-drag removes the translational nullspace. Positive coefficients give a symmetric positive definite matrix. Sparse LU solves the matrix scaled by median self-drag; the residual is checked in the original dimensional equation against 1e-11. Diagnostics include matrix size, nonzeros, solver, residual, solve time and failure. With no active lubrication the original Stokes operation order is used directly for bitwise far-field equality.

## Contact in the resistance metric

Instantaneous touching contacts solve

\[
\min_U\tfrac12(U-U_0)^TA(U-U_0),\qquad G U\ge0,
\]

where pair rows give \((U_j-U_i)\cdot n\), and wall rows give \(-U_i\cdot n_w\). There is no tangential constraint. A SciPy bounded quadratic dual followed by active-set polishing is checked using primal feasibility and complementary slackness; an optimizer success flag alone is insufficient. Dense checkpoints have persistent contacts, so the old diagonal NNLS is not used by this backend.

For RK stage proposals and the final displacement, a **finite-horizon nonpenetration guard** uses the same A metric at the step start, with \(G U\ge-\max(h,0)/\Delta t\) for predicted contacts, and zero lower bound for touching pairs. This permits positive subfloor gaps to close to zero; it never pins them at the floor. The spatial projection then handles residual geometric error. This constrained RK rule is a numerical contact treatment and still requires timestep sensitivity; it is not a film drainage solver or added damping. The legacy path retains its old RK and diagonal contact treatment.

All projection calls, including intermediate RK stages, contribute to count, maximum and RMS correction, and corrections above two radii are separately counted. States distinguish `LUBRICATION_ACTIVE`, `SUBFLOOR`, `HARD_CONTACT` and `PROJECTION_FALLBACK`; subfloor/hard-contact counts can overlap by definition.

## Timestep and predictive neighborhood

For a closing gap h above the coefficient floor, use \(\Delta t\le0.5h/c\). Wall closing speed is \(U\cdot n_w>0\); pair closing speed is \(-(U_j-U_i)\cdot n>0\). Separating and tangent motion do not limit the step. Gaps at or below the floor are excluded; the finite-horizon contact guard handles overlap risk. There is no `(h-floor)/c` formula.

A conservative cKDTree radius includes activation plus `2 max|U-mean(U)| * remaining_time`. Thus fast closing pairs outside activation enter the gap guard before they could cross the whole layer. Common cluster translation does not shrink the pair timestep. Wall limits use the exact cached query. The legacy absolute displacement limiter is retained only on the legacy path.

## Scope and references

Not included: full thin-film drainage; Young–Laplace contact deformation; tangential/rotational lubrication; shell repulsion; adhesion; fusion; coalescence; high-fidelity near-contact acoustic interaction; history/added-mass forces; feedback on CFD/acoustics. Residual tangential rearrangement is possible and does not justify automatically adding friction or locking clusters. Low Reynolds number is an assumption, not a consequence of a converged linear solve.

- Ball & Melrose (1997), [doi:10.1016/S0378-4371(97)00412-3](https://doi.org/10.1016/S0378-4371(97)00412-3): pair-resistance/many-body mobility framework and squeeze mode. Rigid-sphere suspension context; not SonoVue validation.
- Costa et al. (2015), [doi:10.1103/PhysRevE.92.053012](https://doi.org/10.1103/PhysRevE.92.053012), [author preprint, equations 32–34](https://arxiv.org/pdf/1506.01880): normal asymptotic functions and saturation concept. This project uses the explicitly specified relative-speed pair convention and activation parameters, not an IBM resolution calibration.
- Chan, Klaseboer & Manica (2011), [doi:10.1039/C0SM00812E](https://doi.org/10.1039/C0SM00812E): deformable drop/bubble drainage and interface-mobility context. A lipid shell motivates an immobile-interface approximation but these references do not directly calibrate SonoVue or validate this model for it.
