# NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1 — review

**NONAFFINE_LOCALIZATION_VERIFIED** for the frozen intact quadratic axial benchmark. Boundary localization and bulk-model agreement are separate conclusions. This is a surrogate-dependent manufactured-data audit, not a general physical traction law.

Branch `clot-clearance/pfpd-progress-20261005`; HEAD `695e75166aa864809460ea9a1845fe72abce229c`. Frozen configuration: [PREDECLARED_CONFIG.json](PREDECLARED_CONFIG.json).

## 1. Why affine B1 was insufficient for generalization

The previous affine benchmark has zero full-stencil force and zero local body force. Consequently its negative finite-domain residual can be interpreted as boundary conversion on that surrogate. Here the intact bulk force is nonzero. Negating the entire finite-domain force would cancel legitimate interior mechanics as well as truncation effects. [Previous review](../nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md) remains unchanged.

## 2. Manufactured nonaffine field

On the physical N19 cube, $L=1.125\,\mathrm{mm}$, $h=62.5\,\mu\mathrm m$, $\delta=251.25\,\mu\mathrm m$, $\delta/h=4.02$:

$$u_x=u_y=0,\qquad u_z=\epsilon Z+\beta Z^2/L,\qquad \epsilon=0.003,\quad\beta=0.001,$$
$$a=1+\epsilon+2\beta Z/L,\qquad F=\operatorname{diag}(1,1,a),\qquad J=a\in[1.002,1.004].$$

SVK is explicit: $\mu=1000\,\mathrm{Pa}$, $\lambda=0$, $K=666.6666666666666\,\mathrm{Pa}$. Physical boundary-node trapezoid volumes are retained. Phase, AT1, Gc, phase VI, stability, fracture, erosion and CFD are unused. No material defaults are used.

The grid expression is the previous round's exact floating representation $L=18h$. The decimal representation of the same length differs by one floating bit. An initial exact-coordinate storage check stopped before any regression result was evaluated; [IMPLEMENTATION_NOTES.json](IMPLEMENTATION_NOTES.json) preserves this correction and the original logs. No physical parameter or scientific threshold was adjusted.

## 3. Exact local continuum body force and traction

Independent closed forms, implemented without the general constitutive routine:

$$E_{zz}=\tfrac12(a^2-1),\quad S_{zz}=\mu(a^2-1),\quad P_{zz}=\mu(a^3-a),\quad W=\tfrac\mu4(a^2-1)^2,$$
$$b_z=-\mu(3a^2-1)\frac{2\beta}{L},\quad b_x=b_y=0,\qquad t=PN.$$

Top traction is $+P_{zz}(L/2)e_z$, bottom traction is $-P_{zz}(-L/2)e_z$ and side traction is zero. The formulas follow by differentiating the supplied SVK potential and the specified displacement field; their source is the user's [archived task, sections 4–6](TASK_REQUEST.txt), checked independently against the unchanged constitutive implementation. No external theorem about nonlinear BAQP boundary equivalence is asserted.

Eight-point Gauss quadrature is polynomial-exact for the body-force, energy and selected virtual-work integrands. All nine affine variations, three translations and the quadratic axial variation are checked.

Continuum force closure = 4.008214e-14 on the stress-area scale; max work closure = 3.999089e-14 on the stress-volume scale. See [EXACT_NONAFFINE_REFERENCE.json](EXACT_NONAFFINE_REFERENCE.json).

## 4. How the fully closed BAQP bulk oracle is constructed

The nonaffine adapter subclasses the already verified `AffineStencil`; its neighbor, endpoint volume, K, reconstructed F, bond stress, Z and force-assembly methods are inherited unchanged. Only reference translation and manufactured current nodal positions change. Each of the 19 absolute Z levels uses an N27 scaffold centered at that point; the polynomial extends to all required nodes. No global N27 bond state is allocated. D0/D1/D2/D3 contain 1/256/1797/5821 required nodes, with maximum coordinate offsets 0/4/8/12h. All required nodes have volume h³, and every required D2 family is complete. The scaffold has boundary nodes outside the dependency stencil; none participates in the evaluated center force. 

Three selected levels are repeated on translated N29 scaffolds with observation offset (1,-1,1); three reflected fields use beta→−beta, Z→−Z and the reflected force vector; three beta=0 cases recover zero bulk force. The selective result is also compared with full N19 evaluation at three nonaffine points. [BULK_STENCIL_CLOSURE.json](BULK_STENCIL_CLOSURE.json), [BULK_ORACLE_VALIDATION.json](BULK_ORACLE_VALIDATION.json), and [SELECTIVE_VS_FULL_NONAFFINE.json](SELECTIVE_VS_FULL_NONAFFINE.json) retain all checks. The nominal physical/oracle Z levels differ only by the one-ulp L representation described above (max 1.0842e-19 m); [Z_LEVEL_MAPPING_PRECISION.json](Z_LEVEL_MAPPING_PRECISION.json) records it explicitly. No interpolation is used.

## 5. Bulk BAQP force vs local Div P

$$q_{\rm bulk}(Z)=f_{\rm closed}(Z)/h^3,\qquad e_{\rm model}=q_{\rm bulk}+b_{\rm local}.$$

Both densities have units N/m³. Their difference is a finite-horizon bulk-model/discretization difference, not boundary contamination. The predeclared descriptive SMALL threshold is 0.1% of max local body-force density; it is not a general modeling tolerance. The second reference density is the characteristic stress-gradient scale $2\mu|2\beta/L|=3555.555555555556\,\mathrm{N/m^3}$; no small-strain approximation is used in the actual stress or body force.

| Quantity | Value |
| --- | --- |
| max_local_body_N_m3 | 3.598308e+03 |
| max_bulk_q_N_m3 | 3.598308e+03 |
| max_model_error_N_m3 | 7.521586e-05 |
| model_error_relative_to_body | 2.090312e-08 |
| model_error_relative_to_stress_gradient | 2.115446e-08 |
| BULK_MODEL_DIFFERENCE | SMALL |

The tiny variation of the mismatch across Z contains reconstructed-stress floating-point scatter. Its nonzero magnitude is resolved independently of the equilibrium solves. [BULK_BAQP_REFERENCE.csv](BULK_BAQP_REFERENCE.csv) saves all vector components.

## 6. Definition of the localized boundary correction

$$f_{{\rm bulk,target},i}=q_{\rm bulk}(Z_i)V_i,\qquad f_{{\rm BC,localized},i}=f_{{\rm bulk,target},i}-f_{{\rm FD},i}.$$

Physical $V_i$ includes face, edge and corner factors. This converts the physical finite-domain force into the bulk BAQP force at the specified surrogate. It is not a final physical traction law. Section 11 of the task defines the signed expression used here; section 20's reversed verbal subtraction is interpreted using that main definition.

The correction is frozen throughout equilibrium. The common local normalization is
$$P_{\rm ref}=\max_{\rm N19}\|P\|_F,\quad b_{\rm ref}=\max_{\rm N19}\|b\|,\quad F_{\rm ref}=\max(P_{\rm ref}h^2,b_{\rm ref}h^3).$$

F_ref = 3.143775000e-08 N. N27 uses this same physical-N19 reference scale. Beta=0 regression against the frozen old B1 gives relative L2 error 2.298770e-13, below 1e-10.

## 7. Why naive -r_FD is wrong for nonaffine fields

$$f_{\rm BC,naive}=-f_{\rm FD},\qquad f_{\rm BC,naive}-f_{\rm BC,localized}=-q_{\rm bulk}V.$$

The naive quantity is diagnostic only. It removes the legitimate interior force. At the fully closed N27 center its magnitude/F_ref is 2.786077e-02, while the localized correction is 3.196871e-15. This distinguishes interior loading from boundary conversion directly.

## 8. Nodewise residual decomposition

$$f_{\rm FD}+f_{\rm BC,localized}+b_{\rm local}V=e_{\rm model}V.$$

All nodes, six face distances, d/h, d/delta and six force channels are saved in [NODEWISE_RESIDUAL_DECOMPOSITION.npz](NODEWISE_RESIDUAL_DECOMPOSITION.npz). These are nodal forces in N; density arrays are named separately.

| Identity | max nodal norm / F_ref |
| --- | --- |
| Nodewise model decomposition | 4.556740e-17 |
| Naive minus localized | 4.604547e-17 |
| NL target equilibrium | 5.262340e-17 |

All prescribed 1e-10 gates pass. Algebraic closure alone is not evidence of a physical boundary law.

## 9. Global force / moment accounting

The correct force identity retains physical quadrature and assembly terms:
$$F_{\rm BC}-F_{\rm surface}=F_{\rm model}-(F_{\rm body,nodal}-F_{\rm body,exact})-F_{\rm internal}.$$

The user-requested comparison with F_model alone is therefore only exact up to physical body quadrature and floating arithmetic. Quadratic b is not integrated exactly by nodal trapezoid volumes. No volume or old manifest is altered to erase this difference. The composite trapezoid formula independently predicts $\Delta F_{\rm body,z}=-4\mu\beta^3h^2=-1.5625\times10^{-14}\,\mathrm N$ for this cube.

| Axial resultant term | N |
| --- | --- |
| surface_local | 5.108133e-06 |
| body_local_exact | -5.108133e-06 |
| body_local_nodal | -5.108133e-06 |
| BC_localized | 5.108133e-06 |
| model | 1.070725e-13 |
| body_quadrature_error | -1.562500e-14 |
| internal | 3.887749e-23 |
| BC_minus_surface | 1.226979e-13 |
| closure | 4.145559e-19 |

Full force closure/F_ref = 1.318656e-11. Reference-centroid moment accounting retains the same body-quadrature and internal terms; closure/(F_ref L) = 1.142021e-14. BAQP angular momentum refers to current-position moments. A zero reference-position internal moment is not imposed; this axial symmetric case is consistent with zero at floating precision. [FORCE_BALANCE.json](FORCE_BALANCE.json) and [MOMENT_BALANCE.json](MOMENT_BALANCE.json) report actual vectors, including the simple-identity discrepancy.

## 10. Virtual-work accounting

The existing affine/translation variation infrastructure is reused, with $\delta u_z=Z^2/L$ added. For all 13 modes:
$$\sum_i(f_{\rm FD}+f_{\rm BC,localized}+b_{\rm local}V)_i\cdot\delta u_i
=\sum_i(e_{\rm model}V)_i\cdot\delta u_i.$$

[VIRTUAL_WORK_DECOMPOSITION.csv](VIRTUAL_WORK_DECOMPOSITION.csv) records finite-domain internal, converted boundary, local body, corrected residual and bulk-model work. Corrected work is not forced to zero. Continuum work checks are separately retained in the exact reference JSON. Maximum discrete work closure/(F_ref L) is 4.985885e-15.

## 11. Boundary-correction spatial support

| Correction | d50 µm / h / δ | d90 µm / h / δ | d95 µm / h / δ | d99 µm / h / δ |
| --- | --- | --- | --- | --- |
| localized | 62.5 / 1 / 0.248756 | 125.0 / 2 / 0.497512 | 187.5 / 3 / 0.746269 | 250.0 / 4 / 0.995025 |
| naive | 62.5 / 1 / 0.248756 | 187.5 / 3 / 0.746269 | 250.0 / 4 / 0.995025 | 375.0 / 6 / 1.492537 |

These magnitude-weighted quantiles describe **CONVERSION SUPPORT, NOT PHYSICAL BC THICKNESS**. Mean/max profiles and cumulative magnitude use nearest distance to all six physical faces. N19 has no fully closed three-hop physical center stencil, so its near-zero tail cannot establish full bulk closure. The naive correction has a longer tail because it contains the physical bulk force.

## 12. N27 fully closed center result

Only specimen extent changes to 26h; the same physical polynomial retains L=1.125 mm. Predeclared points approach from a lateral x face at Z=0, holding the local field value constant while changing stencil depth. The classical traction on that face is zero; these four samples do not establish every face profile. An independent shifted N29 scaffold supplies their bulk reference.

| Offset | d/δ | Localized norm N | Localized/F_ref | Naive/F_ref | Full graph closed |
| --- | --- | --- | --- | --- | --- |
| [-9, 0, 0] | 9.950249e-01 | 2.887976e-12 | 9.186332e-05 | 2.795263e-02 | False |
| [-5, 0, 0] | 1.990050e+00 | 2.880797e-15 | 9.163496e-08 | 2.786086e-02 | False |
| [-1, 0, 0] | 2.985075e+00 | 9.636654e-23 | 3.065313e-15 | 2.786077e-02 | False |
| [0, 0, 0] | 3.233831e+00 | 1.005024e-22 | 3.196871e-15 | 2.786077e-02 | True |

The 12h point is numerically small but includes nonbulk endpoint volumes in its theoretical dependency graph. Only the 13h center has verified full closure. Numerical cancellation at an earlier depth is not a theorem about minimal sufficient closure or physical BC thickness.

## 13. BAQP-manufactured recovery

$$b_{\rm NL}=-q_{\rm bulk},\qquad f_{\rm ext,NL}=b_{\rm NL}V+f_{\rm BC,localized}.$$

This has a zero target residual but that substitution identity is not the recovery test. The unchanged existing mechanics adapter runs its existing TRF/LSMR solver from a nonzero smooth perturbation. Both z faces are fixed to the manufactured positions; all other particles are free. The prescribed perturbation vanishes on all faces. Frozen external volume forces include fixed nodes; constraint reactions are saved separately.

For nonaffine data, finite-domain reconstructed Fbar at the exact nodal target is not the analytic continuum F near a boundary. The predeclared solver recovery metric is max norm(Fbar(solved)−Fbar(target))/0.004; analytic-gradient error is reported separately and never used to hide reconstruction bias.

| NL metric | Value |
| --- | --- |
| initial_perturbation_max_m | 8.670532e-08 |
| force_residual | 2.030006e-11 |
| F_recovery_error_max | 1.628608e-10 |
| F_analytic_error_max | 2.975839e-02 |
| displacement_error_max_m | 1.893480e-16 |
| energy_relative_deviation | -2.838715e-13 |
| PASS | True |

Exact target Fbar-to-analytic bias = 2.975839e-02 of imposed-gradient scale. The recovery gates remain force<=1e-8 and discrete target-gradient error<=1e-5. This does not claim analytic F is exactly reconstructed at truncated nodes.

## 14. Classical-manufactured result

The second solve uses `b_local V + f_BC_localized` with the same frozen conversion, constraints and initial perturbation. Its target residual is independently `e_model V`. No negative model-error load is added to the boundary correction.

| Classical metric | Value |
| --- | --- |
| force_residual | 2.029948e-11 |
| F_recovery_error_max | 1.113802e-08 |
| F_analytic_error_max | 2.975839e-02 |
| displacement_error_max_m | 8.130976e-15 |
| energy_deviation_J | -2.005678e-20 |
| energy_relative_deviation | -1.509266e-09 |
| PASS | True |

| Cross-control diagnostic | Value |
| --- | --- |
| external_difference_equals_model_scaled | 8.996840e-17 |
| force_response_closure_relative_to_model | 2.462991e-04 |
| NL_to_classical_max_displacement_difference_m | 8.119207e-15 |
| NL_recovery_error_over_classical_deviation | 2.328724e-02 |
| classical_target_residual_existing_solver_scale | 3.402658e-09 |
| classical_target_max_residual_local_scaled | 5.841145e-10 |
| local_continuum_energy_J | 1.333181e-11 |
| BAQP_target_minus_local_continuum_energy_J | -4.271731e-14 |
| classical_minus_local_continuum_energy_J | -4.271733e-14 |

The nonzero target residual establishes that the local manufactured field is not an exact BAQP equilibrium. The load difference between controls is precisely the independently measured model residual. The target itself lies below the broad 1e-8 solver acceptance tolerance, so acceptance alone would not establish exactness; both independent perturbed-start solves continue to their own numerical termination. The cross-control force-response balance leaves 0.025% of the bulk-model load at numerical residual level. The NL recovery displacement error is 2.329% of the classical deviation; this comparison indicates numerical resolution, not a rigorous error bound. A converged classical solve is not proof that every last digit of its displacement deviation measures bulk nonlocality. The classical-field status refers to non-exactness under the measured bulk operator, with response precision quantified above. The primary energy deviation compares BAQP energies at solved and target nodal states. The additional analytic-continuum energy comparison retains finite-domain quadrature and nonlocal-energy differences; it is not solely a solver error.

## 15. What part is boundary physics?

The difference between the finite physical-domain force and the independently closed bulk force is a finite-domain conversion load on this surrogate. It decays to numerical zero at the fully closed center. That spatial and algebraic evidence supports localization. The load also contains discretization features of truncated reconstruction, which remain under the unchanged BAQP operator. Neither this definition nor manufactured recovery identifies a unique experimentally calibrated boundary law.

## 16. What part is finite-horizon bulk-model difference?

`q_bulk + b_local` is evaluated without physical-domain truncation and remains in the classical target residual. It is classified SMALL using the frozen threshold, independently of the boundary-localization verdict. Global differences also include the separately quantified physical volume-quadrature term. At fixed h and delta this round cannot separate nonlocal-horizon error from lattice discretization error, nor establish an asymptotic limit. Small bulk force-density mismatch does not imply identical finite-domain and local continuum energies. No kernel, normalization, horizon, material or force law is adjusted.

## 17. Optional shear result

**NONAFFINE_SHEAR=NOT_RUN**, as predeclared. This round completes the required axial oracle, identities, N27 localization and two equilibrium controls. No mixed traction or additional science campaign is run.

## 18. Implication for future curved clot-wall mechanics

A future surrogate-dependent conversion must retain bulk mechanics and compatible body force explicitly. An affine negative-residual prescription cannot be transferred unqualified. The present Cartesian quadratic axial benchmark does not validate pressure/shear conversion on a curved clot, wall adhesion, delamination, fracture or CFD coupling. No surface operator correction, ghost treatment, perifit weights, optimized weights or variable horizon is implemented. Historical C5/C7 states and claims are unchanged.

## 19. One next task only

After human review: **NONLOCAL_BOUNDARY_NONAFFINE_SHEAR_LOCALIZATION_V1** — apply this same already specified decomposition to the quadratic shear field with analytic traction on every face, retaining the independently measured bulk residual and reporting solver resolution limits. Do not begin it automatically.

| Status | Value |
| --- | --- |
| NONAFFINE_AXIAL_REFERENCE | VERIFIED |
| BULK_BAQP_NONAFFINE_ORACLE | VERIFIED |
| AFFINE_LIMIT_RECOVERY | PASS |
| NODEWISE_DECOMPOSITION | PASS |
| NAIVE_NONAFFINE_CANCELLATION | CONTAINS_INTERIOR_PHYSICS |
| BOUNDARY_LOCALIZATION | VERIFIED |
| N27_CLOSED_CENTER_BOUNDARY_CORRECTION | PASS |
| BULK_MODEL_DIFFERENCE | SMALL |
| NL_MANUFACTURED_RECOVERY | PASS |
| CLASSICAL_MANUFACTURED_FIELD | APPROXIMATE_DUE_TO_BULK_NONLOCALITY |
| NONAFFINE_SHEAR | NOT_RUN |
| BOUNDARY_DATA_CONVERSION_STATUS | NONAFFINE_LOCALIZATION_VERIFIED |
| SURFACE_OPERATOR_CORRECTION | NOT_MODIFIED |
| C5_C7_ACTIVATION | NOT_RUN |
| FRACTURE | NOT_RUN |
| CLOT_GEOMETRY | NOT_RUN |

Protected existing files: 251; changed: 0. Only the nine requested tests are run; [TARGETED_TESTS.log](TARGETED_TESTS.log) records their result. No old assertion or manifest is modified. The initial invocation, exact-storage precondition interruption and prerequisite check are preserved; none accepted a scientific state. [IMPLEMENTATION_NOTES.json](IMPLEMENTATION_NOTES.json) distinguishes those implementation issues from the scientific gates.

Resources: 8 CPU oracle workers; 1 serial equilibrium worker; BLAS 1; GPU_NOT_USED_BY_DESIGN. Measured scientific stage time 1789.74 s, peak process-tree RSS 2.211 GiB (shared pages counted conservatively), minimum available RAM 52.7%. No packages installed. The two large equilibrium solves never overlap.

Five figures: [600 dpi PNG and vector PDF](../../results/nonlocal_boundary_nonaffine_localization_v1/figures). Array units and legacy solver-field conventions are in [DATA_DICTIONARY.md](DATA_DICTIONARY.md). Reproduce numerical stages with `PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 <python> -B scripts/nonaffine_boundary_audit.py --stage <oracle|decomposition|n27|NL|classical>` in that order. Existing outputs are the accepted record; use a separate output copy for reruns. Plotting is `scripts/plot_nonaffine_boundary.py`.

**Stopped for human review. No subsequent scientific development is run.**
