# REDUCED_NEAR_WALL_HYDRODYNAMICS_V0 — PHASE 1R

**STATUS: PASS for reference resolution and standalone implementation/algebra. PHYSICS_REFERENCE_STATUS: PARTIAL for the combined RMBW resistance/shear kernel.** The standalone C++ kernel was implemented, compiled and tested after the reference gates passed. No Phase 2 or engine integration was performed. The previous `20260916_190052` stage remains immutable and BLOCKED.

**Scope and location.** Rigid sphere, infinite rigid no-slip plane, Stokes/overdamped response, gap-dependent excess resistance and analytic wall-tangential simple shear only. Local archive: `/home/lzy/projects/compre_output/reduced_near_wall_hydrodynamics_v0/20260916_192326`. Vast results: `/workspace/microbubble_lammps/results/reduced_near_wall_hydrodynamics_v0_20260916_192326`. CPU build: `/workspace/microbubble_lammps/work/reduced_near_wall_hydrodynamics_v0_20260916_192326/build`. All new code and evidence are confined to this stage and its new work directory. Original request texts, baseline identities, build commands and mirror receipt are retained under provenance/.

**Reference choice.** The new primary physics basis is [Chaoui & Feuillebois 2003](https://doi.org/10.1093/qjmam/56.3.381), Table 17 and Eqs. 6.7/6.8, with the checked [2012 corrigendum](https://doi.org/10.1093/qjmam/hbs012). It supplies the requested high-precision free-sphere shear representation and discusses discrepancies in the older GCB printed translation expression. GCB is now HISTORICAL_ASYMPTOTIC_REFERENCE. Its old printed expressions never enter the kernel or frozen closure. The prior printed-expression/Table3 discrepancy is retained as prior evidence, without changing constants to force agreement.

**Acquisition and correction.** The 2003 publisher PDF endpoint returned 403; HAL's record was a bibliographic notice, not a deposited PDF. The user then supplied `/home/lzy/projects/560381.pdf`, a published library copy with an Oxford/Universitaetsbibliothek Bern watermark. Its identity and required pages were inspected visually. SHA256: `2d572961b10cc7b46ae2a666f8dfc80123ee303b118d964c510ae86a739dc3e5`. Hierarchy B applies to this material. For 2012, the complete single-page publisher GIF was downloaded from the official notice HTML and visually checked. SHA256: `69c6cbb350127dd8b334998b04276641cf376b550e25c34ba2514078932265b8`. This is hierarchy A primary page imagery; the PDF itself was not obtained and no PDF hash is fabricated.

The corrigendum corrects only the opening §2.1 field representation and explicitly leaves the rest unaffected. Impact on Eqs. 6.1, 6.2, 6.7, 6.8, 6.9, Table 17, both §6 normalizations and free-sphere shear results: **DOES_NOT_AFFECT_CLOSURE**. No closure coefficient correction is needed. The full itemized assessment is reference/CF2003_CORRIGENDUM_AUDIT.md.

**Table 17 and normalization.** Primary pages 381, 383, 402–404 and 407–408 were rendered with the existing PyMuPDF 1.24.10 environment and visually reviewed. TRANSCRIPTION_A is manual visual entry from p407; TRANSCRIPTION_B is separate primary-PDF text extraction. All 40 literal strings match, including signs, digits and exponents: 20 coefficients per column, degree 19, ascending j=0..19. The same assistant performed the two extraction methods; these are not separate reviewers or separate theories. The initial incomplete text-only record is retained as raw evidence and is superseded by the fully populated CF2003_TABLE17_VERIFICATION.csv.

Sphere radius a; center height l=d; surface gap h=l-a; epsilon=h/a=l/a-1. FU=Ux/(kappa*l), FOMEGA=Omega_y/(kappa/2), with far-field normalization both tending to 1. For n=+z and g=+kappa*x, positive angular velocity points +y by the right-hand rule. The fitted reciprocal polynomials are evaluated at x=ln(epsilon) using Horner. Their small-gap representation is never extrapolated to the far field.

**Certified reference domain: [0.001,0.2].** The p407 applicability discussion places the small/large approximation crossover at 0.4; the project uses the smaller requested domain, shared with the RMBW lower limit. The paper reports fit errors below 1e-11 for the small-gap representation; this is source-reported fit accuracy, not an independent measurement of exact physical error here. All 10,009 dense/reference points are finite, positive and below one, with increasing FU/FOMEGA. A 70-digit Bernstein coefficient check over 64 log-gap subintervals additionally gives inverse polynomials >1 and negative derivatives; these are numerical bounds, not an interval-arithmetic proof. Horner versus 70-digit evaluation differed by at most 1.904e-14 absolutely.

| epsilon | FU | FOMEGA | Omega*a/U |
|---:|---:|---:|---:|
| 0.001 | 0.396956188956 | 0.427254784676 | 0.537626025988 |
| 0.002 | 0.428551683205 | 0.458862765646 | 0.534295966915 |
| 0.003202 | 0.452922369965 | 0.482986049115 | 0.531486735772 |
| 0.005 | 0.478550714882 | 0.508121190657 | 0.528254591218 |
| 0.01 | 0.524172440352 | 0.552308702291 | 0.521622520759 |
| 0.02 | 0.578022652014 | 0.603645244088 | 0.51192549355 |
| 0.05 | 0.663652448974 | 0.683839503151 | 0.490675291181 |
| 0.08 | 0.713858743604 | 0.730308218375 | 0.47363103651 |
| 0.1 | 0.738891202506 | 0.753397617746 | 0.463469400435 |
| 0.2 | 0.818796424219 | 0.827117666305 | 0.420901155362 |

**Historical diagnostics.** At epsilon=.001, differences from transcribed GCB Table3 (.3966,.4268) are approximately 0.08981% for FU and 0.10656% for FOMEGA. These exceed rounding alone and are consistent with a lower-order historical asymptotic comparison; they are not used as acceptance tolerances. The rolling ratio is 0.537626 at .001 and rises toward the historical contact value 0.5676 as the gap decreases within the domain. No exact equality at finite gap or extrapolation to contact is claimed.

The corrected first-order translation formula has relative differences 4.3202e-4, 6.5616e-5 and 4.6680e-4 at the three required points below .005. On the denser grid its largest difference just below .005 is **1.3437e-3**. Thus the user's approximate order-1e-3 cross-check is consistent, while a strict uniform <=1e-3 diagnostic is explicitly FALSE. An initially coded strict assertion was removed because it was stronger than the requested approximate historical check; neither primary coefficients nor certified domain was changed. The old rotational expression is a low-order diagnostic, with about 6.30% discrepancy by .2; it is not a precision oracle.

**Kernel and resistance interpretation.** src/reduced_wall_types.hpp, reduced_wall_kernel.hpp and reduced_wall_kernel.cpp implement explicit input/result types, an HDF5 read-only table loader, resistance evaluation and analytic shear RHS evaluation. The generated coefficient header comes from the frozen JSON's checked decimal literals and records its hash. Resistance returns full global 6x6 SI excess and total matrices, plus the bulk diagonal. DOF order is [Vx,Vy,Vz,Omega_x,Omega_y,Omega_z]; conjugate external loads are [Fx,Fy,Fz,Tx,Ty,Tz], with torque about the sphere center.

The frozen RMBW file has 8,211 rows and SHA256 `71e62ef6bc22edc66ae0c4d1dcc224378c2f6c5d3848578f3fbedd2cd9a43f98`. Its stored R_wall_excess_scaled already equals R_total_scaled-I: bulk is added exactly once. S=diag(sqrt(6*pi*mu*a) x3,sqrt(8*pi*mu*a^3) x3); local SI excess=S E S. TT, TR and RR scale as mu*a, mu*a^2 and mu*a^3. Q=[t1 t2 n] must be orthonormal and right-handed; T=diag(Q,Q); global=T local T^T. All components use linear interpolation in ln(epsilon). The resistance data-query range is [.001,20], while certified shear is only [.001,.2]. There is no clamping, exterior zero fallback or extrapolation.

The RHS is exactly the contracted architecture: g_t=(I-n*n^T)g; q_free=[d*g_t,0.5*n cross g_t]; q_target=[FU*d*g_t,0.5*FOMEGA*n cross g_t]; b_wall_shear=R_total*q_target-R_bulk*q_free. Tests actually solve R_total*q=R_bulk*q_free+b_wall_shear. The kernel never overwrites solved velocities. No wall-normal force is added. ANALYTIC is the only accepted shear source; GRADIENT_PROXY, PALABOS_VALIDATED and unknown enums reject explicitly. Nonpositive h returns INVALID_GAP; nonfinite input, invalid frames and unsupported domains reject with named statuses.

**Validation.** C++ tests passed for 10,681 valid cases and 27 rejected invalid cases. K covers all 10,009 dense certified reference points. Additional probes cover 96 deterministic random geometries/radii/viscosities/shears, tangent rotations/sign reversal, zero/reversed/doubled shear and radius/viscosity scaling. Input generation documents a <=1 ULP SI-gap fixture adjustment at closed epsilon endpoints; the kernel itself never clamps.

| Gate | C++ worst / result | Independent Python worst / result | Limit |
|---|---:|---:|---:|
| A | 1.74e-16 | 1.74e-16 | 1e-10 |
| B SPD / positive power | PASS | min scaled eigenvalue 1.08249722 | positive |
| C | 5.88e-16 | 5.88e-16 | 1e-10 |
| D | 7.16e-15 | 7.04e-15 | 1e-09 |
| E | 0 | 6.95e-15 | 1e-09 |
| F | 0 | 0 | 1e-12 |
| G | 0 | 0 | 1e-10 |
| H | 0 | 0 | 1e-10 |
| I | 7.25e-14 | 7.25e-14 | 1e-12 |
| J | 8.42e-15 | 3.19e-15 | 1e-08 |
| K | 6.66e-16 | 6.66e-16 | 1e-08 |

The independent validator separately reads the original HDF5, uses NumPy interpolation/tensor rotation, reconstructs bulk and shear RHS, checks scaled eigenvalues, and performs NumPy linear solves of both its own assembly and C++ raw matrices/RHS. Maximum independent excess-matrix relative error: 4.369e-16; shear RHS error relative to characteristic load: 2.15e-14. All 27 C++ rejection statuses matched the input contract and 12 Python reference exterior/nonfinite queries rejected. Python and C++ agree as **independent implementations of one CF2003 physics basis**; K is only assembly consistency.

Build used existing GCC 13.3, CMake 3.28 and HDF5 1.10.10. The installed HDF5 is built with MPI support, so its existing MPI headers/libraries are linked; the executable is run as one ordinary CPU process, without mpirun, GPU work or engine integration. No packages were installed. Build/run logs and compiled-source checks are in provenance/.

**Artifacts and limits.** Frozen reference JSON/CSV/SHA256, the Python evaluator, two coefficient extractions, source page images, A–K validation JSON, all raw CSVs and six static PNGs are archived. FU_CF2003_vs_epsilon.png, FOMEGA_CF2003_vs_epsilon.png, rolling_ratio_vs_epsilon.png and GCB_vs_CF2003_diagnostic.png show reference behavior and the unequal source hierarchy. resistance_vs_epsilon.png and shear_rhs_components_vs_epsilon.png show resistance and actual C++ RHS. At the flat-axis plot configuration, normal force is exactly zero in the saved data; tilted-frame residuals remain below the contracted roundoff gate.

RMBW finite-gap parallel RR accuracy, continuous-gap TR/RT accuracy and its historical rotation-constant discrepancy remain unresolved. The table also has sharp inherited transitions immediately above epsilon .01 and .1; these appear in the RHS and are recorded in validation/RMBW_TRANSITION_DIAGNOSTIC.json. No smoothing or source alteration was performed. Reference certification of FU/FOMEGA and algebraic recovery do not validate every RMBW resistance coefficient or force derivative. Therefore REFERENCE_SOURCE_VALIDATED=YES, ALGEBRA_PASS=YES, INDEPENDENT_IMPLEMENTATION_PASS=YES, free-shear reference status=PASS, and combined PHYSICS_REFERENCE_STATUS=PARTIAL.

**Preservation and stop boundary.** Before/after baseline checks retain the expected Git branch/HEAD and only the two allowed untracked environment reports. All 34 inherited key paths exist, and FrozenFlow/STL identities match. Old-stage file inventories and the final local/remote payload receipts are retained in provenance/. The original user PDF remains unchanged. No LAMMPS/Palabos modifications, real vascular run, WSS, RBC, adhesion, BIE/QBX, fixed-multiblob rerun, GPU rewrite, Git commit or push occurred. Phase 2 is not authorized; the task stops here for review.
