# Candidate data dictionary

All arrays use SI units and zero-based node/component indexing. Force is the actual nodal force in N,
not force density. F axes are node, output component, reference derivative. No phase variables exist.

| File group | Contents |
| --- | --- |
| CANDIDATE_CONFIGURATION.json | Frozen steps, states, scales, gates, material, primary cache hash and resource policy |
| PD2_CACHE_IDENTITY.json | Original N19 cache/config/g identities; new N7 verification-cache provenance |
| PD2_KINEMATIC_REGRESSION.json | Full-node previous-oracle and analytic comparisons for all three quadratic fields |
| N7_REFERENCE_GEOMETRY.npz | Reference coordinates m, volumes m³, topology, original kernel, directed indices |
| N7_PD2_REFERENCE_CACHE.npz | One-time small-grid moment/weight cache from the unchanged RK routines; g in m⁻¹ includes Vj |
| C1_ENERGY_REGRESSION.json | Historical and generic-C1 energies J in the three small-grid states |
| C1_FORCE_REGRESSION.npz | All historical/generic C1 force vectors N, actual x in m, X and physical volumes |
| C1_FORCE_REGRESSION_SUMMARY.json | Vector maximum and V-weighted RMS errors, stable force scales, component differences |
| DELTA_FTILDE_IDENTITY.json | Analytic/direct-perturbation comparison; derivative units m⁻¹ for unit dimensionless v |
| CONSTITUTIVE_STRESS_CHECK.json | Four prescribed F tensors, analytical P and componentwise central FD in Pa |
| *_STATE.npz | Actual target x in m, Fbar, variational force N and energy J; six N19 states plus selected N7 state |
| OBJECTIVITY_RESULTS.json | Rigid F/E/U/f metrics, inherited noise floor and deformed-state rotational covariance |
| GLOBAL_BALANCE_RESULTS.csv | Complete force/current-position moment vectors and normalizing scales for 11 non-rigid validation states |
| FD_*_DIRECTION.npz | Six saved deterministic v arrays; unit Frobenius norm and zero arithmetic mean translation |
| ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv | All 18 values: physical epsilon, U± J, Dfd and −f·v N, stable denominator, normalized/unfloored errors |
| ENERGY_FORCE_FD_SUMMARY.json | Six middle-step gates and full three-step sequences, without step selection |
| SELECTED_DOF_FORCE_ORACLE.csv | Twelve coordinate finite differences; U±, force components and F_ref are included |
| AFFINE_ENERGY_COMPARISON.json | Two N19 affine energy and full-population bondwise Ftilde/energy-density comparisons |
| QUADRATIC_ENERGY_DECOMPOSITION.csv | U_C1, U_PD2, U1, U2 in J; kinematic correction and residual quadrature difference; 32/64 Gauss check |
| CALL_TIMINGS.json / PERFORMANCE_SUMMARY.json | Final implementation timings; sampled RSS and process high-water RSS; thread/BLAS policy |
| TARGETED_TESTS.log / .json | Only the twelve requested focused tests |
| PROTECTED_INPUT_HASHES.json | Pre-run SHA256 of 600 prior files; unchanged historical manifests remain authoritative |
| OUTPUT_MANIFEST.json | SHA256 and bytes of this round's new files only; excludes the manifest itself |

The original N19 geometry and full weight cache remain at
`../higher_order_baqp_rk_oracle_v1/`; they are referenced and checked, not copied or overwritten.
The exact coordinate/volume/family identity and m_b match are checked in the driver.

Reference state -> frozen gradient -> nodal F -> bond projection -> old intact constitutive energy
-> one bond pass + one directed adjoint pass -> nodal force.
No stress-divergence shortcut, old-Z/PD2 control, Hessian, equilibrium or boundary conversion is used.

`run_audit.py` is the numerical gate driver. `plot_figures.py` and `write_report.py` only consume saved
outputs. Run from the project root with PYTHONPATH=., OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1,
MKL_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1 and the existing shared Python environment.
The `tests/test_variational_candidate.py` suite is the only pytest target; disable the pytest cache.

Figures use all requested saved data. The force scatter uses all small-grid force components;
the FD scatter includes every step of the six declared directions. All text is English, PNG is
600 dpi, and PDFs contain vector graphics. No random data subsampling is used.

`assembly_implementation_check/` preserves the first passing implementation's small summaries and
log before replacing block-sized full-node reduction buffers by direct index accumulation.
Those files are explicitly development evidence; the top-level outputs correspond to final code.

Full NPZ shapes/dtypes and top-level CSV schemas are in [DATA_SCHEMA.json](DATA_SCHEMA.json).
