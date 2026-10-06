# Data dictionary and reproduction

All data describe fixed exact target configurations. No equilibrium or time history was generated.
Python IDs are zero-based. Tensor axes are node, displacement component, physical derivative.
The complete array shapes/dtypes and CSV headers/row counts are in [DATA_SCHEMA.json](DATA_SCHEMA.json).

| Files | Meaning / units |
| --- | --- |
| REFERENCE_GEOMETRY.npz | X in m; volume in m³; undirected pairs; directed source/target/ptr; saved energy weights; node/bond depth divided by δ; geometric classes |
| PD2_REFERENCE_CACHE.npz | M in m³, nine descending singular values in m³, rank, condition, coefficients in m⁻⁴, directed g in m⁻¹ including Vj |
| RK2_DIAG_REFERENCE_CACHE.npz | Same conventions for the six-term negative control |
| PD2_NODE_CONDITIONING.csv | One record per node, all singular values, condition flag and family size |
| PD2_CONDITIONING_GROUPS.csv | All-node, four geometric-class, and nearest-face-depth statistics |
| POLYNOMIAL_REPRODUCTION.npz / PD2_POLYNOMIAL_REPRODUCTION.csv | Exact/computed physical gradients and residuals; linear scale 1, quadratic scale L in m; all 6,859 × 9 entries |
| PD2_AFFINE_REGRESSION.npz / .json | Two affine computed/exact nodal gradients and normalized errors |
| *_NODAL_STATE.npz | X, x in m; analytic F and C1/PD2/RK2_DIAG nodal F, dimensionless |
| PD2_AXIAL_QUADRATIC.csv / PD2_SHEAR_QUADRATIC.csv / PD2_CROSS_TERM.csv | All nodes for all three reconstruction spaces; F error divided by 0.004 |
| NODE_GEOMETRY_STATISTICS.csv | Maximum and physical-volume weighted RMS per field/model/geometric class |
| F_DEPTH_PROFILES.csv | All-node/all-bond layer reductions; Fbar uses nodal volume, Ftilde uses original energy weights |
| *_BOND_DIAGNOSTICS.npz | All 690,020 bonds: F error is nonnegative Frobenius norm; Psi error is signed computed-minus-exact density in Pa |
| *_NODE_PSI.npz | Analytic positive-energy density and signed C1/PD2 errors in Pa |
| PD2_FTILDE_DIAGNOSTIC.csv | Bond-layer F errors, absolute and divided by 0.004 |
| PD2_PSI_PLUS_DIAGNOSTIC.csv | Node/bond-layer Psi errors; max and weighted RMS are absolute magnitudes normalized by the single analytic reference per field |
| C1_PD2_SUMMARY.json | Global max and weighted RMS, absolute and normalized; per-field Psi_ref in Pa |
| julia_crosscheck/PYTHON_REPRESENTATIVES.npz | Only 4 representative nodes and their 551 directed weights, not a duplicate of the entire reference cache |
| julia_crosscheck/*_F.csv | Four × nine row-major tensor entries; Python/Julia pairs |
| julia_crosscheck/julia_gradient_weights.csv | No header: original i/j in one-based indexing, then physical gx, gy, gz |
| PERIDYNAMICSJL_CROSSCHECK.json | Four representative tensors per field, exact package provenance and weight comparison |
| PROTECTED_INPUT_HASHES.json | Pre-run hashes of 496 protected historical files; old manifests remain untouched |
| OUTPUT_MANIFEST.json | SHA256 and byte size of every file in the three new namespaces except the manifest itself |

All layer reductions use the full relevant population, without random sampling. Geometry classes
count coordinate axes with nearest opposite-face distance < δ; they do not assert force-stencil closure.
The P2 cache is reference-only and displacement-independent. Cache/config hashes are recorded in
[CACHE_PROVENANCE.json](CACHE_PROVENANCE.json); input-array identity is in [INPUT_IDENTITY.json](INPUT_IDENTITY.json).

Review figures use only saved outputs. To regenerate them from the project root:

```bash
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B \
  verification/higher_order_baqp_rk_oracle_v1/plot_figures.py
```

The scientific computation entry point is `run_audit.py`; the independent installed-package comparison
is `run_julia_crosscheck.py`. These are retained for reproducibility and were not rerun during report
generation. The 12 focused tests are in `tests/test_pd2_oracle.py`; their completed result is recorded
in [TARGETED_TESTS.json](TARGETED_TESTS.json). No full project suite was run.

Reference PDFs are retained with source URLs and SHA256 in [SOURCE_PROVENANCE.json](references/SOURCE_PROVENANCE.json).
The current audit does not replace or edit historical scientific manifests, active context, or results.
