# perifit / BAQP compatibility — review only

Status: `PERIFIT = COMPATIBILITY_REVIEW_ONLY`; `BAQP_SPECIFIC_SURFACE_CORRECTION = NOT_IMPLEMENTED`.

The requested filename uses `PERFIT`; the package name discussed here is **perifit**. Source links and the paper's scope are in [LITERATURE_TO_BAQP_MAPPING.md](LITERATURE_TO_BAQP_MAPPING.md).

Shojaei's Eqs. 45–47 and 54 distinguish BB bond-stretch/energy targets with distance power 3 from OSB dilatation/deviatoric targets with distance power 2. Equation 44 symmetrizes endpoint weights; Eq. 50 fits discrete operators to full-horizon PD polynomial targets. Those targets do not define a BAQP correction. [Full primary PDF](https://tore.tuhh.de/bitstreams/8fda7c78-cb40-401e-b8ce-206cd986fafb/download).

The following is a **proposed compatibility checklist derived from our BAQP energy**, not a method claimed by that paper:

1. Define a reference quadrature and weighted shape tensor \(K_i^w\), retaining invertibility and exact affine reproduction. A scalar factor at a node may cancel from a normalized shape reconstruction; changes to other energy weights then remain. Each occurrence must be derived rather than adjusted independently.
2. Match homogeneous total energy to \(|\Omega|W(F)\) for enough independent normal and shear gradients. Passing a single axial test cannot determine tensor consistency.
3. Derive the force as the first variation of the **same weighted energy**, including all paths through the bond-associated gradient and Z-mediated terms. Numerical energy/force directional derivatives and reciprocal tangent tests are necessary.
4. Preserve translation/rotation invariance, pairwise action–reaction assembly, reference-traction work, and a meaningful full-horizon limit. Symmetric endpoint weights alone cannot guarantee all of these for BAQP.
5. Separate geometric quadrature weights from material parameters and phase-dependent quantities. Keep any optimized weights frozen with respect to deformation, unless their derivatives are explicitly included.
6. Use full-horizon **BAQP-to-BAQP** targets and independent withheld deformation modes. Do not fit weights to C5/C7 activation strains, boundary reaction alone, or an assumed classical local limit.

Could one scalar nodal weight suffice? That is **not established**. Symmetric nonnegative pair weights can be built into a variational energy, but a scalar family may be unable to satisfy all independent tensor moments, energy and Z-mediated force requirements near faces/edges/corners. A post-hoc multiplier on final force generally has no demonstrated common energy; it is not an acceptable shortcut.

Boundary conversion and operator correction solve different design problems. B1 leaves the operator unchanged and changes external data. A hypothetical optimized operator would change the discrete energy and its variation and would need a new scientific audit. This round supplies no evidence that such an operator change is required and performs no optimization.
