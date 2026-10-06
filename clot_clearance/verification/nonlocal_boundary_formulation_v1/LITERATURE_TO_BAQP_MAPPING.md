# Literature-to-BAQP mapping

Full-text reading preceded implementation. The four uploaded PDF files themselves were not accessible in the supplied attachment directory. Public primary full texts were therefore located by title: arXiv versions for D'Elia and Yu, the Scabbia publisher article reproduced in its full-text repository page, and the Shojaei publisher PDF hosted by TU Hamburg. This is not a byte-level authentication of the uploaded copies. Local PDF SHA256 values and access paths are recorded in [SOURCE_PROVENANCE.json](references/SOURCE_PROVENANCE.json). Raw PDFs are cached outside the repository; no PDF-only conclusions are silently presented as new numerical results.

## Source contributions and limits

| Source and sections read | Contribution used | What was not transferred |
|---|---|---|
| [D'Elia, Tian, Yu, 2020; arXiv full text](https://arxiv.org/pdf/1906.04259), §2 interaction domain, §3.1 Eqs. 13–15, §3.2 Eq. 16, Remark 3.1 | A local surrogate supplies nonlocal exterior volume data. Neumann data use a nonlocal flux/operator evaluation; Dirichlet data prescribe the surrogate. | The scalar operator, its exterior-domain equations and convergence proof are not BAQP elasticity theorems. Our B1 evaluates an interior finite-domain residual; B2 constrains interior physical support. Both differences are explicit adaptations. |
| [Yu, You, Trask, 2021; full text](https://arxiv.org/pdf/2101.01434), §3, Eqs. 3.1–3.3, §5.1 patch-test setup | Missing exterior interactions motivate traction treatment; sharp-surface normals/tangents and volumetric interaction regions are distinct. | LPS dilatation, its correction tensor, constants and Eq. 3.3 traction terms are model-specific. None is inserted into BAQP. Published asymptotic compatibility does not establish ours. |
| [Scabbia, Zaccariotto, Galvanetto, 2021; full publisher text](https://www.researchgate.net/publication/352787971_A_novel_and_effective_way_to_impose_boundary_conditions_and_to_mitigate_the_surface_effect_in_state-based_Peridynamics), §2.2, §3 Eqs. 30–34, §4.1–4.2 Eqs. 42–54; [author repository record](https://www.research.unipd.it/handle/11577/3399436) | Fictitious kinematics follow Taylor continuation without independent DOFs; force-flux interpretation helps separate boundary approximation from surface effects. The construction also specifies weighted volume and extrapolated dilatation. | Their OSB formulas do not define BAQP ghost-side shape tensors, bond-associated gradients, Z fields, physical energy or reaction pullback. Affine kinematic extrapolation alone is insufficient. B3 is deferred. |
| [Shojaei et al., 2026; publisher](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.70385), [full PDF](https://tore.tuhh.de/bitstreams/8fda7c78-cb40-401e-b8ce-206cd986fafb/download), §3–4, Eqs. 44–47, 50, 54 | Symmetric endpoint weights fit truncated operators to full-horizon PD targets. The demonstrated BB/OSB targets use model-specific polynomial moments and distance powers. | BAQP/NOSB validation is not supplied by perifit. Multiplying existing BAQP final forces by optimized weights would not establish energy consistency. No weights are fitted here. |

## This repository's discrete reasoning

Let \(\mathcal E_h\) be the already audited BAQP intact energy and \(r_h=-\nabla_u\mathcal E_h\) its unchanged internal nodal force. For an analytic local solution \(u_L\), define a frozen conversion \(f_h=-\chi_{\mathcal S}r_h(u_L)\), where \(\mathcal S\) is the recorded dependency support. This maps a known local solution to a discrete force representation; its sign follows our internal-force convention. It is not a direct transcription of an exterior scalar flux formula.

If the affine BAQP energy equals \(|\Omega|W(F)\), differentiating it along affine variations gives

\[
\sum_i [-r_{h,i}(u_L)]\cdot A(X_i-X_c)=|\Omega|P(F):A.
\]

Therefore affine work consistency is expected from the validated variational mechanics, provided the support retains the required force. It remains useful for finding sign, quadrature, omitted-face and support mistakes; it is not independent proof of a general nonlocal Neumann problem. The continuum target is evaluated independently from \(PN\), rather than inferred from the discrete force.

Translation gives the net-force test. Antisymmetric affine variations relate to the reference-coordinate moment. Nine general affine variations test the full first-moment tensor, including shear entries. These comparisons go beyond the pointwise cancellation identity, but they still span only affine test functions.

The finite-domain adaptation is most transparent when the load retains nonzero support far from the boundary. That case must be reported as surrogate-dependent volume loading, without silently reinterpreting every force as a face traction. A future nonaffine test must distinguish this particular-load construction from a transferable physical boundary law.

## Separation of responsibilities

| Layer | This round |
|---|---|
| Bulk constitutive energy and BAQP variation | Existing sources retained and SHA checked |
| Boundary data conversion | New explicit geometry, analytic surrogate, frozen force conversion and independent checks |
| Surface/operator correction | Separate question; no changed kernel, K/F/Z weights or constitutive law |
| Physical clot-wall interface | Future displacement/traction input only; no adhesion, separation or fracture law |

One-horizon B2 is deliberately a strong **in-domain** reference; no theorem of equivalence to exterior volume Dirichlet is asserted. A wider constraint chosen after seeing the answer would defeat this comparison. B3 would require a derivation of coefficients and virtual-work-consistent ghost elimination before it could be called an independent BAQP control.
