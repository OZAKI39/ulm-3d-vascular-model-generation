# External dependencies and literature provenance

No third-party clone, environment, package cache, build tree, or paper PDF is newly synchronized.

| Resource | Exact identity / paths | Use / provenance |
|---|---|---|
| [Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl) | `b6e86563845a83b5cb06c805af19a817bde6b106`, 0.6.0-DEV; `src/` RKC routines `rkc_weights!`, `rkc_defgrad!`, and `RKCMaterial(monomial=:PD2)` | MIT, copyright Kai Partmann. Installed source remains local under `clot_clearance/verification/pfpd_candidate_v1/external/Peridynamics.jl`. Source-file hashes and specific RKC paths are recorded in the RK audit's `JULIA_PACKAGE_SOURCE_HASHES.json`; package source was unchanged during scientific comparison. |
| Julia | 1.10.12 | Recorded version used by the accepted independent gradient check; no package update or Julia rerun in handoff. |
| PETSc / petsc4py | 3.25.3 / 3.25.3 | Installed binary environment queried during handoff. Historical VI solver provenance remains in `clot_clearance/verification/pfpd_at1_petsc_vi_solver_v1/provenance/`; not a PD2 intact dependency. No unverified source commit is invented. |
| Hillman–Pasetto–Zhou, generalized RK peridynamics | §§4.1–4.2, Eqs. 63–70, 77–81 | Literature identity, URL, original PDF hash/size: [RK source provenance](../clot_clearance/verification/higher_order_baqp_rk_oracle_v1/references/SOURCE_PROVENANCE.json). |
| Behzadinasab–Trask–Bazilevs Part I | arXiv:2004.11477, §2.2, Eqs. 3/5, Example 1 | Complete 3D quadratic derivative space. Part II arXiv:2004.11478 is future reference only. No claim to transplant the entire dynamical formulation. |
| Boundary literature | D'Elia–Tian–Yu; Yu–You–Trask; Scabbia et al.; Shojaei et al. | [Source identity](../clot_clearance/verification/nonlocal_boundary_formulation_v1/references/SOURCE_PROVENANCE.json) and [implementation mapping](../clot_clearance/verification/nonlocal_boundary_formulation_v1/LITERATURE_TO_BAQP_MAPPING.md). BB/OSB/perifit algorithms were not silently transferred into BAQP. |

This file summarizes already accepted provenance; no new literature-derived scientific claim or external theorem is introduced by the handoff.
