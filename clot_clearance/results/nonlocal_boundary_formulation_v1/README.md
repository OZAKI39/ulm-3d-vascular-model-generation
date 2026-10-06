# NONLOCAL_BOUNDARY_FORMULATION_V1

Review: [REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md](../../verification/nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md).

Five figure pairs in `figures/` are generated from the saved benchmark JSON/NPZ/CSV files. Status: `VERIFIED_FOR_INTACT_AFFINE_MECHANICS`.

| Figure | PNG 600 dpi | Vector PDF |
| --- | --- | --- |
| boundary_formulation_schematic | [PNG](figures/boundary_formulation_schematic.png) | [PDF](figures/boundary_formulation_schematic.pdf) |
| normal_traction_candidate_comparison | [PNG](figures/normal_traction_candidate_comparison.png) | [PDF](figures/normal_traction_candidate_comparison.pdf) |
| equivalent_boundary_force_profile | [PNG](figures/equivalent_boundary_force_profile.png) | [PDF](figures/equivalent_boundary_force_profile.pdf) |
| normal_affine_field_profiles | [PNG](figures/normal_affine_field_profiles.png) | [PDF](figures/normal_affine_field_profiles.pdf) |
| shear_traction_consistency | [PNG](figures/shear_traction_consistency.png) | [PDF](figures/shear_traction_consistency.pdf) |

The schematic uses the actual central particle slice and one visual arrow scale. Force-support quantiles describe this conversion, not a physical collar thickness. The geometric core is not a closed bulk stencil. No fracture or CFD result is shown.
