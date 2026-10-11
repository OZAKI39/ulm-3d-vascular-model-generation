from analysis.traction_failure_v1.common import *
assert (B/'GEOMETRY_PREFLIGHT.json').exists()
assert not (B/'MAPPING_DESIGN_FREEZE.json').exists()
save(B/'MAPPING_DESIGN_FREEZE.json',dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),formal_fields=1,new_CFD=0,new_fracture=0,
 selection='Material nodes only, exact distance to exposed triangles <= dx. Face support distance <= delta; nearest-surface normal dot face normal > 0; face lies no more than dx/4 inside local tangent halfspace. No anchors, no dilation/rescaling.',
 node_surface_distance_m=DX,face_support_radius_m=DELTA,normal_dot_min=0.,halfspace_tolerance_m=DX/4,
 initial_kernel='Wendland C2: (1-q)^4(4q+1), q=distance/delta; nonnegative column normalization on actual face centroid. No parameter scan.',
 target_moment='Original degree-2 three-point triangular quadrature, pressure linear and viscous gradient cellwise constant. Centroid moment discrepancy reported separately.',
 correction='One exact six-equality minimum quadratic-norm correction. Per-component D_i=|f_initial_i|^2+(0.01*sum_face|F|*Aeff_i/Atotal)^2. Min sum |df_i|^2/D_i; df=D C^T(C D C^T)^-1 residual. Length-normalized moments. No additional fitted coefficients.',
 gates=dict(all_faces_supported=True,force_and_moment_scaled_error=1e-8,constraint_rank=6,normalized_gram_condition_max=1e8,correction_L2_fraction_max=.10,
 pressure_or_shear_reverse_load_fraction_max=.01,pressure_outward_load_fraction_max=.01,top_one_percent_load_share_change_max=.02,regional_load_share_change_max=.05,single_node_total_force_share_max=.02),
 gate_scope='Engineering admissibility only; pressure/shear/total checked independently. Normal-viscous component conservation and full diagnostics; relative local direction gates exclude this sign-changing tiny component. Node-relative tails reported without hidden epsilon. Regions x quartiles, contact rim and interface; initial-to-corrected regional change gate. >10% is conservatively unresolved, not a material law.',
 geometry_preflight_sha256=sha(B/'GEOMETRY_PREFLIGHT.json')))
print((B/'MAPPING_DESIGN_FREEZE.json').read_text())
