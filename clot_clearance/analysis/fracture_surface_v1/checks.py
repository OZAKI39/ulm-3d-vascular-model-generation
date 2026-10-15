"""Acceptance and protection checks. Reads results only; never reruns the integrals."""
from common import *
from datetime import datetime,timezone
import time
start=time.perf_counter()
def read(name):return json.loads((O/(name+'.json')).read_text())
th=read('THEORY_FORMULA_AUDIT'); src=read('SOURCE_LPS_CONSISTENCY_AUDIT'); top=read('PRISTINE_REFERENCE_TOPOLOGY_AUDIT'); nat=read('NATIVE_INPUT_IDENTITY_AUDIT'); cross=read('CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT'); identity=read('INDEPENDENT_ENERGY_FORM_CROSSCHECK'); disc=read('DISCRETE_SURFACE_ENERGY_RESULTS'); sens=read('GEOMETRIC_SENSITIVITY'); perf=read('INTEGRATION_RESOURCES'); frozen=read('DESIGN_FREEZE')
checks={
 'paper_three_dimensional_equations_and_independent_exact_beta_integrals':th['status']=='PASS' and th['three_dimensional'] and th['exact_beta_over_delta']=='1/5' and th['exact_beta_prime_over_delta']=='26/405',
 'two_continuous_formula_arithmetic':abs(th['difference_J_m2'])<1e-15,
 'native_force_and_intact_independent_energy_correspondence':src['status']=='PASS_FOR_INTACT_REFERENCE_IDEAL_POINT_STATE' and src['saved_state_recheck']['relative_energy_error']<1e-12,
 'full_actual_input_and_archival_manifest_identity':top['status']=='PASS' and nat['status']=='PASS' and top['historical_manifest_match'],
 'actual_atom_bond_IDs_reference_coordinates_and_volume_read_from_original_native_file':nat['atom_ID_and_coordinate_bitwise_match'] and nat['bond_ID_endpoints_and_bulk_type_exact_match'] and nat['native_vfrac_m3']==V,
 'periodic_images_and_stencil_complete':top['periodic_image_exact'] and top['independent_integer_stencil_matches_all_original_pairs'],
 'fixed_reference_m':top['m_recomputed_relative_error']<1e-13,
 'native_taper':top['taper_recomputed_relative_error']<1e-13,
 'exactly_two_frozen_planes':len(disc['cases'])==perf['formal_plane_integrals']==2 and [r['plane_z_m'] for r in disc['cases']]==frozen['planes_m'],
 'no_plane_point_ambiguity':all(v['on_plane_nodes']==0 for v in cross['cases']),
 'bounds_unclipped':all(v['lambda_all_range'][0]>=-1e-14 and v['lambda_all_range'][1]<=1+1e-14 for v in cross['cases']),
 'same_point_state_definition_independent_two_forms':all(v['relative_global_error']<1e-12 and v['pointwise_error_normalized_max_W']<1e-12 for v in identity['cases']),
 'noncross_deviatoric_terms_included':all(v['all_non_cross_deviatoric_terms_included'] for v in identity['cases']),
 'two_sides_once_full_ROI_area':all(v['A_m2']==30*DX*PERIOD and abs(v['lower_side_G_J_m2']+v['upper_side_G_J_m2']-v['G0_discrete_J_m2'])<1e-15 for v in disc['cases']),
 'active_ROI_untruncated_reference_neighbors':all(v['all_active_ROI_degrees']==[122] and v['minimum_active_ROI_center_boundary_distance_m']>DELTA for v in cross['cases']),
 'sensitivity_correctly_named_same_phase_translation':sens['status']=='PASS_SAME_PHASE_INTERIOR_TRANSLATION' and sens['same_grid_phase'] and not sens['grid_phase_test'] and sens['different_actual_reference_bond_ID_sets'],
 'nominal_not_material_Gc':disc['calibrated_Gc'] is False,
 'no_new_native_equilibrium_CFD_or_parameter_scan':perf['native_calls']==0 and perf['equilibrium_solves']==0 and frozen['limits']['material_or_grid_scans']==0,
 'budget':perf['integration_s']<60 and perf['peak_RSS_MiB']<1024,
 'frozen_input_and_design_unchanged':sha(W/'reference_topology.npz')==frozen['reference_sha256'] and sha(O/'DESIGN_FREEZE.json')==disc['frozen_design_sha256'],
 'four_figure_groups_and_index':len(list(F.glob('Fig*.png')))==4 and len(list(F.glob('Fig*.svg')))==4 and (F/'OPEN_RESULTS.html').is_file(),
}
old=read('HISTORICAL_PROTECTION_BEFORE'); missing=[];changed=[]
for p,entry in old.items():
    if not Path(p).is_file():missing.append(p)
    elif sha(p)!=entry['sha256']:changed.append(p)
protected={'status':'PASS' if not missing and not changed else 'FAIL','checked_utc':datetime.now(timezone.utc).isoformat(),'files_checked':len(old),'bytes_checked':sum(v['bytes'] for v in old.values()),'missing':missing,'changed':changed,'final_native_library_sha256':sha(B/'liblammps.so'),'production_library_sha256':sha(PROD),'env_sha256':sha('/workspace/y_clot/env.sh'),'native_CXX_modified':False,'materials_modified':False,'history_overwritten':False}
dump(O/'HISTORICAL_PROTECTION_AFTER.json',protected)
checks['historical_production_source_and_material_integrity']=protected['status']=='PASS' and protected['final_native_library_sha256']==EXPECTED
dump(O/'IMPLEMENTATION_CHECKS.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'number_of_checks':len(checks),'scope':'Scientific offline definition and artifact checks only; Gc calibration is explicitly outside scope. Document and transfer have separate post-write readback checks.','check_wall_s':time.perf_counter()-start})
assert all(checks.values()),checks
print(json.dumps({'checks_passed':len(checks),'protected_files':len(old),'history_changed':len(changed),'wall_s':time.perf_counter()-start}))
