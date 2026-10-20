"""Local quadrature limits and traceable historical shear context; no new cases."""
from common import *
from shear_loading import *
def verify_reference_seal():
 old=json.loads((B/'DELIVERY_MANIFEST.json').read_text());index={q['path']:q for q in old['files']};checks={}
 for name in ['reference.npz','reference.native','body.data']:
  p=B/'provenance'/name;wanted=index[str(p.relative_to(R))]['sha256'];actual=sha(O/'provenance'/name);assert sha(p)==actual==wanted
  checks[name]=dict(sha256=actual,prior_sealed_sha256=wanted)
 dump(O/'REFERENCE_SEAL_MATCH.json',dict(status='PASS',prior_manifest_sha256=sha(B/'DELIVERY_MANIFEST.json'),files=checks))
def main():
 verify_reference_seal()
 ref=load_npz(O/'provenance/reference.npz');load=ShearLoad(ref,'stokes');ids=load.loaded_ids;X=ref['X'];areas=load.W@load.area;dual=np.full(len(ids),DX**2)
 for axis in [0,1]:dual[np.isclose(X[ids,axis],X[ids,axis].min(),rtol=0,atol=1e-14)|np.isclose(X[ids,axis],X[ids,axis].max(),rtol=0,atol=1e-14)]*=.5
 force=load.force(HOLD);density=force[ids,0]/areas[ids];ideal=load.amplitude(HOLD)*np.sin(WAVENUMBER*X[ids,1])
 record=dict(status='QUANTIFIED_DISCRETIZATION_LIMIT',area_ratio_to_structured_control_area_min=float((areas[ids]/dual).min()),area_ratio_to_structured_control_area_max=float((areas[ids]/dual).max()),max_force_per_mapped_area_error_over_peak=float(np.max(abs(density-ideal))/load.physical_peak),sign_reversed_relative_to_nodal_mode=int(np.sum(density*ideal<0)),maximum_convex_area_normalized_shear_over_face_peak=float(np.max(abs(density))/load.peak(HOLD)),interpretation='Nonnegative P1 triangle quadrature has nonuniform lumped nodal area versus the uniform particle control area (alternating diagonals). No signed extrapolation, but local force/stress resolution is not mesh-converged. Frozen protocol unchanged; do not claim pointwise continuum force accuracy.')
 np.savez_compressed(O/'SMALL_MAPPING_LOCAL_DIAGNOSTIC.npz',atom_ids=ids+1,mapped_area_m2=areas[ids],structured_control_area_m2=dual,actual_force_per_mapped_area_Pa=density,analytic_wall_shear_at_node_Pa=ideal);dump(O/'SMALL_MAPPING_LOCAL_DIAGNOSTIC.json',record)
 old=BRACKETING/'cases/surface_shear_1000pa_continuation_2/SUMMARY.json';q=json.loads(old.read_text());orig=load_npz(BRACKETING/'provenance/reference.npz');assert all(np.array_equal(ref[k],orig[k]) for k in ['X','bonds','bond_types','r0','initial_cut'])
 dump(O/'HISTORICAL_SHEAR_CONTEXT.json',dict(source=str(old),source_sha256=sha(old),native_driver_sha256=sha(R/'analysis/clot_mechanical_load_requirement_bracketing_v1/native_case_v2.py'),previous_traction_code_sha256=sha(R/'analysis/clot_mechanical_load_requirement_bracketing_v1/prescribed_traction.py'),identical_reference_geometry_and_original_bonds=True,previous_area_m2=q['loaded_area_m2'],previous_protocol='Uniform +x, cosine ramp 0-2ms, hold to2.5ms, no unload',previous_peak_Pa=1000,previous_history_max_free_nodes=q['history_max_free_nodes'],previous_work_J=q['last_row']['W_J'],current_changes=['opposing transverse continuous shear','unsteady fluid-template startup lag','triangle versus quad quadrature','2.4-2.8ms drive-down and decay through3.6ms'],scope='Read-only contextual basis for 1000 Pa selection; different spatial/time/quadrature protocol, not a universal threshold or single-factor comparator.'))
 print('LOCAL_DISCRETIZATION_AUDIT',record)
if __name__=='__main__':main()
