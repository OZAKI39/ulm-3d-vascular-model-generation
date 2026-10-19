from common import *
def main():
 geometry=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());elastic=json.loads((O/'ELASTIC_PATCH_CHECK.json').read_text());native=json.loads((O/'NATIVE_STATE_SURFACE_CHECK.json').read_text());g=json.loads((O/'ADMISSION_GATES.json').read_text())
 assert geometry['status']=='REFERENCE_GEOMETRY_PASS' and elastic['pass'] and native['status']=='PASS'
 assert native['full_surface_native']['max_relative_R']<g['max_relative_energy_residual'] and native['library_sha256']==sha(LIB)
 dump(O/'LPS_DISCRETIZATION_CHECK.json',dict(status='PASS',geometry='IRREGULAR_REFERENCE_GEOMETRY.json',elastic='ELASTIC_PATCH_CHECK.json',state_and_surface='NATIVE_STATE_SURFACE_CHECK.json',gates='ADMISSION_GATES.json',formal_cases_allowed=1,reference_force_max_N=0.,reference_force_evidence='Native reference smoke check and zero-start full-surface check',elastic_relative_changes=elastic['relative_modulus_changes'],limitations=['Small prescribed affine coupon, not macroscopic fracture calibration','Mapped-control-cell L2 remains a partial surface proxy','Single conditioned random reference; no statistical or microstructural validation']))
 dump(O/'V4_IMPLEMENTATION_FREEZE.json',dict(sources={p.name:sha(p) for p in sorted(A.iterdir()) if p.suffix in ['.py','.cpp','.h','.json','.sh','.conf']},library_sha256=sha(LIB),reference_sha256=sha(O/'provenance/reference.npz'),material_sha256=sha(O/'BIPHASIC_FIELD.npz'),surface_cache_sha256=sha(O/'provenance/mapped_dual_face_crossings.npz')))
 print('V4_ADMITTED_ONE_TRAJECTORY',flush=True)
if __name__=='__main__':main()
