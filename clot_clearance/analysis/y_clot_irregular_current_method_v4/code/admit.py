from common import *
import datetime
def main():
 geometry=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());elastic=json.loads((O/'ELASTIC_CHECK.json').read_text());native=json.loads((O/'KEY_CHECKS.json').read_text())
 assert geometry['status']=='REFERENCE_GEOMETRY_PASS' and elastic['status']=='PASS' and native['status']=='PASS'
 original=json.loads((V4_O/'V4_IMPLEMENTATION_FREEZE.json').read_text());assert sha(LIB)==original['library_sha256']
 names=['common.py','energy.py','native.py','run_case.py','mapped_helpers.py','mapped_geometry.cpp','supported_surface.py','core_surface.py','current_surface.py','shear_field_library.py','admission_gates.json','prepare_reference.py','prepare_boundary.py','checks.py','check_elastic.py','patch_io.py']
 freeze=dict(sources={name:sha(A/name) for name in names},library_sha256=sha(LIB),reference_sha256=sha(O/'provenance/reference.npz'),material_sha256=sha(O/'BIPHASIC_FIELD.npz'),surface_cache_sha256=sha(O/'provenance/mapped_dual_face_crossings.npz'),boundary_sha256=sha(O/'BOUNDARY_CLASSIFICATION.npz'))
 dump(O/'V4_IMPLEMENTATION_FREEZE.json',freeze)
 dump(O/'SOURCE_FREEZE.json',dict(local_base_commit='ba96738',source_method=str(V4_O),source_Y_case=str(Y_OLD),method_library_sha256=sha(LIB),source_Y_reference_sha256=sha(Y_OLD/'provenance/reference.npz'),source_Y_material_sha256=sha(Y_OLD/'BIPHASIC_FIELD.npz'),source_Y_boundary_sha256=sha(Y_OLD/'BOUNDARY_CLASSIFICATION.npz'),frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),formal_cases_allowed=1,protocol='0-2ms admission then same legal restart to3.6ms; no extra physical cases',changes=['V4 irregular interior reference and actual nonuniform cell volumes on original Y discrete mask','Rebuild fixed reference neighborhoods with original continuous-shape void/attachment exclusion','Preserve per-original-ID material phase and integrated strong volume; preserve unnotched Y boundary and anchors','Reuse unchanged V4 native constitutive/contact binary; Y wall and corrected attachment-face rules retained','Wall damping keeps ratio0.5 using each actual rho*Vi mass'],labels=[HET_LABEL,LABEL],CFD_or_acoustics_updated=False))
 dump(O/'LPS_DISCRETIZATION_CHECK.json',dict(status='PASS',geometry='IRREGULAR_REFERENCE_GEOMETRY.json',elastic='ELASTIC_CHECK.json',native='KEY_CHECKS.json',reuse_V4_native_library=True,formal_cases_allowed=1))
 print('Y_V4_ADMITTED_ONE_TRAJECTORY',flush=True)
if __name__=='__main__':main()
