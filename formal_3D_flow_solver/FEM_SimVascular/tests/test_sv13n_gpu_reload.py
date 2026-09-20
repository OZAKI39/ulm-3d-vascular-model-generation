from sv13n_support import *
def test_final_candidate_integrity():
 d=accepted('gpu_steady_candidate')
 assert d['reload']['status']=='PASS' and d['reload']['velocity_finite'] and d['reload']['pressure_finite']
 assert d['measurement']['wall_noslip_pass']
 policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
 assert max(d['measurement']['epsilon_Q'],d['measurement']['epsilon_mass'])<=policy['mass_limit']
 assert d['reload']['sha256']==d['VTU_sha256']
 assert hashlib.sha256((ROOT/d['VTU']).read_bytes()).hexdigest()==d['VTU_sha256']
 p=Path(d['checkpoint']['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==d['checkpoint']['sha256']
 assert d['classification']=='GPU_STEADY_CANDIDATE' and d['scientific_equivalence']=='DEFERRED'
 assert d['solver_sha256']==accepted('svmp_gpu_build')['executable_sha256']
 assert d['PETSc_library_sha256']==accepted('petsc_gpu13_build')['library_sha256']
