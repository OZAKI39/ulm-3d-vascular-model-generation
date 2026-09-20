from sv13o_support import *
def test_final_candidate_hashes_and_native_history():
 d=accepted('optimized_steady_candidate');assert d['reload']['status']=='PASS' and d['scientific_equivalence']=='DEFERRED'
 assert hashlib.sha256(Path(d['VTU']).read_bytes()).hexdigest()==d['VTU_sha256']
 cp=Path(d['checkpoint']['path']);assert hashlib.sha256(cp.read_bytes()).hexdigest()==d['checkpoint']['sha256']
 dt=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())['dt_s'];assert checkpoint_one_rank(cp,d['stop_step'],dt)['complete_integration_history']
 assert d['measurement']['wall_noslip_pass'] and max(d['measurement']['epsilon_Q'],d['measurement']['epsilon_mass'])<=1e-6
 assert not d['production_changed'] and d['CPU_production']=='CPU_EARLY_STOP_PRODUCTION'

