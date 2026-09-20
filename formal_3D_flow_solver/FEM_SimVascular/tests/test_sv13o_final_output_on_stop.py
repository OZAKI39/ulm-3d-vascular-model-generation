from sv13o_support import *
def test_small_official_case_proves_off_cadence_final_output():
 d=accepted('OFFICIAL_OUTPUT_STOP_acceptance')
 assert d['steps']==11 and d['VTU_steps']==[10,11]
 assert d['checkpoint']['step']==11 and d['reload']['status']=='PASS'
 assert d['checkpoint']['nodes']==d['reload']['points'] and d['normal_exit']
 c=accepted('output_consistency')['cases']['OFFICIAL_OUTPUT_STOP']
 assert c['output_consistency']=='PASS'
 assert c['VTU']['time_s']==pytest.approx(c['checkpoint']['time_s'],rel=1e-12,abs=0.)
def test_full_run_final_VTU_matches_restart_time():
 d=accepted('optimized_steady_candidate');r=accepted('REAL_VASCULAR_GPU_PERF_acceptance')
 assert d['stop_step']==d['checkpoint']['step']==r['VTU_steps'][-1]
 assert d['physical_time_s']==pytest.approx(d['checkpoint']['time_s'],rel=1e-12)
 assert d['reload']['time_s']==pytest.approx(d['checkpoint']['time_s'],rel=1e-12,abs=0.)
 c=accepted('output_consistency')['cases']['REAL_VASCULAR_GPU_PERF']
 assert c['output_consistency']=='PASS' and c['solver_status']=='PASS'
