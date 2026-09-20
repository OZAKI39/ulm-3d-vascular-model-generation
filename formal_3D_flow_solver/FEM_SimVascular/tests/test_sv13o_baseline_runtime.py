from sv13o_support import *
def test_observed_stage_N_is_not_a_median():
 d=read('baseline_runtime');n=json.loads((ROOT/'reports/sv1_3n/gpu_steady_candidate.json').read_text())
 assert d['classification']=='DEVELOPMENT_BASELINE_OBSERVED' and d['not_a_repeated_benchmark']
 assert d['wall_time_s']==n['wall_time_s'] and d['total_iterations']==n['total_KSP_iterations']
def test_A_is_actual_low_IO_window():
 d=accepted('PERF_A_acceptance');assert d['steps']==10 and d['VTU_steps']==[70]
 assert d['wall_time_s']>0 and not d['profile_excluded_from_selection']
 for s in d['runtime_semantics']:assert s['PC']=='asm' and s['fill_level']==2 and s['restart']==100

