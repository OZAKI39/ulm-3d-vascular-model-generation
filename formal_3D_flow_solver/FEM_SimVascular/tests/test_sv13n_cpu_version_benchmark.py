from sv13n_support import *
def test_all_cpu_medians_recomputed():
 d=accepted('cpu_version_benchmark');assert set(d['groups'])=={'OLD_CPU1','OLD_CPU4','NEW_CPU1','NEW_CPU4'}
 for g in d['groups'].values():assert timing_gate([read(n+'_acceptance') for n in g['names']])==g['median_s']
