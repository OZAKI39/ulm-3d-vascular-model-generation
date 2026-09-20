from sv13o_support import *

def test_single_permitted_profile_preserves_stage_rows_and_factor_cost():
 d=accepted('profile_summary');assert d['runs']['excluded_from_winner_selection']
 assert d['runs']['runs']==len(d['profiles'])==1
 p=d['profiles']['PROFILE_A'];e=p['events']
 assert p['nested_events_are_not_additive'] and p['GPU_residency']=='NOT_MEASURED'
 assert e['MatAssemblyBegin']['count']==62
 assert len(e['MatAssemblyBegin']['stage_rows'])==3
 assert e['MatAssemblyEnd']['time_s']==pytest.approx(sum(r['time_s'] for r in e['MatAssemblyEnd']['stage_rows']))
 assert e['PCSetUpOnBlocks']['time_s']>e['MatLUFactorNum']['time_s']>e['KSPSolve']['time_s']>e['PCApply']['time_s']
 assert p['stages']['PETSc Solve']['percent_total']>80
