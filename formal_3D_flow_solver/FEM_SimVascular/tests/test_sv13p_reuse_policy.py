from sv13p_support import *
def test_actual_rebuilds_once_per_timestep():
 for name in ('P1_SMOKE','P1_WINDOW'):
  d=accepted(name+'_acceptance');r=d['reuse'];assert r['PC_rebuild_count']==d['steps']<r['KSP_solve_count']
  seen=set()
  for t in r['trace']:
   assert t['reuse']==int(t['step'] in seen);seen.add(t['step'])
def test_cross_timestep_reuse_rejected():
 d=accepted('P1_SMOKE_acceptance');e=read('remote/P1_SMOKE_execution')
 log=(ROOT/'logs/sv1_3p/remote/P1_SMOKE.log').read_text().replace('step=62 equation=0 reuse=0','step=62 equation=0 reuse=1')
 with pytest.raises(GateError):reuse_gate(log,d['profile'],e['history'])
