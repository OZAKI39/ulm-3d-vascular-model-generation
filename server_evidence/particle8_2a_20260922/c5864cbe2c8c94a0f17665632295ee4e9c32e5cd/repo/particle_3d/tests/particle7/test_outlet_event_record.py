import pytest
def test_exit_time(engine):
 engine.step_to(.002)
 assert engine.exits
 for e in engine.exits:
  assert e['outlet_role']=='OUTLET_01'
  assert e['lifetime']==pytest.approx(1e-8/(1e-12/(200e-6)**2),rel=1e-12)
  assert e['geometry_provenance'] and e['trajectory_summary']['exit_position_m'][2]==pytest.approx(1e-8)

def test_certified_mover_early_outlet_time(engine):
 from dataclasses import replace
 import numpy as np
 engine.step_to(.0002)
 assert engine.active
 tag=next(iter(engine.active)); start=engine.time_s; p=engine.active[tag]
 crossing=start+(1e-8-p.position[2])/p.velocity[2]
 def mover(active,t0,t1):
  endpoint=min(t1,min(t0+(1e-8-x.position[2])/x.velocity[2] for x in active.values())) if active else t1
  return {i:replace(x,position=x.position+x.velocity*(endpoint-t0)) for i,x in active.items()},endpoint
 engine.mover=mover; engine.advance(.001)
 e=next(e for e in engine.exits if e['particle_id']==tag)
 assert e['exit_time']==pytest.approx(crossing,abs=1e-18)
