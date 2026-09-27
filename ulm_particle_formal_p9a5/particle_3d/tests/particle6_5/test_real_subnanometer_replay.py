import json,numpy as np
from particle_3d.particle65_cases import real_fixture,records_for,MU
from particle_3d.particle65_motion import Particle65Stepper
from particle_3d.lammps_neighbors import ValidationNeighborPolicy

def test_saved_replay_and_historical(report,real_saved):
 d=json.loads((report/'data/07_historical_p5.json').read_text())
 assert 2.3e-11<d['minimum_h_geom_m']<2.4e-11
 assert d['original_two_mb_initial_rejection']['status']=='INITIAL_STATE_BELOW_CONTINUUM_HANDOFF'
 assert real_saved['summary']['handoff_event_count']>0 and real_saved['summary']['all_accepted_above_lower']
 assert real_saved['initialization']['original_center_unchanged'] and real_saved['initialization']['original_radius_unchanged']

def test_fresh_fem_steps_match_saved(real_saved,repo):
 # This is a regression against historical 0.352841 mm/s evidence, so load
 # the archived input explicitly. Current production never discovers it.
 from particle_3d.particle6_stepper import bind_query_dependency
 from particle_3d.audit import read_frozen
 legacy=repo/'formal_3D_flow_solver/FEM_SimVascular/legacy_inputs/LEGACY_0P352841_MMPS'
 fixture=bind_query_dependency(real_fixture.__wrapped__,{'read_frozen':lambda _:read_frozen(legacy)})
 shapes,provider,wall,classifier,provenance=fixture();query=ValidationNeighborPolicy(20e-6,.5e-6,'PERMANENT_REAL_REPLAY_TEST')
 stepper=Particle65Stepper(records_for(shapes,provider),query,provider,MU,wall=wall,boundary_classifier=classifier)
 dt=real_saved['summary']['dt_s']
 for k in range(1,4):stepper.step_to(k*dt)
 for actual,saved in zip(stepper.accepted_worlds,real_saved['states']):
  assert actual['time_s']==saved['time_s']
  for a,b in zip(actual['particles'],saved['particles']):
   assert set(a)==set(b)
   assert all(np.array_equal(a[k],b[k]) for k in a)
