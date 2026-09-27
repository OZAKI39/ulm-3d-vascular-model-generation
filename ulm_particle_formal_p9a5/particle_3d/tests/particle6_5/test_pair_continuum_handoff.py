import numpy as np

def test_pair_handoff_and_symmetry(pair_run):
 d=pair_run;assert d['summary']['handoff_event_count']==1
 for s in d['states']:
  a,b=s['particles'];assert abs(a['center_m'][0]+b['center_m'][0])<1e-20
  assert a['center_m'][1]==b['center_m'][1]
 assert abs(d['summary']['final_velocity_m_s'][0][0]-d['summary']['final_velocity_m_s'][1][0])<1e-18
