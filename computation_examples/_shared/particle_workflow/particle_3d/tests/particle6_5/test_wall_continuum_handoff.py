import numpy as np

def test_wall_handoff_and_tangent(wall_run):
 d=wall_run;assert d['summary']['handoff_event_count']==1
 assert d['summary']['minimum_h_geom_m']>1.99999999e-9
 assert abs(d['summary']['final_velocity_m_s'][0][2])<1e-18
 assert np.isclose(d['summary']['tangential_displacement_m'],2e-6*.032,rtol=1e-14)
 assert any(s['projection']['handoff_constraints'] for s in d['states'])
 assert all(any(r['coefficient_kg_s']>0 for r in s['projection']['potential_blocks']) for s in d['states'] if s['projection']['handoff_constraints'])
