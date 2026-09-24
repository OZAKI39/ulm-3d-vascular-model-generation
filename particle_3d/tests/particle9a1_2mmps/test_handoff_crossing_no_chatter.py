from particle_3d.particle65_cases import run_case

def test_real_time_event_remainder():
 d=run_case('wall',.032,.032)
 assert d['summary']['final_time_s']==.032
 assert d['summary']['rejected_trials']==0
 assert any(r.get('handoff_event') for r in d['ledger'])
 assert max(r['depth'] for r in d['ledger'])==0
 assert d['summary']['tangential_displacement_m']>6.39e-8
 assert abs(sum(r['dt_s'] for r in d['ledger'])-.032)<1e-16


def test_finite_edge_interior_crossing_with_clear_endpoints():
 import numpy as np
 from particle_3d.nearfield_handoff import first_handoff_event,wall_handoff_certificate
 from particle_3d.nearfield_regularization import NearFieldRegularizationV1
 from particle_3d.wall_geometry import WallGeometry
 from particle_3d.particle_shapes import Sphere
 # Closest feature is the x-directed finite edge; both endpoints are clear.
 wall=WallGeometry([np.array([[-1e-5,0,0],[1e-5,0,0],[0,1e-5,0]])])
 p=NearFieldRegularizationV1();s=Sphere([0,-2e-6,.9e-6],1e-6);v=np.array([0,4e-3,0]);dt=.001
 # End above the wall face is actually penetrating, use a vertex pass instead.
 s=Sphere([1.1e-5,-2e-6,0],1e-6);v=np.array([0,4e-3,0])
 event=first_handoff_event({1:s},{1:v},dt,wall,p)
 assert event is not None and 0<event['time_from_start_s']<dt
 assert event['full_candidate_endpoint_g_nf_m']>0
 end=s.moved(s.center_m+event['time_from_start_s']*v)
 assert wall_handoff_certificate(s,end,wall,2e-9)[0]
