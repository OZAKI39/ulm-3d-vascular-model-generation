"""Exact geometry replacement must preserve clearance/handoff without projection."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'PARTICLE_CFD_INPUT_CONTRACT.json').read_text())
source=c['source_server_root'] if str(ROOT).startswith('/workspace') else c['source_local_root']
sys.path[:0]=[str(Path(source).parent),str(ROOT/'gmsh_reference')]
from exact_cylinder_particle import (ExactCylinder,exact_gap,exact_certificate,exact_first_handoff,exact_constraints)
from particle_3d.particle_shapes import Sphere
from particle_3d.nearfield_regularization import NearFieldRegularizationV1

def test_exact_cylinder_gap_normal_and_contact():
 wall=ExactCylinder(.0006,.0072);s=Sphere([.0006-1e-6-2e-9,0,.0036],1e-6)
 g=exact_gap(s,wall)
 assert abs(g.gap_m-2e-9)<g.roundoff_m
 np.testing.assert_array_equal(g.normal_inward,[-1.,0.,0.])
 contacts,rows=exact_constraints({0:s},wall,NearFieldRegularizationV1())
 assert len(contacts)==1 and rows[0]['contact_type']=='CONTINUUM_HANDOFF_CONTACT'
 assert np.linalg.norm(g.wall_point_m[:2])==wall.radius_m

def test_path_certificate_checks_outward_endpoint_and_safe_chord():
 wall=ExactCylinder(.0006,.0072);a=1e-6;r=wall.radius_m-a-2e-9
 start=Sphere([r,0,.0036],a);end=Sphere([r*np.cos(.1),r*np.sin(.1),.0037],a)
 assert exact_certificate(start,end,wall,2e-9)[0]
 assert not exact_certificate(start,Sphere([r+1e-9,0,.0036],a),wall,2e-9)[0]

def test_first_event_stops_at_original_handoff_floor():
 wall=ExactCylinder(.0006,.0072);s=Sphere([.0006-1e-6-10e-9,0,.0036],1e-6);v=np.array([1e-5,0.,0.])
 event=exact_first_handoff({0:s},{0:v},.001,wall,NearFieldRegularizationV1())
 assert event and abs(event['time_from_start_s']-.0008)<1e-10
 end=s.moved(s.center_m+event['time_from_start_s']*v)
 assert exact_certificate(s,end,wall,2e-9)[0]
 assert abs(exact_gap(end,wall).gap_m-2e-9)<=exact_gap(end,wall).roundoff_m
