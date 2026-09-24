import numpy as np
from particle_3d.nearfield_handoff import first_handoff_event,handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle3_cases import plane_triangle

def test_no_early_activation():
 wall=WallGeometry([plane_triangle()]);s={1:Sphere([0,0,1.01e-6],1e-6)};p=NearFieldRegularizationV1()
 assert not handoff_constraints(s,wall,p)[0]
 assert first_handoff_event(s,{1:np.array([0,0,-1e-6])},1e-5,wall,p) is None


def test_safe_glancing_path_composes_original_certificates():
 import json
 from pathlib import Path
 from particle_3d.nearfield_handoff import wall_handoff_certificate,partitioned_wall_handoff_certificate
 r=json.loads((Path(__file__).resolve().parents[2]/'reports/particle9a1_2mmps/reference/id12_safe_glancing_path.json').read_text())
 w=WallGeometry([r['triangle']]);s=Sphere(r['position'],r['radius']);end=s.moved(s.center_m+r['dt']*np.array(r['velocity']))
 assert not wall_handoff_certificate(s,end,w,2e-9)[0]
 assert first_handoff_event({1:s},{1:np.array(r['velocity'])},r['dt'],w,NearFieldRegularizationV1()) is None
 ok,proof=partitioned_wall_handoff_certificate(s,end,w,2e-9)
 assert ok and proof['proof_partition_only']
 assert len(proof['subcertificates'])==2
 assert proof['subcertificates'][0]['t0_fraction']==0
 assert proof['subcertificates'][-1]['t1_fraction']==1
 # A genuinely intersecting path is not made safe by proof partitioning.
 bad=s.moved(s.center_m+np.array([0,0,-1e-6]))
 assert not partitioned_wall_handoff_certificate(s,bad,w,2e-9)[0]
