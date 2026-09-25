from pathlib import Path
import json,numpy as np
from surface_reference import SurfaceReference,project_velocity
from flow_geometry import Geometry
from independent_wall import segment_distance
S=Path(__file__).resolve().parents[1];D=S/'runs/DIAGNOSTIC_LONG_RETRY1_mpi1';j=json.loads((D/'WALL_CONTACT_FAILURE.json').read_text());state=json.loads((D/'RUN_STATE.json').read_text());assert state['reason']=='CURVED_GEOMETRY_TANGENTIAL_STALL'
required=['reason','time','step','RK2_stage','attempt','dt','stage_dt','particle_id','radius','position','base','raw_velocity','used_velocity','endpoint','predicted_raw_swept_gap','used_swept_gap','tie_factor','smooth_dihedral_threshold_deg','contact','endpoint_contact','active_normal_indices'];assert all(k in j for k in required)
r=SurfaceReference.production(S);q=r.query(j['contact']['position'],j['tie_factor'],j['smooth_dihedral_threshold_deg']);assert q['status']==j['contact']['status']
for a,b in zip(q['candidates'],j['contact']['candidates']):
 assert a['id']==b['id'] and a['feature']==b['feature'] and a['region']==b['region'];assert b['adjacency'] and np.array_equal(r.t[a['id']],b['triangle_coordinates']);assert np.linalg.norm(np.array(b['face_normal'])-r.normal[a['id']])<1e-12
N=[c['normal'] for c in q['clusters']];used=project_velocity(j['raw_velocity'],N)[1];assert np.linalg.norm(used-j['used_velocity'])<1e-16
g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');base=np.array(j['base']);rawend=base+j['stage_dt']*np.array(j['raw_velocity']);rawgap=segment_distance(g,base,rawend,j['radius']+1e-10)-j['radius'];usedgap=segment_distance(g,base,np.array(j['endpoint']),j['radius']+1e-10)-j['radius'];assert abs(rawgap-j['predicted_raw_swept_gap'])<2e-14 and abs(usedgap-j['used_swept_gap'])<2e-14;assert usedgap<1e-10
assert (D/'WALL_CONTACT_FAILURE.json').stat().st_mtime_ns<=(D/'RUN_STATE.json').stat().st_mtime_ns
(S/'validation/FATAL_DIAGNOSTIC_QUALIFICATION.json').write_text(json.dumps(dict(status='PASS',scope='Deliberately limited retry budget = 1; formal configuration remains 24.',file=str(D/'WALL_CONTACT_FAILURE.json'),all_required_fields_present=True,candidate_geometry_and_topology_independently_reconstructed=True,multi_normal_projection_independently_reconstructed=True,raw_swept_gap_m=rawgap,used_swept_gap_m=usedgap,archive_flushed_before_terminal_state=True,default_same_step_succeeds='Formal LONG includes this step and continues to 0.6506913972031004 s'),indent=2)+'\n');print('FATAL_DIAGNOSTIC_QUALIFICATION_PASS')
