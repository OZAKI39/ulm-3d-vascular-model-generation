from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission

def test_wall_and_zero_offset(source):
 s,w,c,v=channel(); a=FiniteSizeAdmission(s,source,wall=w,velocity=v); e=source.next_mb(); e.update(particle_id=1,attempt_count=0)
 p,status,_=a.check(e,[0,0,0],{}); assert p is not None and p.position[2]==0
 p,status,_=a.check(e,[99.999e-6,0,0],{}); assert p is None and 'WALL' in status

def test_whole_plug_wall_path(source):
 from particle_3d.injection_admission import PlugPathAdmission,plug_side_wall_margin,event_shape
 from particle_3d.wall_gap import wall_gap
 s,w,c,v=channel(); a=PlugPathAdmission(s,source,width=200e-6,velocity=v); e=source.next_rbc(); e.update(particle_id=1,attempt_count=0)
 p=a.attempt(e,{})
 assert p is not None and plug_side_wall_margin(p.shape(),200e-6)>0
 for z in [0,5e-9,1e-8]:
  assert plug_side_wall_margin(p.shape().moved([p.position[0],p.position[1],z]),200e-6)>0

def test_actual_remaining_time_after_birth(report):
 import json,numpy as np
 d=json.loads((report/'data/13_isolated_mb_transport.json').read_text())
 assert d['status']=='PASS' and d['artificial_inlet_offset_m']==0
 states=d['states']; first=states[0]['particle']; last=states[-1]['particle']
 assert first['position']==d['birth_position_m']==d['candidates'][-1]['position_m']
 assert first['radius']==last['radius']==d['geometry']['radius_m']
 assert np.linalg.norm(np.array(last['position'])-first['position'])>0
 assert d['model']=='UNCHANGED_P65_SPHERE_NORMAL_NEARFIELD'
