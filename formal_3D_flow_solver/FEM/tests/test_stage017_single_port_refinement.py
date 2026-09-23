from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance,feedback

def test_two_bad_ports_change_only_largest_one():
 m=measured_good();m['quality']['cap_adjacent_below_0_1']=80;m['quality']['total_below_0_1']=80
 m['quality']['low_quality_nearest_boundary_counts'].update(INLET=50,OUTLET_01=30)
 d=feedback(0,m,volume_acceptance(m,BASE,POLICY,True),POLICY)
 H={p:1. for p in CONTRACT['ports']};new=dict(H);new[d['changed_port']]/=d['H_divisor']
 assert [p for p in H if H[p]!=new[p]]==['inlet']

def test_ties_are_deterministic():
 m=measured_good();m['quality']['cap_adjacent_below_0_1']=80;m['quality']['total_below_0_1']=80
 m['quality']['low_quality_nearest_boundary_counts'].update(INLET=40,OUTLET_01=40)
 assert feedback(0,m,volume_acceptance(m,BASE,POLICY,True),POLICY)['changed_port']=='inlet'


def test_combining_a_finer_tested_port_reuses_three_other_real_cap_geometries():
 from fem3d.adaptive_surface import combine_ports
 r=result();chosen={name:next(t for t in s['trials'] if t['trial']==s['selected_trial']) for name,s in r['initial_surface_searches'].items()}
 caps={name:dict(np.load(OUT/t['path']/'cap_mesh.npz')) for name,t in chosen.items()}
 before,_=combine_ports(source(),CONTRACT,caps,POLICY)
 changed=dict(caps);changed['inlet']=dict(np.load(OUT/'surface_trials/inlet/trial_00/cap_mesh.npz'))
 after,qc=combine_ports(source(),CONTRACT,changed,POLICY)
 assert qc['status']=='PASS'
 def geometric_triangles(data,tag):
  triangles=data['points_m'][data['triangles'][data['facet_tags']==tag]]
  return {tuple(sorted(tuple(point) for point in tri)) for tri in triangles}
 for name,port in CONTRACT['ports'].items():
  same=geometric_triangles(before,port['entity_id'])==geometric_triangles(after,port['entity_id'])
  assert same==(name!='inlet')
