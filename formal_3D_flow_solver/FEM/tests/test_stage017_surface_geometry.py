from stage017_helpers import *
import pytest
from fem3d.cap_remesh import check_rim,rim_loop,coverage_check,check_wall
from fem3d.planar_port import verify_plane,verify_metrics

def test_selected_closed_geometry_and_wall_are_unchanged():
 r=selected_recomputed();d=dict(np.load(OUT/'selected/surface/tagged_surface_si.npz'))
 assert check_wall(source(),d)['maximum_wall_displacement_m']==0
 assert r['boundary_fidelity']['maximum_boundary_displacement_m']==0
 for p in r['ports'].values():
  assert p['rim']['maximum_rim_displacement_m']==0
  assert p['relative_projected_area_error']<=1e-12 and p['relative_vector_area_error']<=1e-12

@pytest.mark.parametrize('attack',['rim_moved','rim_split','normal','region','hole','overlap'])
def test_malicious_geometry_rejected(attack):
 p,t=square();_,edges=rim_loop(t)
 with pytest.raises(ValueError):
  if attack=='rim_moved':
   changed=p.copy();changed[0,2]+=1e-15;check_rim(p,edges,changed,t)
  elif attack=='rim_split':
   changed=np.vstack([p,[.5,0,0]]);tri=np.vstack([[[0,5,4],[5,1,4]],t[1:]]);check_rim(p,edges,changed,tri)
  elif attack=='normal':
   port=CONTRACT['ports']['inlet'];verify_plane(port,port['plane_origin_m'],-np.array(port['outward_normal']),port['basis'])
  elif attack=='region':
   port=CONTRACT['ports']['inlet'];other=deepcopy(port);other['formal_projected_area_m2']*=1.01;verify_metrics(port,other,POLICY['geometry'])
  elif attack=='hole':coverage_check(p[:,:2],t[1:],np.arange(4))
  else:coverage_check(p[:,:2],np.vstack([t,t[:1]]),np.arange(4))
