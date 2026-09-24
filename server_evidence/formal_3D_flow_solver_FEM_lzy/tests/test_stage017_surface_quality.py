from stage017_helpers import *
import pytest
from fem3d.adaptive_port import surface_quality
from fem3d.cap_remesh import triangle_quality

def test_selected_every_port_quality_recomputed():
 d=dict(np.load(OUT/'selected/surface/tagged_surface_si.npz'))
 for name,p in CONTRACT['ports'].items():
  r=surface_quality(d['points_m'],d['triangles'][d['facet_tags']==p['entity_id']],POLICY)
  assert r['status']=='PASS' and r['low_count']==0 and r['q_tri']['P5']>=.45 and r['q_tri']['median']>=.70

def test_degenerate_and_bad_aspect_rejected():
 for points in (np.array([[0.,0.,0.],[1,0,0],[2,0,0]]),np.array([[0.,0.,0.],[1,0,0],[0,1e-5,0]])):
  assert surface_quality(points,np.array([[0,1,2]]),POLICY)['status']=='FAIL'

def test_surface_triangle_count_alone_does_not_reject():
 n=24;x,y=np.meshgrid(np.arange(n+1),np.arange(n+1));points=np.c_[x.ravel(),y.ravel(),np.zeros(x.size)]
 tri=[]
 for j in range(n):
  for i in range(n):
   a=j*(n+1)+i;tri.extend([[a,a+1,a+n+2],[a,a+n+2,a+n+1]])
 assert len(tri)>800
 assert surface_quality(points,np.array(tri),POLICY)['status']=='PASS'
