#!/usr/bin/env python3
"""Actual Gmsh searches on scaled ellipses and noncircular convex polygons."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from fem3d.audit import write_json,sha256,timestamp
from fem3d.adaptive_port import surface_search,surface_quality
from fem3d.adaptive_surface import triangulate_polygon
from fem3d.cap_remesh import coverage_check,signed_area,cross2
p=ROOT/'inputs/stage01_7/acceptance_policy.json';policy=json.loads(p.read_text());lock=json.loads((p.parent/'freeze_lock.json').read_text());assert sha256(p)==lock['policy_sha256']
out=ROOT/'outputs/stage01_7/synthetic';out.mkdir(exist_ok=True,parents=True)
assert not (out/'results.json').exists(),'Synthetic evidence already exists'
results={}
for shape in policy['synthetic']['shapes']:
 theta=np.arange(48)*2*np.pi/48
 if shape=='ellipse':rim=np.c_[np.cos(theta),.65*np.sin(theta)]
 else:
  radius=1+.02*np.cos(5*theta)+.015*np.sin(3*theta);rim=np.c_[radius*np.cos(theta),radius*np.sin(theta)]
 edge=np.roll(rim,-1,axis=0)-rim
 assert np.all(cross2(edge,np.roll(edge,-1,axis=0))>0)
 cases=[]
 for scale in policy['synthetic']['scales']:
  actual=rim*scale;hrim=float(np.median(np.linalg.norm(np.roll(actual,-1,axis=0)-actual,axis=1)));req=np.sqrt(signed_area(actual)/np.pi)
  folder=out/shape/str(scale);folder.mkdir(parents=True)
  def evaluate(H,index,reason):
   xy,tri,meta=triangulate_polygon(actual,hrim,H,policy);coverage=coverage_check(xy,tri,np.arange(len(actual)),policy['geometry']['coverage_relative_roundoff_tolerance'])
   quality=surface_quality(np.c_[xy,np.zeros(len(xy))],tri,policy)
   np.savez_compressed(folder/f'trial_{index:02d}.npz',xy=xy,triangles=tri)
   return {'trial':f'trial_{index:02d}','geometry_status':coverage['status'],'triangle_count':len(tri),'quality':quality,'metadata':meta}
  search=surface_search(evaluate,3*hrim,hrim,req,policy);write_json(folder/'search.json',search)
  assert search['status']=='PASS';chosen=np.load(folder/(search['selected_trial']+'.npz'))
  cases.append({'scale':scale,'search':search,'normalized_H':search['selected_H']/scale,'normalized_points':(chosen['xy']/scale).tolist(),'triangles':chosen['triangles'].tolist()})
 reference=cases[0]
 for case in cases[1:]:
  assert case['search']['selected_trial']==reference['search']['selected_trial']
  assert case['triangles']==reference['triangles']
  np.testing.assert_allclose(case['normalized_points'],reference['normalized_points'],rtol=0,atol=policy['synthetic']['normalized_quality_tolerance'])
  for a,b in zip(case['search']['trials'],reference['search']['trials']):
   assert a['triangle_count']==b['triangle_count']
   np.testing.assert_allclose(list(a['quality']['q_tri'].values()),list(b['quality']['q_tri'].values()),atol=policy['synthetic']['normalized_quality_tolerance'],rtol=0)
 results[shape]={'status':'PASS','cases':cases,'actual_rim_used':True,'circle_fit':False}
write_json(out/'results.json',{'status':'PASS','timestamp':timestamp(),'policy_sha256':sha256(p),'shapes':results})
print(json.dumps({k:[{'scale':v['scale'],'selected_trial':v['search']['selected_trial'],'H':v['search']['selected_H'],'triangles':len(v['triangles'])} for v in d['cases']] for k,d in results.items()},indent=2))
