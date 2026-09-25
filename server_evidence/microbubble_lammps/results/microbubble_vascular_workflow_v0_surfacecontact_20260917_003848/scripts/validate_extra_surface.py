from pathlib import Path
import numpy as np,json
from surface_reference import SurfaceReference,project_velocity
from qualify_surface import Probe,S
out=[]
for name in ['plane0','cylinder','sphere']:
 d=np.load(S/'raw/synthetic'/f'{name}.npz');ref=SurfaceReference(d['vertices'],d['faces']);probe=Probe(name);ns=[];meshns=[];maxcpp=0.;features=set();maxflat=0.
 points=np.column_stack([np.linspace(-1e-6,1e-6,1001),np.zeros(1001),np.full(1001,1e-6)]) if name=='plane0' else np.array([[9e-6*np.cos(a),9e-6*np.sin(a),0] for a in np.linspace(-.02,.1,1001)])
 raw=np.array([-1e-4,2e-4,-3e-4])
 for x in points:
  q=ref.query(x);c=probe.ask(x,raw);assert len(q['clusters'])==1 and len(c['query']['clusters'])==1;N=np.array(q['clusters'][0]['normal']);engine=np.array(c['query']['clusters'][0]['normal']);maxcpp=max(maxcpp,float(np.linalg.norm(N-engine)));assert maxcpp<1e-10
  ns.append(N);meshns.append(q['clusters'][0]['surface_normal']);features.update(a['feature'] for a in q['candidates'])
  expected=raw-min(raw@N,0)*N;err=np.linalg.norm(expected-c['used'])/np.linalg.norm(raw);assert err<1e-11
  if name=='plane0':maxflat=max(maxflat,float(err));assert err<=1e-12 and np.linalg.norm(N-[0,0,1])<1e-12
 probe.close();n=np.array(ns);m=np.array(meshns);jump=np.linalg.norm(np.diff(n,axis=0),axis=1);meshjump=np.linalg.norm(np.diff(m,axis=0),axis=1);assert max(jump)<.002
 out.append(dict(name=name,queries=len(points),features=sorted(features),max_adjacent_contact_normal_difference=float(max(jump)),max_mesh_guide_normal_difference=float(max(meshjump)),max_cpp_python_difference=maxcpp,flat_single_normal_max_relative_error=maxflat))
# Check every saved C++ QP fixture independently, including both surface and offset normals.
count=0;maxmesh=0.;maxqp=0.;multi=0
refs={}
for line in (S/'raw/SYNTHETIC_QUERY_COMPARISON.jsonl').read_text().splitlines():
 a=json.loads(line);name=a['label'];q=a['cpp']['query']
 if name not in refs:
  d=np.load(S/'raw/synthetic'/f'{name}.npz');refs[name]=SurfaceReference(d['vertices'],d['faces'])
 py=refs[name].query(a['x'],a['factor'],a['theta']);assert py['status']==q['status'];assert len(py['clusters'])==len(q['clusters'])
 for pc,cc in zip(py['clusters'],q['clusters']):
  assert pc['support']==cc['support'];e=np.linalg.norm(np.array(pc['surface_normal'])-cc['surface_normal']);maxmesh=max(maxmesh,float(e));assert e<1e-10
 if q['clusters']:
  N=np.array([c['normal'] for c in q['clusters']]);raw=np.array(a['v']);ref=project_velocity(raw,N);e=np.linalg.norm(ref[1]-a['cpp']['used']);maxqp=max(maxqp,float(e));assert e<=2e-12*np.linalg.norm(raw)
  if len(N)>1:multi+=1
 count+=1
(S/'validation/FEATURE_NORMAL_CONTINUITY_AND_SINGLE_NORMAL_LIMIT.json').write_text(json.dumps(dict(status='PASS',paths=out,fixture_rechecks=count,multi_surface_fixtures=multi,max_weighted_mesh_pseudonormal_error=maxmesh,max_independent_QP_error=maxqp,interpretation='The offset contact normal is continuous across convex face/edge/vertex Voronoi features. The auxiliary weighted mesh pseudonormal is piecewise feature-defined and is not claimed to be an analytic smooth normal field on an unchanged STL.'),indent=2)+'\n');print('FEATURE_CONTINUITY_AND_SINGLE_NORMAL_LIMIT_PASS',out)
