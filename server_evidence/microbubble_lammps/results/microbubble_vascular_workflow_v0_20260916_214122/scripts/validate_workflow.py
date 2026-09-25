"""Independent preflight validator; runtime checks are NOT_RUN after mandatory STOP.
This file must never promote a preflight-only stage to WORKFLOW PASS.
"""
from pathlib import Path
import json,csv,hashlib,subprocess,os
import numpy as np,vtk
from scipy.spatial import ConvexHull
from vtk.util.numpy_support import vtk_to_numpy
S=Path(__file__).resolve().parents[1];G=dict(np.load(S/'geometry/GEOMETRY_ARRAYS.npz'));m=json.loads((S/'geometry/BOUNDARY_MANIFEST.json').read_text());r=json.loads((S/'validation/INITIALIZATION_PREFLIGHT.json').read_text());checks={}
def read(p):
 v=vtk.vtkXMLPolyDataReader();v.SetFileName(str(p));v.Update();a=v.GetOutput();assert a.GetNumberOfCells()>0;return a
wall=read(S/'geometry/WALL_ONLY.vtp');wallids=vtk_to_numpy(wall.GetCellData().GetArray('production_face_id'));expected=np.flatnonzero(np.isin(G['classes'],[0,1,2]));assert np.array_equal(wallids,expected);allids=list(wallids)
for port in m['ports']:
 mesh=read(S/port['surface']);ids=vtk_to_numpy(mesh.GetCellData().GetArray('production_face_id'));assert not np.intersect1d(ids,wallids).size;allids.extend(ids)
assert sorted(allids)==list(range(len(G['faces'])));checks['wall_and_ports_disjoint_complete_partition']=True;checks['wall_triangle_count']=len(wallids);checks['inlet_count']=1;checks['outlet_count']=3
# Verify every derived STL record against the source reference, preserving all vertex bits.
raw=(S/'geometry/CLOSED_REFERENCE.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('x','<f4',(3,3)),('att','<u2')]);source=np.frombuffer(raw,offset=84,dtype=dt)
for name,ids in [('WALL_ONLY',wallids)]+[(p['name'],vtk_to_numpy(read(S/p['surface']).GetCellData().GetArray('production_face_id'))) for p in m['ports']]:
 b=(S/'geometry'/f'{name}.stl').read_bytes();assert b[84:]==source[ids].tobytes()
checks['derived_stl_vertices_bitwise_original_subsets']=True
proof=[]
for c in r['cases']:
 p=S/'cases'/c['case'];data=np.load(p/'INJECTION_SECTION.npz');xyz=data['points'];normal=data['normal'];area3=.5*abs(np.cross(xyz,np.roll(xyz,-1,axis=0)).sum(axis=0)@normal);assert abs(area3-c['section_area_m2'])<1e-19
 xy=(xyz-data['plane_center'])@np.c_[data['axis1'],data['axis2']];hull=ConvexHull(xy);generous_area=float(hull.volume);rows=list(csv.DictReader((p/'POPULATION_DRAWN.csv').open()));radii=np.array([float(row['radius_m']) for row in rows]);draws=np.array([float(row['sample_uniform_0_1']) for row in rows]);assert np.array_equal(draws,np.random.default_rng(c['seed']).random(c['count']));assert all(abs(float(row['diameter_m'])-2*float(row['radius_m']))<1e-20 for row in rows)
 if c['case']!='W1':assert np.pi*np.sum(radii*radii)>generous_area
 proof.append({'case':c['case'],'cross_section_area_3d_m2':float(area3),'convex_hull_area_generous_upper_bound_m2':generous_area,'demand_m2':float(np.pi*np.sum(radii*radii)),'necessary_area_condition_even_with_convex_hull':bool(np.pi*np.sum(radii*radii)<=generous_area)})
checks['independent_3d_area_and_convex_hull_proof']=proof
loc=vtk.vtkStaticCellLocator();loc.SetDataSet(wall);loc.BuildLocator();closed=read(S/'geometry/CLOSED_REFERENCE.vtp');signed=vtk.vtkImplicitPolyDataDistance();signed.SetInput(closed)
row=next(csv.DictReader((S/'cases/W1/BUBBLES_INITIAL.csv').open()));x=[float(row['initial_'+a]) for a in 'xyz'];hit=[0.,0.,0.];cid=vtk.reference(0);sid=vtk.reference(0);d2=vtk.reference(0.);loc.FindClosestPoint(x,hit,cid,sid,d2);gap=np.sqrt(float(d2))-float(row['radius_m']);assert signed.EvaluateFunction(x)<0;assert gap>=r['cases'][0]['initial_clearance_m'];assert abs(gap-float(row['initial_gap_m']))<1e-15;checks['W1_independent_sphere_inside_and_wall_gap_m']=float(gap)
# Verify repeatability with exact byte comparison; captures finite position rejection deterministically.
files=[f'cases/{name}/POPULATION_DRAWN.csv' for name in ['W1','W2','W3']]+['cases/W1/BUBBLES_INITIAL.csv'];before={p:hashlib.sha256((S/p).read_bytes()).hexdigest() for p in files};result=subprocess.run(['/usr/bin/python3','-B',str(S/'scripts/initialization_preflight.py')],capture_output=True,text=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'));(S/'validation/REPEAT_INITIALIZATION.log').write_text(result.stdout+result.stderr);assert result.returncode==0
assert all(hashlib.sha256((S/p).read_bytes()).hexdigest()==h for p,h in before.items());checks['deterministic_inputs_byte_identical_on_repeat']=True
# Native dynamics were deliberately not started after user stop condition.
assert not list((S/'raw').glob('*TRAJECT*'));checks['no_fabricated_trajectories']=True
inputs=json.loads((S/'provenance/INPUT_HASHES_BEFORE.json').read_text());changed=[p for p,h in inputs.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h['sha256']];assert not changed;checks['local_recorded_input_hashes_unchanged']=len(inputs)
remote=json.loads((S/'provenance/REMOTE_PREFLIGHT.json').read_text());matches=[]
for rp,info in remote['inputs'].items():
 candidates=[(p,v) for p,v in inputs.items() if p.split('/')[-1]==rp.split('/')[-1] and ('stable_rigid_sphere_near_field_v0' in rp)==('stable_rigid_sphere_near_field_v0' in p) and ('/src/' in rp)==('/src/' in p)]
 if candidates:assert any(v['sha256']==info['sha256'] for p,v in candidates);matches.append(rp)
checks['remote_local_frozen_input_matches']=len(matches)
beforegit=json.loads((S/'provenance/GIT_BEFORE.json').read_text());after={};env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
for name,record in beforegit.items():
 after[name]={'root':record['root']}
 for key,args in [('branch',['branch','--show-current']),('commit',['rev-parse','HEAD']),('status',['status','--short'])]:after[name][key]=subprocess.run(['git','-C',record['root'],*args],capture_output=True,text=True,check=True,env=env).stdout
 assert after[name]==record
(S/'provenance/GIT_AFTER.json').write_text(json.dumps(after,indent=2)+'\n');checks['both_git_repositories_unchanged']=True
for p in (S/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
result={'PREFLIGHT_VALIDATION_STATUS':'PASS','WORKFLOW_TECHNICAL_STATUS':'BLOCKED','STOP_REASON':'STOP_INLET_BATCH_INITIALIZATION_INFEASIBLE','checks':checks,'runtime_validation_not_run':['native FlowFieldSampler execution','G1-G5 C++ geometry/swept crossing suite','RK2 trajectory continuity','dynamic bubble non-overlap','swept wall exclusion','outlet crossing lifecycle','near-wall residence recomputation','MPI1/MPI4 trajectory consistency'],'no_runtime_PASS_inferred_from_preflight':True};(S/'validation/PREFLIGHT_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
