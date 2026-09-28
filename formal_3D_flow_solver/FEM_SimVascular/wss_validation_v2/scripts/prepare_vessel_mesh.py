"""Parameterised isolated remesh export/QC; never overwrites formal mesh."""
from pathlib import Path
import argparse,json,sys,importlib.util,subprocess,os,time
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from case_common import *
from vascular_validation.geometry import triangles_of,surface_stats,enclosed_volume,closest_surface_distances,surface_samples
from vascular_validation.mesh_diagnostics import volume_validity,quality_summary
spec=importlib.util.spec_from_file_location('original_mesh_audit',C/'mesh_generate/scripts/audit_mesh.py');original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)

def prepare(case):
 raw=case/'raw_mesh';g=pv.read(raw/'volume.vtu');assert set(g.cells_dict)=={pv.CellType.TETRA};x=np.array(g.points,dtype=float);t=g.cells_dict[pv.CellType.TETRA].astype(np.int64)
 validity=volume_validity(x,t);b,own=boundary(x,t);lookup={tuple(sorted(row)):i for i,row in enumerate(b)};tags=np.zeros(len(b),np.int32);tree=cKDTree(x);maxdist=0.
 for tag in ROLES:
  face=pv.read(raw/('face_%d.vtp'%tag));dist,ids=tree.query(face.points);maxdist=max(maxdist,float(dist.max()));assert dist.max()<1e-12
  ft=ids[triangles_of(face)]
  if tag!=1:original.one_port(ft)
  for row in ft:
   i=lookup[tuple(sorted(row))];assert tags[i]==0;tags[i]=tag
 assert np.all(tags>0);write_mesh(case,x,t,b,tags,own)
 source=np.load(C/'inputs/fem_reference/exterior_surface.npz');sx=source['points_m'];sb=source['triangles'];st=source['facet_tags'];oldwall=wss.surface(sx,sb[st==1])[0];newwall=wss.surface(x,b[tags==1])[0]
 d1=closest_surface_distances(surface_samples(oldwall),newwall);d2=closest_surface_distances(surface_samples(newwall),oldwall);dd=np.r_[d1,d2];ports={}
 for tag,role in ROLES.items():
  if tag==1:continue
  before=surface_stats(sx,sb[st==tag]);after=surface_stats(x,b[tags==tag]);normaldot=np.dot(before['normal'],after['normal']);assert normaldot>0
  ports[role]=dict(source=before,actual=after,area_relative_change=after['area_m2']/before['area_m2']-1,centroid_shift_m=float(np.linalg.norm(np.array(before['centroid_m'])-after['centroid_m'])),normal_dot=float(normaldot))
 dump(case/'reports/geometry_qc.json',dict(validity=validity,source_surface_sha256=sha(C/'outputs/mesh_and_flow/model/imported.vtp'),wall_bidirectional_distance_m=dict(max=float(dd.max()),p95=float(np.quantile(dd,.95)),rms=float(np.sqrt(np.mean(dd**2)))),source_volume_m3=abs(enclosed_volume(sx,sb)),actual_volume_m3=validity['volume_m3'],ports=ports,face_vertex_max_distance_m=maxdist))
 hcase=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1')
 assert json.loads((V/'data/fine_backend_step5.json').read_text())['backend_equivalence_work_gate']
 policy=vascular_policy(case.name,'vascular_mesh_CFD');inlet=ports['INLET']['actual'];policy.update(A_in_m2=inlet['area_m2'],Umean_m_s=policy['Q_target_m3_s']/inlet['area_m2'],Dh_m=inlet['hydraulic_diameter_m']);policy['Re']=policy['Umean_m_s']*policy['Dh_m']/policy['nu_m2_s']
 dump(case/'policy.json',policy);xml_case(case,policy['dt_s'],steps=500)
 opts=case/'run/PETSC_OPTIONS.txt';opts.write_text(opts.read_text().replace('-mat_type aijcusparse -vec_type cuda','-mat_type aij -vec_type standard'));lock_case(case)
 print(json.dumps({'case':case.name,'nodes':len(x),'tetra':len(t),'ports':{k:v['area_relative_change'] for k,v in ports.items()},'wall_distance_max_m':float(dd.max())},indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--edge-factor',type=float,required=True);a=p.parse_args();case=a.case.resolve();assert case.is_relative_to(V) and not (case/'raw_mesh/volume.vtu').exists();case.mkdir(parents=True,exist_ok=True)
 href=json.loads((C/'configs/mesh_policy.json').read_text())['h_ref_m'];request=dict(source_surface=str(C/'outputs/mesh_and_flow/model/imported.vtp'),output=str(case/'raw_mesh'),global_edge_size_m=href*a.edge_factor,factor_relative_to_formal=a.edge_factor)
 dump(case/'mesh_request.json',request)
 distribution=json.loads((C/'reports/mesh_and_flow/simvascular_distribution.json').read_text());launcher=C/distribution['path']/'simvascular'
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',QT_QPA_PLATFORM='offscreen',OMP_NUM_THREADS='1',XDG_CONFIG_HOME=str(V/'runtime/app_config'),XDG_CACHE_HOME=str(V/'runtime/app_cache'),WSS_MESH_REQUEST=str(case/'mesh_request.json'));env.pop('PYTHONPATH',None)
 cmd=[str(launcher),'-python','--',str(V/'scripts/remesh_vessel_embedded.py')]
 start=time.time()
 with (case/'mesh_generation.log').open('w') as f:r=subprocess.run(cmd,env=env,cwd=case,stdout=f,stderr=subprocess.STDOUT)
 dump(case/'mesh_launch.json',dict(command=cmd,returncode=r.returncode,elapsed_s=time.time()-start,request_sha256=sha(case/'mesh_request.json'),launcher_sha256=sha(launcher)))
 assert r.returncode==0,'Official mesher failed; retain log, no silent fallback'
 prepare(case)
