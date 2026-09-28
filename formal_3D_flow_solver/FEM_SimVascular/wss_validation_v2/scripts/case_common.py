from pathlib import Path
import sys,json,hashlib,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
V=Path(__file__).resolve().parents[1];C=V.parent
sys.path[:0]=[str(C/'solver_support/src'),str(C/'mesh_generate/src')]
from flow_solver_support import wss
from vascular_validation.mesh_diagnostics import tetra_sicn
ROLES={1:'WALL',2:'OUTLET_03',3:'OUTLET_01',4:'INLET',5:'OUTLET_02'}
MU=.00345312;RHO=1056.;DT=1.5040600916882215e-7

def vascular_policy(name,kind):
 source=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/policy.json')
 old=json.loads(source.read_text())
 return dict(case=name,kind=kind,dt_s=old['dt_s'],Q_target_m3_s=old['Q_target_m3_s'],mu_Pa_s=MU,nu_m2_s=old['nu_m2_s'],rho_kg_m3=RHO,A_in_m2=old['A_in_m2'],Umean_m_s=old['Umean_m_s'],Dh_m=old['Dh_m'],Re=old['Re'],minimum_steps=20,steady_change_limit=1e-7,steady_consecutive_intervals=3,mass_limit=1e-6,maximum_wall_time_s=43200,maximum_total_steps=500,save_interval_steps=5,MPI_ranks=8,initial_state='zero; no restart or interpolated solution',linear_algebra_backend='CPU PETSc aij/standard; ASM overlap2, subdomain ILU2',dt_rule='Fixed original H0 dt on every mesh; not recalculated from new cap area or mesh size',baseline_policy_sha256=sha(source),baseline_time_step_derivation_reference=old)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,a):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
def boundary(x,t):
 f=t[:,[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]].reshape(-1,3)
 _,first,count=np.unique(wss.face_keys(f),return_index=True,return_counts=True)
 assert count.max()==2
 b=f[first[count==1]].copy();own=first[count==1]//4
 g=x[b];n=np.cross(g[:,1]-g[:,0],g[:,2]-g[:,0]);direction=g.mean(axis=1)-x[t[own]].mean(axis=1)
 flip=np.einsum('ij,ij->i',n,direction)<0;b[flip]=b[flip][:,[0,2,1]]
 return b,own

def write_mesh(case,x,t,b,tags,own):
 case=Path(case);out=case/'SV_MESH';(out/'mesh-surfaces').mkdir(parents=True,exist_ok=True)
 det=np.linalg.det(x[t[:,1:]]-x[t[:,:1]]);assert np.all(det>0)
 q=tetra_sicn(x,t);assert q.min()>0
 grid=pv.UnstructuredGrid(np.column_stack([np.full(len(t),4),t]).ravel(),np.full(len(t),pv.CellType.TETRA,np.uint8),x)
 grid.point_data['GlobalNodeID']=np.arange(1,len(x)+1,dtype=np.int32);grid.cell_data['GlobalElementID']=np.arange(1,len(t)+1,dtype=np.int32);grid.cell_data['minSICN']=q
 grid.save(out/'mesh-complete.mesh.vtu')
 for tag in [None,*np.unique(tags)]:
  sel=np.ones(len(b),bool) if tag is None else tags==tag
  face,ids=wss.surface(x,b[sel]);del face.point_data['GlobalNodeID_zero_based'];face.point_data['GlobalNodeID']=(ids+1).astype(np.int32)
  face.cell_data['GlobalElementID']=(own[sel]+1).astype(np.int32);face.cell_data['ModelFaceID']=tags[sel]
  path=out/'mesh-complete.exterior.vtp' if tag is None else out/'mesh-surfaces'/(ROLES[int(tag)]+'.vtp')
  face.save(path)
 np.savez_compressed(out/'mesh_arrays.npz',points_m=x,tetra=t,boundary_triangles=b,facet_tags=tags,adjacent_tetra=own,min_sicn=q)
 return grid

def xml_case(case,dt,steps=500,pipe=False):
 case=Path(case);run=case/'run';run.mkdir(parents=True,exist_ok=True)
 tree=ET.parse(C/'wss_audit/inputs/H0/solver.xml');root=tree.getroot();g=root.find('GeneralSimulationParameters')
 for key,val in {'Time_step_size':dt,'Number_of_time_steps':steps,'Increment_in_saving_VTK_files':5,'Increment_in_saving_restart_files':5}.items():g.find(key).text=str(val)
 if pipe:
  mesh=root.find('Add_mesh')
  for e in list(mesh.findall('Add_face')):
   if e.get('name') not in ['WALL','INLET','OUTLET_03']:mesh.remove(e)
  eq=root.find('Add_equation')
  for e in list(eq.findall('Add_BC')):
   if e.get('name') not in ['WALL','INLET','OUTLET_03']:eq.remove(e)
  inlet=eq.find('Add_BC[@name="INLET"]')
  inlet.find('Profile').text='User_defined';inlet.find('Impose_flux').text='false';inlet.find('Value').text='-1.0'
  ET.SubElement(inlet,'Spatial_profile_file_path').text='inlet_profile.txt'
  outlet=eq.find('Add_BC[@name="OUTLET_03"]')
  for e in list(outlet):outlet.remove(e)
  ET.SubElement(outlet,'Type').text='Trac';ET.SubElement(outlet,'Traction_values_file_path').text='outlet_traction.vtp'
 ET.indent(tree,space='  ');tree.write(run/'solver.xml',encoding='utf-8',xml_declaration=True)
 # Same discretization and KSP/PC as current case; omit redundant debug matrix printing only.
 opts='-ksp_type gmres -ksp_pc_side right -ksp_norm_type unpreconditioned -ksp_rtol 1e-10 -ksp_atol 1e-24 -ksp_max_it 2000 -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_monitor_true_residual -ksp_converged_reason -ksp_view -options_view -options_left -use_gpu_aware_mpi 0 -mat_type aijcusparse -vec_type cuda -ksp_gmres_restart 100 -pc_type asm -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type ilu -sub_pc_factor_levels 2 -sv_pc_adaptive_rebuild true -log_view :petsc_profile.txt -log_view_gpu_time'
 (run/'PETSC_OPTIONS.txt').write_text(opts+'\n')
 if not pipe:
  source=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml')
  def flattened(node,path=''):
   key=path+'/'+node.tag+str(sorted(node.attrib.items()))
   result={key:(node.text or '').strip()} if not len(node) else {}
   for child in node:result.update(flattened(child,key))
   return result
  before=flattened(ET.parse(source).getroot());after=flattened(root)
  changes={k:dict(before=before.get(k),after=after.get(k)) for k in set(before)|set(after) if before.get(k)!=after.get(k)}
  allowed=['Time_step_size','Number_of_time_steps','Increment_in_saving_VTK_files','Increment_in_saving_restart_files']
  assert all('GeneralSimulationParameters' in k and any(k.endswith('/'+s+'[]') for s in allowed) for k in changes),changes
  dump(case/'reports/H0_physical_configuration_identity.json',dict(reference_xml=str(source),reference_sha256=sha(source),new_sha256=sha(run/'solver.xml'),only_time_or_output_changes=True,semantic_changes=changes,material_boundary_equation_and_solver_parameters_identical=True))
 return tree

def lock_case(case):
 case=Path(case);paths=[p for root in [case/'SV_MESH',case/'run'] for p in root.rglob('*') if p.is_file()]+[case/'policy.json'];dump(case/'input_hashes.json',{str(p.relative_to(case)):sha(p) for p in paths})
