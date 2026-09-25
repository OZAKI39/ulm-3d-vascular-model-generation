from pathlib import Path
import json,hashlib,shutil,subprocess,os
import numpy as np, yaml, vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk,numpy_to_vtkIdTypeArray
S=Path(__file__).resolve().parents[1];cfg=yaml.safe_load((S/'configs/source_paths.yaml').read_text());V=Path(cfg['vascular_project']);O=V/'outputs';A=Path(cfg['aligned_run']);C=Path(cfg['cfd_surface_run']);U=Path(cfg['ultraliser_run']);B=Path(cfg['stable_rigid_stage']);W=Path(cfg['wall_geometry_donor_stage']);D=Path(cfg['sonovue_root']);inputs={}
def reg(p,role):
 p=Path(p);b=p.read_bytes();h=hashlib.sha256(b).hexdigest();inputs[str(p)]={'sha256':h,'bytes':len(b),'role':role};return h
def j(p,role):reg(p,role);return json.loads(p.read_text())
def save(name,obj):(S/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def loadmesh(p):
 reg(p,'tagged or source geometry');r=vtk.vtkXMLPolyDataReader();r.SetFileName(str(p));r.Update();m=r.GetOutput();return m,vtk_to_numpy(m.GetPoints().GetData()),vtk_to_numpy(m.GetPolys().GetData()).reshape(-1,4)[:,1:]
prod=A/'geometry/cfd_surface_axis_aligned_inlet_m.stl';prodsha=reg(prod,'production closed STL');assert prodsha=='840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb';assert reg(W/'provenance/closed_geometry_m.stl','historical production geometry copy')==prodsha
m,ap,af=loadmesh(A/'geometry/cfd_surface_axis_aligned_inlet_um.vtp');cm,cp,cf=loadmesh(C/'geometry/cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp');um,up,uf=loadmesh(U/'geometry/lumen_surface_um.vtp');tf=j(A/'transform/anatomical_to_cfd_transform.json','rigid transform');R=np.array(tf['rotation_matrix_3x3']);t=np.array(tf['forward_homogeneous_transform_4x4_um'])[:3,3];assert np.array_equal(af,cf);err=np.linalg.norm(ap-cp@R.T-t,axis=1).max();assert err<1e-9
raw=prod.read_bytes();n=int(np.frombuffer(raw[80:84],'<u4')[0]);dtype=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attribute','<u2')]);records=np.frombuffer(raw,offset=84,dtype=dtype);tri=records['vertices'].astype(float);assert len(tri)==len(af);quant=float(np.linalg.norm(tri*1e6-ap[af],axis=2).max());assert quant<2e-5
points,back=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);faces=back.reshape(-1,3);cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)/2
names=[m.GetCellData().GetAbstractArray(i).GetName() for i in range(m.GetCellData().GetNumberOfArrays())];tag={name:vtk_to_numpy(m.GetCellData().GetArray(name)) for name in ['CellEntityIds','SurfaceRegionId','RemeshEntityId','boundary_index','boundary_type_code']}
for name,v in tag.items():assert np.array_equal(v,vtk_to_numpy(cm.GetCellData().GetArray(name)))
# Workflow wall includes artificial extension sidewalls and the local surgery collar; caps are ports.
classes=np.full(n,9,int);wall=tag['CellEntityIds']==1;classes[wall&(tag['SurfaceRegionId']==0)&(tag['RemeshEntityId']==1)]=0;classes[wall&(tag['SurfaceRegionId']==1)]=1;classes[wall&(tag['SurfaceRegionId']==0)&(tag['RemeshEntityId']==2)]=2;classes[(tag['boundary_type_code']==1)&~wall]=3;classes[(tag['boundary_type_code']==2)&~wall]=4;assert not np.any(classes==9),'STOP_PORT_CLASSIFICATION'
from scipy.spatial import cKDTree
dist,match=cKDTree(up.astype(np.float32).astype(float)).query(cp);keys={tuple(sorted(f)) for f in uf};bio=classes==0;assert np.max(dist[cf[bio]])==0 and all(tuple(sorted(f)) in keys for f in match[cf[bio]])
meta=j(U/'input/metadata.json','Ultraliser ROI/run metadata');RO=O/'sampling/20260825_133201_radius_plus_structure_k5/roi_library'/f"{meta['roi_id']}.npz";reg(RO,'saved actual ROI');roi=np.load(RO);swc=np.loadtxt(U/'input/roi_core.swc');assert reg(U/'input/roi_core.swc','canonical ROI SWC')==meta['canonical_swc_sha256'];sd={int(r[0]):r for r in swc};mp=meta['swc_node_id_by_local_node_id'];assert all(np.array_equal(sd[int(mp[str(i)])][2:5],p) for i,p in enumerate(roi['local_node_positions_um']))
for p in [A/'input/source_provenance.json',C/'input/frozen_open_geometry_reference.json',C/'qc/far_core_exact_preservation_qc.json',O/'cfd_preprocess/global_to_roi_anchor003274_20260825_183628/input/geometry_reference.json',U/'input/source_swc_stl_model_generate.yaml']:
 reg(p,'upstream lineage evidence');dest=S/'provenance/upstream'/p.relative_to(O);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
# Read requested files without restoring deleted worktree entrypoints.
for rel in ['swc_roi_generate.py','swc_stl_model_generate.py','cfd_surface_prepare.py','configs/swc_roi_generate.yaml','configs/swc_stl_model_generate.yaml','configs/cfd_surface_prepare.yaml','docs/CFD_SURFACE_PREPARE.md','docs/ULTRALISER_PIPELINE.md','docs/CFD_PREPROCESS.md','utils/cfd_surface_prepare/vmtk_qc.py','utils/cfd_surface_prepare/guarded_remesh.py']:
 p=V/rel;dest=S/'provenance/vascular_source_readback'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():reg(p,'vascular source read-only');dest.write_bytes(p.read_bytes())
 else:
  b=subprocess.run(['git','-C',str(V),'show','HEAD:'+rel],check=True,capture_output=True).stdout;dest.write_bytes(b)
# The baseline is archived as donor source, not mislabelled as an implemented workflow executable.
for p in list((B/'src').glob('*'))+[B/'CMakeLists.txt',B/'provenance/PROJECT_BUILD_PROVENANCE.json']:
 if p.is_file():reg(p,'stable rigid baseline donor');dest=S/'provenance/stable_rigid_source'/p.name;dest.parent.mkdir(exist_ok=True);shutil.copyfile(p,dest)
for name in ['wall_distance.cpp','wall_distance.hpp']:
 p=W/'src'/name;reg(p,'pure geometry donor audit');dest=S/'provenance/wall_geometry_donor'/name;dest.parent.mkdir(exist_ok=True);shutil.copyfile(p,dest)
for name in ['frozen_flow.cpp','frozen_flow.hpp']:shutil.copyfile(B/'src'/name,S/'src'/name)
field=W/'fields/FROZEN_FLOW_FIELD_V0.h5';h=reg(field,'frozen engineering transient flow');assert h=='7828ed32b0cdd1e0e85836a564b6737128e7dbe30181ddd43110233b37b64a7f';fc=j(W/'provenance/FROZEN_FLOW_FIELD_CONTRACT.json','field geometry/units contract');assert fc['geometry_sha256']==prodsha;shutil.copyfile(field,S/'fields/FROZEN_FLOW_FIELD_V0.h5');save('contracts/FROZEN_FLOW_CONTRACT.json',fc)
for p in [D/'src/sonovue_sampler.py',D/'input/FROZEN_SONOVUE_HISTOGRAM.csv',D/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json']:
 reg(p,'frozen SonoVue distribution');shutil.copyfile(p,S/'inputs'/p.name)
def poly(ids):
 mesh=vtk.vtkPolyData();v=vtk.vtkPoints();v.SetData(numpy_to_vtk(points,deep=True));mesh.SetPoints(v);cells=vtk.vtkCellArray();cells.SetCells(len(ids),numpy_to_vtkIdTypeArray(np.c_[np.full(len(ids),3),faces[ids]].astype(np.int64).ravel(),deep=True));mesh.SetPolys(cells)
 for name,v in {'production_face_id':ids,'region_class':classes[ids],**{key:val[ids] for key,val in tag.items()}}.items():a=numpy_to_vtk(v,deep=True);a.SetName(name);mesh.GetCellData().AddArray(a)
 return mesh
def export(name,ids):
 ids=np.array(ids,int);mesh=poly(ids);writer=vtk.vtkXMLPolyDataWriter();writer.SetDataModeToBinary();writer.SetFileName(str(S/'geometry'/f'{name}.vtp'));writer.SetInputData(mesh);assert writer.Write()==1
 # Exact subset of original binary STL records; no geometry round-trip quantization.
 with (S/'geometry'/f'{name}.stl').open('wb') as f:f.write(raw[:80]);f.write(np.array([len(ids)],'<u4').tobytes());f.write(records[ids].tobytes())
 return mesh
export('WALL_ONLY',np.flatnonzero(wall));export('CLOSED_REFERENCE',np.arange(n));shutil.copyfile(prod,S/'geometry/CLOSED_REFERENCE.stl')
ports=[];outlets=sorted(np.unique(tag['boundary_index'][classes==4]))
for typ in [3,4]:
 for b in sorted(np.unique(tag['boundary_index'][classes==typ])):
  ids=np.flatnonzero((classes==typ)&(tag['boundary_index']==b));name='INLET' if typ==3 else f'OUTLET_{outlets.index(b)}';export(name,ids);normal=cross[ids].sum(axis=0);normal/=np.linalg.norm(normal);center=(tri[ids].mean(axis=1)*area[ids,None]).sum(axis=0)/area[ids].sum();ports.append({'name':name,'role':'INLET_PORT' if typ==3 else 'OUTLET_PORT','workflow_outlet_id':None if typ==3 else outlets.index(b),'source_boundary_index':int(b),'source_cell_entity_id':int(tag['CellEntityIds'][ids[0]]),'triangle_count':len(ids),'area_m2':float(area[ids].sum()),'center_m':center.tolist(),'outward_unit_normal':normal.tolist(),'surface':f'geometry/{name}.vtp','stl':f'geometry/{name}.stl'})
labels={0:'BIOLOGICAL_CORE_WALL',1:'ARTIFICIAL_EXTENSION_WALL',2:'ARTIFICIAL_BOUNDARY_SURGERY_WALL',3:'INLET_PORT',4:'OUTLET_PORT',9:'UNKNOWN'}
manifest={'status':'PASS','units':'m','production_stl_sha256':prodsha,'actual_vtp_arrays':names,'wall_cell_selection':'CellEntityIds==1; no inlet or outlet cap included','biological_event_selection':'region_class==0 only','ports':ports,'regions':[{'code':k,'label':v,'triangle_count':int(np.sum(classes==k)),'area_m2':float(area[classes==k].sum())} for k,v in labels.items()],'no_geometry_modification':'STL assets are exact record subsets; VTP uses same coordinates; only identical coordinate indices unified for adjacency.'};save('geometry/BOUNDARY_MANIFEST.json',manifest)
np.savez_compressed(S/'geometry/GEOMETRY_ARRAYS.npz',points=points,faces=faces,classes=classes,area=area)
provenance={'status':'PASS','production_geometry_sha256':prodsha,'ROI':meta['roi_id'],'Ultraliser_run':U.name,'CFD_preprocess_run':'global_to_roi_anchor003274_20260825_183628','CFD_surface_prepare_run':C.name,'axis_alignment_run':A.name,'rigid_transform_vertex_error_um':float(err),'stl_vtp_max_quantization_um':quant,'exact_original_Ultraliser_biological_faces':int(bio.sum()),'source_cell_tags_preserved':True,'frozen_field_geometry_contract_matches':True,'default_anchor_not_used_as_identity_evidence':True,'input_records':'provenance/INPUT_HASHES_BEFORE.json'};save('GEOMETRY_PROVENANCE.json',provenance)
config={'description':'Frozen-flow vascular microbubble transport with bubble-bubble hydrodynamics and geometric hard-wall exclusion. Preflight; runtime not yet implemented.','paths':{'frozen_flow':'fields/FROZEN_FLOW_FIELD_V0.h5','wall_surface':'geometry/WALL_ONLY.stl','closed_reference':'geometry/CLOSED_REFERENCE.stl','inlet_surface':'geometry/INLET.vtp','outlet_surfaces':[p['surface'] for p in ports if p['role']=='OUTLET_PORT']},'bubbles':{'counts':[1,8,16],'random_seed':42,'sonovue_distribution':'inputs/FROZEN_SONOVUE_HISTOGRAM.csv','input_quantity':'diameter_um; convert to radius_m by diameter*0.5e-6','resample_sizes':False},'injection':{'mode':'BATCH_AT_T0_SINGLE_PLANE','continuous_injection':False,'position_sampling':'uniform by inlet cap triangle area, translated inward, followed by exact finite-radius and fluid rejection','offset_rule':'max(4*flow_dx_m, maximum_radius_in_case_m + initial_clearance_m)','initial_clearance_rule':'one frozen-flow voxel, flow_dx_m','maximum_position_attempts_per_bubble':5000},'simulation':{'maximum_time_s':.02,'maximum_steps':200000,'dt_max_s':5.791749420578671e-5,'C_adv':.25,'C_gap':.4,'minimum_dt_s':1e-12,'maximum_retry_count':20,'integrator':'inherit stable rigid RK2, only after initialization preflight passes'},'near_wall':{'thresholds':[1.,.5,.2,.1],'event_type':'GEOMETRIC_ONLY','biological_region_code':0},'output':{'trajectory_stride_steps':1,'visualization_stride_steps':1},'physics':{'R_bulk':True,'R_bubble_bubble_excess':True,'R_wall':False,'ambient_wall_shear':False,'adhesion':False,'RBC':False,'buoyancy':False,'lift':False,'particle_inertia':False,'TWIST':'PENDING'},'FLOW_PHYSICS_STATUS':'ENGINEERING_TRANSIENT_FIELD_ONLY'}
(S/'configs/workflow_v0.yaml').write_text(yaml.safe_dump(config,sort_keys=False));save('provenance/INPUT_HASHES_BEFORE.json',inputs);print(json.dumps({'geometry':provenance,'ports':ports},indent=2))
