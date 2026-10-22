from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pyvista as pv
import trimesh
ROOT=Path(__file__).resolve().parents[1]
C=Path('/home/lzy/projects/vascular_printing/outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi')
P=C/'print_fixture_design_all_ports_aligned'
source=P/'core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl'
a=json.loads((P/'all_port_attachment_summary.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest()=='b9fb2d8a57b58b877f940f6c2486f923d93ac71b42793593f4092180dcb7e35e'
t=trimesh.load_mesh(source,process=True)
assert t.is_watertight and t.is_winding_consistent and len(t.split())==1
# Source in mm. Translate to a centered canonical simulation frame, then scale to SI.
origin=t.bounds.mean(axis=0); points=(t.vertices-origin)*1e-3
tri=np.asarray(t.faces); tags=np.ones(len(tri),np.int32)
ports={}; role_ids={'I1':4,'O1':3,'O2':5,'O3':2}
for name,tag in role_ids.items():
 v=a['per_port'][name]; center=np.asarray(v['outside_end_mm']);wall=np.asarray(v['wall_target_mm'])
 normal=(center-wall)/np.linalg.norm(center-wall)
 plane=np.max(np.abs((t.vertices[tri]-center)@normal),axis=1)
 near=np.max(np.linalg.norm(t.vertices[tri]-center,axis=2),axis=1)
 mask=(plane<2e-5)&(near<1.2*v['radius_mm'])
 assert mask.sum()>3,(name,plane.min(),mask.sum())
 assert (tags[mask]==1).all()
 tags[mask]=tag
 area=t.area_faces[mask].sum()*1e-6
 normals=t.face_normals[mask]
 assert np.min(normals@normal)>.999
 ports[name]=dict(face_id=tag,role={'I1':'INLET','O1':'OUTLET_01','O2':'OUTLET_02','O3':'OUTLET_03'}[name],
  center_m=((center-origin)*1e-3).tolist(),center_source_mm=center.tolist(),outward_normal=normal.tolist(),radius_m=v['radius_mm']*1e-3,
  area_m2=area,triangles=int(mask.sum()),mean_velocity_m_s=(3e-7 if name=='I1' else 1e-7)/area,
  Reynolds_number=(3e-7 if name=='I1' else 1e-7)/area*2*v['radius_mm']*1e-3/(.00345312/1056))
surf=pv.PolyData(points,np.c_[np.full(len(tri),3),tri].ravel());surf.cell_data['ModelFaceID']=tags
surf.save(ROOT/'inputs/labelled_core_m.vtp')
np.savez_compressed(ROOT/'inputs/source_surface.npz',points_m=points,triangles=tri,facet_tags=tags)
shutil.copy2(source,ROOT/'inputs/source_core_mm.stl')
transforms={}
for n,d in [('0','candidate_00_upright'),('15','candidate_15_side')]:
 src=C/'final_abs_casting_mold/candidate_0_and_15'/d
 data=json.loads((src/'print_transform.json').read_text());T=np.asarray(data['transform_4x4']);R=T[:3,:3]
 assert np.allclose(R@R.T,np.eye(3),atol=1e-12) and np.isclose(np.linalg.det(R),1)
 transforms[n]={'source_mm_to_candidate_mm':T.tolist(),'canonical_origin_source_mm':origin.tolist(),'candidate_id':int(n)}
 for f in ['manifest.json','print_transform.json','casting_restore_transform.json']:
  shutil.copy2(src/f,ROOT/'inputs'/('candidate_'+n+'_'+f))
for f in ['print_transform.json','compact_manifest.json']:
 shutil.copy2(C/f,ROOT/'inputs'/f)
for f in ['compensated.swc','original_radius.swc']:
 shutil.copy2(C/'candidates/BALANCED'/f,ROOT/'inputs'/f)
shutil.copy2(C/'candidates/BALANCED/VascularMD/compensated_vmd_smooth.swc',ROOT/'inputs/fitted_centerline.swc')
# Keep original source summary: port identity derives from manufacturing provenance, never boundary tag sorting.
shutil.copy2(P/'all_port_attachment_summary.json',ROOT/'inputs/port_attachment_summary.json')
report=dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),geometry='four_port_core_including_engineered_tubes_no_box',
 units={'source':'mm','simulation':'m'},origin_source_mm=origin.tolist(),points=len(points),triangles=len(tri),closed=t.is_watertight,
 volume_m3=t.volume*1e-9,ports=ports,pose_transforms=transforms,material={'density_kg_m3':1056,'dynamic_viscosity_Pa_s':.00345312},
 Q_in_m3_s=3e-7,Q_in_uL_min=18000,outlet_target_m3_s=1e-7,gravity=False,
 pose_equivalence='Rigid stationary Newtonian model without gravity: one physical solution, rotated fields for two print poses. Printed box excluded.')
(ROOT/'inputs/geometry.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
