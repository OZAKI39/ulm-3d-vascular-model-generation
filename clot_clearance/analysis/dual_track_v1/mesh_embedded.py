# set_walls is a mesher inclusion mask here: all PLC faces must survive for
# volume-only meshing. Physical no-slip is assigned ONLY to tags 1 and 5 in solver.xml.
# Runs only with the existing SimVascular embedded Python 3.5.
import json,os,time,sv,vtk
with open(os.environ['DUAL_TRACK_MESH_REQUEST']) as f:r=json.load(f)
out=r['output'];os.makedirs(out,exist_ok=True)
assert not os.path.exists(os.path.join(out,'volume.vtu'))
print('ENTER load_model',flush=True);m=sv.meshing.TetGen();m.load_model(r['surface']);print('LOADED',m.get_model_face_ids(),flush=True);m.set_walls([1,2,3,4,5]);assert sorted(m.get_model_face_ids())==[1,2,3,4,5]
o=sv.meshing.TetGenOptions(global_edge_size=r['global_edge_m'],surface_mesh_flag=False,volume_mesh_flag=True)
o.use_mmg=False;o.no_bisect=True;o.no_merge=True
o.sphere_refinement=[dict(edge_size=r['local_edge_m'],radius=r['sphere_radius_m'],center=r['sphere_center_m'])];o.sphere_refinement_on=True
with open(os.path.join(out,'OPTIONS.json'),'w') as f:json.dump(o.get_values(),f,indent=2)
t=time.time();m.generate_mesh(o);grid=m.get_mesh();assert grid.GetNumberOfCells()>0;m.write_mesh(os.path.join(out,'volume.vtu'))
for fid in m.get_model_face_ids():
 w=vtk.vtkXMLPolyDataWriter();w.SetFileName(os.path.join(out,'face_'+str(fid)+'.vtp'));w.SetInputData(m.get_face_polydata(fid));assert w.Write()==1
with open(os.path.join(out,'GENERATION.json'),'w') as f:json.dump(dict(points=grid.GetNumberOfPoints(),cells=grid.GetNumberOfCells(),wall_s=time.time()-t),f,indent=2)
