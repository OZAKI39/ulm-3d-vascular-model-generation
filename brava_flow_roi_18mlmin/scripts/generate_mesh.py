# Executed with official SimVascular embedded Python 3.5.
import os,json,time,sv,vtk
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out=os.path.join(ROOT,'mesh','generation')
os.makedirs(out,exist_ok=True)
assert not os.path.exists(os.path.join(out,'volume.vtu'))
mesher=sv.meshing.TetGen();mesher.load_model(os.path.join(ROOT,'inputs','labelled_core_m.vtp'))
mesher.set_walls([1]);assert sorted(mesher.get_model_face_ids())==[1,2,3,4,5]
options=sv.meshing.TetGenOptions(global_edge_size=0.0002,surface_mesh_flag=True,volume_mesh_flag=True)
options.use_mmg=False
start=time.time();mesher.generate_mesh(options)
mesh=mesher.get_mesh();assert mesh.GetNumberOfCells()>0
mesher.write_mesh(os.path.join(out,'volume.vtu'))
def write_surface(data,name):
 w=vtk.vtkXMLPolyDataWriter();w.SetFileName(os.path.join(out,name));w.SetInputData(data);assert w.Write()==1
write_surface(mesher.get_surface(),'surface.vtp')
for fid in mesher.get_model_face_ids():write_surface(mesher.get_face_polydata(fid),'face_'+str(fid)+'.vtp')
r=dict(status='GENERATED_PENDING_AUDIT',global_edge_size_m=.0002,points=mesh.GetNumberOfPoints(),cells=mesh.GetNumberOfCells(),elapsed_s=time.time()-start,generator='SimVascular 2023-05-31 sv.meshing.TetGen',boundary_layers=False)
with open(os.path.join(out,'generation.json'),'w') as f:json.dump(r,f,indent=2)
print('MESH_GENERATION='+json.dumps(r))
