# Official SimVascular embedded Python 3.5 only; no legacy fallback semantics.
import os,json,time,hashlib
import sv,vtk
p=os.environ['WSS_MESH_REQUEST']
with open(p) as f:request=json.load(f)
source=request['source_surface'];out=request['output'];h=float(request['global_edge_size_m'])
assert h>0 and 'wss_validation_v2' in out and not os.path.exists(os.path.join(out,'volume.vtu'))
os.makedirs(out,exist_ok=True)
mesher=sv.meshing.TetGen();mesher.load_model(source);mesher.set_walls([1])
assert sorted(mesher.get_model_face_ids())==[1,2,3,4,5]
options=sv.meshing.TetGenOptions(global_edge_size=h,surface_mesh_flag=True,volume_mesh_flag=True);options.use_mmg=False
start=time.time();mesher.generate_mesh(options);grid=mesher.get_mesh();assert grid.GetNumberOfCells()>0
mesher.write_mesh(os.path.join(out,'volume.vtu'))
def save(data,name):
 w=vtk.vtkXMLPolyDataWriter();w.SetFileName(os.path.join(out,name));w.SetInputData(data);assert w.Write()==1
save(mesher.get_surface(),'surface.vtp')
for tag in mesher.get_model_face_ids():save(mesher.get_face_polydata(tag),'face_'+str(tag)+'.vtp')
record=dict(status='GENERATED',generator='sv.meshing.TetGen',request=request,source_sha256=hashlib.sha256(open(source,'rb').read()).hexdigest(),global_edge_size_m=h,nodes=grid.GetNumberOfPoints(),tetra=grid.GetNumberOfCells(),elapsed_s=time.time()-start,boundary_layers=False,anisotropy=False,mmg=False)
with open(os.path.join(out,'generation.json'),'w') as f:json.dump(record,f,indent=2)
print(json.dumps(record,indent=2))
