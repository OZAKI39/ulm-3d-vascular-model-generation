"""Run only with SimVascular's embedded Python. Input and new output are explicit."""
import os, sys, json, time

def main():
    import sv, vtk
    source, output=sys.argv[1:]
    if os.path.exists(output):raise ValueError('Output already exists')
    os.makedirs(output)
    mesher=sv.meshing.TetGen();mesher.load_model(source)
    mesher.set_walls([1]);assert sorted(mesher.get_model_face_ids())==[1,2,3,4,5]
    options=sv.meshing.TetGenOptions(global_edge_size=0.0002,surface_mesh_flag=True,volume_mesh_flag=True)
    options.use_mmg=False
    start=time.time();mesher.generate_mesh(options)
    mesh=mesher.get_mesh();assert mesh.GetNumberOfCells()>0
    mesher.write_mesh(os.path.join(output,'volume.vtu'))
    def write_surface(data,name):
        writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(os.path.join(output,name));writer.SetInputData(data);assert writer.Write()==1
    write_surface(mesher.get_surface(),'surface.vtp')
    for fid in mesher.get_model_face_ids():write_surface(mesher.get_face_polydata(fid),'face_'+str(fid)+'.vtp')
    result=dict(status='GENERATED_PENDING_AUDIT',global_edge_size_m=.0002,points=mesh.GetNumberOfPoints(),cells=mesh.GetNumberOfCells(),elapsed_s=time.time()-start,generator='SimVascular 2023-05-31 sv.meshing.TetGen',boundary_layers=False)
    with open(os.path.join(output,'generation.json'),'w') as f:json.dump(result,f,indent=2)
    print(json.dumps(result))

if __name__=='__main__':main()

