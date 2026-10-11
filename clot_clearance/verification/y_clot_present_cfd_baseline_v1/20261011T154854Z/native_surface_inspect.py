import sv,vtk,os
m=sv.meshing.TetGen();m.load_model('/workspace/y_clot/project/verification/y_clot_present_cfd_baseline_v1/20261011T154854Z/clot_present_lumen_int32.vtp')
print("methods",dir(m),flush=True)
s=m.get_model_polydata();print("points",s.GetNumberOfPoints(),"cells",s.GetNumberOfCells(),flush=True)
w=vtk.vtkXMLPolyDataWriter();w.SetFileName('/workspace/y_clot/project/verification/y_clot_present_cfd_baseline_v1/20261011T154854Z/native_loaded_surface.vtp');w.SetInputData(s);w.Write()
