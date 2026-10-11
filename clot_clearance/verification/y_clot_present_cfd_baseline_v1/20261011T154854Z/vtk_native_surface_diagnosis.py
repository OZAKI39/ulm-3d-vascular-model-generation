import vtk,json
r=vtk.vtkXMLPolyDataReader();r.SetFileName('/workspace/y_clot/project/verification/y_clot_present_cfd_baseline_v1/20261011T154854Z/clot_present_lumen_int32.vtp');r.Update();s=r.GetOutput();print("vtk",vtk.vtkVersion.GetVTKVersion(),"read",s.GetNumberOfPoints(),s.GetNumberOfCells(),flush=True)
c=vtk.vtkCleanPolyData();c.SetInputData(s);c.Update();s=c.GetOutput();s.BuildLinks();print("clean",s.GetNumberOfPoints(),s.GetNumberOfCells(),flush=True)
ids=vtk.vtkIdList();bad=[]
for i in range(s.GetNumberOfCells()):
 t=s.GetCell(i).GetPointIds()
 for j in range(t.GetNumberOfIds()):
  a=t.GetId(j);b=t.GetId((j+1)%t.GetNumberOfIds());s.GetCellEdgeNeighbors(i,a,b,ids)
  if ids.GetNumberOfIds()!=1:bad.append([i,a,b,ids.GetNumberOfIds()])
print("bad",len(bad),bad[:6],flush=True)
w=vtk.vtkXMLPolyDataWriter();w.SetFileName('/workspace/y_clot/project/verification/y_clot_present_cfd_baseline_v1/20261011T154854Z/vtk_native_cleaned.vtp');w.SetInputData(s);w.Write()
with open('/workspace/y_clot/project/verification/y_clot_present_cfd_baseline_v1/20261011T154854Z/NATIVE_CLEAN_DIAG.json',"w") as f:json.dump(dict(vtk=vtk.vtkVersion.GetVTKVersion(),points=s.GetNumberOfPoints(),cells=s.GetNumberOfCells(),bad_edges=bad),f)
