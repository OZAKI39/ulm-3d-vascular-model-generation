"""WSL post-download check against the unchanged closed STL, every Case C record."""
from pathlib import Path
import sys,json
import vtk,numpy as np
from vtk.util.numpy_support import numpy_to_vtk,vtk_to_numpy
root=Path(sys.argv[1]);p=root/'cases/case_c_real';r=np.genfromtxt(p/'trajectory_rank0.csv',delimiter=',',names=True);xyz=np.column_stack([r[k] for k in ['x','y','z']])
reader=vtk.vtkSTLReader();reader.SetFileName(str(root/'provenance/closed_geometry_m.stl'));reader.Update();surface=reader.GetOutput()
distance=vtk.vtkImplicitPolyDataDistance();distance.SetInput(surface);d=np.array([abs(distance.EvaluateFunction(x)) for x in xyz])
points=vtk.vtkPoints();points.SetData(numpy_to_vtk(xyz,deep=True));cloud=vtk.vtkPolyData();cloud.SetPoints(points)
enclosed=vtk.vtkSelectEnclosedPoints();enclosed.SetInputData(cloud);enclosed.SetSurfaceData(surface);enclosed.SetTolerance(1e-9);enclosed.Update();inside=vtk_to_numpy(enclosed.GetOutput().GetPointData().GetArray('SelectedPoints'))!=0
ct=json.loads((root/'contracts/CASE_C_SAFE_INTERIOR_CONTRACT.json').read_text());fc=json.loads((root/'contracts/FROZEN_FLOW_FIELD_CONTRACT.json').read_text());radius=r['diameter_m']/2
checks={'all_centers_inside_closed_STL':bool(inside.all()),'all_center_distances_gt_radius_plus_3dx':bool(np.all(d>radius+3*fc['dx_m'])),'recorded_Lipschitz_bound_conservative':bool(np.all(d>=r['center_wall_lower_bound_m']-1e-15))}
a={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'number_of_trajectory_points':len(xyz),'min_exact_center_wall_distance_m':float(d.min()),'min_exact_surface_clearance_m':float(np.min(d-radius)),'required_surface_margin_m':3*fc['dx_m'],'method':'VTK triangle distance and closed-surface enclosure at every saved step; surface clearance = center-wall distance - radius'}
np.savetxt(p/'EXACT_GEOMETRY_TRAJECTORY_AUDIT.csv',np.column_stack([r['step'],d,d-radius,inside]),delimiter=',',header='step,center_wall_m,surface_clearance_m,inside',comments='',fmt='%.17g')
(root/'LOCAL_EXACT_GEOMETRY_AUDIT.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a));raise SystemExit(0 if a['status']=='PASS' else 1)
