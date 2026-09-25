from pathlib import Path
import json,hashlib,math
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk
C=Path(__file__).resolve().parent;G=C/'geometry';fit=json.loads((G/'RBC_GEOMETRY_FIT_AUDIT.json').read_text());meta=json.loads((G/'NATIVE_MESH_GENERATION.json').read_text());dx=meta['dx_m'];dxum=dx*1e6
v=np.loadtxt(G/'RBC_REFERENCE_VERTICES.csv',delimiter=',',skiprows=1)[:,1:];f=np.loadtxt(G/'RBC_REFERENCE_TRIANGLES.csv',delimiter=',',skiprows=1,dtype=int)[:,1:]
reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(G/'BEST_RIGID_PLACEMENT.vtp'));reader.Update();rbc=reader.GetOutput();rv=vtk_to_numpy(rbc.GetPoints().GetData()).copy()
reader2=vtk.vtkXMLPolyDataReader();reader2.SetFileName(str(G/'FROZEN_LUMEN_LU.vtp'));reader2.Update();wall=reader2.GetOutput();wv=vtk_to_numpy(wall.GetPoints().GetData());wf=vtk_to_numpy(wall.GetPolys().GetData()).reshape(-1,4)[:,1:];tri=wv[wf]
best=fit['best_pose'];computed=(v-v.mean(axis=0))@np.array(best['rotation']).T+np.array(best['center']);assert np.max(np.abs(computed-rv))<1e-10
# Independent solid-angle winding classification, not vtkImplicitPolyDataDistance sign.
winding=[]
for i in range(0,len(rv),8):
 p=rv[i:i+8,None,:];a=tri[None,:,0,:]-p;b=tri[None,:,1,:]-p;c=tri[None,:,2,:]-p;la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
 numerator=np.einsum('ijk,ijk->ij',a,np.cross(b,c));denominator=la*lb*lc+np.einsum('ijk,ijk->ij',a,b)*lc+np.einsum('ijk,ijk->ij',b,c)*la+np.einsum('ijk,ijk->ij',c,a)*lb
 winding.extend((np.sum(2*np.arctan2(numerator,denominator),axis=1)/(4*math.pi)).tolist())
winding=np.array(winding);inside=np.abs(winding)>.5;assert np.all(np.minimum(np.abs(winding),np.abs(np.abs(winding)-1))<1e-6),'Winding ambiguous'
locator=vtk.vtkStaticCellLocator();locator.SetDataSet(wall);locator.BuildLocator();dist=[]
for p in rv:
 closest=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.);locator.FindClosestPoint(p,closest,cell,sub,d2);dist.append(math.sqrt(float(d2)))
clearance=np.array(dist)*np.where(inside,1.,-1.)
# Independent VTK mass properties for native reference mesh (before rigid translation).
props=vtk.vtkMassProperties();props.SetInputData(rbc);props.Update();volume=props.GetVolume()*dxum**3;area=props.GetSurfaceArea()*dxum**2
checks=dict(native_mesh_finite=bool(np.isfinite(v).all()),mesh_size=len(v)==642 and len(f)==1280,rigid_transform_identity=True,volume_identity=abs(volume-fit['mesh']['actual_volume_um3'])<1e-7,area_identity=abs(area-fit['mesh']['actual_area_um2'])<1e-7,signed_distance_identity=abs(float(clearance.min())-best['min_distance'])<1e-8,outside_vertex_count_identity=int(np.sum(~inside))==best['outside_vertices'],failure_gate_honored=fit['status']=='FAIL' and fit['continuous_verified_safe_placements']==0)
result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,scope='Independent STATIC geometry audit; no coupled solver ran, so SOLVER_FINALIZER_IDENTITY remains NOT_RUN',volume_um3=volume,area_um2=area,target_volume_um3=50,relative_volume_difference_from_nominal_target=(volume-50)/50,minimum_wall_clearance_lu=float(clearance.min()),minimum_wall_clearance_um=float(clearance.min()*dxum),outside_vertices=int(np.sum(~inside)),winding_min=float(winding.min()),winding_max=float(winding.max()),wall_penetration_metric='VERIFIED_STATIC_GEOMETRY_ONLY',RBC_runtime_timesteps=0)
(G/'INDEPENDENT_GEOMETRY_FINAL_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['status']=='PASS'
