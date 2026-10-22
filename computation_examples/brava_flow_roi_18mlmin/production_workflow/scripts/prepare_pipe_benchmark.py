"""Three predeclared mature-mesher policies, two measured radii, one analytic reference."""
from pathlib import Path
import json,hashlib
import numpy as np
import pyvista as pv
R=Path(__file__).resolve().parents[1];O=R/'straight_pipe';O.mkdir(exist_ok=True)
scale=json.loads((R/'vessel_scale/SCALE.json').read_text());radii=scale['radius_for_pipe'];rho=1056.;mu=.00345312
policies={
 'P0':dict(description='Ordinary uniform tetra',h_bulk_over_R=.4,angular_vertices=128,axial_spacing_over_R=.4,surface_remesh=True,layers=0),
 'P1':dict(description='Radius-scaled bulk tetra',h_bulk_over_R=1/6,angular_vertices=128,axial_spacing_over_R=1/6,surface_remesh=True,layers=0),
 'P2':dict(description='Radius-scaled bulk plus native boundary layers; actual surface remeshing geometry error measured',h_bulk_over_R=1/6,angular_vertices=1024,axial_spacing_over_R=1/3,surface_remesh=True,layers=4,edge_size_fraction=.5,layer_decreasing_ratio=.8,constant_thickness=True)
}
config=dict(radii=radii,policies=policies,length_over_R=12.,analysis_z_over_L=[.4,.5,.6],mean_velocity_m_s=.1,rho_kg_m3=rho,mu_Pa_s=mu,wall_distances_m=[1e-7,1e-6,1e-5],relative_wall_distances_R=[.05,.1,.25,.5],azimuth_count=32,azimuth_offset=0.371,velocity_relative_L2_limit=.05,gradient_WSS_relative_L2_limit=.08,physical_samples='Same exact analytical cylinder coordinates for every policy; no wall-relative repositioning or missing-probe deletion. Geometry coverage and PDE error reported separately.',maximum_tetra_per_pipe=750000,maximum_meshing_s_per_case=300,maximum_CFD_s_per_case=900,maximum_CFD_cases=6,automatic_BraVa_launch=False,poiseuille=dict(u='2*Umean*(1-r^2/R^2)',du_dr='-4*Umean*r/R^2',WSS='4*mu*Umean/R',Q='pi*R^2*Umean',dp_dz='-8*mu*Umean/R^2'))
(O/'BENCHMARK_PROTOCOL.json').write_text(json.dumps(config,indent=2)+'\n')
for rn,rad in radii.items():
 for pn,p in policies.items():
  name=rn.split('_')[0]+'_'+pn;d=O/name;d.mkdir(exist_ok=False)
  length=12*rad;nt=p['angular_vertices'];nz=int(np.ceil(12/p['axial_spacing_over_R']));theta=np.arange(nt)*2*np.pi/nt;z=np.linspace(0,length,nz+1)
  pts=np.array([[rad*np.cos(t),rad*np.sin(t),zz] for zz in z for t in theta]);pts=np.vstack([pts,[0,0,0],[0,0,length]])
  faces=[];tags=[]
  for k in range(nz):
   for j in range(nt):
    a=k*nt+j;b=k*nt+(j+1)%nt;c=(k+1)*nt+(j+1)%nt;e=(k+1)*nt+j;faces.extend([[a,b,c],[a,c,e]]);tags.extend([1,1])
  for j in range(nt):
   faces.extend([[len(pts)-2,(j+1)%nt,j],[len(pts)-1,nz*nt+j,nz*nt+(j+1)%nt]]);tags.extend([2,3])
  surf=pv.PolyData(pts,np.c_[np.full(len(faces),3),faces].ravel());surf.cell_data['ModelFaceID']=np.array(tags,dtype=np.int32);surf.save(d/'model.vtp')
  spec=dict(name=name,radius_m=rad,length_m=length,h_bulk_m=rad*p['h_bulk_over_R'],policy=p,policy_name=pn,expected_max_PLC_radial_chord_error_m=rad*(1-np.cos(np.pi/nt)),analytic=dict(Q_m3_s=np.pi*rad**2*.1,WSS_Pa=4*mu*.1/rad,pressure_gradient_Pa_m=-8*mu*.1/rad**2),model_sha256=hashlib.sha256((d/'model.vtp').read_bytes()).hexdigest())
  (d/'INPUT.json').write_text(json.dumps(spec,indent=2)+'\n')
print(json.dumps(config,indent=2))
