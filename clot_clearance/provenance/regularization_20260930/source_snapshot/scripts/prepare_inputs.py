"""Generate canonical straight-pipe SI VTK and local traction import fixtures."""
import json,sys
from pathlib import Path
import numpy as np
import pyvista as pv
from scipy.spatial import Delaunay
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pd_clot.geometry import make_cloud
from pd_clot.streaming import pipe_field,make_provider,FileStreaming
from pd_clot.output import write_json


def main():
    root=Path(__file__).resolve().parents[1];out=root/'inputs';out.mkdir(exist_ok=True)
    c=json.loads((root/'configs/straight_pipe.json').read_text());p=c['pipe'];R=p['radius_m'];L=p['length_m']
    yz=[[0.,0.]]
    for radius in np.linspace(R/10,R,10):
        yz.extend([[radius*np.cos(a),radius*np.sin(a)] for a in np.linspace(0,2*np.pi,64,endpoint=False)])
    yz=np.asarray(yz);tri=np.sort(Delaunay(yz).simplices,axis=1);nc=len(yz)
    points=np.vstack([np.column_stack((np.full(nc,x),yz)) for x in np.linspace(-L/2,L/2,21)])
    tetra=[]
    for level in range(20):
        a,b,d=(tri+level*nc).T;A,B,D=a+nc,b+nc,d+nc
        tetra.extend(np.column_stack((a,b,d,D)));tetra.extend(np.column_stack((a,b,B,D)));tetra.extend(np.column_stack((a,A,B,D)))
    tetra=np.asarray(tetra,dtype=np.int64)
    det=np.linalg.det(points[tetra[:,1:]]-points[tetra[:,:1]])
    flip=det<0;tetra[flip,:2]=tetra[flip,:2][:,::-1]
    grid=pv.UnstructuredGrid(np.column_stack((np.full(len(tetra),4),tetra)).ravel(),np.full(len(tetra),10,np.uint8),points)
    u,pressure,grad=pipe_field(points,p);grid['Velocity']=u;grid['Pressure']=pressure
    grid.save(out/'straight_pipe_poiseuille.vtu')
    tube=pv.Cylinder(center=(0,0,0),direction=(1,0,0),radius=R,height=L,resolution=96,capping=False)
    tube.save(out/'straight_pipe_wall.vtp')
    g=make_cloud(c['clot']);static=json.loads(json.dumps(c))
    static['streaming']['oscillatory_amplitude']=0
    provider=make_provider(static);traction=provider.traction(g.face_position,g.face_normal,0)
    np.savetxt(out/'synthetic_traction.csv',np.column_stack((g.face_position,traction)),delimiter=',',header='x,y,z,tx,ty,tz',comments='')
    imported=FileStreaming(out/'synthetic_traction.csv',p['viscosity_Pa_s'])
    error=np.max(np.abs(imported.traction(g.face_position,g.face_normal,0)-traction))
    static['name']='analytic_static_reference';static['simulation']['number_of_macro_steps']=1;static['damage']['enabled']=False
    (root/'configs/analytic_static.json').write_text(json.dumps(static,indent=2)+'\n')
    static['name']='file_traction_demo';static['streaming']['provider']='file';static['streaming']['path']='../inputs/synthetic_traction.csv'
    (root/'configs/file_traction.json').write_text(json.dumps(static,indent=2)+'\n')
    write_json(out/'PIPE_IDENTITY.json',dict(kind='analytical_Poiseuille_no_CFD_solve',units='SI',points=len(points),tetrahedra=len(tetra),
        mu_Pa_s=p['viscosity_Pa_s'],rho_kg_m3=p['density_kg_m3'],Q_m3_s=p['flow_rate_m3_s'],
        mean_velocity_m_s=p['flow_rate_m3_s']/(np.pi*R*R),Re=2*p['density_kg_m3']*p['flow_rate_m3_s']/(np.pi*R*p['viscosity_Pa_s']),
        maximum_velocity_m_s=2*p['flow_rate_m3_s']/(np.pi*R*R),pressure_drop_Pa=8*p['viscosity_Pa_s']*p['flow_rate_m3_s']*L/(np.pi*R**4),
        minimum_tetra_volume_m3=float(np.abs(det).min()/6),file_traction_roundtrip_error_Pa=float(error),
        limitation='Empty-pipe analytical solution; clot does not alter flow; polygonal boundary approximates circle'))
    print('Generated straight pipe and traction importer fixture:',out)


if __name__=='__main__':main()
