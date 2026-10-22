from pathlib import Path
import sys,json,shutil
import numpy as np
from scipy.spatial import Delaunay
from case_common import *

def build(nr,nt,nz,dt=1e-6,suffix=''):
 case=V/'stage2'/('pipe_nr%d'%nr+suffix);assert not case.exists(),case
 R=4e-6;L=40e-6;U=.002;xy=[[0.,0.]]
 for r in np.linspace(R/nr,R,nr):xy.extend(np.column_stack([r*np.cos(np.arange(nt)*2*np.pi/nt),r*np.sin(np.arange(nt)*2*np.pi/nt)]))
 xy=np.array(xy);tri=np.sort(Delaunay(xy/R).simplices,axis=1);n=len(xy)
 x=np.vstack([np.column_stack([xy,np.full(n,z)]) for z in np.linspace(0,L,nz+1)])
 a,b,c=(tri.T);slab=np.vstack([np.column_stack([a,b,c,c+n]),np.column_stack([a,b,b+n,c+n]),np.column_stack([a,a+n,b+n,c+n])])
 t=np.vstack([slab+k*n for k in range(nz)]);negative=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])<0;t[negative,:2]=t[negative,:2][:,::-1]
 b,own=boundary(x,t);g=x[b];tags=np.ones(len(b),np.int32);tags[np.all(abs(g[:,:,2])<1e-15,axis=1)]=4;tags[np.all(abs(g[:,:,2]-L)<1e-15,axis=1)]=2
 assert set(tags)=={1,2,4};write_mesh(case,x,t,b,tags,own)
 u=np.zeros_like(x);u[:,2]=2*U*(1-np.sum(x[:,:2]**2,axis=1)/R**2);u[np.unique(b[tags==1]),2]=0
 pp=8*MU*U*(L-x[:,2])/R**2
 av=.5*np.cross(g[:,1]-g[:,0],g[:,2]-g[:,0]);Qh=-np.einsum('ij,ij->i',av[tags==4],u[b[tags==4]].mean(axis=1)).sum()
 xml_case(case,dt,pipe=True)
 ids=np.unique(b[tags==4]);np.savetxt(case/'run/inlet_profile.txt',np.column_stack([ids,u[ids,2]]),fmt=['%d','%.17g'])
 face=pv.read(case/'SV_MESH/mesh-surfaces/OUTLET_03.vtp');tr=np.zeros((face.n_points,3));tr[:,:2]=-4*MU*U*face.points[:,:2]/R**2;face.point_data['Traction']=tr;face.save(case/'run/outlet_traction.vtp')
 dump(case/'policy.json',dict(case=case.name,kind='pipe_CFD',dt_s=dt,Q_target_m3_s=float(Qh),Q_analytical_m3_s=float(np.pi*R**2*U),Umean_m_s=U,mu_Pa_s=MU,nu_m2_s=MU/RHO,R_m=R,L_m=L,radial_intervals=nr,circumferential_intervals=nt,axial_intervals=nz,minimum_steps=20,steady_change_limit=1e-7,steady_consecutive_intervals=3,mass_limit=1e-6,maximum_wall_time_s=21600))
 np.savez_compressed(case/'analytical_nodal.npz',velocity_m_s=u,pressure_pa=pp)
 lock_case(case);print(case.name,len(x),len(t),'Qh/Q',Qh/(np.pi*R**2*U),flush=True)

if __name__=='__main__':
 for nr,nt,nz in [(4,24,10),(8,48,20),(16,96,40)]:build(nr,nt,nz)
 # Same mesh, halved time step: no geometrical variation.
 src=V/'stage2/pipe_nr16';dst=V/'stage2/pipe_nr16_halfdt';shutil.copytree(src,dst)
 pol=json.loads((dst/'policy.json').read_text());pol.update(case=dst.name,dt_s=5e-7);dump(dst/'policy.json',pol)
 tree=ET.parse(dst/'run/solver.xml');tree.find('.//Time_step_size').text='5e-7';tree.write(dst/'run/solver.xml',encoding='utf-8',xml_declaration=True);lock_case(dst)
