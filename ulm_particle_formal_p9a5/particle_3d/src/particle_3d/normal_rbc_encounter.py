"""Undeformed biconcave RBC geometry and a restricted equatorial contact demo.

Shape: Evans–Fung profile (Dao et al. 2003, Eq. 7), initialized at the
project's RBC diameter/volume. The model does not solve a cell membrane.
Equatorial contact is exact: every body point has cylindrical radius <= R,
and the rim (R, z=0) belongs to the surface. A sphere centered in that plane,
outside the rim, therefore has gap |relative_center|-R-a.
"""
import numpy as np
from scipy.special import beta
from .particle_shapes import Sphere,roundoff_length
from .particle4_motion import isolated_sphere_contact_path

RBC_RADIUS=3.2875785225588278e-6
RBC_VOLUME=4.320087977386719e-17
MB_RADIUS=.588015722765549e-6
COEFFICIENTS=np.array([.207161,2.002558,-1.122762])
SHAPE_SOURCE='https://www.mit.edu/~mingdao/papers/JMPS_2003_Red_Blood_Cell.pdf'

def volume_scale(radius=RBC_RADIUS,volume=RBC_VOLUME):
    # V = 4 pi integral_0^R r z_top(r) dr, evaluated with beta integrals.
    integral=.5*sum(c*beta(k+1,1.5) for k,c in enumerate(COEFFICIENTS))
    return volume/(2*np.pi*radius**3*integral)

def surface_points(theta,phi):
    theta=np.asarray(theta);phi=np.asarray(phi);s=np.sin(theta);c=np.cos(theta)
    z=.5*RBC_RADIUS*volume_scale()*c*(COEFFICIENTS[0]+COEFFICIENTS[1]*s*s+COEFFICIENTS[2]*s**4)
    return np.stack(np.broadcast_arrays(RBC_RADIUS*s*np.cos(phi),RBC_RADIUS*s*np.sin(phi),z),axis=-1)

def mesh_arrays(nt=128,np_=192):
    # Single poles and shared periodic ring vertices form a closed manifold.
    rings=np.linspace(0,np.pi,nt+1)[1:-1];phis=np.arange(np_)*2*np.pi/np_
    points=np.vstack([surface_points(0.,0.),surface_points(rings[:,None],phis[None,:]).reshape(-1,3),surface_points(np.pi,0.)])
    faces=[]
    for j in range(np_):faces.append([0,1+j,1+(j+1)%np_])
    for k in range(nt-2):
        for j in range(np_):
            a=1+k*np_+j;b=1+k*np_+(j+1)%np_;c=a+np_;d=b+np_
            faces.extend([[a,c,b],[b,c,d]])
    end=len(points)-1;start=1+(nt-2)*np_
    for j in range(np_):faces.append([start+j,end,start+(j+1)%np_])
    return points,np.array(faces,dtype=int)

def advance_equatorial(x,y,vx,vy,dt):
    x,y,vx,vy=map(lambda v:np.asarray(v,dtype=float),(x,y,vx,vy))
    if not np.isfinite(dt) or dt<0:raise ValueError('Nonnegative finite time required')
    if abs(x[2]-y[2])>1e-18 or abs(vx[2])+abs(vy[2])>1e-18:
        raise ValueError('This demonstration only supports fixed-orientation equatorial contact')
    d=x-y;w=vx-vy;radius=RBC_RADIUS+MB_RADIUS;pad=roundoff_length(x,y,radius)
    gap=np.linalg.norm(d)-radius
    if gap< -pad:raise ValueError('Initial overlap')
    if dt==0:
        velocities=np.array([vx,vy])
        if gap<=pad and d@w<0:
            normal=d/np.linalg.norm(d);correction=.5*(w@normal)*normal
            velocities[0]-=correction;velocities[1]+=correction
        return x.copy(),y.copy(),velocities,dict(hit_time_s=None,contact_duration_s=0.,gap_m=float(gap),
            free_velocities_m_s=[vx.tolist(),vy.tolist()],position_projection=False,dt_s=0.)
    hit=None
    if gap<=pad and d@w<0:hit=0.
    elif d@w<0 and w@w>0:
        c=d@d-radius**2;dw=d@w;disc=dw*dw-(w@w)*c
        if disc>=0:
            root=c/(-dw+np.sqrt(disc))
            if 0<=root<=dt:hit=float(root)
    v={101:vx.copy(),203:vy.copy()};contact=0.
    if hit is None:a=x+dt*vx;b=y+dt*vy
    else:
        atx=x+hit*vx;aty=y+hit*vy
        exact=isolated_sphere_contact_path({101:Sphere(atx,RBC_RADIUS),203:Sphere(aty,MB_RADIUS)},
                                          {101:vx,203:vy},dt-hit)
        if exact is None:a=atx+(dt-hit)*vx;b=aty+(dt-hit)*vy
        else:
            shapes,v,rec=exact;a=shapes[101].center_m;b=shapes[203].center_m;contact=rec['contact_duration_s']
    gap=float(np.linalg.norm(a-b)-radius)
    if gap< -roundoff_length(a,b,radius):raise ValueError('Nonpenetration failure')
    return a,b,np.array([v[101],v[203]]),dict(hit_time_s=hit,contact_duration_s=contact,gap_m=gap,
        free_velocities_m_s=[vx.tolist(),vy.tolist()],position_projection=False,dt_s=dt)

def poiseuille(position,radius,mean_speed=.002):
    p=np.asarray(position);v=2*mean_speed*(1-(p[1]**2+p[2]**2)/radius**2)
    if v<0:raise ValueError('Outside idealized vessel')
    return np.array([v,0.,0.])
