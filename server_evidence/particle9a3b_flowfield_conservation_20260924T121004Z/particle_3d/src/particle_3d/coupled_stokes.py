"""Regularized-Stokeslet mobility of normal rigid RBC + shell-immobilized MB.

All numerical coordinates are nondimensional: length 1 um, speed 1 mm/s,
viscosity 0.00345312 Pa s. Surface forces are integrated nodal forces.
The fluid boundary solve includes the stationary cylindrical wall, both body
surfaces, and independent zero total force / zero total torque constraints.
There is no collision impulse, prescribed particle velocity, or repulsion.
"""
import numpy as np
from scipy.linalg import lu_factor,lu_solve,solve
from scipy.spatial.transform import Rotation
from .normal_rbc_encounter import RBC_RADIUS,MB_RADIUS,surface_points

LENGTH=1e-6;SPEED=1e-3;TIME=LENGTH/SPEED;MU=.00345312
CONFIGS={
    'coarse':dict(rbc_theta=10,rbc_phi=16,mb_n=42,wall_nx=29,wall_np=24,wall_half=24.,epsilon_factor=.30),
    'medium':dict(rbc_theta=14,rbc_phi=24,mb_n=82,wall_nx=37,wall_np=32,wall_half=24.,epsilon_factor=.30),
    'fine':dict(rbc_theta=18,rbc_phi=32,mb_n=162,wall_nx=49,wall_np=40,wall_half=24.,epsilon_factor=.30),
}

def rbc_nodes(nt,np_):
    theta=np.linspace(0,np.pi,nt+1)[1:-1];phi=np.arange(np_)*2*np.pi/np_
    return np.vstack([surface_points(0.,0.),surface_points(theta[:,None],phi[None,:]).reshape(-1,3),surface_points(np.pi,0.)])/LENGTH

def sphere_nodes(n,radius):
    # Exact antipodal pairs: no artificial center or inversion asymmetry.
    k=np.arange(n//2);z=(k+.5)/(n//2);phi=k*np.pi*(3-np.sqrt(5))
    p=np.c_[np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z]
    return np.vstack([p,-p])*radius

def cylinder_nodes(nx,np_,half,radius,shift=False):
    x=np.linspace(-half,half,nx);phi=(np.arange(np_)+(.5 if shift else 0))*2*np.pi/np_
    if shift:x=(x[:-1]+x[1:])/2
    xx,pp=np.meshgrid(x,phi,indexing='ij')
    return np.c_[xx.ravel(),radius*np.cos(pp.ravel()),radius*np.sin(pp.ravel())]

def kernel(target,source,epsilon):
    """Cortez 2001 three-dimensional regularized Stokeslet, mu*=1.

Source-dependent regularization resolves different surface length scales.
Consequently the collocation matrix is not assumed symmetric; use pivoted LU.
"""
    r=np.asarray(target)[:,None,:]-np.asarray(source)[None,:,:]
    r2=np.einsum('ijk,ijk->ij',r,r);e2=np.asarray(epsilon)**2
    scale=1/(8*np.pi*(r2+e2)**1.5)
    tensor=r[:,:,:,None]*r[:,:,None,:]
    tensor+=(r2+2*e2)[:,:,None,None]*np.eye(3)
    tensor*=scale[:,:,None,None]
    return tensor.transpose(0,2,1,3).reshape(3*len(target),3*len(source))

def background(points,radius=7.,mean=2.):
    points=np.asarray(points);v=np.zeros_like(points);v[:,0]=2*mean*(1-np.sum(points[:,1:]**2,axis=1)/radius**2)
    return v

def rigid_map(relative):
    n=len(relative);k=np.zeros((n,3,6));k[:,:,:3]=np.eye(3)
    x,y,z=np.asarray(relative).T
    k[:,0,4]=z;k[:,0,5]=-y;k[:,1,3]=-z;k[:,1,5]=x;k[:,2,3]=y;k[:,2,4]=-x
    return k.reshape(3*n,6)

class Mobility:
    def __init__(self,level='coarse',radius=7.,config=None):
        self.config=dict(CONFIGS[level] if config is None else config);c=self.config;self.radius=radius
        self.local=[rbc_nodes(c['rbc_theta'],c['rbc_phi']),sphere_nodes(c['mb_n'],MB_RADIUS/LENGTH)]
        self.counts=[len(p) for p in self.local];factor=c['epsilon_factor']
        self.eps_body=np.r_[np.full(self.counts[0],c.get('rbc_epsilon_factor',factor)*np.sqrt(79./self.counts[0])),
                            np.full(self.counts[1],c.get('mb_epsilon_factor',factor)*np.sqrt(4*np.pi*(MB_RADIUS/LENGTH)**2/self.counts[1]))]
        self.wall=cylinder_nodes(c['wall_nx'],c['wall_np'],c['wall_half'],radius)
        spacing=np.sqrt((2*c['wall_half']/(c['wall_nx']-1))*(2*np.pi*radius/c['wall_np']))
        self.eps_wall=np.full(len(self.wall),c.get('wall_epsilon_factor',factor)*spacing)
        self.wall_lu=lu_factor(kernel(self.wall,self.wall,self.eps_wall),check_finite=False)

    def solve(self,centers,rotations,*,active=(True,True),audit=True):
        """Return both 6DOF velocities and fluid forces in a translating query frame.

Only the finite wall quadrature window follows mean x; wall velocity is ZERO.
The original idealized cylinder is uniform and stationary.
"""
        centers=np.asarray(centers);rotations=np.asarray(rotations)
        origin=np.array([centers[:,0].mean(),0.,0.]);relative=[];parts=[];ids=[];eps=[]
        for i in range(2):
            if not active[i]:continue
            p=self.local[i]@rotations[i].T;relative.append(p);parts.append(p+centers[i]-origin);ids.append(i)
            start=0 if i==0 else self.counts[0];eps.extend(self.eps_body[start:start+self.counts[i]])
        points=np.vstack(parts);eps=np.array(eps);n=len(points);dofs=6*len(parts)
        k=np.zeros((3*n,dofs));offset=0
        for i,p in enumerate(relative):k[offset:offset+3*len(p),6*i:6*i+6]=rigid_map(p);offset+=3*len(p)
        pw=kernel(points,self.wall,self.eps_wall);wp=kernel(self.wall,points,eps)
        elimination=lu_solve(self.wall_lu,wp,check_finite=False)
        effective=kernel(points,points,eps)-pw@elimination
        rhs=background(points,self.radius).reshape(-1)
        inv=solve(effective,np.c_[k,rhs],assume_a='gen',check_finite=False)
        resistance=k.T@inv[:,:dofs];velocity=solve(resistance,k.T@inv[:,-1],assume_a='gen',check_finite=False)
        forces=inv[:,:dofs]@velocity-inv[:,-1];wallforces=-elimination@forces
        surface_error=effective@forces-k@velocity+rhs;balance=k.T@forces
        fullvelocity=np.full((2,6),np.nan)
        for j,i in enumerate(ids):fullvelocity[i]=velocity[6*j:6*j+6]
        record=dict(max_surface_collocation_error=float(np.max(abs(surface_error))),
            max_force_torque_balance_error=float(np.max(abs(balance))),
            resistance_reciprocity_relative_error=float(np.linalg.norm(resistance-resistance.T)/np.linalg.norm(resistance)),
            all_finite=bool(np.isfinite(velocity).all() and np.isfinite(forces).all()),
            body_nodes=[len(p) for p in parts],wall_nodes=len(self.wall),
            epsilon_body_um=[float(self.eps_body[0]),float(self.eps_body[-1])],epsilon_wall_um=float(self.eps_wall[0]))
        return dict(velocity=fullvelocity,body_points=points,body_forces=forces.reshape(-1,3),
            wall_forces=wallforces.reshape(-1,3),body_epsilon=eps,origin=origin,centers=centers,
            rotations=rotations,active=active,audit=record)

    def sample(self,points,solution,batch=256):
        p=np.asarray(points)-solution['origin'];output=background(p,self.radius)
        for start in range(0,len(p),batch):
            q=p[start:start+batch]
            v=kernel(q,solution['body_points'],solution['body_epsilon'])@solution['body_forces'].ravel()
            v+=kernel(q,self.wall,self.eps_wall)@solution['wall_forces'].ravel()
            output[start:start+batch]+=v.reshape(-1,3)
        return output

    def probe_residuals(self,solution):
        errors=[]
        probes=rbc_nodes(self.config['rbc_theta']+1,self.config['rbc_phi']+3)
        for i in range(2):
            if not solution['active'][i]:continue
            local=probes if i==0 else sphere_nodes(self.config['mb_n']+20,MB_RADIUS/LENGTH)
            r=local@solution['rotations'][i].T;points=r+solution['centers'][i]
            expected=solution['velocity'][i,:3]+np.cross(solution['velocity'][i,3:],r)
            err=np.linalg.norm(self.sample(points,solution)-expected,axis=1)
            errors.append(dict(body=i,max_slip_mm_s=float(err.max()),rms_slip_mm_s=float(np.sqrt(np.mean(err*err)))))
        wall=cylinder_nodes(self.config['wall_nx'],self.config['wall_np'],self.config['wall_half'],self.radius,True)
        wall=wall[np.abs(wall[:,0])<12]+solution['origin']
        err=np.linalg.norm(self.sample(wall,solution),axis=1)
        return dict(body=errors,wall_max_slip_mm_s=float(err.max()),wall_rms_slip_mm_s=float(np.sqrt(np.mean(err*err))))

def increment(centers,rotations,velocity,dt):
    x=np.asarray(centers)+dt*np.asarray(velocity)[:,:3]
    q=np.array([Rotation.from_rotvec(dt*v[3:]).as_matrix()@r for v,r in zip(velocity,rotations)])
    return x,q
