"""Frictionless kinematic velocity constraint only; no forces or position edits."""
from itertools import combinations
import numpy as np
from .particle_shapes import Sphere,Capsule,EPS


class ContactConstraintError(RuntimeError):pass


def contact_velocity(free_velocity,free_omega,center,gap,*,translation_only=False):
    v=np.asarray(free_velocity,dtype=float);omega=np.asarray(free_omega,dtype=float)
    if gap.state!='TOUCHING':return v.copy(),omega.copy()
    arm=gap.particle_point_m-np.asarray(center)
    normal=gap.normal_inward
    vn=float(normal@(v if translation_only else v+np.cross(omega,arm)))
    return v-min(vn,0.)*normal,omega.copy()


def constrain_contacts(shape,free_velocity,free_omega,contacts):
    """Minimum necessary center-velocity correction for simultaneous contacts.

For one normal this is exactly V_free-min(v_contact·n,0)n. At an edge,
several independent unilateral normals can be active; enumerate <=3 active
constraints of the three-dimensional Euclidean projection QP. No tangential
friction or Omega correction is introduced.
"""
    v=np.asarray(free_velocity,dtype=float);omega=np.asarray(free_omega,dtype=float)
    if not contacts:return v.copy(),omega.copy(),0.
    ns=[];bs=[]
    for _,gap in contacts:
        n=gap.normal_inward
        b=0. if isinstance(shape,(Sphere,Capsule)) else -float(n@np.cross(omega,gap.particle_point_m-shape.center_m))
        if any(np.linalg.norm(n-old)<256*EPS and abs(b-ob)<=256*EPS*max(abs(b),abs(ob),np.linalg.norm(v),np.finfo(float).tiny) for old,ob in zip(ns,bs)):continue
        ns.append(n);bs.append(b)
    n=np.array(ns);b=np.array(bs);tol=1024*EPS*max(np.linalg.norm(v),np.max(np.abs(b)),np.finfo(float).tiny)
    if np.all(n@v>=b-tol):return v.copy(),omega.copy(),0.
    candidates=[]
    for count in range(1,min(3,len(n))+1):
        for ids in combinations(range(len(n)),count):
            a=n[list(ids)];rhs=b[list(ids)]-a@v
            multipliers=np.linalg.lstsq(a@a.T,rhs,rcond=128*EPS)[0]
            corrected=v+a.T@multipliers
            if np.min(multipliers)>=-tol and np.all(n@corrected>=b-tol):
                candidates.append((float(np.linalg.norm(corrected-v)),tuple(ids),corrected))
    if not candidates:raise ContactConstraintError('Simultaneous normal constraints are infeasible with fixed free Omega')
    _,_,corrected=min(candidates,key=lambda x:(x[0],x[1]))
    error=float(max(0.,np.max(b-n@corrected)))
    return corrected,omega.copy(),error
