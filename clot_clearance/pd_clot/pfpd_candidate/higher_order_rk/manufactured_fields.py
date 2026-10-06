"""Only one new quadratic field: u_x=0.004 X Z/L. Prior fields are reused."""
import numpy as np
from ..nonlocal_boundary.nonaffine_manufactured import QuadraticAxial,QuadraticShear

def exact(field,X,L):
    if field in ('axial','shear'):
        f=QuadraticAxial() if field=='axial' else QuadraticShear()
        return f.position(X),f.response(X[:,2])['F']
    u=np.zeros_like(X);F=np.broadcast_to(np.eye(3),(len(X),3,3)).copy()
    if field=='cross':
        u[:,0]=.004*X[:,0]*X[:,2]/L;F[:,0,0]+=.004*X[:,2]/L;F[:,0,2]=.004*X[:,0]/L
    elif field in ('affine_axial','affine_shear'):
        a=2 if field=='affine_axial' else 0
        u[:,a]=.003*X[:,2];F[:,a,2]+=.003
    else:raise ValueError(field)
    return X+u,F
