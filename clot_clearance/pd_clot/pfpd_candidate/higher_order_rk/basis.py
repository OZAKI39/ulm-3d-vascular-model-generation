"""Part I, Section 2.2, Example 1: nine monomials, no constant term."""
import numpy as np

NAMES=('X','Y','Z','X2','Y2','Z2','XY','XZ','YZ')

def basis(x,dimension=9):
    x=np.asarray(x);a,b,c=x[...,0],x[...,1],x[...,2]
    return np.stack((a,b,c,a*a,b*b,c*c,a*b,a*c,b*c),axis=-1)[...,:dimension]

def selector(delta,dimension=9):
    if delta<=0 or dimension not in (3,6,9):raise ValueError('Invalid selector')
    result=np.zeros((dimension,3));result[:3]=np.eye(3)/delta
    return result

def physical_polynomials(X):
    """Values and exact physical gradients of all nine global polynomials."""
    p=basis(X);grad=np.zeros((len(X),9,3));grad[:,:3,:]=np.eye(3)
    for c in range(3):grad[:,3+c,c]=2*X[:,c]
    for k,(a,b) in enumerate([(0,1),(0,2),(1,2)],6):
        grad[:,k,a]=X[:,b];grad[:,k,b]=X[:,a]
    return p,grad
