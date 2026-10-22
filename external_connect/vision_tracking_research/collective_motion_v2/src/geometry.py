"""Project-specific six-parameter and coordinate adapters, using SciPy/F10U."""
import numpy as np
from scipy.linalg import expm
from legacy import generator, parameters

def matrix_from_six(p):
    vx,vy,d,w,s1,s2=np.asarray(p,float)
    return np.array([vx,vy]),np.array([[d/2+s1,s2-w/2],[s2+w/2,d/2-s1]])

def transform(p,dt,center):
    v,A=matrix_from_six(p)
    return expm(generator(v,A,np.asarray(center))*dt)

def endpoint_six(H,dt,shape):
    h,w=shape
    return parameters(H[:2,:2],H[:2,2],dt,np.array([w/2,h/2]))

def resize_map(original_shape,new_shape):
    """OpenCV/PyTorch align_corners=False pixel-center convention."""
    sy,sx=np.asarray(new_shape,float)/np.asarray(original_shape,float)
    return np.array([[sx,0,(sx-1)/2],[0,sy,(sy-1)/2],[0,0,1.]])

def original_transform(H_resized,original_shape,new_shape):
    S=resize_map(original_shape,new_shape)
    return np.linalg.solve(S,H_resized@S)

def apply_transform(H,points):
    points=np.asarray(points)
    return (np.c_[points,np.ones(len(points))]@H.T)[:,:2]

def affine_residual(points,displacement):
    X=np.c_[np.ones(len(points)),np.asarray(points)]
    B=np.linalg.lstsq(X,displacement,rcond=None)[0]
    return float(np.sqrt(np.mean(np.sum((X@B-displacement)**2,axis=1))))
