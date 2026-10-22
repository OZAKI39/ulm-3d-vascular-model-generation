"""Unchanged RBC palette/camera extracted from the retired flow renderer."""
import numpy as np
import matplotlib
from matplotlib.colors import ListedColormap

CMAP=ListedColormap(matplotlib.colormaps['turbo'](np.linspace(.10,.96,256)))

def camera(p,center,radius,angle,scale):
    theta=np.deg2rad(angle);phi=np.deg2rad(22)
    direction=np.array([np.cos(theta)*np.cos(phi),np.sin(theta)*np.cos(phi),np.sin(phi)])
    p.camera.position=center+4.5*radius*direction
    p.camera.focal_point=center
    p.camera.up=(0,0,1)
    p.camera.parallel_projection=True
    p.camera.parallel_scale=scale
    p.camera.clipping_range=(.01,1000)
