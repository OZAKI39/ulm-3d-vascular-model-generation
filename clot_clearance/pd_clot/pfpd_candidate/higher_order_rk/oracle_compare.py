"""Comparison primitives for the optional unmodified Julia package oracle."""
import numpy as np

def tensor_difference(python,julia,scale=.004):
    a=np.asarray(python)-np.asarray(julia)
    return float(np.linalg.norm(a,axis=(-2,-1)).max()/scale)

def representative_nodes(X,h):
    # Center, face, edge, corner; exact physical boundary nodes stress truncation.
    targets=np.array([[0,0,0],[0,0,-9*h],[-9*h,0,-9*h],[-9*h,-9*h,-9*h]])
    return np.array([np.linalg.norm(X-p,axis=1).argmin() for p in targets])
