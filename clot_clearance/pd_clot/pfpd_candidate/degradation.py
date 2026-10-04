"""Unambiguous g/h factors for prescribed s only. NO evolution/history API."""
import numpy as np


def bounded(s):
    a=np.asarray(s,float)
    if not np.isfinite(a).all() or np.any((a<0)|(a>1)):raise ValueError('s must be finite in [0,1]')
    return a


def energetic(s):return (1-bounded(s))**2


def kinematic(s,sc=.95):
    a=bounded(s)
    if not np.isfinite(sc) or not 0<=sc<=1:raise ValueError('sc must be in [0,1]')
    if sc==1:return np.ones_like(a) # Disabled, including at s=1.
    return np.where(a<=sc,1.,((1-a)/(1-sc))**2)
