"""Exact circular-pipe reference in SI; never substituted for a production CFD field."""
import numpy as np

def reference(points_m,radius_m,length_m,mean_velocity_m_s,mu_Pa_s):
    x=np.asarray(points_m,dtype=float)
    if x.shape[-1]!=3 or min(radius_m,length_m,mu_Pa_s)<=0 or mean_velocity_m_s<0:raise ValueError('Invalid pipe reference parameters')
    radial2=np.sum(x[...,:2]**2,axis=-1)
    if np.any(radial2>radius_m**2*(1+1e-12)):raise ValueError('Reference requested outside the analytical cylinder')
    u=np.zeros_like(x);u[...,2]=2*mean_velocity_m_s*(1-radial2/radius_m**2)
    G=np.zeros(x.shape[:-1]+(3,3));G[...,2,:2]=-4*mean_velocity_m_s*x[...,:2]/radius_m**2
    dp_dz=-8*mu_Pa_s*mean_velocity_m_s/radius_m**2
    p=dp_dz*(x[...,2]-length_m)
    return dict(u_m_s=u,G_s_inv=G,p_Pa=p,du_dr_s_inv=-4*mean_velocity_m_s*np.sqrt(radial2)/radius_m**2,WSS_Pa=4*mu_Pa_s*mean_velocity_m_s/radius_m,Q_m3_s=np.pi*radius_m**2*mean_velocity_m_s,pressure_gradient_Pa_m=dp_dz)
