"""SI Poiseuille reference and pipe frame, independent of solver assembly."""
import numpy as np


def pipe_frame(origin,axis):
    origin=np.asarray(origin,dtype=float)
    axis=np.asarray(axis,dtype=float)
    if not np.isfinite(origin).all() or not np.isfinite(axis).all() or np.linalg.norm(axis)==0:
        raise ValueError("Invalid pipe frame")
    axis=axis/np.linalg.norm(axis)
    candidate=np.eye(3)[np.argmin(abs(axis))]
    first=np.cross(axis,candidate); first/=np.linalg.norm(first)
    second=np.cross(axis,first)
    return origin,axis,first,second


def pipe_theory(radius,length,mu,signed_flow,rho):
    if radius<=0 or length<=0 or mu<=0 or rho<=0:
        raise ValueError("Positive SI parameters required")
    mean=signed_flow/(np.pi*radius**2)
    return {"mean_velocity_m_s":float(mean),"centerline_velocity_m_s":float(2*mean),
        "delta_pressure_pa":float(8*mu*length*signed_flow/(np.pi*radius**4)),
        "reynolds_number":float(rho*abs(mean)*2*radius/mu)}


def poiseuille_at(points,origin,axis,radius,length,mu,signed_flow):
    origin,axis,_,_=pipe_frame(origin,axis)
    offset=np.asarray(points)-origin
    axial=offset@axis
    radial=offset-axial[...,None]*axis
    mean=signed_flow/(np.pi*radius**2)
    speed=2*mean*(1-np.sum(radial**2,axis=-1)/radius**2)
    pressure=8*mu*length*signed_flow/(np.pi*radius**4)*(1-axial/length)
    return speed[...,None]*axis,pressure
