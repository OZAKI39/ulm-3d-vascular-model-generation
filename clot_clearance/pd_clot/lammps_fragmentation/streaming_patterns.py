"""Replaceable provider emits physical per-surface-node force vectors in N."""
import numpy as np
from .schema import validate_force_table

def build_streaming_load(request, geometry):
    r=request;X=geometry['X'];n=geometry['normals'];mask=geometry['surface'];R=r.get('radius_m',np.ptp(X,axis=0).max()/2)
    center=np.asarray(r.get('load_center_m',[0.,0.,R]),float);width=float(r.get('width_m',.4*R));direction=np.asarray(r.get('direction',[1.,0.,0.]),float);direction/=np.linalg.norm(direction)
    axis=np.asarray(r.get('rotation_axis',[0.,0.,1.]),float);axis/=np.linalg.norm(axis)
    if r['pattern_id']=='LOCAL_GAUSSIAN_SHEAR':
        raw=np.exp(-np.sum((X-center)**2,axis=1)/(2*width**2))[:,None]*direction
    elif r['pattern_id']=='OPPOSING_SHEAR':
        # Opposing x-directed patches on either side of the crown.
        c1=center+[-.50*R,0.,-.15*R];c2=center+[.50*R,0.,-.15*R]
        w1=np.exp(-np.sum((X-c1)**2,axis=1)/(2*width**2));w2=np.exp(-np.sum((X-c2)**2,axis=1)/(2*width**2))
        raw=(w2-w1)[:,None]*direction
    elif r['pattern_id']=='VORTEX_TANGENTIAL':
        raw=np.cross(axis,X-np.asarray(r.get('center_m',[0.,0.,0.]))) / R
    else: raise ValueError('Only the three declared analytic patterns are available')
    tangent=raw-np.sum(raw*n,axis=1)[:,None]*n;tangent[~mask]=0.
    F=tangent*geometry['areas'][:,None];scale=r['load_budget_N']/np.linalg.norm(F,axis=1).sum();F*=scale
    ids=np.flatnonzero(mask)+1;validate_force_table(ids,F[mask],len(X))
    traction=np.linalg.norm(F[mask],axis=1)/geometry['areas'][mask]
    check=float(np.max(np.abs(np.sum(F[mask]*n[mask],axis=1))))
    return dict(force_N=F,info=dict(provider_id=r['load_provider'],pattern_id=r['pattern_id'],surface_region='FREE_CURVED_REFERENCE_SURFACE',
        center_m=center,width_m=width,direction=direction,rotation_axis=axis,
        amplitude_N=r['load_budget_N'],integrated_force_budget_N=float(np.linalg.norm(F,axis=1).sum()),peak_local_traction_Pa=float(traction.max()),
        resultant_force_N=F.sum(0),resultant_torque_Nm=np.cross(X,F).sum(0),tangency_absolute_error_N=check,
        ramp_time_s=r['ramp_time_s'],hold_time_s=r['hold_time_s'],ramp_rule='Half-cosine 0 to 1, then hold',normalization=r['normalization'],
        source=r['source'],scope='EXPLORATORY_ONLY',load_identity=r['load_identity'],
        frame='Reference-surface material-node force vectors; continue acting on those nodes after detachment; not a solved Eulerian fluid field'))
