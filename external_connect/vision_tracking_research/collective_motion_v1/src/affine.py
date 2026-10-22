"""Interval-equivalent six parameters: log of a homogeneous endpoint transform."""
import time
import numpy as np
from scipy.linalg import logm
from scipy.optimize import least_squares
from scipy.spatial import ConvexHull
from common import NAMES,P

def parameters(F, b, dt, center):
    if dt<=0 or not np.isfinite(dt): raise ValueError('NONPOSITIVE_TIME')
    if np.linalg.det(F)<=0 or np.linalg.cond(F)>20: raise ValueError('NONINVERTIBLE_OR_ILL_CONDITIONED_TRANSFORM')
    H=np.eye(3);H[:2,:2]=F;H[:2,2]=b
    L=logm(H)/dt
    if np.max(np.abs(np.imag(L)))>1e-7: raise ValueError('NO_REAL_LOG_GENERATOR')
    L=L.real; A=L[:2,:2];v=L[:2,2]+A@center
    return np.array([v[0],v[1],np.trace(A),A[1,0]-A[0,1],(A[0,0]-A[1,1])/2,(A[0,1]+A[1,0])/2])

def empty(reason, method='UNKNOWN'):
    return {'method':method,'values':{k:None for k in NAMES},'valid':{k:False for k in NAMES},
            'reasons':{k:reason for k in NAMES},'observability_status':'NOT_OBSERVABLE',
            'uncertainty':None,'sample_count':0,'coverage_fraction':0.,'condition_number':None,'residual':None}

def fit_affine(xy, displacement, dt, shape, method, bootstrap=0):
    start=time.perf_counter();h,w=shape[:2];center=np.array([w/2,h/2]);scale=max(h,w)/2
    xy=np.asarray(xy,float).reshape(-1,2);d=np.asarray(displacement,float).reshape(-1,2)
    good=np.isfinite(xy).all(1)&np.isfinite(d).all(1);xy=xy[good];d=d[good]
    out=empty('INSUFFICIENT_OR_DEGENERATE_SUPPORT',method);out['sample_count']=len(xy)
    if len(xy)<P['piv']['min_vectors']: return out
    X=np.c_[np.ones(len(xy)),(xy-center)/scale]
    rank=np.linalg.matrix_rank(X);cond=np.linalg.cond(X)
    try: coverage=ConvexHull(xy).volume/(h*w)
    except Exception: coverage=0.
    out.update(design_rank=int(rank),condition_number=float(cond),coverage_fraction=float(coverage))
    if rank<3 or cond>P['piv']['condition_max'] or coverage<P['piv']['coverage_min']: return out
    initial=np.linalg.lstsq(X,d,rcond=None)[0]
    fit=least_squares(lambda a:(X@a.reshape(3,2)-d).ravel(),initial.ravel(),loss='huber',f_scale=.4,max_nfev=50)
    B=fit.x.reshape(3,2);pred=X@B;res=np.linalg.norm(pred-d,axis=1)
    F=np.eye(2)+B[1:].T/scale;b=B[0]-(F-np.eye(2))@center
    quadrants=[]
    for right in [False,True]:
        for bottom in [False,True]:
            mask=((xy[:,0]>=center[0])==right)&((xy[:,1]>=center[1])==bottom)
            quadrants.append(float(np.median(res[mask])) if mask.any() else None)
    out.update(F=F.tolist(),b=b.tolist(),residual=float(np.sqrt(np.mean(res**2))),optimizer_success=bool(fit.success),dt_s=dt,reference_center_xy=center.tolist(),residual_quadrants_px=quadrants)
    try: value=parameters(F,b,dt,center)
    except ValueError as exc: return {**out,'failure':str(exc)}
    out['raw_estimate']=value.tolist()
    if not fit.success or out['residual']>P['piv']['residual_max_px']:
        out['reasons']={k:'NON_AFFINE_OR_INCONSISTENT_VECTORS' for k in NAMES};return out
    # Spatial block bootstrap (overlapping PIV vectors are NOT independent).
    # These are empirical, nominal 95% intervals, whose coverage is measured separately.
    samples=[];rng=np.random.default_rng(1329)
    blocks=((xy[:,0]>=center[0]).astype(int)+2*(xy[:,1]>=center[1]).astype(int))
    for _ in range(bootstrap):
        selection=np.concatenate([np.flatnonzero(blocks==j) for j in rng.integers(0,4,4)])
        if len(selection)<6 or np.linalg.matrix_rank(X[selection])<3:continue
        C=np.linalg.lstsq(X[selection],d[selection],rcond=None)[0];FF=np.eye(2)+C[1:].T/scale;bb=C[0]-(FF-np.eye(2))@center
        try:samples.append(parameters(FF,bb,dt,center))
        except ValueError:pass
    if len(samples)>=8:
        ci=np.percentile(samples,[2.5,97.5],axis=0).T
        out['uncertainty']={'method':'spatial_quadrant_block_bootstrap_nominal_95','calibrated':False,'n':len(samples),'intervals':dict(zip(NAMES,ci.tolist()))}
    out.update(values=dict(zip(NAMES,value.tolist())),valid={k:True for k in NAMES},reasons={k:None for k in NAMES},observability_status='INTERVAL_AFFINE_ESTIMATE',fit_ms=(time.perf_counter()-start)*1000)
    return out
