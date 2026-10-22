"""Local fixed-ID rotation and paired motion diagnostics, without new physics."""
import time
import numpy as np,h5py
from scipy.spatial.distance import cdist,pdist
from scipy.sparse.csgraph import connected_components
from .common import *
from .metrics import rotation

def first_contact(path,lf=False):
    if lf:
        records=read(path.parent/'STEP_RECORDS.json')
        return next((r['time_s'] for r in records if r['min_gap_m']<=1e-9),None)
    with h5py.File(path) as h:
        x=h['positions_m'];t=h['time_s'];n=x.shape[1];ii,jj=np.triu_indices(n,1)
        for lo in range(0,len(x),128):
            q=x[lo:lo+128];m=(np.linalg.norm(q[:,ii]-q[:,jj],axis=2)-2*R).min(1);z=np.flatnonzero(m<=1e-9)
            if len(z):return float(t[lo+int(z[0])])
    return None

def interval_rotation(t,x):
    q=.5*(x[1:]+x[:-1]);v=np.diff(x,axis=0)/np.diff(t)[:,None,None];q-=q.mean(1)[:,None,:];v-=v.mean(1)[:,None,:];om=np.cross(q,v)[:,:,2].sum(1)/np.sum(q*q,axis=(1,2));return om,float(np.degrees(np.sum(om*np.diff(t))))

def main():
    begin=time.perf_counter();sci=read(BASE/'SCIENTIFIC_RESULTS.json');runs=read(BASE/'ALL_HF_RUNS.json');by={r['id']:r for r in runs};local=[];paired=[]
    for r in runs:
        if r['status']!='PASS':continue
        path=BASE/r['trajectory_path']
        with h5py.File(path) as h:
            t=h['time_s'][:];indices=np.unique(np.rint(np.linspace(0,len(t)-1,101)).astype(int));xx=h['positions_m'][indices];ix=np.minimum(indices,len(t)-2);vel=np.asarray([h['force_samples/velocity_m_s'][int(k)] for k in ix]);n,labels=connected_components(cdist(xx[-1],xx[-1])<=20e-6,directed=False);groups=[]
            for label in range(n):
                ids=np.flatnonzero(labels==label)
                if len(ids)<3:continue
                rr=[rotation(q[ids],v[ids]) for q,v in zip(xx,vel)];om=np.array([z['Omega_z_s'] for z in rr]);groups.append(dict(particle_ids=ids,N=len(ids),definition='Endpoint 20 um component; retrospectively selected fixed IDs, not preregistered success',time_s=t[indices],Omega_z_s=om,cumulative_angle_deg=float(np.degrees(np.trapezoid(om,t[indices]))),mean_rigid_fraction=float(np.mean([z['rigid_fraction'] for z in rr])),scope='Coarse sampled local diagnostic; no claim of sustained controllable mode or bubble spin'))
            local.append(dict(id=r['id'],fixed_ID_endpoint_groups=groups))
        lp=BASE/'lf'/r['id']/'SUMMARY.json'
        if not lp.exists():continue
        l=read(lp)
        if l['status']!='PASS':continue
        with h5py.File(BASE/l['trajectory_path']) as a,h5py.File(path) as b:
            t=a['sample_time_s'][:];xa=a['sample_positions_m'][:];ht=b['time_s'][:]
            # A 5 ms / 400 display grid includes half-microsecond times.
            # Compare only exact shared saved times; never interpolate HF truth.
            ix=np.minimum(np.searchsorted(ht,t-1e-13),len(ht)-1);keep=abs(ht[ix]-t)<1e-12
            discarded=int((~keep).sum());t=t[keep];xa=xa[keep];ix=ix[keep]
            assert len(t)>=3 and t[0]==0 and abs(t[-1]-ht[-1])<1e-12
            xb=b['positions_m'][ix]
        oa,aa=interval_rotation(t,xa);ob,ab=interval_rotation(t,xb);qa=xa[-1]-xa[-1].mean(0);qb=xb[-1]-xb[-1].mean(0);ca=qa.T@qa/len(qa);cb=qb.T@qb/len(qb);tha=first_contact(BASE/l['trajectory_path'],True);thb=first_contact(path)
        paired.append(dict(id=r['id'],shared_samples=len(t),nonshared_display_samples_omitted=discarded,rotation_measure='Same-interval midpoint angular rate derived from stored coordinates, not instantaneous force velocity',cumulative_LF_deg=aa,cumulative_HF_deg=ab,cumulative_difference_deg=abs(aa-ab),Omega_RMS_difference_s=float(np.sqrt(np.mean((oa-ob)**2))),endpoint_covariance_relative_error=float(np.linalg.norm(ca-cb)/np.linalg.norm(cb)),endpoint_nearest_neighbor_max_error_m=float(np.max(abs(np.sort(cdist(xa[-1],xa[-1]),axis=1)[:,1]-np.sort(cdist(xb[-1],xb[-1]),axis=1)[:,1]))),first_geometric_contact_LF_s=tha,first_geometric_contact_HF_s=thb,contact_presence_agrees=(tha is None)==(thb is None),contact_time_difference_s=abs(tha-thb) if tha is not None and thb is not None else None,contact_definition='surface gap <= 1 nm at accepted nodes; not a fusion event'))
    write(BASE/'LOCAL_ROTATION_CHECKS.json',dict(rows=local,wall_s=time.perf_counter()-begin,identity=identity()));write(BASE/'PAIRED_MOTION_CHECKS.json',dict(rows=paired,identity=identity()))
if __name__=='__main__':
    main();write(BASE/'POST_ANALYSIS_COMPLETE.json',dict(identity=identity(),diagnostic_time_alignment='Exact shared accepted HF times only'))
