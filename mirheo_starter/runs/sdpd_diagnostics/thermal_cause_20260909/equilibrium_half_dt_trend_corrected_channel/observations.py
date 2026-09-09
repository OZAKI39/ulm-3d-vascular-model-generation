"""Small, GPU-import-free observation helpers; copied into each diagnostic run."""
import numpy as np


def cuda_layout(interface):
    shape=tuple(int(n) for n in interface['shape']);dtype=np.dtype(interface['typestr'])
    if not shape or any(n<0 for n in shape) or dtype.hasobject:raise ValueError('INVALID_CUDA_LAYOUT')
    strides=interface.get('strides')
    if strides is None:
        strides=[];stride=dtype.itemsize
        for n in reversed(shape):strides.insert(0,stride);stride*=n
    strides=tuple(int(s) for s in strides)
    if len(strides)!=len(shape) or any(s<0 for s in strides):raise ValueError('UNSUPPORTED_CUDA_STRIDES')
    size=0 if any(n==0 for n in shape) else sum((n-1)*s for n,s in zip(shape,strides))+dtype.itemsize
    if size and not interface['data'][0]:raise ValueError('NULL_NONEMPTY_CUDA_POINTER')
    return shape,dtype,strides,size


def phase_audit(pre_x,pre_v,force,post_x,post_v,dt,mass,domain):
    arrays=[np.asarray(a,float) for a in (pre_x,pre_v,force,post_x,post_v)]
    if len({a.shape for a in arrays})!=1 or arrays[0].ndim!=2 or arrays[0].shape[1]!=3:
        raise ValueError('PHASE_SHAPE_MISMATCH')
    if not all(np.isfinite(a).all() for a in arrays):raise ValueError('NONFINITE_PHASE_STATE')
    x,v,f,x1,v1=arrays
    dv=v1-v-dt*f/mass;dx=x1-x-dt*v1;L=np.asarray(domain)
    dx-=L*np.rint(dx/L)
    return {'kick_max_abs':float(np.max(abs(dv),initial=0)),
            'drift_PBC_max_abs':float(np.max(abs(dx),initial=0)),
            'kick_RMS':float(np.sqrt(np.mean(dv**2))) if dv.size else None,
            'drift_PBC_RMS':float(np.sqrt(np.mean(dx**2))) if dx.size else None}


def velocity_statistics(vel,mass):
    vel=np.asarray(vel,float);N=len(vel)
    if N<=1:raise ValueError('INSUFFICIENT_THERMAL_DOF')
    centered=vel-vel.mean(axis=0);energy=np.sum(centered**2,axis=1)*mass
    order=np.sort(energy)
    return {'N':N,'COM_velocity':vel.mean(axis=0).tolist(),
            'component_kBT':(mass*np.sum(centered**2,axis=0)/(N-1)).tolist(),
            'COM_kBT':float(energy.sum()/(3*(N-1))),
            'max_speed':float(np.max(np.linalg.norm(vel,axis=1))),
            'top_one_percent_thermal_energy_fraction':float(order[-max(1,int(np.ceil(.01*N))):].sum()/energy.sum()) if energy.sum() else None}
