"""Independent CPU energy-gradient check against archived native force probes."""
import numpy as np
from py_scripts.single_rbc_benchmark.physics import read_off
from py_scripts.single_rbc_benchmark.analysis import read_csv
from py_scripts.fluid_physics.common import write_json
from .forensics import OLD_RUNS


def audit(output):
    v,f=read_off(output/'common_reference.off')
    incident={}
    for i,face in enumerate(f):
        for k in range(3):
            a,b,c=int(face[k]),int(face[(k+1)%3]),int(face[(k+2)%3])
            incident.setdefault(tuple(sorted((a,b))),[]).append((i,c,a,b))
    pairs=np.array([[a[0],b[0],a[1],b[1]] for a,b in incident.values()])
    def energy(x,sign,theta):
        tri=x[f];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.linalg.norm(norm,axis=1)[:,None]
        na,nb=norm[pairs[:,0]],norm[pairs[:,1]]
        cosine=np.clip(np.sum(na*nb,axis=1),-1,1)
        orientation=np.sign(np.sum((na-nb)*(x[pairs[:,3]]-x[pairs[:,2]]),axis=1))*sign
        return 8*np.sum(1-np.cos(theta)*cosine-np.sin(theta)*orientation*np.sqrt(np.maximum(0,1-cosine*cosine)))
    a=read_csv(OLD_RUNS/'membrane_probe_candidate/native_response.csv');a=a[(a['mode']=='shear')&(a['epsilon']==0)]
    native=np.column_stack([a[k] for k in ('fx','fy','fz')])
    rows=[];best=None
    for sign in (-1,1):
        force=np.zeros_like(v);h=1e-5
        for i in range(len(v)):
            for j in range(3):
                p=v.copy();m=v.copy();p[i,j]+=h;m[i,j]-=h
                force[i,j]=-(energy(p,sign,np.deg2rad(6.97))-energy(m,sign,np.deg2rad(6.97)))/(2*h)
        error=float(np.linalg.norm(force-native)/np.linalg.norm(native))
        rows.append(dict(orientation_convention=sign,relative_force_l2_error=error,
                         force_rms=float(np.sqrt(np.mean(np.sum(force**2,axis=1)))),finite_difference_h=h))
        if best is None or error<best[0]:best=(error,force)
    np.savetxt(output/'initial_bending_cpu_vs_native.csv',np.column_stack([np.arange(len(v)),native,best[1]]),delimiter=',',
               header='vertex,native_fx,native_fy,native_fz,cpu_bend_fx,cpu_bend_fy,cpu_bend_fz',comments='')
    result=dict(method='central finite difference of constant-angle Kantor energy on the archived reference mesh; both face orientation conventions tested, no trajectory fit',
                energy='8 * sum_edges(1 - cos(theta_signed - theta0)), theta0=6.97 degrees',comparisons=rows,
                interpretation='Agreement supports bending reference stress as the source of the initial residual force, not as proof of the subsequent native crash cause.',
                native_force_contains_all_terms=True,precision='CPU float64 energy; archived native float32 force',new_gpu_seconds=0)
    write_json(output/'initial_bending_force_audit.json',result)
    return result
