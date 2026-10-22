"""Force-equilibrium stability, winding, and independent null calibration."""
import time
import numpy as np,h5py
from scipy.optimize import root
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs
from .common import *
from .field import response_context,direct_force
from .metrics import spatial

def main():
    begin=time.perf_counter();verify();rows=[];winding=[];failed_roots=[]
    for f in [2500000.,2400000.,2600000.]:
        c=response_context(f);model=WaterDiscs(c['discs'].elements,f)
        for name in ['X_PAIR','Y_PAIR','XY_INPHASE','XY_QUADRATURE','ROT_PLUS','ROT_MINUS']:
            a=read(BASE/'ACTIONS.json')[name];sol=[]
            def force(x):return direct_force(model,c['response'],np.array([x]),a['amplitudes'],a['phases'])[2][0]
            def fun(uv):return force(np.r_[uv*1e-3,0])[:2]/1e-12
            for x in np.linspace(-.45,.45,5):
                for y in np.linspace(-.4,.4,5):
                    try:rr=root(fun,[x,y])
                    except ValueError as ex:
                        failed_roots.append(dict(frequency=f,action=name,seed_mm=[x,y],reason=str(ex)));continue
                    q=np.r_[rr.x*1e-3,0]
                    if not rr.success or max(abs(q[:2]))>.5e-3 or np.linalg.norm(fun(rr.x))>1e-5 or any(np.linalg.norm(q-z['position_m'])<1e-6 for z in sol):continue
                    eps=1e-7;J=np.stack([(force(q+np.eye(3)[j]*eps)-force(q-np.eye(3)[j]*eps))/(2*eps) for j in range(3)],axis=1)
                    if not np.isfinite(J).all():continue
                    ev=np.linalg.eigvals(J);sol.append(dict(position_m=q,force_N=force(q),jacobian_N_m=J,eigenvalues_real=ev.real,eigenvalues_imag=ev.imag,stable_3D=bool(np.all(ev.real< -1e-12)),stable_XY=bool(np.all(np.linalg.eigvals(J[:2,:2]).real< -1e-12))))
            rows.append(dict(frequency=f,action=name,roots=sol,scope='Local F1 linearization at z=0; no F2, no proof of nonlinear particle pattern'))
            if name.startswith('ROT'):
                for rad in [50e-6,150e-6,300e-6,450e-6]:
                    ang=np.linspace(0,2*np.pi,721);q=np.c_[rad*np.cos(ang),rad*np.sin(ang),np.zeros(len(ang))];P,G,F,v=direct_force(model,c['response'],q,a['amplitudes'],a['phases']);wind=(np.unwrap(np.angle(P))[-1]-np.unwrap(np.angle(P))[0])/(2*np.pi);tan=np.c_[-np.sin(ang),np.cos(ang),np.zeros(len(ang))];tau=np.sum(F*tan,axis=1);winding.append(dict(frequency=f,action=name,radius_m=rad,phase_winding=wind,pressure_min_Pa=float(abs(P).min()),force_valid=bool(v.all()),mean_tangential_F1_N=float(np.nanmean(tau)),rms_tangential_F1_N=float(np.sqrt(np.nanmean(tau**2)))))
    write(BASE/'STATIC_EXTRA.json',dict(equilibria=rows,winding=winding,failed_root_searches=failed_roots,wall_s=time.perf_counter()-begin,identity=identity()))
    begin=time.perf_counter();rng=np.random.default_rng(462099);records=[]
    for n in [32,64,128,200]:
        peaks=[];false_stripe=0;false_lattice=0
        for j in range(128):
            x=rng.uniform([-.0005,-.0004,-.00002],[.0005,.0004,.00002],(n,3));m=spatial(x);peaks.append(m['structure_peak_S']);false_stripe+=m['stripe_descriptor_pass'];false_lattice+=m['lattice_descriptor_pass']
        records.append(dict(N=n,samples=128,null_peak_S=peaks,empirical_99=float(np.quantile(peaks,.99)),admission_peak=max(12.,float(np.quantile(peaks,.99))),raw_false_stripe=false_stripe,raw_false_lattice=false_lattice))
    write(BASE/'METRIC_NULL_CALIBRATION.json',dict(rows=records,scope='Analysis-only random geometry, not acoustic trajectories; calibration independent of physical discoveries',wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
