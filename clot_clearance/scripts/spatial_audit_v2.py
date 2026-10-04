"""Isolated V2 audit adapters. All physical energies/forces and branch solves reuse verified code."""
from pathlib import Path
from unittest.mock import patch
from functools import lru_cache
import json, os, platform, shutil, time, resource
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator
from pd_clot.geometry import Cloud
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.energy import Material, constitutive
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate import coupling_stability_audit as ca
from pd_clot.pfpd_candidate import branch_stability_tracking as tracker
from pd_clot.pfpd_candidate.normalization import c0

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/pfpd_at1_spatial_discretization_audit_v2'
RESULTS=ROOT/'results/pfpd_at1_spatial_discretization_audit_v2'
L=0.00025
DELTA=0.00025125
LEVELS=(3,5,7)
EPS=(1e-5,3e-6,1e-6)
save=ca.save_json
read=ca.read_json

def declare():
    cfg=read(ROOT/'verification/pfpd_at1_branch_stability_tracking_v1/provenance/PREDECLARED_CONFIG.json')
    cfg.update(task='PFPD_AT1_SPATIAL_DISCRETIZATION_AUDIT_V2',
        family='CORRECTED_BOUNDARY_NODE_FAMILY',physical_L_m=L,physical_delta_m=DELTA,
        volume_relative_tolerance=2e-15,surface_relative_tolerance=2e-15,
        affine_F_absolute_tolerance=2e-12,affine_energy_relative_tolerance=2e-11,
        affine_balance_relative_tolerance=2e-12,external_relative_tolerance=1e-9,
        analytic_strain_absolute_tolerance=1e-13,proposal_relative_tolerance=1e-10,
        matrix_free_Hv_relative_tolerance=5e-5,matrix_free_sign_relative_tolerance=.02,
        stationary_jacobian_phase_step=.001,
        stationary_jacobian_step_rationale='Second-order one-sided derivatives of the quadratic phase-only gradient; separately normalized split vectors avoid cancellation under residual row scaling. Scientific HVP retains the original three epsilons.',
        spatial_interval_midpoint_tolerance=2e-5,onset_reaction_relative_tolerance=.02,
        stability_length_scale_m=L*ca.physics_config()['phase']['threshold_engineering_strain'],
        stability_energy_scale_J=ca.physics_config()['phase']['Gc_demo_J_per_m2']*L**2,
        matrix_free_scheme='central if phase bounds allow; split active-bound components and use legal second-order one-sided derivatives; algebraic mechanical probes are not physical branch states',
        initial_window=[.009,.00925,.0095,.00975,.009875,.0099375,.00996875,.009984375,.0099921875,.01,.010125,.01025,.0105,.011],
        schedules={'A':list(np.arange(.009,.01100001,.0005)), 'B':list(np.arange(.009,.01100001,.00025))},
        resource_limits={'reserve_total_RAM_fraction':.2,'dense_fraction_available_limit':.4,'serial_heavy_jobs':True})
    save(OUT/'provenance/PREDECLARED_CONFIG.json',cfg)
    import scipy,numba
    save(OUT/'provenance/ENVIRONMENT.json',dict(cpu=Path('/proc/cpuinfo').read_text().split('model name')[1].splitlines()[0],
        logical_cores=os.cpu_count(),memory=memory(),disk=shutil.disk_usage(ROOT)._asdict(),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,
        BLAS={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}))

def memory():
    m={k:int(v.split()[0])*1024 for k,v in (line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())}
    return dict(total=m['MemTotal'],available=m['MemAvailable'])

@lru_cache(None)
def fixture(n):
    h=L/(n-1);axis=np.linspace(-L/2,L/2,n);ijk=np.indices((n,n,n)).reshape(3,-1).T
    X=axis[ijk];w=np.full(n,h);w[[0,-1]]=h/2;V=np.prod(w[ijk],axis=1)
    top=ijk[:,2]==n-1;bottom=ijk[:,2]==0;surface=np.any((ijk==0)|(ijk==n-1),axis=1)
    fp=[];fn=[];fa=[]
    for a in range(3):
        for side,idx in [(-1,0),(1,n-1)]:
            ids=np.flatnonzero(ijk[:,a]==idx);fp.extend(ids)
            fn.extend(np.tile(np.eye(3)[a]*side,(len(ids),1)))
            fa.extend(np.prod(w[ijk[ids][:,[b for b in range(3) if b!=a]]],axis=1))
    fp=np.array(fp);pairs=cKDTree(X).query_pairs(DELTA,output_type='ndarray')
    f=Family.from_arrays(X,V,pairs,DELTA)
    cloud=Cloud(X,V,top|bottom,surface,fp,X[fp].copy(),np.array(fn),np.array(fa),f.pairs,f.xi,f.length,f.omega,h,DELTA)
    cfg=ca.physics_config()['material'];mat=Material(cfg['model'],cfg['mu'],cfg['bulk'])
    return cloud,f,mat,top

class Problem(ca.AuditProblem):
    def __init__(self,n,strain,lower=None,common_scale=False):
        data=fixture(n)
        if lower is None:lower=np.zeros(len(data[1].pairs))
        with patch.object(ca,'fixture',lambda config:data):
            super().__init__(ca.physics_config(),strain,lower)
        self.n=n
        if common_scale:
            cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json')
            self.length_scale=cfg['stability_length_scale_m'];self.energy_scale=cfg['stability_energy_scale_J']

def affine(p,strain=None):
    e=p.strain if strain is None else strain
    x=p.family.X.copy();x[:,2]+=e*(x[:,2]+L/2)
    return p.prescribe(x)

def aggregate(f,bond_values):
    a=np.zeros(len(f.X));den=np.zeros(len(f.X));w=2*f.wbar*f.volume[f.pairs[:,0]]*f.volume[f.pairs[:,1]]
    for ids in f.pairs.T:
        np.add.at(a,ids,w*np.abs(bond_values));np.add.at(den,ids,w)
    return a/den

def participation(a):
    a=np.abs(np.asarray(a));return dict(PR=float(np.sum(a*a)**2/(len(a)*np.sum(a**4))) if np.any(a) else None,
        top10_fraction=float(np.sort(a)[-int(np.ceil(.1*len(a))):].sum()/a.sum()) if np.any(a) else None)

def hessian_operator(p,z,eps,phase_eps=None):
    """Linear combination of legal bound-aware FD derivatives; no dense allocation."""
    z=np.asarray(z);g0=p.objective_gradient(z)[1];calls=[0]
    def grad(t):calls[0]+=1;return p.objective_gradient(t)[1]
    def forward(v):
        norm=np.linalg.norm(v);u=v/norm;step=eps if phase_eps is None else phase_eps
        q=u[p.nx:];s=z[p.nx:]
        if np.any(q>0):step=min(step,.49*np.min((p.sc-s[q>0])/q[q>0]))
        if np.any(q<0):step=min(step,.49*np.min((s[q<0]-p.lower[q<0])/(-q[q<0])))
        return norm*(-3*g0+4*grad(z+step*u)-grad(z+2*step*u))/(2*step)
    def mv(v):
        v=np.asarray(v).ravel();norm=np.linalg.norm(v)
        if norm==0:return np.zeros_like(v)
        v=v/norm;s=z[p.nx:];q=v[p.nx:]
        central=(s-p.lower>=2*eps)&(p.sc-s>=2*eps)
        if np.all(central):return norm*(grad(z+eps*v)-grad(z-eps*v))/(2*eps)
        base=v.copy();base[p.nx:][~central]=0
        bn=np.linalg.norm(base)
        out=bn*(grad(z+eps*base/bn)-grad(z-eps*base/bn))/(2*eps) if bn else np.zeros_like(v)
        for sign in [1,-1]:
            u=np.zeros_like(v);u[p.nx:]=np.where(~central,np.maximum(sign*q,0),0)
            # Lower-bound opening is forward; upper-bound closing is backward.
            upper=(p.sc-s<2*eps);lo=u.copy();lo[p.nx:][upper]=0
            if np.any(lo):out+=sign*forward(lo)
            hi=u.copy();hi[p.nx:][~upper]=0
            if np.any(hi):out-=sign*forward(-hi)
        return norm*out
    op=LinearOperator((len(z),len(z)),matvec=mv,rmatvec=mv,dtype=float);op.calls=calls
    return op

def matrix_free_residual_jacobian(p,z):
    """Same chain rule as tracker.residual_jacobian, expressed as operators."""
    cfg=tracker.settings()
    H=hessian_operator(p,z,cfg['residual_jacobian_epsilon'],phase_eps=cfg['stationary_jacobian_phase_step'])
    x,s=p.unpack(z);state=p.state(x,s);g=p.objective_gradient(z)[1]
    q=g[p.nx:]*p.energy_scale/p.weights;Y=state.W0;D=np.maximum(Y,p.yc);above=Y>p.yc
    a=p.energy_scale/p.weights*(1/D+above*q/(2*(1-s)*D**2))
    b=-q*above*Y/((1-s)*D**2);arg=s-q/D;clipped=(arg<p.lower)|(arg>p.sc)
    a[clipped]=0;b[clipped]=1
    scale=np.r_[np.full(p.nx,-p.energy_scale/(p.length_scale*p.config['mechanical_equilibrium']['force_floor_N'])),a]
    def diagonal(v):
        v=np.asarray(v).ravel();return np.r_[np.zeros(p.nx),b*v[p.nx:]]
    return LinearOperator(H.shape,matvec=lambda v:scale*(H@np.asarray(v).ravel())+diagonal(v),
        rmatvec=lambda v:H@(scale*np.asarray(v).ravel())+diagonal(v),dtype=float)

def stationary(p,x,s):
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json')
    if p.n>3:cfg['least_squares']['tr_solver']='lsmr'
    with patch.object(tracker,'settings',lambda:cfg):
        if p.n>3:
            with patch.object(tracker,'residual_jacobian',matrix_free_residual_jacobian):
                return tracker.solve_stationary_branch_state(p,x,s)
        return tracker.solve_stationary_branch_state(p,x,s)

def geometry_mechanics():
    rows=[];reactions=[]
    for n in LEVELS:
        start=time.time();p=Problem(n,0);c,f,mat,top=fixture(n);N=len(c.X);B=len(f.pairs);d=p.nx+B;mem=memory()
        estimate=dict(Nx=p.nx,Ns=B,total_DOF=d,dense_Hessian_bytes=8*d*d,temporary_fd_bytes=8*d*d*4,
            memory=mem,dense_permitted=8*d*d<=.4*mem['available'] and 32*d*d<mem['available']-.2*mem['total'])
        save(OUT/f'corrected_C{n}/MEMORY_ESTIMATE.json',estimate)
        state=evaluate(f,f.X,mat);eig=np.linalg.eigvalsh(state.K);rank=(eig>eig[:,-1,None]*1e-12).sum(axis=1)
        cond=eig[:,-1]/eig[:,0];counts=np.bincount(f.pairs.ravel(),minlength=N);dist=np.min(L/2-abs(f.X),axis=1)
        def stats(v):return dict(min=float(np.min(v)),mean=float(np.mean(v)),median=float(np.median(v)),max=float(np.max(v)))
        row=dict(level=f'C{n}',name=f'corrected_boundary_node_C{n}',N=n,particles=N,bonds=B,h_m=c.spacing,delta_m=f.horizon,
            delta_over_h=f.horizon/c.spacing,volume_m3=float(c.volume.sum()),area_m2=L**2,
            top_surface_sum_m2=float(c.face_area[c.face_normal[:,2]==1].sum()),bottom_surface_sum_m2=float(c.face_area[c.face_normal[:,2]==-1].sum()),
            family=stats(counts),surface_family=stats(counts[c.exposed]),interior_family=stats(counts[~c.exposed]),
            surface_nodes=int(c.exposed.sum()),interior_nodes=int((~c.exposed).sum()),rank_histogram={str(k):int((rank==k).sum()) for k in np.unique(rank)},
            K_condition=stats(cond),all_families_truncated=bool(np.all(dist<DELTA)),full_horizon_families=int((dist>=DELTA).sum()),
            Gc= p.config['phase']['Gc_demo_J_per_m2'],c0=c0(p.config['phase']['kernel']),Yc=p.yc)
        assert np.all(rank==3)
        assert abs(c.volume.sum()/L**3-1)<2e-15
        assert abs(row['top_surface_sum_m2']/L**2-1)<2e-15 and abs(row['bottom_surface_sum_m2']/L**2-1)<2e-15
        np.savez_compressed(OUT/f'geometry/C{n}.npz',X=f.X,volume=f.volume,pairs=f.pairs,family_size=counts,boundary_distance=dist,top=top,
            surface=c.exposed,K_condition=cond,bond_midpoints=f.X[f.pairs].mean(axis=1),face_particle=c.face_particle,face_normal=c.face_normal,face_area=c.face_area)
        theta=.37;rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
        cases=dict(translation=np.eye(3),rotation=rot,tension=np.diag([1,1,1.01]),compression=np.diag([1,1,.99]),
            shear=np.array([[1,.013,0],[0,1,0],[0,0,1]]),isotropic=np.eye(3)*1.01)
        checks=[]
        for name,F in cases.items():
            x=f.X@F.T+np.array([.00003,-.00002,.00004]);s=evaluate(f,x,mat);W,_=constitutive(F,mat)
            totalF=np.linalg.norm(s.force.sum(axis=0));force_scale=max(np.linalg.norm(s.force),1e-12)
            moment=np.linalg.norm(np.cross(x,s.force).sum(axis=0));energy_err=abs(s.energy-W*L**3)/max(W*L**3,1e-24)
            check=dict(case=name,Fbar_error=float(abs(s.Fbar-F).max()),Ftilde_error=float(abs(s.Ftilde-F).max()),
                energy_density_Pa=float(s.energy/L**3),analytic_density_Pa=float(W),energy_relative_error=float(energy_err),
                total_energy_J=s.energy,force_balance_relative=float(totalF/force_scale),moment_balance_relative=float(moment/(L*force_scale)))
            check['passed']=bool(check['Fbar_error']<2e-12 and check['Ftilde_error']<2e-12 and energy_err<2e-11 and
                check['force_balance_relative']<2e-12 and check['moment_balance_relative']<2e-12)
            checks.append(check)
        save(OUT/f'mechanics/C{n}_AFFINE.json',checks);row['intact_passed']=all(t['passed'] for t in checks)
        for e in [.005,.01,.015]:
            pp=Problem(n,e);x=affine(pp);s=pp.state(x,np.zeros(B));R=float(-s.force[top,2].sum())
            reactions.append(dict(level=f'C{n}',strain=e,reaction_N=R,nominal_stress_Pa=R/L**2,energy_J=s.energy,**{'metrics':pp.metrics(x,np.zeros(B))}))
        def driver(e):
            F=np.diag([1,1,1+e]);s=evaluate(f,f.X@F.T,mat);return float(np.average(s.W0,weights=p.weights))
        threshold=brentq(lambda e:driver(e)-p.yc,.001,.02,xtol=5e-16)
        row['analytic_threshold']=threshold;row['threshold_driver_spread_Pa']=float(np.ptp(evaluate(f,f.X@np.diag([1,1,1+threshold]),mat).W0))
        row['wall_seconds']=time.time()-start;row['process_peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        rows.append(row);save(OUT/f'corrected_C{n}/GEOMETRY_MECHANICS.json',row)
        print(f'C{n}: particles={N} bonds={B} affine_pass={row["intact_passed"]} threshold={threshold:.16g} free_force={reactions[-2]["metrics"]["force_residual"]:.6g}',flush=True)
    save(OUT/'geometry/LEVELS.json',rows);save(OUT/'mechanics/AFFINE_REACTIONS.json',reactions)
    save(OUT/'geometry/FIXED_HORIZON_CHECK.json',dict(passed=all(r['delta_m']==DELTA for r in rows),values=[r['delta_m'] for r in rows]))
    proposal=read(ROOT/'verification/pfpd_at1_spatial_discretization_audit_v1/comparisons/ISOLATED_COARSE_QUADRATURE_COMPARISON.json')
    current=next(r for r in reactions if r['level']=='C3' and r['strain']==.01)
    expected=dict(volume_m3=1.5625e-11,energy_J=1.5781640625000852e-12,reaction_N=1.2688125000000346e-6)
    errors=dict(volume=abs(rows[0]['volume_m3']/expected['volume_m3']-1),energy=abs(current['energy_J']/expected['energy_J']-1),reaction=abs(current['reaction_N']/expected['reaction_N']-1))
    save(OUT/'comparisons/CORRECTED_C3_PROPOSAL_CHECK.json',dict(passed=max(errors.values())<1e-10,expected=expected,recomputed=current,relative_errors=errors,prior_proposal=proposal))
    assert max(errors.values())<1e-10
    assert all(r['intact_passed'] for r in rows)
    assert np.ptp([r['analytic_threshold'] for r in rows])<1e-13
    return rows

if __name__=='__main__':
    declare();geometry_mechanics()
