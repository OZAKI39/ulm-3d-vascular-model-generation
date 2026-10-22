"""Extended FNM caches and force-based pattern feasibility, not particle motion."""
import time,json,resource
import numpy as np,h5py
from scipy.optimize import brentq
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs,CachedWaterBasis
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context as old_context
from scripts.straight_tube_exp1.f10ae.response import ExtendedResponse
from .common import *
def response_context(f=2500000.):
    c=old_context(f)
    if f==2500000.:
        assert read(AE/'LOW_PRESSURE_ADMISSION.json')['production_extension_admitted']
        c['response']=ExtendedResponse(c['response'],c['physics'])
    return c
def direct_force(model,response,x,a,p,order=96):
    h,g=model.evaluate(x,order);u=np.asarray(a)*np.exp(1j*np.asarray(p));P=u@h;G=np.einsum('i,inj->nj',u,g)
    valid=(abs(P)==0)|((abs(P)>=response.knots[0])&(abs(P)<=response.knots[-1]))
    F=np.full((len(P),3),np.nan);F[valid]=-.5*np.real(response.query(P[valid])[:,0,None]*G[valid].conj())
    return P,G,F,valid
def build(f):
    start=time.perf_counter();folder=BASE/'fields'/str(int(f));folder.mkdir(parents=True,exist_ok=False)
    cfg=read(CONFIG)['field'];elements=read(OLD/'CONFIG_SNAPSHOT.json')['geometry']['pzt'];model=WaterDiscs(elements,f)
    rng=np.random.default_rng(462001);lo=np.array(cfg['bounds_m'][0]);hi=np.array(cfg['bounds_m'][1]);q=rng.uniform(lo+30e-6,hi-30e-6,(128,3))
    p,g=model.evaluate(q,96);p48,g48=model.evaluate(q,48);pr,gr=model.rayleigh_reference(q[:12],96,256);eps=1e-9
    gd=np.stack([(model.evaluate(q+np.eye(3)[j]*eps,96)[0]-model.evaluate(q-np.eye(3)[j]*eps,96)[0])/(2*eps) for j in range(3)],-1)
    norm=lambda a,b:float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-99))
    checks=dict(quadrature_pressure=norm(p48,p),quadrature_gradient=norm(g48,g),rayleigh_pressure=norm(pr,p[:,:12]),rayleigh_gradient=norm(gr,g[:,:12]),analytic_gradient=norm(gd,g));cache=[]
    for name,step in [('coarse',cfg['coarse_step_m']),('fine',cfg['fine_step_m'])]:
        axes=[np.linspace(l,h,round((h-l)/step)+1) for l,h in zip(lo,hi)];shape=tuple(map(len,axes));P=np.empty((4,)+shape,complex);flat=P.reshape(4,-1)
        # Bounded memory: only a chunk of coordinates and quadrature work arrays.
        for begin in range(0,np.prod(shape),2048):
            ii=np.arange(begin,min(begin+2048,np.prod(shape)));ix=np.unravel_index(ii,shape);x=np.column_stack([axes[j][ix[j]] for j in range(3)]);flat[:,begin:begin+len(ii)]=model.evaluate(x,48)[0]
        meta=dict(model='HOMOGENEOUS_WATER_CONTROL_PROXY',frequency_Hz=f,sound_speed_m_s=1500.,elements=elements,quadrature_order=48,harmonic_convention='Re[P exp(+i omega t)]',identity=identity())
        path=folder/(name+'.h5')
        with h5py.File(path,'x') as h:
            for s,a in zip('xyz',axes):h['axes/'+s+'_m']=a
            h.create_dataset('pressure_basis',data=P,compression='gzip');h.attrs['metadata_json']=json.dumps(meta)
        cp,cg=CachedWaterBasis(path).query_basis(q);checks[name+'_pressure']=norm(cp,p);checks[name+'_gradient']=norm(cg,g);cache.append(dict(path=str(path.relative_to(BASE)),sha256=sha(path),bytes=path.stat().st_size,shape=shape))
    limits={k:cfg['quadrature_limit'] for k in checks}
    limits.update(analytic_gradient=1e-6,coarse_pressure=1e-4,fine_pressure=1e-5,coarse_gradient=.005,fine_gradient=.0013)
    passed=all(checks[k]<=limits[k] for k in checks)
    out=dict(frequency_Hz=f,checks=checks,limits=limits,PASS=passed,caches=cache,wall_s=time.perf_counter()-start,rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,identity=identity())
    write(folder/'FIELD_GATE.json',out);print('FIELD',f,passed,checks,flush=True);return out
def diagnostic(f):
    c=response_context(f);model=WaterDiscs(c['discs'].elements,f);actions=read(BASE/'ACTIONS.json');folder=BASE/'fields'/str(int(f));axis=np.linspace(-.5e-3,.5e-3,101);xx,yy=np.meshgrid(axis,axis,indexing='ij');xy=np.c_[xx.ravel(),yy.ravel(),np.zeros(xx.size)];rows=[]
    with h5py.File(folder/'DIAGNOSTICS.h5','x') as h:
        h['xy_points_m']=xy;h['xy_axis_m']=axis
        for name,a in actions.items():
            P,G,F,v=direct_force(model,c['response'],xy,a['amplitudes'],a['phases']);I=-np.imag(P.conj()[:,None]*G)/(2*c['physics'].density_kg_m3*2*np.pi*f);r=xy-xy.mean(0);relative=F-np.nanmean(F,0);tangent=np.cross(r,relative)[:,2];curl=np.gradient(F[:,1].reshape(xx.shape),axis,axis=0)-np.gradient(F[:,0].reshape(xx.shape),axis,axis=1)
            for k,z in dict(pressure_Pa=P,gradient_Pa_m=G,primary_force_N=F,valid=v,intensity_W_m2=I,curl_z_N_m=curl).items():h[name+'/'+k]=z
            rows.append(dict(action=name,frequency=f,valid_fraction=float(v.mean()),P_max_Pa=float(abs(P).max()),mean_torque_z_N_m=float(np.nanmean(tangent)),relative_tangential_force_rms_N=float(np.sqrt(np.nanmean((tangent/np.maximum(np.linalg.norm(r,axis=1),1e-15))**2))),curl_rms_N_m=float(np.sqrt(np.nanmean(curl**2))),scope='STATIC_PRIMARY_FORCE_ONLY; F2 requires a particle configuration'))
        # Actual F1 equilibria; pressure extrema alone are not labeled traps.
        line=np.linspace(-.55e-3,.55e-3,2201);x=np.c_[line,np.zeros((len(line),2))];traps=[]
        for deg in [0,45,90,180]:
            a=[0,0,1000,1000];p=[0,0,0,np.radians(deg)];P,G,F,v=direct_force(model,c['response'],x,a,p)
            h[f'line_{deg}/x_m']=line;h[f'line_{deg}/pressure_Pa']=P;h[f'line_{deg}/force_N']=F
            roots=[]
            for j in range(len(line)-1):
                if not(v[j] and v[j+1]) or not(F[j,0]>0>=F[j+1,0]):continue
                fun=lambda t:direct_force(model,c['response'],np.array([[t,0,0]]),a,p)[2][0,0]
                try:
                    z=brentq(fun,line[j],line[j+1],xtol=1e-13);der=(fun(z+1e-8)-fun(z-1e-8))/2e-8;roots.append(dict(x_m=z,dFdx_N_m=der,stable_X=bool(der<0)))
                except (ValueError,RuntimeError):pass
            traps.append(dict(phase_deg=deg,stable_X_roots=roots,spacing_m=np.diff([r['x_m'] for r in roots]),ideal_half_wavelength_m=1500/f/2,scope='X-direction restoring force; full 3D and dynamical stability not implied'))
    write(folder/'STATIC_FEASIBILITY.json',dict(rows=rows,stripes=traps,identity=identity()))
def main():
    verify();start=time.perf_counter();rr=[]
    for f in [2500000.,2400000.,2600000.]:
        r=build(f);rr.append(r);diagnostic(f)
    write(BASE/'FIELD_BUILD.json',dict(rows=rr,all_pass=all(r['PASS'] for r in rr),wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
