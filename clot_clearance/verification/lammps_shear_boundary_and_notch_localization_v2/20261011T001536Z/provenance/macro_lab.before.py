"""Bounded macroscopic LPS tests, free interiors and independent read-only energy.
Run with the private migrated environment; one MPI rank per case. No production edits.
"""
from pathlib import Path
import argparse,ctypes,csv,json,math,os,time,traceback
import numpy as np
from scipy.spatial import cKDTree
from .interface_model import attach_wall,write_typed_data
from .lammps_writer import force_conversion

K=100000/3;G=940.;RHO=1000.;DELTA=150.5e-6

def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False))
def geometry(dx,kind,periodic=False):
    dims=np.array([.0006,.0006,.001]);counts=np.rint(dims/dx).astype(int)
    axes=[(np.arange(n)+.5)*dx-d/2 for n,d in zip(counts,dims)]
    X=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);N=len(X)
    b=(cKDTree(X+[.0003,.0003,.0005],boxsize=[.0006,.0006,10.]) if periodic else cKDTree(X)).query_pairs(DELTA*(1+1e-12),output_type='ndarray')
    b=b[np.lexsort((b[:,1],b[:,0]))]
    g=dict(X=X,bonds=b,bond_types=np.ones(len(b),int),clot_count=N,anchors=np.zeros(N,bool),base=X[:,2]<X[:,2].min()+dx/2,surface=np.zeros(N,bool),areas=np.zeros(N),normals=np.zeros_like(X),info={})
    if kind=='M6':g=attach_wall(g,dict(dx_m=dx,delta_m=DELTA,rho_kg_m3=RHO,anchor_layers=3))
    removed=np.empty((0,2),int)
    if kind=='M5':
        a,c=X[b[:,0]],X[b[:,1]];cross=((a[:,2]<-1e-12)!=(c[:,2]<-1e-12))
        frac=np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(b)),where=c[:,2]!=a[:,2])
        xx=a[:,0]+frac*(c[:,0]-a[:,0]);cut=cross&(xx<-.0001)
        removed=b[cut];g['bonds']=b[~cut];g['bond_types']=g['bond_types'][~cut]
    return g,removed

class Energy:
    def __init__(self,g,dx,lib,periodic=False):
        self.X=g['X'];self.b=np.ascontiguousarray(g['bonds'],dtype=np.int32);self.V=dx**3
        self.periodic=periodic
        delta=self.X[self.b[:,1]]-self.X[self.b[:,0]]
        if periodic:delta[:,:2]-=.0006*np.rint(delta[:,:2]/.0006)
        self.r=np.linalg.norm(delta,axis=1)
        self.h=np.where(abs(self.r-DELTA)<=dx/2,(DELTA+dx/2-self.r)/dx,1.)
        self.m=np.bincount(self.b.ravel(),weights=np.repeat(self.r*self.V*self.h,2),minlength=len(self.X))
        self.lib=ctypes.CDLL(str(lib));self.fun=self.lib.evaluate_periodic if periodic else self.lib.evaluate;self.fun.restype=ctypes.c_double
        self.fun.argtypes=[ctypes.c_int,ctypes.c_int]+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p]+([ctypes.c_double] if periodic else [])
    def eval(self,x,alive,force=False):
        f=np.zeros_like(x) if force else None
        arr=[x,self.b,alive,self.r,self.h,self.m]
        u=self.fun(len(x),len(self.b),*[a.ctypes.data for a in arr],self.V,K,G,f.ctypes.data if force else None,*([.0006] if self.periodic else []))
        return (u,f) if force else u

def run(out,case,dx=50e-6,dt=.5e-6,slow=1.,static=False):
    start=time.perf_counter();out.mkdir(parents=True,exist_ok=False);os.chdir(out)
    kind=case[:2] if case.startswith('M') else 'M5'
    periodic=case=='M3P'
    g,removed=geometry(dx,kind,periodic);X=g['X'];n=len(X);nc=g['clot_count'];b=g['bonds'];types=g['bond_types']
    r=dict(dx_m=dx,delta_m=DELTA,rho_kg_m3=RHO,K_Pa=K,G_Pa=G,s00_bulk=.15,s00_interface=.140624859375)
    write_typed_data(out/'body.data',g,r)
    if periodic:
        lines=(out/'body.data').read_text().splitlines();lines=[('-0.0003 0.0003 '+line.split()[-2]+' '+line.split()[-1]) if line.endswith(('xlo xhi','ylo yhi')) else line for line in lines];(out/'body.data').write_text('\n'.join(lines)+'\n')
    # Physical 150 um grips; axial tangents remain free except M2 prescribed shear.
    mask=np.zeros((n,3),bool);sign=np.zeros((n,3));grip=150e-6
    low=X[:,2]<-.0005+grip-1e-12;high=(X[:,2]>.0005-grip+1e-12)&~g['anchors']
    if kind=='M6':low=g['anchors'];mask[low]=True;mask[high,2]=True;sign[high,2]=1
    else:
        mask[low|high,2]=True;sign[low,2]=-.5;sign[high,2]=.5
    if kind=='M2':
        mask[low|high,0]=True;sign[:,0]=sign[:,2];sign[:,2]=0
    if kind=='M3' and not periodic:
        mask[np.abs(X[:,0])>.0003-grip+1e-12,0]=True
        mask[np.abs(X[:,1])>.0003-grip+1e-12,1]=True
    H=float(X[high,2].mean()-(0 if kind=='M6' else X[low,2].mean()))
    # Interface origin is the anchor/clot midplane, z=-0.5 mm.
    if kind=='M6':H=float(X[high,2].mean()+.0005)
    small=kind in ['M1','M2','M3'];amp=(-.008 if kind in ['M1','M3'] else .008 if kind=='M2' else .8)*(.00085 if kind in ['M4','M5'] else H)
    ramp=(.010 if small else .024)*slow;duration=2*ramp if small else ramp;steps=0 if static else round(duration/dt)
    roi=(np.abs(X[:,0])<.0003-DELTA)&(np.abs(X[:,1])<.0003-DELTA)&(np.abs(X[:,2])<.0005-grip-DELTA)&~g['anchors']
    centered=X[roi]-X[roi].mean(0);roi_pinv=np.linalg.pinv(centered)
    plan=dict(case=case,kind=kind,material=r,dt_s=dt,steps=steps,duration_s=duration,ramp_s=ramp,amplitude_m=amp,nodes=n,clot_nodes=nc,bonds=len(b),removed_notch_bonds=len(removed),size_m=[.0006,.0006,.001],grip_m=grip,grip_center_distance_m=H,area_m2=.0006**2,roi_nodes=int(roi.sum()),roi_definition='|x|,|y|<150um approximately; |z|<199.5um; centroid least squares full deformation gradient',notch='z=0 from x=-300um to -100um, through y; projected single plane area',interface='Inherited attach_wall: only basal clot layer connects to 3 anchor layers; not a full-horizon bonded half-space',loading='Raised cosine; load-unload for M1-M3, monotone for M4-M6; fixed physical time across dt/dx',bc='fix move NULL components integrate naturally; no overlapping nve group; only selected boundary DOFs prescribed',damping_s_inv=500,stop='Fixed duration, no outcome-based stopping',uncertainty='Free-interior macro response/independent energy; '+('time-step sensitivity' if case=='T1' else 'two-level spacing sensitivity, fixed horizon; contact distance changes with dx' if case=='D1' else kind),decision='Assess mapping/identifiability; never tune production',output='Every accepted step scalar force, ROI strain and energy; independent native topology events; initial/first/peak/end double states',thresholds=dict(modulus_relative=.15,nu_absolute=.02,KE_over_U=.05,dt_relative=.10,dx_relative=.15,energy_residual=.05),fit='Loading 20%-70% imposed small amplitude; no outcome-selected fit interval')
    plan['transverse_periodic_fixed_cell']=periodic
    if periodic:plan['bc']='Periodic x/y with fixed600um cell; only top/bottom axial grips prescribed; all other DOFs free; not a rerun with density/damping changes'
    save(out/'test_plan.json',plan);np.savez_compressed(out/'geometry.npz',**{k:g[k] for k in ['X','bonds','bond_types','anchors']},prescribed_mask=mask,roi=roi,removed_notch_bonds=removed)
    root=Path(__file__).resolve().parents[2];lib=root/('analysis/material_macro_audit_v1/energy_periodic.so' if periodic else 'analysis/material_macro_audit_v1/energy.so');energy=Energy(g,dx,lib,periodic)
    from lammps import lammps
    L=lammps(cmdargs=['-log',str(out/'log.lammps'),'-screen','none'])
    text=f'''units si\ndimension 3\nboundary m m m\natom_style bond\natom_modify map array sort 0 0\nnewton on off\nspecial_bonds lj 0 1 1 coul 1 1 1\nlattice sc {dx:.17g}\nread_data body.data extra/bond/per/atom 320 extra/special/per/atom 360\nneighbor {2*dx:.17g} bin\nneigh_modify every 1 delay 0 check yes\ncomm_modify cutoff {3*DELTA:.17g}\nfix props all property/atom d_vfrac d_theta d_damage ghost yes\nset group all d_vfrac {dx**3:.17g}\npair_style bpm/peri contact_distance_factor 1.0\npair_coeff * * {18*K/(math.pi*DELTA**4):.17g} {DELTA:.17g}\nbond_style bpm/peri\nbond_coeff 1 lps {K:.17g} {G} {DELTA} .15 0\nbond_coeff 2 lps {K:.17g} {G} {DELTA} .140624859375 0\nbond_coeff 3 lps {K:.17g} {G} {DELTA} 1 0\nvelocity all set 0 0 0 units box\n{force_conversion()}\ntimestep {dt:.17g}\nvariable zero equal 0\nvariable disp equal {amp:.17g}*0.5*(1-cos(PI*time/{ramp:.17g}))\nvariable vel equal {amp:.17g}*0.5*PI/{ramp:.17g}*sin(PI*time/{ramp:.17g})\n'''
    # Disjoint integrator groups. Sign values are constant within each group.
    keys=np.column_stack([mask.astype(int),sign]);unique,inverse=np.unique(keys,axis=0,return_inverse=True)
    for j,key in enumerate(unique):
        ids=np.flatnonzero(inverse==j)+1;text+=f'group g{j} id '+ ' '.join(map(str,ids))+'\n'
        if not mask[ids[0]-1].any():text+=f'fix integrate{j} g{j} nve\n';continue
        uv=[];vv=[]
        for d in range(3):
            if not key[d]:uv.append('NULL');vv.append('NULL')
            else:
                text+=f'variable u{j}_{d} equal {key[d+3]:.17g}*v_disp\nvariable v{j}_{d} equal {key[d+3]:.17g}*v_vel\n';uv.append(f'v_u{j}_{d}');vv.append(f'v_v{j}_{d}')
        text+=f'fix integrate{j} g{j} move variable '+ ' '.join(uv+vv)+'\n'
    text+='''compute contactPE all pe pair
compute bondPE all pe bond
fix pe_record all ave/time 1 1 1 c_contactPE c_bondPE file native_energy_density.dat
fix audit all external pf/callback 1 1
thermo 1000
thermo_style custom step time atoms bonds ke c_contactPE c_bondPE
thermo_modify lost error lost/bond error flush yes
'''
    if periodic:text=text.replace('boundary m m m','boundary p p m')
    (out/'in.lammps').write_text(text+'# Python registers callback then run '+str(steps)+'\n');L.commands_string(text)
    alive=np.ones(len(b),np.uint8);keyb=b[:,0]*n+b[:,1];prev=None;W=Dn=Db=0.;events=[];peak=-np.inf;first=False;calls=0;forcechecks=[];negative_releases=[]
    stream=(out/'history.csv').open('w');writer=csv.writer(stream);writer.writerow(['step','time_s','displacement_m','boundary_strain','stress_Pa','reaction_N','eps_x','eps_y','eps_z','gamma_xz','K_free_J','U_bulk_J','U_contact_J','D_num_J','D_break_independent_J','W_ext_J','R_energy_J','native_bond_pe_times_V_J','broken_bulk','broken_interface','broken_wall','KE_over_U','max_speed_m_s','roi_affine_residual_m'])
    def state(label,step,x,v):np.savez_compressed(out/(label+'.npz'),step=step,time_s=step*dt,x_m=x,v_full_m_s=v,alive=alive)
    def callback(caller,step,nlocal,tags,x,fext):
      nonlocal alive,prev,W,Dn,Db,peak,first,calls
      try:
        assert nlocal==n and np.array_equal(tags,np.arange(1,n+1)),'One-rank ordering required'
        step=int(step);calls+=1;x=np.ascontiguousarray(x);fm=L.numpy.extract_atom('f')[:n].copy();vh=L.numpy.extract_atom('v')[:n].copy()
        fd=-500*RHO*dx**3*vh;fd[mask]=0;fext[:]=fd
        nb=int(L.extract_global('nbonds'));changed=nb!=int(alive.sum())
        if changed:
            _,raw=L.gather_bonds();a=np.ctypeslib.as_array(raw).reshape(-1,3);ends=np.sort(a[:,1:].astype(np.int64)-1,axis=1);idx=np.searchsorted(keyb,ends[:,0]*n+ends[:,1]);new=np.zeros_like(alive);new[idx]=1
            assert np.all(new<=alive);broken=np.flatnonzero(alive&~new);before=energy.eval(x,alive);after=energy.eval(x,new);release=before-after;Db+=release
            if release<0:negative_releases.append([step,release])
            events.extend((step,int(k),int(types[k])) for k in broken);alive=new
        U,fp=energy.eval(x,alive,True)
        Uc=float(L.extract_compute('contactPE',0,0))*dx**3;native=float(L.extract_compute('bondPE',0,0))*dx**3
        vf=vh.copy();vf[~mask]+=.5*dt*(fm+fd)[~mask]/(RHO*dx**3)
        kinetic=.5*RHO*dx**3*float(np.sum(vf[~mask]**2))
        if prev is not None and step>prev[0]:
            dxstep=x-prev[1];W-=float(np.sum(.5*(fm+prev[2])[mask]*dxstep[mask]));Dn-=float(np.sum(.5*(fd+prev[3])*dxstep))
        if not changed and Uc<1e-20 and step in [100,1000,5000]:
            err=float(np.linalg.norm(fp-fm)/max(np.linalg.norm(fm),1e-30));q=int(np.argmax(abs(fp).ravel()));node,dim=divmod(q,3);hh=1e-10;xx=x.copy();xx[node,dim]+=hh;up=energy.eval(xx,alive);xx[node,dim]-=2*hh;um=energy.eval(xx,alive);fdcheck=-(up-um)/(2*hh)
            forcechecks.append(dict(step=step,native_force_relative_error=err,finite_difference_relative_error=abs(fdcheck-fp[node,dim])/max(abs(fp[node,dim]),1e-30),native_energy_relative_error=(native-U)/max(abs(U),1e-30)))
        direction=0 if kind=='M2' else 2;reaction=-fm[high,direction].sum();stress=reaction/plan['area_m2']
        disp=float((x[high,direction]-X[high,direction]).mean()-(0 if kind=='M6' else (x[low,direction]-X[low,direction]).mean()))
        F=(roi_pinv@(x[roi]-x[roi].mean(0))).T;strain=F-np.eye(3);res=float(np.sqrt(np.mean(((x[roi]-x[roi].mean(0))-centered@F.T)**2)))
        counts=[int(np.sum((1-alive)[types==t])) for t in [1,2,3]]
        writer.writerow([step,step*dt,disp,disp/H,stress,reaction,strain[0,0],strain[1,1],strain[2,2],strain[0,2]+strain[2,0],kinetic,U,Uc,Dn,Db,W,W-kinetic-U-Uc-Dn-Db,native,*counts,kinetic/U if U>1e-30 else 0,float(np.linalg.norm(vf,axis=1).max()),res])
        if step==0:state('initial',step,x,vf)
        if changed and not first:state('first_break',step,x,vf);first=True
        # Keep peak in memory; only write one accepted peak state at the end.
        if abs(reaction)>peak:peak=abs(reaction);callback.peak=(step,x.copy(),vf.copy(),alive.copy())
        if step==steps:state('end',step,x,vf)
        prev=(step,x.copy(),fm,fd)
        if step%1000==0:stream.flush()
      except BaseException:
        (out/'CALLBACK_ERROR.txt').write_text(traceback.format_exc());raise
    L.set_fix_external_callback('audit',callback)
    try:
        L.command(f'run {steps}')
        p=callback.peak;np.savez_compressed(out/'peak.npz',step=p[0],time_s=p[0]*dt,x_m=p[1],v_full_m_s=p[2],alive=p[3])
        np.savez_compressed(out/'break_events.npz',events=np.asarray(events,dtype=np.int64).reshape(-1,3))
        save(out/'DIAGNOSTIC_CHECKS.json',dict(force_checks=forcechecks,negative_topology_releases=negative_releases,callbacks=calls,expected=steps+1,energy_definition='U=.5*(K-5G/3)*V*sum(theta^2)+7.5G*V^2*sum_alive(h*e^2/r0*(1/m_i+1/m_j)); theta recomputed on topology, m fixed initial; native contact PE*V',energy_release='U(x,old_topology)-U(x,new_topology); independent of work balance; native force uses pre-break theta on break step, residual retains lag effects',wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time(),mpi_ranks=1,gpu_compute=False))
        save(out/'EXECUTION.json',dict(status='COMPLETED',steps=steps,wall_seconds=time.perf_counter()-start))
    finally:stream.close();L.close()

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--case',required=True);a.add_argument('--dx',type=float,default=50e-6);a.add_argument('--dt',type=float,default=.5e-6);a.add_argument('--slow',type=float,default=1);a.add_argument('--static',action='store_true')
    z=a.parse_args();run(z.output.resolve(),z.case,z.dx,z.dt,z.slow,z.static)
