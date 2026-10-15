from common import *
import ctypes as ct,numpy as np,time,resource,os
from lammps import lammps
q=O/'static_preflight_verified';q.mkdir(exist_ok=False);t=time.perf_counter()
ref=np.load(PRIOR/'final_sort/reference.npz');priorstate=np.load(PRIOR/'final_sort/T0.npz');X=ref['X'];b=ref['bonds'];m=ref['m'];V=float(ref['V'][0]);n=len(X)
L=lammps(cmdargs=['-log',str(q/'log.lammps'),'-screen','none']);ident=loaded(L);dump(q/'LOADED_LIBRARY.json',ident)
probe=ct.CDLL(str(O/'observer.so'));probe.observe.argtypes=[ct.c_void_p,ct.c_int];probe.observe.restype=ct.POINTER(ct.c_double);cbtype=ct.CFUNCTYPE(None,ct.c_longlong,ct.c_int);calls=[]
@cbtype
def cb(step,phase):calls.append([step,phase,float(L.numpy.extract_atom('v')[:L.extract_global('nlocal')].max())])
probe.register_observer.argtypes=[cbtype];probe.register_observer(cb)
deck=(PRIOR/'final_sort/in.lammps').read_text();deck=deck[:deck.index('run 0')];deck=deck.replace('sort 2 5.0000000000000002e-05','sort 0 0') # static only, not V3
L.commands_string(deck);L.command('fix audit all audit/readonly/finalv3');L.command('run 0')
ids=L.numpy.extract_atom('id')[:n].copy();L.numpy.extract_atom('x')[:n]=priorstate['x'][ids-1];L.command('run 1 pre no post no')
ids=L.numpy.extract_atom('id')[:n].copy();order=np.argsort(ids);x=L.numpy.extract_atom('x')[:n].copy()[order];f=L.numpy.extract_atom('f')[:n].copy()[order]*V;nm=np.ctypeslib.as_array(probe.observe(L.lmp,0),shape=(n,)).copy()[order]
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
energy=Energy({'X':X,'bonds':b},50e-6,R/'analysis/material_macro_audit_v1/energy.so');U,Fref=energy.eval(x,np.ones(len(b),np.uint8),True)
eps=50e-6*1e-5;xp=x.copy();xm=x.copy();xp[0,0]+=eps;xm[0,0]-=eps;fd=-(energy.eval(xp,np.ones(len(b),np.uint8),False)-energy.eval(xm,np.ones(len(b),np.uint8),False))/(2*eps)
errF=float(np.linalg.norm(f-Fref,axis=1).max());scale=float(np.linalg.norm(Fref,axis=1).max());errm=float(np.max(abs(nm-m)/m));PE=L.get_thermo('pe')*V
res=dict(status='PASS',loaded=ident,m_max_relative=errm,m_max_absolute=float(abs(nm-m).max()),force_error_N=errF,force_scale_N=scale,force_relative=errF/scale,independent_energy_J=float(U),native_PE_times_V_J=PE,native_PE_relative_mismatch=(PE-U)/U,energy_gradient_FD_N=float(fd),reference_Fx_atom1_N=float(Fref[0,0]),FD_relative_to_global_force_scale=abs(fd-Fref[0,0])/scale,observer_calls=calls,static_steps=1,no_integrator=True,wall_s=time.perf_counter()-t,production_sha256=sha(PROD))
assert errm<=2e-12 and errF/scale<1e-10 and abs(fd-Fref[0,0])/scale<1e-7 and len(calls)==2
np.savez_compressed(q/'static_state.npz',X=X,x=x,bonds=b,force_N=f,reference_force_N=Fref,m_native=nm,m_reference=m)
L.close();dump(O/'STATIC_PREFLIGHT.json',res);print(json.dumps(res))
