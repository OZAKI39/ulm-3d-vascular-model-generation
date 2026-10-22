"""Frozen-state decomposition using original force and resistance implementations."""
import time
import numpy as np
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.nearfield import pair_reference
from .common import *

def run(c):
 start=time.perf_counter();rows=[]
 for initial,x in c.states.items():
  c.x0=x.copy();sim=c.engine('coarse');op=resistance(x,R);zeta=sim.zeta
  for action in ['OFF','XM','XP','PAIR','YM','YP']:
   ev=c.schedules[action].segments[0].event();sim.control(ev)
   try:
    f=sim.fields(x);pb,gb=c.bases['coarse'].query_basis(x);u=c.bases['coarse'].vector(ev['amplitude_Pa'],ev['phase_rad'],ev['enabled']);p=u@pb;g=np.einsum('i,inj->nj',u,gb);v=c.response.query(p);ref1=-.5*np.real(v[:,0,None]*g.conj());ref2=pair_reference(x,v,p,c.physics,sim.sim.chi_zero)
    finite=[];h=1e-8
    for axis in range(3):
     e=np.eye(3)[axis]*h;pp,_=c.bases['coarse'].query_basis(x+e);pm,_=c.bases['coarse'].query_basis(x-e);finite.append((u@pp-u@pm)/(2*h))
    dg=np.stack(finite,-1);errors={key:float(np.max(abs(f[key]-ref))) for key,ref in [('pressure_Pa',p),('primary_force_N',ref1),('secondary_force_N',ref2)]}
    valid=all(errors[k]<=max(np.max(abs(f[k]))*1e-9,1e-24) for k in errors)
    vg=np.linalg.solve(op['A'],f['total_force_N'].ravel()/zeta).reshape(32,3);v1=np.linalg.solve(op['A'],f['primary_force_N'].ravel()/zeta).reshape(32,3);v2=np.linalg.solve(op['A'],f['secondary_force_N'].ravel()/zeta).reshape(32,3)
    rows.append(dict(initial=initial,action=action,status='PASS' if valid else 'FAIL',position_m=x,pressure_Pa=f['pressure_Pa'],gradient_Pa_m=g,**{k:v for k,v in f.items() if k!='pressure_Pa'},unconstrained_velocity_m_s=vg,primary_only_frozen_velocity_m_s=v1,secondary_only_frozen_velocity_m_s=v2,net_secondary_N=f['secondary_force_N'].sum(0),net_primary_N=f['primary_force_N'].sum(0),GPU_CPU_max_errors=errors,gradient_finite_difference_relative=float(np.linalg.norm(dg-g)/max(np.linalg.norm(g),1e-30)),wall_check=sim.wall_check(x),scope='Frozen force/mobility decomposition only; not an ablated physical trajectory. Y proxy-specific.'))
   except Exception as exc:rows.append(dict(initial=initial,action=action,status='FAIL',reason=repr(exc)));sim.sim.error.zero_()
 write(c.out/'FORCE_DIAGNOSTICS.json',dict(rows=rows,wall_s=time.perf_counter()-start,identity=c.ident))
 return rows
