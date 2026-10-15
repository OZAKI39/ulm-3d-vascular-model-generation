from common import *
from model import SofteningLPS
import time,resource,ctypes,sys
from analysis.timelayer_v1.model import LPS
start=time.perf_counter();ref=np.load(O/'provenance/reference.npz');m=SofteningLPS(ref);old=LPS(ref);X=m.X;alive=np.ones(m.nb,bool)
rng=np.random.default_rng(20261013)
states={'uniform_small':X@np.diag([1.002,1.001,1.0005]),'shear_nonaffine':X@np.array([[1,.002,0],[0,1,.001],[0,0,1.]]).T+rng.normal(size=X.shape)*DX*.0002}
records=[]
# Load ONLY the existing final double-fix build (scoped dynamic linker path).
from lammps import lammps
l=lammps(cmdargs=['-log',str(O/'native_baseline.log'),'-screen','none'])
maps=sorted(set(line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in line))
assert maps and all(sha(p)==EXPECTED for p in maps)
l.commands_string(f'''units si
dimension 3
boundary f f f
atom_style bond
atom_modify map array sort 0 0
newton on off
read_data {O/'provenance/body.data'} extra/bond/per/atom 150 extra/special/per/atom 200
lattice sc {DX:.17g}
neighbor {DX:.17g} bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*DELTA:.17g}
fix props all property/atom d_vfrac d_theta ghost yes
set group all d_vfrac {m.V:.17g}
pair_style zero {DELTA:.17g}
pair_coeff * *
bond_style bpm/peri
bond_coeff 1 lps {K:.17g} {G} {DELTA} .15 0
special_bonds lj 0 1 1 coul 1 1 1
thermo 1
thermo_style custom step atoms bonds pe
run 0
''')
for name,x in states.items():
 ids=l.numpy.extract_atom('id')[:m.n].copy();l.numpy.extract_atom('x')[:m.n]=x[ids-1];l.command('run 0')
 ids=l.numpy.extract_atom('id')[:m.n].copy();fn=l.numpy.extract_atom('f')[:m.n].copy()[np.argsort(ids)]*m.V
 th,fo,uo=old.eval(x,alive);a=m.evaluate(x,np.zeros(m.nb));scale=float(np.linalg.norm(fo));q=dict(name=name,U_J=a['U'],native_force_relative_L2=float(np.linalg.norm(fn-a['force'])/scale),old_force_relative_L2=float(np.linalg.norm(fo-a['force'])/scale),old_energy_relative=float(abs(uo-a['U'])/uo),theta_max_abs=float(abs(th-a['theta']).max()),net_force_N=np.sum(a['force'],0),net_torque_Nm=np.cross(x,a['force']).sum(0),native_PE_times_V_diagnostic_only_J=l.get_thermo('pe')*m.V,edges_remaining=int(l.get_thermo('bonds')))
 records.append(q);np.savez_compressed(O/f'baseline_{name}.npz',x=x,force_native_N=fn,force_candidate_N=a['force'],force_existing_N=fo,theta_candidate=a['theta'],theta_existing=th,U_candidate_J=a['U'],U_existing_J=uo)
l.close()
# Recover historical native post-break force under its unchanged binary topology.
z=np.load(O/'provenance/postbreak_same_coordinates.npz');a=m.baseline(z['x'],z['alive']);historical=dict(relative_force_L2=float(np.linalg.norm(a['force']-z['force_N'])/np.linalg.norm(z['force_N'])),theta_max_abs=float(abs(a['theta']-z['theta']).max()),broken_bonds=int((~z['alive']).sum()))
passed=all(z['native_force_relative_L2']<2e-9 and z['old_force_relative_L2']<2e-12 and z['old_energy_relative']<2e-12 and z['edges_remaining']==m.nb for z in records) and historical['relative_force_L2']<2e-9
baseline=dict(status='PASS' if passed else 'FAIL',records=records,historical_postbreak=historical,loaded_library_sha256={p:sha(p) for p in maps},positive_energy_condition=dict(K_minus_5G_over3_Pa=K-5*G/3,min_diagonal_N_per_m=float(m.D.min()),m_fixed=True),native_PE_not_energy_reference=True)
dump(O/'BASELINE_LPS_ENERGY_FORCE_CONSISTENCY.json',baseline);dump(O/'INTACT_ELASTIC_LIMIT_AUDIT.json',baseline)
assert passed,'baseline mismatch: stop'
x=X@np.array([[1.011,.003,0],[0,.994,.002],[0,0,1.003]]).T
# Diagnostic-only damage, held fixed for BOTH signs of each coordinate perturbation.
d=.1+.65*rng.random(m.nb);a=m.evaluate(x,d);dofs=np.r_[np.argsort(abs(a['force']).ravel())[-8:],rng.choice(x.size,8,replace=False)]
rows=[];Fscale=float(abs(a['force']).max())
for h in [DX*1e-5,DX*2e-6]:
 for dof in dofs:
  xp=x.copy();xm=x.copy();xp.flat[dof]+=h;xm.flat[dof]-=h
  fd=-(m.evaluate(xp,d)['U']-m.evaluate(xm,d)['U'])/(2*h);f=float(a['force'].flat[dof]);rows.append(dict(atom_id=int(dof//3+1),axis=int(dof%3),h_m=h,FD_N=fd,force_N=f,absolute_error_N=abs(fd-f),relative_to_global_scale=abs(fd-f)/Fscale,relative_to_component=abs(fd-f)/max(abs(f),1e-30)))
force_error=max(z['relative_to_global_scale'] for z in rows)
net=np.linalg.norm(a['force'].sum(0));torque=np.linalg.norm(np.cross(x,a['force']).sum(0));force_norm=float(np.linalg.norm(a['force']));length=float(np.linalg.norm(x,axis=1).max())
# Damage derivative: genuine signed central difference and signed finite drop.
ychecks=[]
for k in rng.choice(m.nb,8,replace=False):
 dp=d.copy();dm=d.copy();dp[k]+=1e-5;dm[k]-=1e-5
 yn=-(m.evaluate(x,dp)['U']-m.evaluate(x,dm)['U'])/2e-5
 ychecks.append(dict(edge=int(k),analytic_Y_J=float(a['Y'][k]),FD_Y_J=yn,relative_error=abs(yn-a['Y'][k])/max(a['Y'][k],1e-30)))
monotone=[]
for label,xx in [('mixed',x),('compression',X*.98),('tension',X*1.17)]:
 dnew=np.minimum(d+.08,1);u0=m.evaluate(xx,d);u1=m.evaluate(xx,dnew);monotone.append(dict(state=label,U_before_J=u0['U'],U_after_J=u1['U'],signed_dissipation_J=u0['U']-u1['U'],minimum_Y_J=float(u0['Y'].min())))
# Deliberate complete cut with nonuniform in-half strain: all cross-cut edge forces zero.
cross=(X[m.i,0]<-30e-6)!=(X[m.j,0]<-30e-6);dc=np.zeros(m.nb);dc[cross]=1;ac=m.evaluate(x,dc)
sep=dict(cross_bonds=int(cross.sum()),max_cross_force_N=float(abs(ac['t'][cross]).max()),net_half_force_N=np.sum(ac['force'][X[:,0]<-30e-6],0),internal_force_norm_N=float(np.linalg.norm(ac['force'])))
mask=X[:,0]<-30e-6;xx=X.copy();xx[mask,1]*=1.002
fa=m.evaluate(xx,dc)['force'];fb=m.evaluate(X,dc)['force'];leak=float(np.linalg.norm(fa[~mask]-fb[~mask]))
sep['disconnected_component_force_leak_N']=leak
checks=dict(disconnected_component_locality=leak<1e-22,gradient=force_error<2e-6,force_balance=net/(force_norm+1e-30)<2e-13,moment_balance=torque/(force_norm*length+1e-30)<2e-13,damage_gradient=max(z['relative_error'] for z in ychecks)<2e-5,monotone_energy=all(z['signed_dissipation_J']>=0 for z in monotone),zero_cross_force=sep['max_cross_force_N']==0,positive_energy=a['U']>=0,dual_energy_identity=abs(a['dual_energy_error_J'])<1e-12*a['U'])
dump(O/'FIXED_DAMAGE_FORCE_GRADIENT_AUDIT.json',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,force_scale_N=Fscale,rows=rows,max_error_record=max(rows,key=lambda z:z['absolute_error_N']),force_balance_N=net,moment_balance_Nm=torque,damage_derivative=ychecks,signed_damage_increments=monotone,full_cut_diagnostic=sep,damage_fixed_between_FD_signs=True))
np.savez_compressed(O/'fixed_damage_diagnostic.npz',x=x,d=d,force=a['force'],Y=a['Y'],theta=a['theta'],v=a['v'],cut_damage=dc,cut_force=ac['force'],cut_edge_force=ac['t'])
dump(O/'CHECK_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,native_instances=1,native_intact_states=2,production_modifications=0))
assert all(checks.values()),checks
print(json.dumps(dict(baseline='PASS',gradient='PASS',max_relative_error=force_error,checks={k:bool(v) for k,v in checks.items()})),flush=True)
