from checks import *
from spatial_model import SpatialLPS
from scipy.sparse.linalg import LinearOperator,eigsh
import resource
with (O/'preflight/STARTED.json').open('x') as f:json.dump(dict(purpose='one targeted new-geometry and mixed-DOF check; no fracture scan'),f)
start=time.perf_counter();cpu0=time.process_time();ids=np.load(O/'provenance/boundary_ids.npz');ref=np.load(O/'provenance/reference.npz');b=ref['bonds'];cut=ref['initial_cut'];free=ref['free_dof'];model=SpatialLPS(ref)
ns=Native('preflight/native',dt=BASE);ns.dynamics(ids['bottom_ids'],ids['top_ids'],ids['free_ids']);s=ns.state();first=s;nativecodes,special=topology(ns);expected=b[~cut]+1
assert np.array_equal(s['atom_ids'],np.arange(1,ns.N+1)) and np.array_equal(ns.field(9),np.arange(1,ns.NB+1))
assert np.array_equal(s['d'],cut.astype(float)) and np.array_equal(s['kappa'],SC*cut)
assert np.array_equal(ns.field(4),ref['m']) and ns.l.extract_atom('mass')[1]==RHO*V
assert np.array_equal(nativecodes,np.sort(expected[:,0]*ns.N+expected[:,1]))
assert s['scalars'][1]==s['scalars'][3]==s['scalars'][4]==0 and s['step']==0 and s['time']==0
assert np.array_equal(s['x'],ref['X']) and s['failed_tension_max_N']==0
assert abs(s['scalars'][0])<1e-20 and np.max(abs(s['force']))<2e-17
save_state(O/'preflight/initial.npz',s)
steps=[];elapsedsteps=time.perf_counter()
for number,h in enumerate([BASE,BASE/2,BASE]):
 prev=s;xp=predict_x(prev,h,ref,ids);ns.l.command(f'timestep {h:.17g}');ns.dt=h;ns.l.command('run 1 pre no post no');s=ns.state();e=model.evaluate(s['x'],s['d']);kk,dd=model.update(s['x'],prev['kappa'],s_complete=SC)
 err=dict(x_VV_m=float(abs(s['x']-xp).max()),v_VV_m_s=float(abs(s['velocity'][free]-prev['velocity'][free]-h/(2*RHO*V)*(s['force'][free]+prev['force'][free])).max()),force_N=float(abs(e['force']-s['material_force']).max()),U_J=float(abs(e['U']-s['scalars'][0])),star_energy_sum_J=float(abs(e['node_energy'].sum()-e['U'])),kappa=float(abs(kk-s['kappa']).max()),damage=float(abs(dd-s['d']).max()))
 assert err['x_VV_m']<2e-18 and err['v_VV_m_s']<2e-12 and err['force_N']<2e-17 and err['U_J']<1e-20 and err['star_energy_sum_J']<1e-20 and err['kappa']<2e-13 and err['damage']<2e-13,err
 assert np.array_equal(s['d'],cut.astype(float)) and s['scalars'][3]==0
 steps.append(dict(step=s['step'],time_s=s['time'],h_s=h,errors=err))
 if number==1:
  cp=O/'preflight/mixed_dof.restart';ns.l.command(f'write_restart {cp}');tc,ec=topology(ns);ns.close();ns=Native('preflight/restart_readback',restart=cp,dt=h);ns.dynamics(ids['bottom_ids'],ids['top_ids'],ids['free_ids']);restored=ns.state();rest=compare(s,restored);tc2,ec2=topology(ns);assert rest['status']=='PASS' and np.array_equal(tc,tc2) and ec==ec2;s=restored
actual_step_wall=time.perf_counter()-elapsedsteps;assert np.max(abs(s['velocity'][ids['free_ids']-1]))>0 and np.max(abs(s['velocity'][np.r_[ids['bottom_ids'],ids['top_ids']]-1,:2]))>0
save_state(O/'preflight/last.npz',s);ns.close()
# Initial precracked fixed-damage tangent; exact linear extension operator, no physical path.
di=(ref['X'][b[:,1]]-ref['X'][b[:,0]])/ref['r0'][:,None];d0=cut.astype(float);dd=d0[model.edge];q=1-dd;den=model.Ddir*dd+model.krefdir*q;w=model.krefdir*q/den;hh=dd/den;denom=1/model.A+np.bincount(model.center,weights=model.l**2*hh,minlength=model.n)
def tangent(y):
 u=np.zeros(ref['X'].shape);u[free]=y;edge=np.sum((u[b[:,1]]-u[b[:,0]])*di,axis=1);ee=edge[model.edge];pp=np.bincount(model.center,weights=model.l*w*ee,minlength=model.n)/denom;td=w*(model.Ddir*ee+model.l*pp[model.center]);ff=model.assemble(td[:model.nb]+td[model.nb:],di)
 return -ff[free]/(RHO*V)
n=int(free.sum());op=LinearOperator((n,n),matvec=tangent,dtype=float);v0=np.sin(np.arange(n)+.123);eigval,eigvec=eigsh(op,k=1,which='LA',tol=1e-6,maxiter=400,v0=v0);lam=float(eigval[0]);resid=float(np.linalg.norm(tangent(eigvec[:,0])-lam*eigvec[:,0])/lam);critical=2/np.sqrt(lam);assert BASE<critical/4
audit=dict(status='PASS',nodes=len(ref['X']),original_bonds=len(b),initial_cut_count=int(cut.sum()),initial_live_bonds=int((~cut).sum()),original_ids_m_mass_exact=True,initial_cut_tension_zero=True,initial_natural_event_count=0,initial_natural_D_J=0,initial_force_max_N=float(abs(first['force']).max()),initial_U_J=float(first['scalars'][0]),diagnostic_physical_steps=3,diagnostic_final_time_s=float(s['time']),diagnostic_has_no_natural_damage=True,steps=steps,mixed_DOF_restart_comparison=rest,restart_native_topology_and_special_exact=True,free_internal_nodes_moved=True,grip_transverse_DOFs_moved=True,no_DOF_integrated_twice=True,initial_tangent_lambda_max_s2=lam,eigensolver_relative_residual=resid,initial_tangent_Verlet_limit_s=critical,base_dt_over_initial_limit=BASE/critical,stability_scope='Initial precracked tangent only; not a nonlinear global time convergence certificate. Runtime damage/energy/velocity gates stay active.',three_step_check_wall_including_restart_s=actual_step_wall,wall_s=time.perf_counter()-start,CPU_s=time.process_time()-cpu0,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
dump(O/'MINIMAL_IMPLEMENTATION_CHECK.json',audit);geo=json.loads((O/'GEOMETRY_AND_INITIAL_CUT_AUDIT.json').read_text());geo['status']='PASS';dump(O/'GEOMETRY_AND_INITIAL_CUT_AUDIT.json',geo)
print('P1_PASS',json.dumps({k:audit[k] for k in ['nodes','original_bonds','initial_cut_count','initial_tangent_Verlet_limit_s','base_dt_over_initial_limit','three_step_check_wall_including_restart_s','wall_s']}))
