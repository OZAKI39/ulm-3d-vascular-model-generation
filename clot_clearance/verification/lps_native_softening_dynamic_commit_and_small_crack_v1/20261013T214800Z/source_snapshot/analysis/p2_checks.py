from native import *
import resource
start=time.perf_counter();ref=np.load(O/'provenance/reference.npz');X=ref['X'];mdl=SofteningLPS(ref);records=[]
for label,sc,stretch in [('intact',.18,.03),('partial018',.18,.164),('partial024',.24,.18),('complete',.18,.2)]:
 x=X*(1+stretch);k,d=mdl.update(x,np.zeros(mdl.nb),s_complete=sc);ns=Native('p2_'+label,sc=sc);st=ns.static(x,k,d);q=mdl.evaluate(x,d);ferr=float(abs(q['force']-st['material_force']).max());uerr=abs(q['U']-st['scalars'][0]);assert ferr<2e-17 and uerr<1e-20
 before=st['scalars'][:6].copy();ns.l.command('run 0');after=ns.state();assert np.array_equal(before,after['scalars'][:6]);records.append(dict(label=label,force_abs_N=ferr,energy_abs_J=uerr,run0_readonly=True,full_zero_force_N=float(abs(st['material_force']).max()) if label=='complete' else None));ns.close()
dump(O/'P2_FIXED_STATE_REGRESSION.json',dict(status='PASS',records=records,wall_s=time.perf_counter()-start,scope='Four frozen static states only; no preliminary dynamic loading family'))
print('P2_STATIC_PASS',records,flush=True)
