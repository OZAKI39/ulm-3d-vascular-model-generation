"""Only new multi-source force/work accounting and real-topology objective checks."""
from common import *
from native import native
from distributed_multisite_loading import *
from multifragment_objective import metrics,shape
import time

start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);results={}
for protocol in ['simultaneous','staggered']:
    ns,C=native(O/'checks'/protocol,protocol);initial=ns.dynamics();contact=C(ref)
    assert not initial['indenter_force'].any() and not initial['site_energies'].any() and not initial['site_reactions'].any()
    x=ref['X'].copy();ids=[]
    for site in ns.tool.sites:
        idx=int(np.argmin(np.linalg.norm(x[:n]-site.config['surface_point_m'],axis=1)));ids.append(idx)
        x[idx]=site.c0+site.direction*(site.radius-.1*DX)
    both=ns.static(x);uc,fc,_=contact.evaluate(x,both['d']);actual=both['force']-both['material_force']-fc
    force_error=float(np.max(abs(actual-both['site_forces'].sum(0))))
    reaction_error=float(np.max(abs(both['site_reactions']+both['site_forces'].sum(1))))
    energy_error=float(np.max(abs(both['site_energies']-both['site_energy_oracle'])))
    assert force_error<1e-12 and reaction_error<1e-15 and energy_error<1e-20
    assert np.all(np.sum(both['site_forces'][:,ids,:].sum(1)*np.array([q.direction for q in ns.tool.sites]),axis=1)>0)
    independent=[];gradient_errors=[]
    for j,site in enumerate(ns.tool.sites):
        isolated=ref['X'].copy();isolated[ids[j]]=x[ids[j]];state=ns.static(isolated)
        assert state['site_energies'][1-j]==0 and not state['site_reactions'][1-j].any()
        independent.append(float(state['site_energies'][j]));epsilon=1e-9;energies=[]
        for sign in [-1,1]:
            shifted=x.copy();shifted[ids[j]]+=sign*epsilon*site.direction;energies.append(float(ns.static(shifted)['indenter_U']))
        finite_difference=-(energies[1]-energies[0])/(2*epsilon);analytic=float(both['site_forces'][j,ids[j]]@site.direction);gradient_errors.append(abs(finite_difference-analytic))
    assert abs(sum(independent)-both['indenter_U'])<1e-20 and max(gradient_errors)<1e-12
    released=ns.static(ref['X']);assert not released['indenter_force'].any() and not released['site_energies'].any()
    work=[];du=[]
    for site in ns.tool.sites:
        t0=site.delay+APPROACH;t1=t0+.0001
        work.append(site.work(ref['X'],ref['X'],t0,t1,32));du.append(site.evaluate(ref['X'],t1)['U']-site.evaluate(ref['X'],t0)['U'])
    work=np.array(work);du=np.array(du)
    assert np.max(abs(work-du))<max(1e-18,float(abs(du).max())*1e-8)
    assert np.array_equal(ns.field(13),S0[ref['bond_types']]) and np.array_equal(ns.field(3),ref['initial_cut'].astype(float))
    assert ns.state()['native_bonds']==int(np.sum(~ref['initial_cut']))
    results[protocol]=dict(status='PASS',native_contact_force_error_N=force_error,per_site_reaction_error_N=reaction_error,per_site_energy_error_J=energy_error,contact_gradient_errors_N=gradient_errors,stationary_material_per_site_work_error_J=abs(work-du),independent_site_additivity=True,unloaded_and_release_zero=True,uniform_onsets_and_original_precuts_unchanged=True)
    ns.close()
graph=module('objective_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py');actual={}
for name in ['uniform','heterogeneous']:
    state=load_npz(HISTORY/'cases'/name/'keys/final.npz');labels,sizes,attached=graph.topology(ref,state['d']);free=sizes[~attached];q=metrics(free,n)
    assert free.sum()==(18 if name=='uniform' else 14) and q['fragment_count_target_band']==0 and q['fine_debris_mass_fraction_detached']==1
    for label in np.flatnonzero(~attached):
        ids=np.flatnonzero(labels==label);boundary,_=graph.boundary(ref,ids);assert not np.any(state['d'][boundary]<1);m=shape(state['x'][ids]);assert m['shape_class']=='UNRESOLVED_FINE_DEBRIS'
    actual[name]=q
dump(O/'KEY_CHECKS.json',dict(status='PASS',multisite=results,historical_real_topology_objective=actual,wall_s=time.perf_counter()-start,library_sha256=sha(LIB),scope='No constitutive/MPI/restart change or repeat matrix; new contact composition and actual-component metrics only'))
print('KEY_MULTISITE_CHECKS_PASS',results)
