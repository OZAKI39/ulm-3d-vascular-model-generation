from common import *
from native import native,predict
from indenter import *
import ctypes,time
start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);model=energy_model(ref);results={}
u,C=native(O/'checks/uniform');h,_=native(O/'checks/heterogeneous',True);expected=load_npz(O/'HETEROGENEITY_FIELD.npz')['s0_by_original_bond_id'];assert np.array_equal(h.field(13),expected)
x=ref['X'].copy();x[:n,0]+=.01*(x[:n,2]-x[:n,2].min());su=u.static(x);sh=h.static(x)
assert np.array_equal(su['force'],sh['force']) and su['scalars'][0]==sh['scalars'][0]
results['initial_elastic_parity']=dict(force_max_diff_N=float(np.max(abs(su['force']-sh['force']))),energy_diff_J=float(su['scalars'][0]-sh['scalars'][0]),onset_original_id_mapping_exact=True)
# A single material consistency fixture, not a new physical loading campaign.
x=ref['X']*1.16;h.static(x);assert h.lib.lpsvar_action(h.l.lmp,2,1,None,None)==0;sh=h.static();oracle=model.evaluate(x,sh['d']);assert np.max(abs(oracle['force']-sh['material_force']))<1e-12 and abs(oracle['U']-sh['scalars'][0])<1e-17
kk,dd=predict(model,ref,x,SC[ref['bond_types']]*ref['initial_cut'],expected);assert np.max(abs(sh['d']-dd))<1e-12
direction=np.sin(np.arange(x.size)).reshape(x.shape);direction/=np.linalg.norm(direction);eps=1e-9;up=model.evaluate(x+eps*direction,sh['d'])['U'];um=model.evaluate(x-eps*direction,sh['d'])['U'];fd=-(up-um)/(2*eps);analytic=float(np.sum(sh['material_force']*direction));assert abs(fd-analytic)<1e-10
# Sorting and restart directly cover the only newly stored field. No MPI sweep.
before_tags=h.l.numpy.extract_atom('id')[:h.N].copy();h.l.commands_string('atom_modify sort 1 0.0001\nrun 1');after_tags=h.l.numpy.extract_atom('id')[:h.N].copy();ss=h.state();assert np.array_equal(ss['d'],sh['d']) and np.array_equal(h.field(13),expected)
restart=O/'checks/onset.restart';h.l.command(f'write_restart {restart}');h.close();v,_=native(O/'checks/restart',True,restart);sr=v.static()
for key in ['x','kappa','d']:assert np.array_equal(ss[key],sr[key]),key
assert np.array_equal(v.field(13),expected) and ss['scalars'][1]==sr['scalars'][1]
prior=sr['scalars'][1];v.static(ref['X']);assert v.lib.lpsvar_action(v.l.lmp,2,2,None,None)==0;unloaded=v.static();assert np.all(unloaded['d']>=sr['d']) and unloaded['scalars'][1]==prior
results['heterogeneous_material']=dict(status='PASS',force_energy_gradient_error_N=abs(fd-analytic),oracle_force_max_N=float(np.max(abs(oracle['force']-sh['material_force']))),sorting_changed_atom_slots=int(np.sum(before_tags!=after_tags)),sorting_and_restart_local_parameters_exact=True,restart_history_exact=True,unloading_no_double_dissipation=True,new_damage_bonds=int(np.sum((sh['d']>0)&~ref['initial_cut'])))
u.close();v.close()
i,C=native(O/'checks/indenter');s0=i.dynamics();contact=C(ref);assert not s0['indenter_force'].any() and s0['indenter_U']==0
x=ref['X'].copy();idx=int(np.argmin(np.linalg.norm(x[:n,:2]-i.tool.c0[:2],axis=1)+abs(x[:n,2]-x[:n,2].max())));x[idx]=i.tool.c0-[0,0,i.tool.radius-.1*DX];s=i.static(x);uc,fc,_=contact.evaluate(x,s['d']);native_force=s['force']-s['material_force']-fc
assert np.max(abs(native_force-s['indenter_force']))<1e-12;assert np.max(abs(s['indenter_reaction']+s['indenter_force'].sum(0)))<1e-15;assert abs(s['indenter_U']-s['indenter_U_oracle'])<1e-20;assert s['indenter_force'][idx,2]<0
energies=[];eps=1e-9
for sign in [-1,1]:
 xx=x.copy();xx[idx,2]+=sign*eps;energies.append(float(i.static(xx)['indenter_U']))
fd=-(energies[1]-energies[0])/(2*eps);assert abs(fd-s['indenter_force'][idx,2])<1e-12
released=i.static(ref['X']);assert not released['indenter_force'].any() and not released['indenter_reaction'].any() and released['indenter_U']==0
wa=i.tool.work(ref['X'],ref['X'],APPROACH,PEAK,32);du=i.tool.evaluate(ref['X'],PEAK)['U']-i.tool.evaluate(ref['X'],APPROACH)['U'];assert abs(wa-du)<max(1e-14,abs(du)*1e-8)
results['indenter']=dict(status='PASS',contact_force_native_oracle_error_N=float(np.max(abs(native_force-s['indenter_force']))),action_reaction_error_N=float(np.max(abs(s['indenter_reaction']+s['indenter_force'].sum(0)))),energy_gradient_error_N=abs(fd-s['indenter_force'][idx,2]),no_contact_and_release_zero=True,stationary_material_moving_tool_work_error_J=abs(wa-du),native_style='fix indent sphere after physical SI conversion')
i.close();dump(O/'KEY_CONSISTENCY_CHECKS.json',dict(status='PASS',checks=results,wall_s=time.perf_counter()-start,library_sha256=sha(LIB),scope='Only new indenter and original-ID onset storage; one combined sorting/restart fixture; no historical benchmark matrix'));print('KEY_CONSISTENCY_PASS',results,flush=True)
