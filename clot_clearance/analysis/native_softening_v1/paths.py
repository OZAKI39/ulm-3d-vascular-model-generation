from native import *
import resource
assert all(json.loads((O/p).read_text())['status']=='PASS' for p in ['UNDAMAGED_NATIVE_PARITY.json','PARTIAL_DAMAGE_FORCE_ENERGY_PARITY.json','FIXED_DAMAGE_GRADIENT_AUDIT.json'])
start=time.perf_counter();summaries={};reference=np.load(PREV/'provenance/reference.npz');X=reference['X'];vec=X/DX
for name,sc in [('soft_sc018',.18),('soft_sc024',.24)]:
 t0=time.perf_counter();old=np.load(PREV/(name+'.npz'));ns=Native(name+'_path',sc);prevx=X.copy();W=0.;records=[];stored={k:[] for k in ['x','force','d','kappa','theta','scalars']};maxforce=0.;maxenergy=0.;maxd=0.;minincrement=0.;maxhealing=0.;dprev=np.zeros(12047);kprev=np.zeros(12047);idempotence=[]
 for it,x in enumerate(old['x']):
  dx=x-prevx;dw=0.
  for g in [.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)]:
   st=ns.eval(prevx+g*dx);dw-=.5*float(np.sum(st['force']*dx))
  ns.eval(x);st=ns.commit(it);U,Diss,drop,events,count,token=st['scalars'][:6];W+=dw;Q=-np.sum(st['force']*vec);res=W-U-Diss
  if it in [0,90,160,353]:
   before=st['scalars'].copy();repeat=ns.commit(it);idempotence.append(dict(frame=it,unchanged_ledger=bool(np.array_equal(before[:6],repeat['scalars'][:6])),unchanged_damage=bool(np.array_equal(st['d'],repeat['d']))));assert idempotence[-1]['unchanged_ledger'] and idempotence[-1]['unchanged_damage']
  maxforce=max(maxforce,float(abs(st['force']-old['force'][it]).max()));maxenergy=max(maxenergy,abs(U-old['series'][it,3]));maxd=max(maxd,float(abs(st['d']-old['d'][it]).max()));minincrement=min(minincrement,drop);maxhealing=max(maxhealing,float(np.max(dprev-st['d'])),float(np.max(kprev-st['kappa'])));
  records.append([old['series'][it,0],Q,U,W,Diss,res,drop,events,count,st['bonds']]);
  for k in stored:stored[k].append(st[k].copy())
  prevx=x;dprev=st['d'];kprev=st['kappa']
  if maxforce>2e-17 or maxenergy>1e-20 or maxd>2e-13 or minincrement < -1e-20 or maxhealing>2e-13:
   np.savez_compressed(O/(name+'_FAILURE.npz'),**st,frame=it);raise RuntimeError('native path stop gate; preserved actual state')
 arr=np.array(records);np.savez_compressed(O/(name+'_native.npz'),**{k:np.array(v) for k,v in stored.items()},series=arr,columns=['s','reaction_N','stored_J','external_work_J','damage_dissipation_J','balance_residual_J','signed_damage_increment_J','failed_events','accepted_commits','native_bonds'],phase=old['phase'],bond_ids=st['bond_ids'],reference_X=X,bonds=reference['bonds'])
 z=dict(status='PASS',states=len(arr),sc=sc,total_native_D_J=float(arr[-1,4]),offline_D_J=float(old['series'][-1,5]),native_offline_D_abs_J=float(abs(arr[-1,4]-old['series'][-1,5])),max_force_abs_N=maxforce,max_energy_abs_J=maxenergy,max_damage_abs=maxd,min_signed_damage_increment_J=minincrement,max_history_decrease=maxhealing,work_balance_max_relative=float(abs(arr[:,5]).max()/max(arr[:,3])),final_native_bonds=int(arr[-1,9]),failure_events=int(arr[-1,7]),accepted_commits=int(arr[-1,8]),final_residual_force_N=float(abs(st['force']).max()),final_U_J=float(arr[-1,2]),peak_generalized_reaction_N=float(arr[:,1].max()),wall_s=time.perf_counter()-t0,native_evaluations=int(st['scalars'][6]),native_internal_solver_wall_s=float(st['scalars'][7]),endpoint_star_solves_per_rank=int(st['scalars'][6])*341,repeat_commit_checks=idempotence)
 assert z['work_balance_max_relative']<2e-8 and z['final_native_bonds']==0 and z['failure_events']==12047 and z['final_residual_force_N']==0 and z['native_offline_D_abs_J']<1e-20
 summaries[name]=z;ns.close();dump(O/'NATIVE_LOADING_UNLOADING_RELOADING.json',dict(status='PASS' if len(summaries)==2 else 'IN_PROGRESS',records=summaries));print(name,z,flush=True)
dump(O/'IRREVERSIBLE_DAMAGE_AND_DISSIPATION.json',dict(status='PASS',records=summaries,raw_signs_retained=True,damage_commit_time_layer='fixed new coordinates after Gauss work with old damage'))
dump(O/'PATH_RESOURCES.json',dict(wall_s=time.perf_counter()-start,CPU_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,native_paths=2,old_instant_reference_reused=True))
