from common import *
from bench_native import native_factory,boundary,update,BASE,RAMP,END,AMPLITUDE,SOURCE
from graph_ops import topology,boundary as edge_boundary
from scipy.spatial import cKDTree
import time,resource,csv,traceback

def prepare():
 root=O/'benchmark';freeze=root/'DESIGN_FREEZE.json';assert not freeze.exists(),'Frozen experiment exists; inspect results, do not repeat automatically'
 counts=(14,8,16);Xc=np.stack(np.meshgrid(*[(np.arange(n)-(n-1)/2)*DX for n in counts],indexing='ij'),axis=-1).reshape(-1,3)
 Xa=np.stack(np.meshgrid(np.unique(Xc[:,0]),np.unique(Xc[:,1]),Xc[:,2].min()-DX*np.arange(1,4),indexing='ij'),axis=-1).reshape(-1,3);X=np.r_[Xc,Xa];n=len(Xc)
 b=cKDTree(X).query_pairs(HORIZON,output_type='ndarray');b.sort(axis=1);b=b[np.lexsort((b[:,1],b[:,0]))];r=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1);h=np.minimum(1,(HORIZON+DX/2-r)/DX);m=np.bincount(b.ravel(),weights=np.repeat(h*r*V,2),minlength=len(X));bt=np.where(b[:,1]<n,1,np.where(b[:,0]<n,2,3))
 cross=(X[b[:,0],2]*X[b[:,1],2]<0)&(bt==1);ci=np.flatnonzero(cross);frac=-X[b[ci,0],2]/(X[b[ci,1],2]-X[b[ci,0],2]);hit=X[b[ci,0]]+frac[:,None]*(X[b[ci,1]]-X[b[ci,0]]);cut=np.zeros(len(b),bool);cut[ci]=hit[:,0]<=-DX+1e-15;top=np.flatnonzero(Xc[:,2]==Xc[:,2].max())+1
 ref=dict(X=X,bonds=b,bond_types=bt,r0=r,taper=h,m=m,V=np.full(len(X),V),clot_count=np.array(n),initial_cut=cut,top_ids=top,cross=cross,cross_ids=ci,cross_hits=hit)
 for name in ['opening_z','shear_x']:
  case=root/name;(case/'provenance').mkdir(parents=True);(case/'keys').mkdir();np.savez_compressed(case/'provenance/reference.npz',**ref)
  text=['LPSVAR_REFERENCE_V1',f'{len(X)} {len(b)} {V:.17g}']+[f'{i+1} {v:.17g}' for i,v in enumerate(m)]+[f'{q+1} {i+1} {j+1} {bt[q]} {r[q]:.17g} {h[q]:.17g} 0 0 0' for q,(i,j) in enumerate(b)]
  (case/'provenance/reference.native').write_text('\n'.join(text)+'\n')
  text=['Synthetic 3D single-edge-precracked adhesive specimen','',f'{len(X)} atoms',f'{len(b)} bonds','','2 atom types','3 bond types','']+['-0.005 0.005 '+a+'lo '+a+'hi' for a in 'xyz']+['','Masses','',f'1 {RHO*V:.17g}',f'2 {RHO*V:.17g}','','Atoms # bond','']+[f'{i+1} 1 {1 if i<n else 2} '+' '.join(f'{v:.17g}' for v in x) for i,x in enumerate(X)]+['','Bonds','']+[f'{q+1} {bt[q]} {i+1} {j+1}' for q,(i,j) in enumerate(b)]
  (case/'provenance/body.data').write_text('\n'.join(text)+'\n')
 dump(freeze,dict(status='FROZEN_BEFORE_NATIVE_RUN',hypothesis='With identical native axial-history softening, the direction of applied displacement relative to a precrack changes notch localization, cohesive/adhesive competition and retained internal connectivity. This tests geometry/loading effects before introducing mixed-mode failure.',cases=['opening_z','shear_x'],only_changed_factor='top prescribed displacement direction z versus x',material=dict(G_Pa=G,K_Pa=K,rho=RHO,dx_m=DX,horizon_m=HORIZON,s0=S0,sc=SC,interface_uncalibrated=True),geometry=dict(material_counts=counts,material_nodes=n,anchors=len(Xa),original_bonds=len(b),typed_bonds=np.bincount(bt)[1:],precut_bonds=int(cut.sum()),notch='reference z=0 segment intersections with x<=-50 um; initial geometry only',m='immutable full original star before precut',dimensions_m=np.array(counts)*DX),loading=dict(displacement_m=AMPLITUDE,ramp_s=RAMP,end_s=END,waveform='half cosine; same 240 um net separation used by the prior +/-120 um small specimen',release='At 2 ms remove displacement constraint on top and integrate its existing x,v. No material bond is cut at release. Grip-connected components before release are not free fragments.',damping='none in either case',flow='none; synthetic isolated mechanical benchmark'),integration=dict(base_dt_s=BASE,min_dt_s=BASE/8,max_dd=.1,method='Exact VV coordinate forecast with typed damage; halve BEFORE native commit, verify forecast against accepted x/d. No rejected native state is reused.',energy_relative_gate=.01),limitations='Finite 3D specimen, thickness=2.66 delta, ligament=2.66 delta. Not continuum-converged, not calibrated Gc, not a real V3 rerun. Only two pre-frozen directions; no amplitude or parameter search.',library_sha256=sha(LIB),frozen_native_bridge_sha256=sha(SOURCE),driver_sha256={p.name:sha(p) for p in A.glob('*.py')},expected_cost='Two ~10000-step cases, measured wall/RSS retained; no arbitrary RSS or wall-clock stop.'))
 return ref

def run_case(name,axis):
 start=time.perf_counter();case=O/'benchmark'/name;Native,Contact=native_factory(case,axis);ns=Native();ref=ns.ref;n=ns.n;b=ref['bonds'];bt=ref['bond_types'];top=ref['top_ids']-1;cut=ref['initial_cut'];model=energy_model(ref);contact=Contact(ref);s=ns.dynamics();initial=s;prev=s;W=0.;Jc=0.;scale=1e-15;err={};rows=[];events=[];substeps=0;last=model.evaluate(s['x'],s['d']);prevfc=contact.evaluate(s['x'],s['d'])[1];defects=np.zeros(4);stop=[];maxfree=0;maxkey=None;first={};lasttop=None
 def error(k,v):err[k]=max(err.get(k,0.),float(v))
 def save(s,reason):np.savez_compressed(case/'keys'/f"step_{s['step']:06d}.npz",**s);events.append(dict(kind='KEY',step=s['step'],time_s=s['time'],reason=reason))
 assert s['time']==0 and s['scalars'][1]==0 and np.array_equal(s['d'],cut.astype(float));save(s,'initial')
 stream=(case/'history.csv').open('w');writer=None;changes=(case/'damage_events.jsonl').open('w')
 try:
  while prev['time']<END-1e-15:
   if not ns.released and prev['time']>=RAMP-1e-15:
    ns.release();prev=ns.state();events.append(dict(kind='GRIP_RELEASE',step=prev['step'],time_s=prev['time']));save(prev,'release')
   mask=np.zeros(ref['X'].shape,bool);mask[:n]=True
   if not ns.released:mask[top,axis]=False
   h=min(BASE,END-prev['time'],(RAMP-prev['time']) if not ns.released else END-prev['time']);depth=0
   while True:
    xp=prev['x'].copy();xp[mask]+=h*prev['velocity'][mask]+h*h/(2*RHO*V)*prev['force'][mask]
    if not ns.released:xp[top,axis]=ref['X'][top,axis]+boundary(prev['time']+h)[0]
    kk,dp=update(model,ref,xp,prev['kappa']);ddp=float(np.max(dp-prev['d']))
    if ddp<=.1:break
    assert h>BASE/8*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1
   if depth:substeps+=1
   if h!=ns.dt:ns.l.command(f'timestep {h:.17g}');ns.dt=h
   ns.l.command('run 1 pre no post no');s=ns.state();p=prev;dx=s['x']-p['x'];dv=s['velocity']-p['velocity'];dd=s['d']-p['d'];sc=s['scalars'];new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],p['d']) if np.any(dd) else new;Uc,fc,pairs=contact.evaluate(s['x'],s['d']);nf=np.flatnonzero((s['d']==1)&(p['d']<1));soft=np.flatnonzero((s['d']>0)&(p['d']==0))
   error('predictor_x_m',np.max(abs(xp-s['x'])));error('material_force_N',np.max(abs(new['force']-s['material_force'])));error('contact_force_N',np.max(abs(s['force']-s['material_force']-fc)));error('energy_J',abs(new['U']-sc[0]));error('dissipation_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));error('kappa',np.max(abs(kk-s['kappa'])));error('damage',np.max(abs(dp-s['d'])));error('irreversibility',max(np.max(-dd),np.max(p['kappa']-s['kappa'])));error('contact_U_J',abs(Uc-s['Uc']));error('KE_J',abs(.5*RHO*V*np.sum(s['velocity'][:n]**2)-s['K']));error('x_VV_m',np.max(abs(dx[mask]-h*p['velocity'][mask]-h*h/(2*RHO*V)*p['force'][mask])));error('v_VV_m_s',np.max(abs(dv[mask]-h/(2*RHO*V)*(p['force'][mask]+s['force'][mask]))))
   if not ns.released:W+=float(-.5*np.sum((p['force'][top,axis]+s['force'][top,axis])*dx[top,axis])+.5*RHO*V*np.sum(s['velocity'][top,axis]**2-p['velocity'][top,axis]**2))
   dJ=0.
   if len(nf):
    rr=np.linalg.norm(s['x'][b[nf,1]]-s['x'][b[nf,0]],axis=1);dJ=float(.5*contact.c*np.sum(np.minimum(rr-DX,0)**2))
   Jc+=dJ;res=W-(sc[0]-initial['scalars'][0]+s['K']-initial['K']+sc[1]+Uc-initial['Uc']-Jc);scale=max(scale,abs(W),abs(sc[0]-initial['scalars'][0])+s['K']+sc[1]+Uc);rel=abs(res)/scale
   dKfree=.5*RHO*V*np.sum(s['velocity'][mask]**2-p['velocity'][mask]**2);defects+=np.array([dKfree-.5*np.sum((p['force'][mask]+s['force'][mask])*dx[mask]),old['U']-last['U']+.5*np.sum((p['material_force']+old['force'])*dx),.5*np.sum((s['material_force']-old['force'])*dx),Uc-p['Uc']-dJ+.5*np.sum((prevfc+fc)*dx)]);error('ledger_decomposition_J',abs(res+defects.sum()))
   assert err['predictor_x_m']<2e-18 and err['x_VV_m']<2e-18 and err['v_VV_m_s']<1e-10,err
   assert max(err['material_force_N'],err['contact_force_N'])<1e-12 and max(err['energy_J'],err['dissipation_increment_J'],err['contact_U_J'],err['ledger_decomposition_J'],err['KE_J'])<1e-17,err
   assert max(err['kappa'],err['damage'],err['irreversibility'])<1e-10 and np.max(dd)<=.1+1e-10,err
   assert rel<=.01,('ENERGY_RESIDUAL',rel)
   assert np.array_equal(s['atom_ids'],np.arange(1,ns.N+1)) and np.array_equal(ns.field(9),np.arange(1,len(b)+1)) and np.array_equal(ns.field(12),bt)
   assert np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any() and not s['d'][bt==3].any()
   assert s['native_bonds']==int(np.sum(s['d']<1)) and sc[3]==int(np.sum((s['d']==1)&~cut)) and sc[4]==s['step'] and sc[5]==s['step']
   assert np.isfinite(s['x']).all() and np.isfinite(s['force']).all()
   if lasttop is None or len(nf):lasttop=topology(ref,s['d'])
   labels,sizes,attached=lasttop;gripped=np.zeros(len(sizes),bool)
   if not ns.released:gripped[labels[top]]=True
   free=sizes[~attached&~gripped];largest=int(free.max(initial=0));anchorfree=sizes[~attached];natural=(s['d']==1)&~cut;tipdist=np.hypot(.5*(ref['X'][b[:,0],0]+ref['X'][b[:,1],0])+DX,.5*(ref['X'][b[:,0],2]+ref['X'][b[:,1],2]));tipfailed=int(np.sum(natural&(bt==1)&(tipdist<=HORIZON)));bulkfailed=int(np.sum(natural&(bt==1)));intfailed=int(np.sum(natural&(bt==2)))
   row=dict(step=s['step'],time_s=s['time'],dt_s=h,displacement_m=boundary(s['time'])[0],reaction_N=-float(np.sum(s['force'][top,axis])),anchor_reaction_N=-float(np.sum(s['force'][n:,axis])),released=int(ns.released),bulk_failed=bulkfailed,interface_failed=intfailed,tip_bulk_failed=tipfailed,detached_nodes=int(free.sum()),detached_components=len(free),largest_detached_nodes=largest,largest_anchor_free_nodes=int(anchorfree.max(initial=0)),N_ge200=int(np.sum(free>=34)),N_ge300=int(np.sum(free>=114)),W_J=W,U_J=float(sc[0]),K_J=s['K'],D_J=float(sc[1]),Uc_J=Uc,Jc_J=Jc,R_J=res,relative_R=rel,max_dd=float(np.max(dd)),max_speed_m_s=float(np.max(np.linalg.norm(s['velocity'],axis=1))),contact_pairs=len(pairs))
   if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
   writer.writerow(row);rows.append(row);changed=np.flatnonzero(dd>0)
   if len(changed):changes.write(json.dumps(dict(step=s['step'],time_s=s['time'],bond_ids=(changed+1).tolist(),d=s['d'][changed].tolist(),kappa=s['kappa'][changed].tolist()))+'\n')
   for kind,ids in [('bulk_soft',soft[bt[soft]==1]),('interface_soft',soft[bt[soft]==2]),('bulk_failure',nf[bt[nf]==1]),('interface_failure',nf[bt[nf]==2])]:
    if len(ids) and kind not in first:first[kind]=dict(step=s['step'],time_s=s['time'],bond_ids=ids+1,reference_midpoints_m=.5*(ref['X'][b[ids,0]]+ref['X'][b[ids,1]]));save(s,kind)
   if largest>maxfree:maxfree=largest;maxkey=s['step'];save(s,'new_maximum_free_component')
   if s['step']%1000==0:save(s,'periodic');stream.flush();changes.flush();print(name,'STEP',s['step'],'ms',s['time']*1e3,'bulk',bulkfailed,'interface',intfailed,'largest_free',largest,'R',rel,'wall',time.perf_counter()-start,flush=True)
   last=new;prevfc=fc;prev=s
 except Exception as exc:
  stop.append(repr(exc));(case/'FAILURE_TRACEBACK.txt').write_text(traceback.format_exc());save(s,'failed_candidate_not_accepted');s=prev
 finally:
  stream.close();changes.close();save(s,'last_accepted');labels,sizes,attached=topology(ref,s['d']);comps=[]
  for q in range(len(sizes)):
   ids=np.flatnonzero(labels==q);out,inside=edge_boundary(ref,ids);grip=not ns.released and bool(np.intersect1d(ids,top).size)
   comps.append(dict(atom_ids=ids+1,nodes=len(ids),status='ATTACHED' if attached[q] else ('GRIPPED' if grip else 'ACTUAL DETACHED'),internal_alive_bulk=int(np.sum((bt[inside]==1)&(s['d'][inside]<1))),internal_reference_bulk=int(np.sum(bt[inside]==1)),external_alive=int(np.sum(s['d'][out]<1)),mass_kg=len(ids)*RHO*V))
  count,raw=ns.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);wanted=b[s['d']<1]+1;topologyok=np.array_equal(np.sort(actual[:,0]*ns.N+actual[:,1]),np.sort(wanted[:,0]*ns.N+wanted[:,1]));assert topologyok
  summary=dict(status='COMPLETED' if not stop else 'FAILED_ADMISSION',stop_reasons=stop,case=name,axis=axis,time_s=float(s['time']),accepted_steps=len(rows),predictive_substeps=substeps,first_events=first,last_row=rows[-1] if rows else None,max_relative_energy_residual=max([r['relative_R'] for r in rows],default=0),max_independent_errors=err,maximum_free_nodes=maxfree,maximum_free_key=maxkey,final_components=comps,native_final_topology_matches=topologyok,library_sha256=sha(LIB),wall_s=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
  dump(case/'SUMMARY.json',summary);dump(case/'EVENTS.json',events);ns.close()
 print('CASE_DONE',name,summary['status'],summary['wall_s'],flush=True);return summary

def main():
 prepare();out=[]
 for name,axis in [('opening_z',2),('shear_x',0)]:
  z=run_case(name,axis);out.append(z)
  if z['status']!='COMPLETED':break
 dump(O/'SMALL_SPECIMEN_VALIDATION.json',dict(status='COMPLETED' if len(out)==2 and all(x['status']=='COMPLETED' for x in out) else 'INCOMPLETE',label='SYNTHETIC BENCHMARK',design_freeze='benchmark/DESIGN_FREEZE.json',cases=out,material_law_changed=False,real_V3_rerun=False))
if __name__=='__main__':main()
