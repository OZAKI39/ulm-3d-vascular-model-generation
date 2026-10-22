"""All-step independent audit and fixed-ID capture/rotation/morphology statistics."""
import time,json,concurrent.futures,multiprocessing
import numpy as np,h5py
from scipy.spatial.distance import pdist,cdist
from scripts.straight_tube_exp1.f10af.metrics import spatial,rotation
from .common import *
from .capture import capture
from .paired import aligned

def roots_for(t):
 phase=round(np.degrees(t['phases'][3]-t['phases'][2]))%360
 rows=read(BASE/'STATIC_RESULTS.json')['rows'];r=next(v for v in rows if v['frequency_Hz']==t['frequency'] and v['phase_deg']==phase) if phase in [0,90,180,270] else next(v for v in rows if v['frequency_Hz']==2500000 and v['phase_deg']==0)
 return np.array([q['position_m'][0] for q in r['roots']])
def audit_one(r):
 _,_,audit,_=bind();return audit.one(r)
def analyze_one(job):
 kind,r=job;begin=time.perf_counter();path=BASE/r['trajectory_path'];roots=roots_for(r['task']);out=BASE/'series';out.mkdir(exist_ok=True)
 prior=out/(kind+'_'+r['id']+'.json')
 if prior.exists():
  saved=read(prior)
  if saved.get('analysis_source_sha256')==sha(Path(__file__)) and saved.get('raw_sha256')==r['trajectory_sha256']:return saved
 with h5py.File(path) as h:
  alltime=h['time_s'][:];allx=h['positions_m'][:];n=allx.shape[1];ix=np.unique(np.rint(np.linspace(0,len(alltime)-1,min(401,len(alltime)))).astype(int))
  if kind=='HF':t=alltime[ix];x=allx[ix];v=h['force_samples/velocity_m_s'][:len(alltime)-1];vv=np.asarray([v[min(k,len(v)-1)] for k in ix]);ft=h['force_samples/total_force_N'][:len(alltime)-1];F1=h['force_samples/primary_force_N'][:len(alltime)-1];F2=h['force_samples/secondary_force_N'][:len(alltime)-1]
  else:t=h['sample_time_s'][:];x=h['sample_positions_m'][:];vv=np.gradient(x,t,axis=0);v=np.diff(allx,axis=0)/np.diff(alltime)[:,None,None];ft=F1=F2=None
  if kind=='LF':
   # Recover original acoustic diagnostics only at stored accepted LF states.
   from scripts.straight_tube_exp1.f10af.field import response_context
   from scripts.straight_tube_exp1.f10ae.schedules import schedule
   from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
   from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
   c=response_context(r['task']['frequency']);bs=CachedWaterBasis(BASE/'fields'/str(int(r['task']['frequency']))/('fine.h5' if r['task'].get('variant')=='FINE' else 'coarse.h5'));sim=WaterAccelerated(allx[0],c['physics'],c['response'],bs,c['cavity'],policy=policy());sim.control(schedule(r['task']).segments[0].event());sampled=[sim.fields(allx[k]) for k in ix];dest=BASE/'lf_force_diagnostics';dest.mkdir(exist_ok=True)
   with h5py.File(dest/(r['id']+'.h5'),'w') as z:
    z['time_s']=alltime[ix];z['positions_m']=allx[ix];z['particle_id']=np.arange(n)
    for key in ['pressure_Pa','primary_force_N','secondary_force_N','total_force_N']:z[key]=np.array([v[key] for v in sampled])
    z.attrs.update(scope='Original acoustic F1/F2 at recorded LF accepted states; no new dynamics',identity_json=json.dumps(identity()))
  # Same fixed-ID COM removal, all accepted steps for cumulative HF angle.
  q=allx[:-1]-allx[:-1].mean(1)[:,None,:];w=v-v.mean(1)[:,None,:];cross=np.cross(q,w)[:,:,2].sum(1);den=np.sum(q*q,axis=(1,2));denxy=np.sum(q[:,:,:2]**2,axis=(1,2));omega=cross/np.maximum(den,1e-99);omxy=cross/np.maximum(denxy,1e-99);rig=np.stack([-omxy[:,None]*q[:,:,1],omxy[:,None]*q[:,:,0],np.zeros(q.shape[:2])],-1);energy=np.sum(w*w,axis=(1,2));fraction=np.where(energy>1e-40,1-np.sum((w-rig)**2,axis=(1,2))/np.maximum(energy,1e-99),0);dt=np.diff(alltime);angle=np.r_[0,np.cumsum(omega*dt)]*180/np.pi
  metrics=[spatial(z) for z in x];S=np.array([m['structure_peak_S'] for m in metrics]);order=np.array([m['stripe_order'] for m in metrics]);occupied=np.array([m['occupied_bands'] for m in metrics]);rg=np.array([m['Rg_3d_m'] for m in metrics]);components=np.array([m['proximity_components'] for m in metrics]);spacing=np.array([m['stripe_spacing_m'] for m in metrics]);orientation=np.array([m['stripe_normal_deg'] for m in metrics]);band=np.abs(x[:,:,0,None]-roots).argmin(2)
  cap={str(int(round(th*1e6))):capture(alltime,allx,roots,th) for th in read(CONFIG)['capture']['thresholds_m']};cap3=capture(alltime,allx,roots,20e-6,full_3d=True)
  physical_period=float(np.mean(np.diff(roots)));amp=np.exp(2j*np.pi*x[:,:,0]/physical_period).mean(1);Sx=n*abs(amp-np.sinc(.001/physical_period))**2
  gain=order-order[0];descriptor=np.array([m['stripe_descriptor_pass'] for m in metrics]);pattern=descriptor & (gain>=.15) & (not descriptor[0]);last_10=t>=max(0,t[-1]-.01);persist=bool(pattern[last_10].all() and t[-1]>=.01)
  # First true near-contact event and final neighbor graph.
  ii,jj=np.triu_indices(n,1);gaps=np.linalg.norm(allx[:,jj]-allx[:,ii],axis=2)-2*R;hit=np.flatnonzero((gaps<=1e-9).any(1));first_contact=float(alltime[hit[0]]) if len(hit) else None
  np.savez_compressed(out/(kind+'_'+r['id']+'.npz'),time_s=t,positions_m=x,velocity_m_s=vv,Rg_m=rg,S=S,Sx=Sx,order=order,order_gain=gain,occupied=occupied,components=components,spacing_m=spacing,orientation_deg=orientation,band=band,omega_s=omega,omega_time_s=alltime[:-1],angle_deg=angle,angle_time_s=alltime,rigid_fraction=fraction,roots_m=roots)
  pr=None
  if kind=='HF':P=h['force_samples/pressure_Pa'][:len(alltime)-1];pr=[float(abs(P).min()),float(abs(P).max())]
  rotmean=float(np.sum(fraction*dt)/max(alltime[-1],1e-99));sign=float(np.sum(dt[np.sign(omega)==np.sign(angle[-1])])/max(alltime[-1],1e-99)) if angle[-1]!=0 else 0.;rotation_pass=abs(angle[-1])>=1 and rotmean>=.8 and sign>=.8
  row=dict(id=r['id'],kind=kind,status=r['status'],state=r['task']['state'],N=n,action=r['task']['action'],task=r['task'],completed_s=float(alltime[-1]),requested_s=r['task']['duration'],capture=cap,full_3d_point_capture_20um=cap3,initial_rg_m=rg[0],final_rg_m=rg[-1],delta_rg_m=rg[-1]-rg[0],COM_displacement_m=float(np.linalg.norm(x[-1].mean(0)-x[0].mean(0))),max_particle_displacement_m=float(np.linalg.norm(x[-1]-x[0],axis=1).max()),structure_initial=S[0],structure_final=S[-1],structure_max=float(S.max()),structure_X_final=float(Sx[-1]),order_gain_final=gain[-1],occupied_final=int(occupied[-1]),stripe_spacing_final_m=spacing[-1],stripe_normal_final_deg=orientation[-1],stripe_descriptor_initial=bool(descriptor[0]),stripe_descriptor_final=bool(descriptor[-1]),stripe_dynamics_candidate=bool(pattern[-1]),stripe_persistence=persist,final_components=int(components[-1]),first_contact_s=first_contact,min_gap_m=float(gaps.min()),final_band=band[-1],angle_deg=float(angle[-1]),mean_rigid_fraction=rotmean,rotation_sign_persistence=sign,rotation_candidate=rotation_pass,pressure_range_Pa=pr,HF_fallback_steps=r.get('HF_fallback_steps'),wall_s=r.get('task_wall_s',r.get('wall_s')),analysis_wall_s=time.perf_counter()-begin,reason=r.get('reasons',r.get('reason')),raw_sha256=r['trajectory_sha256'],analysis_source_sha256=sha(Path(__file__)),identity=identity())
  if ft is not None:row['force_summary']=dict(mean_F1_magnitude_N=float(np.linalg.norm(F1,axis=2).mean()),max_F1_N=float(np.linalg.norm(F1,axis=2).max()),mean_F2_magnitude_N=float(np.linalg.norm(F2,axis=2).mean()),max_F2_N=float(np.linalg.norm(F2,axis=2).max()))
 write(out/(kind+'_'+r['id']+'.json'),row);print('ANALYZED',kind,r['id'],flush=True);return clean(row)

def compare_pair(a,b):
 # Exact HF timestamps matched to LF dense output, neither padded nor interpolated HF.
 from .common import bind
 roots=roots_for(a['task']);path=BASE/a['trajectory_path'];bp=BASE/b['trajectory_path'];r=dict(id=a['id'],HF_status=a['status'],LF_status=b['status'],PASS=False)
 if a['status']!='PASS' or b['status']!='PASS':return r
 with h5py.File(path) as h,h5py.File(bp) as z:
  tt,hx,lx,alignment=aligned(h,z);diff=hx-lx;err=np.linalg.norm(diff,axis=2);cent=np.linalg.norm(diff.mean(1),axis=1);rel=np.linalg.norm(diff-diff.mean(1)[:,None,:],axis=2);hb=abs(hx[:,:,0,None]-roots).argmin(2);lb=abs(lx[:,:,0,None]-roots).argmin(2)
  network_dif=0;pairerr=0;cov_error=0
  for xx,yy in zip(hx,lx):
   dx=pdist(xx);dy=pdist(yy);pairerr=max(pairerr,float(abs(dx-dy).max()));network_dif+=int(np.count_nonzero((dx<=2.1*R)!=(dy<=2.1*R)));cov_error=max(cov_error,float(np.linalg.norm(np.cov(xx.T,bias=True)-np.cov(yy.T,bias=True))))
  ah=read(BASE/'series'/('HF_'+a['id']+'.json'));al=read(BASE/'series'/('LF_'+b['id']+'.json'));caperr=[];capdis=0
  for th in ['10','20','30']:
   for u,v in zip(ah['capture'][th],al['capture'][th]):
    capdis+=int(u['right_censored']!=v['right_censored'])
    if u['first_held_entry_s'] is not None and v['first_held_entry_s'] is not None:caperr.append(abs(u['first_held_entry_s']-v['first_held_entry_s']))
  r.update(max_fixed_ID_m=float(err.max()),COM_error_m=float(cent.max()),relative_ID_error_m=float(rel.max()),pair_distance_error_m=pairerr,covariance_error_m2=cov_error,Rg_error_m=abs(ah['final_rg_m']-al['final_rg_m']),band_disagreements=int((hb!=lb).sum()),neighbor_disagreements=network_dif,capture_outcome_disagreements=capdis,max_capture_time_difference_s=max(caperr,default=0),capture_time_resolution_s=float(np.diff(z['time_s'][:]).max()),structure_peak_difference=abs(ah['structure_final']-al['structure_final']),stripe_spacing_difference_m=abs(ah['stripe_spacing_final_m']-al['stripe_spacing_final_m']),orientation_difference_deg=min(abs(ah['stripe_normal_final_deg']-al['stripe_normal_final_deg']),180-abs(ah['stripe_normal_final_deg']-al['stripe_normal_final_deg'])),cumulative_angle_difference_deg=abs(ah['angle_deg']-al['angle_deg']),first_contact_HF_s=ah['first_contact_s'],first_contact_LF_s=al['first_contact_s'],paired_duration_s=a['completed_s'],alignment=alignment,PASS=bool(err.max()<=50e-9 and (hb==lb).all() and capdis==0))
 return clean(r)
def main():
 begin=time.perf_counter();verify();hf=[read(p) for p in sorted((BASE/'hf').glob('*/SUMMARY.json'))];lf=[read(p) for p in sorted((BASE/'lf').glob('*/SUMMARY.json'))];hf=[r for r in hf if r.get('trajectory_path')];lf=[r for r in lf if r.get('trajectory_path')]
 known={v['id']:v for v in read(BASE/'AUDIT.json')['rows']} if (BASE/'AUDIT.json').exists() else {v['audit']['id']:v['audit'] for v in read(BASE/'BASELINE_RESULTS.json')['rows']}
 todo=[r for r in hf if r['id'] not in known or known[r['id']]['raw_sha256']!=r['trajectory_sha256']]
 with concurrent.futures.ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as ex:
  for v in ex.map(audit_one,todo):known[v['id']]=v
 audits=[known[r['id']] for r in hf]
 write(BASE/'AUDIT.json',dict(rows=audits,all_pass=all(q['audit']['PASS'] for q in audits),checked_steps=sum(q['audit']['checked_steps'] for q in audits),identity=identity()))
 with concurrent.futures.ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as ex:rows=list(ex.map(analyze_one,[('HF',r) for r in hf]+[('LF',r) for r in lf]))
 paired=[compare_pair(a,next(b for b in lf if b['id']==a['id'])) for a in hf if any(b['id']==a['id'] for b in lf)];conv=[]
 for r in hf:
  if r['task'].get('variant') not in ['DT05','FINE']:continue
  parent=r['id'].rsplit('_',1)[0];a=next((a for a in hf if a['id']==parent),None)
  if not a or a['status']!=r['status'] or r['status']!='PASS':continue
  with h5py.File(BASE/a['trajectory_path']) as h,h5py.File(BASE/r['trajectory_path']) as z:
   ht=h['time_s'][:];zt=z['time_s'][:];ix=np.unique(np.rint(np.linspace(0,len(ht)-1,401)).astype(int));zi=np.rint(ht[ix]/r['dt_s']).astype(int);assert abs(ht[ix]-zt[zi]).max()<1e-12;err=float(np.linalg.norm(h['positions_m'][ix]-z['positions_m'][zi],axis=2).max())
  conv.append(dict(id=r['id'],parent=parent,error_m=err,PASS=err<=50e-9,completed_s=r['completed_s']))
 write(BASE/'SCIENTIFIC_RESULTS.json',dict(rows=rows,paired=paired,convergence=conv,audit_pass=all(q['audit']['PASS'] for q in audits),max_HF_time_s=max(r['completed_s'] for r in hf if r['status']=='PASS'),max_LF_time_s=max(r['completed_s'] for r in lf if r['status']=='PASS'),max_paired_admitted_s=max([r['paired_duration_s'] for r in paired if r['PASS']]+[0]),HF_attempts=len(hf),HF_complete=sum(r['status']=='PASS' for r in hf),LF_attempts=len(lf),LF_complete=sum(r['status']=='PASS' for r in lf),complete_campaign=bool((BASE/'LF_LADDER_RESULTS.json').exists()),wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
