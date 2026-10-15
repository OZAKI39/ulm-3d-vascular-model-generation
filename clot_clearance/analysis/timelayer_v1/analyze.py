"""Read-only accepted-state momentum and signed energy ledger.
No counterfactual dynamics or additional LAMMPS evaluation.
"""
from common import *
from model import LPS
import numpy as np,csv,time,sys
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0,str(R))
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy
start=time.perf_counter();freeze=json.loads((O/'TIMESTEP_DESIGN_FREEZE.json').read_text());scale=freeze['reference_scales'];g=np.load(O/'historical_static/reference.npz');model=LPS(g);b=g['bonds'];mass=freeze['mass_per_node_kg'];N=len(g['X']);free=np.ones(N,bool);free[89]=False
eng=Energy({'X':g['X'],'bonds':b},50e-6,R/'analysis/material_macro_audit_v1/energy.so')
def normmax(a):return float(np.max(np.linalg.norm(a,axis=-1)))
def csvout(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def components(alive):
 a,j=b[alive].T;mat=coo_matrix((np.ones(2*len(a)),(np.r_[a,j],np.r_[j,a])),shape=(N,N)).tocsr();nc,lab=connected_components(mat,directed=False)
 return int(nc),lab,np.bincount(lab)
allseries={};summaries={};event_recs={};forces={};energy_details={};impulses={};raw={}
for case in freeze['formal_cases']:
 name=case['name'];dt=case['dt_s'];d=O/name;s=np.load(d/'native_states.npz');raw[name]=s;x=s['x'];v=s['v_full'];ff=s['force_N'];t=s['time_s'];alive=np.unpackbits(s['alive_packed'],axis=1)[:,:len(b)].astype(bool);nn=len(t);fpost=np.zeros_like(ff);lag=np.zeros_like(ff);rows=[];evs=[];fd=[];uvals=[];J=np.zeros((N,3));Jconsistent=J.copy();Jtop=J.copy();Jtotal=J.copy();Jhistory=[];actualhistory=[];maxp=0.;maxKidentity=0.;decomp_err=0.;W=DU=WLAG=WTOP=Q=KICK=0.;merr=float(np.max(abs(s['m']-g['m'])/g['m']))
 U0=model.eval(x[0],alive[0])[2];K0=.5*mass*np.sum(v[0,free]**2);Kprev=K0;previous_U=U0
 for k in range(nn):
  before=alive[k-1] if k else alive[0];broken=before&~alive[k]
  thpre,fpre,upre=model.eval(x[k],before);thpost,fpost[k],u=model.eval(x[k],alive[k]);_,mixed,_=model.eval(x[k],alive[k],thpre);lag[k]=mixed-fpost[k]
  assert normmax(ff[k]-mixed)<2e-9*scale['force_N']
  thetaerr=float(np.max(abs(s['theta_native'][k]-thpre)));assert thetaerr<1e-12
  if broken.any():
   uc,fc=eng.eval(np.ascontiguousarray(x[k]),alive[k].astype(np.uint8),True);cpp_err=normmax(fc-fpost[k]);assert cpp_err<2e-9*scale['force_N']
   direction=np.random.default_rng(123+k).normal(size=x[k].shape);direction/=np.linalg.norm(direction);eps=5e-10
   du=(eng.eval(x[k]+eps*direction,alive[k].astype(np.uint8))-eng.eval(x[k]-eps*direction,alive[k].astype(np.uint8)))/(2*eps);fpred=-float(np.sum(fpost[k]*direction));fderr=abs(du-fpred)/abs(fpred);assert fderr<2e-6
   rec={'step':k,'time_s':float(t[k]),'bond_ids':(np.flatnonzero(broken)+1).tolist(),'atom_id_pairs':(b[broken]+1).tolist(),'reference_midpoints_m':g['X'][b[broken]].mean(1).tolist(),'max_theta_difference':float(abs(thpre-thpost).max()),'DeltaF_topology_max_N':normmax(fpost[k]-fpre),'DeltaF_timelayer_max_all_N':normmax(lag[k]),'DeltaF_timelayer_max_free_N':normmax(lag[k,free]),'max_lag_all_atom_id':int(np.argmax(np.linalg.norm(lag[k],axis=1))+1),'max_lag_free_atom_id':int(np.flatnonzero(free)[np.argmax(np.linalg.norm(lag[k,free],axis=1))]+1),'current_final_half_kick_proxy_max_free_Ns':.5*dt*normmax(lag[k,free]),'full_two_half_kick_proxy_max_free_Ns':dt*normmax(lag[k,free]),'U_pre_J':upre,'U_post_J':u,'DeltaU_topology_J':u-upre,'native_mixed_error_N':normmax(mixed-ff[k]),'independent_CPP_post_force_error_N':cpp_err,'gradient_FD_relative':float(fderr)}
   np.savez_compressed(d/f'decomposition_{k:04d}.npz',x=x[k],X=g['X'],bonds=b,alive_pre=before,alive_post=alive[k],theta_pre=thpre,theta_post=thpost,F_pre_N=fpre,F_mixed_N=mixed,F_post_N=fpost[k],DeltaF_topology_N=fpost[k]-fpre,DeltaF_timelayer_N=lag[k],v_half=s['v_half'][k],v_full=v[k])
   # Same-half-velocity final-kick comparison: algebraic diagnostic, not an alternate run.
   vc=s['v_half'][k,free]+.5*dt*fpost[k,free]/mass
   rec['same_half_velocity_final_kick_DeltaK_J']=float(.5*mass*(np.sum(v[k,free]**2)-np.sum(vc**2)))
   evs.append(rec)
  deltaU=u-upre;DU+=deltaU;K=.5*mass*float(np.sum(v[k,free]**2));dW=dLag=dTop=dQ=dKick=0.;perr=0.
  if k:
   dx=x[k]-x[k-1];avg=.5*(ff[k-1]+ff[k]);dW=-float(np.sum(avg[~free]*dx[~free]));W+=dW
   dLag=float(np.sum(.5*(lag[k-1]+lag[k])*dx));WLAG+=dLag
   dTop=float(np.sum(.5*(fpost[k]-fpre)*dx));WTOP+=dTop
   dQ=upre-previous_U+float(np.sum(.5*(fpost[k-1]+fpre)*dx));Q+=dQ
   dKick=dt**2/(8*mass)*float(np.sum(ff[k,free]**2)-np.sum(ff[k-1,free]**2));KICK+=dKick
   actualK=K-Kprev-float(np.sum(avg[free]*dx[free]));maxKidentity=max(maxKidentity,abs(actualK-dKick))
   J+=.5*dt*(lag[k-1]+lag[k]);Jtop+=.5*dt*(fpost[k]-fpre);Jconsistent+=.5*dt*(fpost[k-1]+fpre);Jtotal+=dt*avg
   actual=mass*(v[k]-v[0]);perr=normmax(actual[free]-Jtotal[free]);maxp=max(maxp,perr)
   assert perr<freeze['gates']['momentum_balance_abs_Ns']
  else:dx=np.zeros_like(x[k]);actual=mass*(v[k]-v[0])
  Renergy=K+u-K0-U0-W-DU;accounted=Q+WLAG+WTOP+KICK;decomp_err=max(decomp_err,abs(Renergy-accounted));previous_U=u;Kprev=K
  a,j=b[before].T;stretch=(np.linalg.norm(x[k,j]-x[k,a],axis=1)-model.r[before])/model.r[before];nc,labels,counts=components(alive[k])
  pred_u=freeze['loading']['amplitude_m']*.5*(1-np.cos(np.pi*min(t[k]/freeze['loading']['ramp_s'],1)))
  load_error=abs(x[k,89,0]-x[0,89,0]-pred_u)
  row=dict(step=k,time_s=float(t[k]),alive_bonds=int(alive[k].sum()),new_breaks=int(broken.sum()),cumulative_breaks=int((~alive[k]).sum()),max_prebreak_stretch_utilization=float(stretch.max()/.15),vmax_all_m_s=normmax(v[k]),vmax_free_m_s=normmax(v[k,free]),driven_vx_m_s=float(v[k,89,0]),neighbor4_vx_m_s=float(v[k,3,0]),max_step_displacement_m=normmax(dx),K_free_J=K,K_all_J=.5*mass*float(np.sum(v[k]**2)),U_LPS_J=u,U_pre_same_x_J=upre,DeltaU_topology_J=deltaU,cumulative_DeltaU_topology_J=DU,W_boundary_J=W,damping_work_J=0.,contact_work_J=0.,energy_residual_J=Renergy,lag_force_max_all_N=normmax(lag[k]),lag_force_max_free_N=normmax(lag[k,free]),J_lag_proxy_max_free_Ns=normmax(J[free]),DeltaP_actual_max_free_Ns=normmax(actual[free]),DeltaP_node4_x_Ns=float(actual[3,0]),J_lag_node4_x_Ns=float(J[3,0]),momentum_reconstruction_error_Ns=perr,W_lag_on_actual_path_J=WLAG,W_topology_switch_J=WTOP,conservative_trapezoid_error_J=Q,VV_endpoint_energy_J=KICK,energy_decomposition_residual_J=Renergy-accounted,components=nc,largest_component_nodes=int(counts.max()),prescribed_position_error_m=float(load_error))
  rows.append(row);Jhistory.append(J.copy());actualhistory.append(actual.copy())
 assert max(abs(r['prescribed_position_error_m']) for r in rows)<2e-18
 assert decomp_err<1e-20 and maxKidentity<1e-20
 csvout(d/'analysis_scalars.csv',rows);allseries[name]=rows;event_recs[name]=evs
 np.savez_compressed(d/'impulse_on_actual_path.npz',time_s=t,J_lag_proxy_Ns=np.asarray(Jhistory),DeltaP_actual_Ns=np.asarray(actualhistory),J_total_force_Ns=Jtotal,J_consistent_pre_Ns=Jconsistent,J_topology_switch_Ns=Jtop,free_node_mask=free)
 nc,labels,counts=components(alive[-1]);first=evs[0] if evs else None
 vmax=np.linalg.norm(v,axis=2);imax=np.unravel_index(np.argmax(vmax),vmax.shape);vfree=vmax.copy();vfree[:,~free]=-np.inf;imaxf=np.unravel_index(np.argmax(vfree),vfree.shape)
 summary={'dt_s':dt,'accepted_steps':case['steps'],'end_time_s':float(t[-1]),'first_event':first,'final_broken_bond_ids':(np.flatnonzero(~alive[-1])+1).tolist(),'final_broken_atom_id_pairs':(b[~alive[-1]]+1).tolist(),'event_count':len(evs),'broken_bonds':int((~alive[-1]).sum()),'final_components':nc,'final_component_sizes':sorted(counts.tolist(),reverse=True),'peak_K_free_J':max(r['K_free_J'] for r in rows),'peak_K_all_J':max(r['K_all_J'] for r in rows),'peak_v_all_m_s':float(vmax[imax]),'peak_v_all_atom_id':int(imax[1]+1),'peak_v_all_time_s':float(t[imax[0]]),'peak_v_free_m_s':float(vmax[imaxf]),'peak_v_free_atom_id':int(imaxf[1]+1),'peak_v_free_time_s':float(t[imaxf[0]]),'max_step_displacement_m':max(r['max_step_displacement_m'] for r in rows),'max_m_relative':merr,'U_initial_J':U0,'U_final_J':rows[-1]['U_LPS_J'],'W_boundary_final_J':W,'signed_topology_sum_J':DU,'signed_topology_positive_J':sum(max(r['DeltaU_topology_J'],0) for r in rows),'signed_topology_negative_J':sum(min(r['DeltaU_topology_J'],0) for r in rows),'max_abs_energy_residual_J':max(abs(r['energy_residual_J']) for r in rows),'max_abs_energy_residual_over_frozen_Uref':max(abs(r['energy_residual_J']) for r in rows)/scale['energy_J'],'final_energy_residual_J':rows[-1]['energy_residual_J'],'lag_work_final_J':WLAG,'topology_switch_work_final_J':WTOP,'conservative_quadrature_final_J':Q,'VV_endpoint_final_J':KICK,'max_energy_decomposition_closure_J':decomp_err,'max_VV_kinetic_identity_error_J':maxKidentity,'max_momentum_balance_abs_Ns':maxp,'max_J_lag_proxy_free_Ns':max(r['J_lag_proxy_max_free_Ns'] for r in rows),'final_J_lag_proxy_free_Ns':rows[-1]['J_lag_proxy_max_free_Ns'],'final_actual_momentum_free_Ns':rows[-1]['DeltaP_actual_max_free_Ns'],'max_force_lag_free_N':max(r['lag_force_max_free_N'] for r in rows),'max_force_lag_all_N':max(r['lag_force_max_all_N'] for r in rows)}
 summaries[name]=summary;dump(d/'CASE_ANALYSIS.json',summary);dump(d/'EVENT_DECOMPOSITION.json',evs)
 # Actual force and accepted-velocity integration identity, not a causal trajectory subtraction.
 impulses[name]={'max_actual_momentum_balance_error_Ns':maxp,'final_J_lag_proxy_free_Ns':summary['final_J_lag_proxy_free_Ns'],'max_J_lag_proxy_free_Ns':summary['max_J_lag_proxy_free_Ns'],'final_actual_momentum_free_Ns':summary['final_actual_momentum_free_Ns'],'loaded_atom_excluded_from_NVE_impulse':90,'decomposition_free_max_error_Ns':normmax((mass*(v[-1]-v[0])-Jconsistent-Jtop-J)[free]),'event_proxies':[{k:e[k] for k in ['step','time_s','current_final_half_kick_proxy_max_free_Ns','full_two_half_kick_proxy_max_free_Ns','same_half_velocity_final_kick_DeltaK_J']} for e in evs]}
 forces[name]={'events':evs,'same_coordinate_force_interpretation':'F_post-F_pre is normal topology force change; F_mixed-F_post is lagged-theta difference; prescribed atom 90 has no NVE half kick'}
 energy_details[name]={k:v for k,v in summary.items() if '_J' in k or 'energy' in k};energy_details[name]['damping_contact_external_body']='All absent in original small fixture and adapted cases; numerical zero is not an inferred or fitted loss'
# Physical-time comparison: exact common accepted times only (fine even indices).
c=raw['dt_0p5us'];h=raw['dt_0p25us'];np.testing.assert_array_equal(c['time_s'],h['time_s'][::2]);np.testing.assert_array_equal(c['x'][0],h['x'][0]);np.testing.assert_array_equal(c['v_full'][0],h['v_full'][0])
cs=allseries['dt_0p5us'];hs=allseries['dt_0p25us'];firsttime=min(summaries[n]['first_event']['time_s'] for n in summaries if summaries[n]['first_event']);cmp=[]
ca=np.unpackbits(c['alive_packed'],axis=1)[:,:len(b)];ha=np.unpackbits(h['alive_packed'],axis=1)[:,:len(b)]
for k,r in enumerate(cs):
 q=hs[2*k];cmp.append(dict(time_s=r['time_s'],coarse_step=k,fine_step=2*k,prebreak=r['time_s']<firsttime,max_x_difference_m=normmax(c['x'][k]-h['x'][2*k]),max_v_difference_m_s=normmax(c['v_full'][k]-h['v_full'][2*k]),max_force_difference_N=normmax(c['force_N'][k]-h['force_N'][2*k]),U_difference_J=r['U_LPS_J']-q['U_LPS_J'],reaction_x_difference_N=float(-c['force_N'][k,89,0]+h['force_N'][2*k,89,0]),theta_max_difference=float(abs(c['theta_native'][k]-h['theta_native'][2*k]).max()),broken_coarse=r['cumulative_breaks'],broken_fine=q['cumulative_breaks'],topology_symmetric_difference=int(np.count_nonzero(ca[k]!=ha[2*k])),K_free_coarse_J=r['K_free_J'],K_free_fine_J=q['K_free_J'],J_lag_max_coarse_Ns=r['J_lag_proxy_max_free_Ns'],J_lag_max_fine_Ns=q['J_lag_proxy_max_free_Ns']))
csvout(O/'TIMESTEP_COMPARISON.csv',cmp)
pre=[r for r in cmp if r['prebreak']];pretests={}
for key,den,gat in [('max_x_difference_m','displacement_m','prebreak_x_difference_over_dx'),('max_v_difference_m_s','velocity_m_s','prebreak_v_difference_over_prescribed_peak'),('max_force_difference_N','force_N','prebreak_force_difference_over_historical_force'),('U_difference_J','energy_J','prebreak_energy_difference_over_historical_energy')]:
 val=max(abs(q[key]) for q in pre);pretests[key]={'max_absolute':val,'reference_scale':scale[den],'ratio':val/scale[den],'frozen_gate':freeze['gates'][gat],'pass':val/scale[den]<=freeze['gates'][gat]}
ec=event_recs['dt_0p5us'];eh=event_recs['dt_0p25us'];bybond={};samed=[]
for name,events in event_recs.items():
 bybond[name]={str(bid):e['time_s'] for e in events for bid in e['bond_ids']}
for bid,tc in bybond['dt_0p5us'].items():
 if bid in bybond['dt_0p25us']:samed.append({'bond_id':int(bid),'coarse_time_s':tc,'fine_time_s':bybond['dt_0p25us'][bid],'difference_s':bybond['dt_0p25us'][bid]-tc})
firstsame=set(ec[0]['bond_ids'])==set(eh[0]['bond_ids']);ratio=eh[0]['full_two_half_kick_proxy_max_free_Ns']/ec[0]['full_two_half_kick_proxy_max_free_Ns'] if firstsame else None
alignment={'status':'SAME_FIRST_EVENT_BONDS' if firstsame else 'FIRST_EVENT_SETS_DIFFER','cases':summaries,'per_bond_true_times':samed,'first_event_time_difference_s':eh[0]['time_s']-ec[0]['time_s'],'first_event_bond_set_equal':firstsame,'final_bond_set_equal':summaries['dt_0p5us']['final_broken_bond_ids']==summaries['dt_0p25us']['final_broken_bond_ids'],'first_common_time_topology_difference':next((r for r in cmp if r['topology_symmetric_difference']),None),'event_interpolation_used':False,'first_event_full_duration_free_impulse_proxy_ratio_half_over_base':ratio}
dump(O/'FRACTURE_EVENT_ALIGNMENT.json',alignment)
dump(O/'PREBREAK_TIMESTEP_RESPONSE.json',{'status':'PASS_FROZEN_ENGINEERING_SCREEN' if all(q['pass'] for q in pretests.values()) else 'SENSITIVE_BEFORE_FRACTURE','physical_time_window_s':[0,max(r['time_s'] for r in pre)],'checks':pretests,'note':'Initial T0 is prestrained and free interior is not equilibrated; this is transient matching, not quasi-static calibration'})
dump(O/'THETA_TIMELAYER_FORCE_AUDIT.json',{'status':'PASS_DECLARED_MIXED_AND_POST_GRADIENT','historical_static':json.loads((O/'HISTORICAL_STATIC_DECOMPOSITION_AUDIT.json').read_text()),'dynamic_cases':forces})
dump(O/'IMPULSE_AND_MOMENTUM_AUDIT.json',{'status':'ACTUAL_NATIVE_VV_BALANCE_PASS_PROXY_DISTINCT_FROM_MEASURED_JUMP','definition':'J_lag_proxy[N]=sum_n dt/2*(DeltaF[n-1]+DeltaF[n]) on actual native trajectory; only free DOFs integrate it. Actual momentum is mass*(v_full[N]-v_full[0]).','force_duration':'Event force enters current final half kick and next initial half kick, if next step executed; terminal event only half weight','cases':impulses,'first_event_proxy_ratio':ratio,'identifiability':'Algebraic decomposition of actual accepted momentum is exact within roundoff. It is not the velocity difference from an unexecuted synchronous-theta trajectory; coordinates and later forces would change in that trajectory.'})
dump(O/'TOPOLOGY_ENERGY_RELEASE_AUDIT.json',{'status':'SIGNED_SAME_COORDINATE_TOPOLOGY_CHANGE_IDENTIFIED','cases':energy_details,'convention':'DeltaU=U_post-U_pre at same actual coordinates; negative is reduction in recoverable potential, not automatically irreversible fracture energy; Gc unavailable','energy_balance':'R=K_free+U-(K0+U0)-W_boundary-sum(DeltaU_topology), with zero damping/contact/body terms in this fixture','residual_exact_decomposition':'R = sum[fixed-topology trapezoid quadrature error + .5*(lag_prev+lag_now).dx + .5*(F_post-F_pre).dx] + dt^2/(8m)*(sum_free(F_final^2)-sum_free(F_initial^2))','interpretation':'Terms use actual trajectory and prescribed-boundary trapezoid work. Time placement of topology change, lagged force work and smooth integration quadrature are separated algebraically; this is not a causal counterfactual or Gc.'})
dump(O/'ANALYSIS_SUMMARY.json',{'cases':summaries,'alignment':{k:v for k,v in alignment.items() if k!='cases'},'prebreak':pretests,'analysis_wall_s':time.perf_counter()-start})
print(json.dumps({'first_coarse':ec[0],'first_fine':eh[0],'prebreak':pretests,'summaries':summaries,'analysis_wall_s':time.perf_counter()-start}))
