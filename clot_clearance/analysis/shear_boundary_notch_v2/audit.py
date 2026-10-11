"""Read-only actual-state diagnostics; reuse the V1 conservative LPS evaluator."""
import sys,json,csv
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import label
from pd_clot.lammps_fragmentation.macro_lab_v1 import Energy,geometry,DELTA,G,K
sys.path.insert(0,str(Path(__file__).parents[1]/'material_macro_audit_v1'))
from postprocess import save,js,fracture_area
ROOT=Path(__file__).resolve().parents[2]
DX=50e-6;A=.0006**2

def regions(X,b):
 m=X[b].mean(1);tip=np.hypot(m[:,0]+.0001,m[:,2]);grip=np.minimum(abs(m[:,2]-.00035),abs(m[:,2]+.00035));surf=np.min(np.array([.0003,.0003,.0005])-abs(m),axis=1)
 reg=np.full(len(b),3,np.int8);reg[surf<=DELTA]=2;reg[grip<=DELTA]=1;reg[tip<=DELTA]=0
 return reg,tip,grip,(abs(m[:,2])<=DELTA/2)&(m[:,0]>=-.0001-DELTA/2)

def fitF(X,x,ids):
 a=X[ids]-X[ids].mean(0);z=x[ids]-x[ids].mean(0);F=np.linalg.lstsq(a,z,rcond=None)[0].T
 return F,float(np.sqrt(np.mean((a@F.T-z)**2)))

def local_F(X,x,periodic):
 tree=cKDTree(X+[.0003,.0003,.0005],boxsize=[.0006,.0006,10.]) if periodic else cKDTree(X)
 neighbors=tree.query_ball_point(tree.data,DELTA*(1+1e-12));out=[]
 for i,ids in enumerate(neighbors):
  a=X[ids]-X[i];z=x[ids]-x[i]
  if periodic:a[:,:2]-=.0006*np.rint(a[:,:2]/.0006);z[:,:2]-=.0006*np.rint(z[:,:2]/.0006)
  out.append(np.linalg.lstsq(a,z,rcond=None)[0].T)
 return np.array(out)

def pair_forces(e,x,alive):
 d=x[e.b[:,1]]-x[e.b[:,0]]
 if e.periodic:d[:,:2]-=.0006*np.rint(d[:,:2]/.0006)
 length=np.linalg.norm(d,axis=1);extension=length-e.r;v=e.V
 theta=3*np.bincount(e.b.ravel(),weights=np.repeat(alive*v*e.h*extension,2),minlength=len(x))/e.m
 coeff=(3*K-5*G)*(theta[e.b[:,0]]/e.m[e.b[:,0]]+theta[e.b[:,1]]/e.m[e.b[:,1]])+15*G*extension/e.r*(1/e.m[e.b[:,0]]+1/e.m[e.b[:,1]])
 return (alive*v*v*e.h*coeff)[:,None]*d/length[:,None]

def shear(p,out):
 out.mkdir(parents=True,exist_ok=True);g=dict(np.load(p/'geometry.npz'));X=g['X'];b=g['bonds'];roi=g['roi'];plan=js(p/'test_plan.json');periodic=plan.get('transverse_periodic_fixed_cell',False)
 e=Energy(g,DX,ROOT/('analysis/material_macro_audit_v1/energy_periodic.so' if periodic else 'analysis/material_macro_audit_v1/energy.so'),periodic)
 a=np.genfromtxt(p/'history.csv',names=True,delimiter=',');frac=abs(a['boundary_strain'])/.008;fit=(a['time_s']<=.010)&(frac>=.2)&(frac<=.7)
 Gnom=np.polyfit(a['boundary_strain'][fit],a['stress_Pa'][fit],1)[0];Groi=np.polyfit(a['gamma_xz'][fit],a['stress_Pa'][fit],1)[0];res={}
 low=X[:,2]<-.00035-1e-12;high=X[:,2]>.00035+1e-12
 for key in ['initial','fit20','fit50','fit70','peak','end']:
  if not (p/(key+'.npz')).exists():continue
  state=dict(np.load(p/(key+'.npz')));x=state['x_m'];alive=state['alive'];U,force=e.eval(np.ascontiguousarray(x),alive,True);F,rms=fitF(X,x,roi);i=int(state['step']);loc=local_F(X,x,periodic);gamma=loc[:,0,2]+loc[:,2,0];pf=pair_forces(e,x,alive)
  assembled=np.zeros_like(x);np.add.at(assembled,b[:,0],pf);np.add.at(assembled,b[:,1],-pf);assert np.linalg.norm(assembled-force)<1e-12
  upper=-force[high].sum(0);lower=-force[low].sum(0);balance=upper+lower;gauge=float((x[high,0]-X[high,0]).mean()-(x[low,0]-X[low,0]).mean())/.00085
  plane=X[:,2]==X[:,2].max();F2=np.linalg.lstsq(X[plane,:2]-X[plane,:2].mean(0),x[plane]-x[plane].mean(0),rcond=None)[0].T;area_ratio=float(np.linalg.norm(np.cross(F2[:,0],F2[:,1])))
  cuts=[]
  for z in np.arange(-.00030,.000301,.00005):
   side=X[:,2]>z;cross=side[b[:,0]]!=side[b[:,1]];traction=np.sum(pf[cross]*(side[b[cross,1]].astype(int)-side[b[cross,0]].astype(int))[:,None],axis=0)/A;cuts.append([z,*traction])
  # Bond virial reference stress proxy, explicitly distinguished from nominal traction.
  ref=X[b[:,1]]-X[b[:,0]]
  if periodic:ref[:,:2]-=.0006*np.rint(ref[:,:2]/.0006)
  pxy=np.zeros(len(X));q=pf[:,0]*ref[:,2]/(2*DX**3);np.add.at(pxy,b[:,0],q);np.add.at(pxy,b[:,1],q)
  profile=[]
  for z in np.unique(X[:,2]):
   ids=abs(X[:,2]-z)<1e-12;profile.append([z,np.mean(x[ids,0]-X[ids,0]),np.mean(gamma[ids]),np.std(gamma[ids]),np.mean(loc[ids,0,2]),np.mean(loc[ids,2,0])])
  res[key]=dict(step=i,time_s=float(state['time_s']),F_ROI=F,ROI_residual_m=rms,gamma_ROI=F[0,2]+F[2,0],gamma_nominal=gauge,rotation_xz=.5*(F[0,2]-F[2,0]),upper_external_reaction_N=upper,lower_external_reaction_N=lower,grip_net_force_N=balance,grip_x_balance_relative=float(abs(balance[0])/max(abs(upper[0]),1e-30)),current_face_area_from_affine_fit_m2=A*area_ratio,reference_cell_face_area_m2=A,point_center_hull_area_without_half_cells_m2=(.0006-DX)**2,current_area_scope='Reference nodal cell face includes half-dx margins; fitted deformed face uses cross product of two in-plane basis vectors, no fictitious hull-area correction',reaction_vs_original_csv_N=float(upper[0]-a['reaction_N'][i]),ROI_local_gamma_std=float(gamma[roi].std()),ROI_local_gamma_CV=float(gamma[roi].std()/max(abs(gamma[roi].mean()),1e-30)),ROI_mean_bond_virial_Pxz_Pa=float(pxy[roi].mean()),virial_scope='Reference bond virial proxy assigned half to each endpoint; spatial stress assignment is nonunique in inhomogeneous nonlocal material',cuts=cuts)
  np.savez_compressed(out/(key+'_fields.npz'),X=X,x=x,F=loc,ROI=roi,prescribed_mask=g['prescribed_mask'],gamma=gamma,force_N=force,profile=np.array(profile),cut_traction=np.array(cuts),bond_virial_Pxz=pxy,step=i,time_s=state['time_s'])
  assert abs(upper[0]-a['reaction_N'][i])<1e-12
 result=dict(case=plan['case'],source=str(p),G_nominal_Pa=Gnom,G_ROI_Pa=Groi,ROI_vs_nominal_strain_slope=np.polyfit(a['boundary_strain'][fit],a['gamma_xz'][fit],1)[0],relative_error_to940=Groi/G-1,fit_time_s=[a['time_s'][fit][0],a['time_s'][fit][-1]],fit_max_KE_U=a['KE_over_U'][fit].max(),fit_end_Dnum_over_U=a['D_num_J'][fit][-1]/a['U_bulk_J'][fit][-1],broken=int(a['broken_bulk'][-1]),max_contact_U_J=a['U_contact_J'].max(),states=res,available_full_fields=list(res),missing_full_fields=('Full fields only at declared accepted snapshots; native boundary reactions sampled every20 steps.' if periodic else 'No original loading-fit-window field states or every-step lower reaction. Reconstruction at saved accepted states only; no invented spatial interpolation.'))
 save(out/'SHEAR_AUDIT.json',result);return result

def notch(p,out):
 out.mkdir(parents=True,exist_ok=True);g=dict(np.load(p/'geometry.npz'));X=g['X'];b=g['bonds'];mid=X[b].mean(1);e=Energy(g,DX,ROOT/'analysis/material_macro_audit_v1/energy.so');reg,dtip,dgrip,band=regions(X,b);events=np.load(p/'break_events.npz')['events'];hist=np.genfromtxt(p/'history.csv',names=True,delimiter=',');plan=js(p/'test_plan.json');times=np.full(len(b),np.inf);times[events[:,1]]=events[:,0]*plan['dt_s'];names=['notch_tip','grip_transition','other_surface','free_interior'];snapshots={}
 cross=((X[b[:,0],2]<-1e-12)!=(X[b[:,1],2]<-1e-12));firststep=int(events[0,0]);peakstep=int(np.argmax(hist['stress_Pa']));firstids=events[events[:,0]==firststep,1]
 for key in ['initial','pre_first_break','first_break','peak','end']:
  if not (p/(key+'.npz')).exists():continue
  z=dict(np.load(p/(key+'.npz')));x=z['x_m'];alive=z['alive'].astype(bool);util=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)/e.r-1)/.15
  assert np.array_equal(~alive,times<=float(z['time_s'])+1e-15)
  snapshots[key]=dict(step=int(z['step']),time_s=float(z['time_s']),region_max_utilization={n:float(util[(reg==j)&(alive|(abs(times-float(z['time_s']))<1e-15))].max()) if np.any((reg==j)&(alive|(abs(times-float(z['time_s']))<1e-15))) else None for j,n in enumerate(names)},utilization_scope='Alive bonds and bonds newly broken at this exact accepted state; previously broken bonds excluded',broken_by_region={n:int(np.sum(~alive&(reg==j))) for j,n in enumerate(names)},broken_total=int((~alive).sum()),target_band_fraction=float(band[~alive].mean()) if (~alive).any() else None,off_projected_plane_fraction=float((~cross)[~alive].mean()) if (~alive).any() else None)
  np.savez_compressed(out/(key+'_risk.npz'),x=x,X=X,bonds=b,alive=alive,utilization=util,region=reg,midpoint=mid,tip_distance=dtip,grip_distance=dgrip,neighbor_count=np.bincount(b.ravel(),minlength=len(X)),weighted_volume=e.m,first_bond_indices=firstids,step=z['step'],time_s=z['time_s'])
 full,_=geometry(DX,'M4');bb=full['bonds'];a,c=X[bb[:,0]],X[bb[:,1]];cc=((a[:,2]<-1e-12)!=(c[:,2]<-1e-12));frac=np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(bb)),where=c[:,2]!=a[:,2]);xx=a[:,0]+frac*(c[:,0]-a[:,0]);cut=cc&(xx<-.0001);assert np.array_equal(bb[~cut],b);assert np.array_equal(bb[cut],g['removed_notch_bonds'])
 offevents=events[~cross[events[:,1]]];node_degree=np.bincount(b.ravel(),minlength=len(X));support=[]
 for k in firstids:
  support.append(dict(bond_index=int(k),atom_ids=(b[k]+1).tolist(),reference_midpoint_um=(mid[k]*1e6).tolist(),region=names[reg[k]],tip_distance_um=dtip[k]*1e6,grip_distance_um=dgrip[k]*1e6,endpoint_neighbors=node_degree[b[k]],endpoint_m_over_max=e.m[b[k]]/e.m.max()))
 # Event voxels are localization diagnostics, never a crack area estimator.
 voxel=np.floor((mid+np.array([.0003,.0003,.0005]))/DX+1e-8).astype(int);shape=(13,13,21);first_coherent=None
 for step in np.unique(events[:,0]):
  ids=events[events[:,0]<=step,1];grid=np.zeros(shape,bool);grid[tuple(voxel[ids].T)]=True;lab,nlab=label(grid);sizes=np.bincount(lab.ravel());eligible=np.flatnonzero(sizes[1:]>=10)+1
  if len(eligible):
   kept=ids[np.isin(lab[tuple(voxel[ids].T)],eligible)];nearfront=(abs(mid[kept,0]+.0001)<=DX)&(abs(mid[kept,2])<=DELTA/2)
   first_coherent=dict(step=int(step),time_s=step*plan['dt_s'],voxel_count_max=int(sizes[eligible].max()),intersects_tip_strip=bool(nearfront.any()),regions={n:int(np.sum(reg[kept]==j)) for j,n in enumerate(names)});break
 bytime=[]
 for step in np.unique(events[:,0]):
  ids=events[events[:,0]<=step,1];bytime.append([step*plan['dt_s'],len(ids),sum(band[ids]),sum(~cross[ids]),*[sum(reg[ids]==j) for j in range(4)]])
 np.savez_compressed(out/'EVENT_LOCALIZATION.npz',events=events,bond_first_break_time_s=times,region=reg,tip_distance_m=dtip,grip_distance_m=dgrip,midpoint_reference_m=mid,cross_plane=cross,target_band=band,cumulative=np.array(bytime),first_bond_indices=firstids)
 area=fracture_area(p,g,np.load(p/'end.npz'),events);save(out/'AREA_ATTEMPT.json',area)
 peak=snapshots['peak'];den=peak['broken_total'];tipfirst=bool(np.any(reg[firstids]==0));coherent=bool(first_coherent and first_coherent['intersects_tip_strip']);bandpass=bool(peak['target_band_fraction'] is not None and peak['target_band_fraction']>=.7);grippass=bool(den and peak['broken_by_region']['grip_transition']/den<=.2)
 result=dict(case=plan['case'],source=str(p),region_names=names,region_initial_bond_counts={n:int(sum(reg==j)) for j,n in enumerate(names)},removed_initial_notch_bonds=int(cut.sum()),unexpected_initial_cross_notch_bonds=0,first_event= support,first_step=firststep,first_time_s=firststep*plan['dt_s'],first_off_plane_time_s=float(offevents[0,0]*plan['dt_s']),notch_front_to_grip_transition_m=.00035,notch_front_to_grip_over_horizon=.00035/DELTA,snapshots=snapshots,first_coherent_voxel_cluster=first_coherent,notch_dominance_checks=dict(first_tip_or_coherent_tip=tipfirst or coherent,prepeak_target_band_ge70percent=bandpass,prepeak_grip_le20percent=grippass),area=area,fracture_work_status='GC_UNAVAILABLE',energy_residual_work_fraction=float(max(abs(hist['R_energy_J']))/max(abs(hist['W_ext_J']))),pre_first_break_field_status='ACTUAL_SAVED' if 'pre_first_break' in snapshots else 'UNAVAILABLE: no coordinates at step(first_break-1); first-break accepted positions reported without claiming pre-break history',first_break_KE_U_before=float(hist['KE_over_U'][firststep-1]))
 save(out/'NOTCH_AUDIT.json',result);return result

if __name__=='__main__':
 old=Path(sys.argv[1]);out=Path(sys.argv[2]);assert (out/'ANALYSIS_PLAN.json').exists()
 s=shear(old/'M2',out/'baseline/M2');n=notch(old/'M5',out/'baseline/M5');m=js(old/'ELASTIC_MAPPING.json')['measured'];s['M1_inferred_G_Pa']=m['E_Pa']/(2*(1+m['nu']));s['M1_inferred_vs_M2_relative']=s['G_ROI_Pa']/s['M1_inferred_G_Pa']-1
 save(out/'SHEAR_BOUNDARY_AUDIT.json',dict(baseline=s));save(out/'NOTCH_FAILURE_DIAGNOSTIC.json',dict(baseline=n));print(json.dumps(dict(Gnom=s['G_nominal_Pa'],Groi=s['G_ROI_Pa'],M1_G=s['M1_inferred_G_Pa'],M2_peak=s['states']['peak'],M5_first=n['first_event'],M5_risk=n['snapshots']['first_break'],M5_peak=n['snapshots']['peak'],first_coherent=n['first_coherent_voxel_cluster']),default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else x,indent=2))
