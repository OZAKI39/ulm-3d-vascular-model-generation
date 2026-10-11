"""Audit completed accepted states; no solver reruns and no retrospective fit tuning."""
import sys,json,csv,math,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
K=100000/3;G=940.;E=9*K*G/(3*K+G);NU=(3*K-2*G)/(2*(3*K+G))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,'axes.labelsize':11,'legend.fontsize':9,'axes.grid':True,'grid.color':'#9099a0','grid.alpha':.55,'figure.facecolor':'white','savefig.facecolor':'white'})
COLORS={'M1':'#236c9d','M2':'#cf8a28','M3':'#6c4c9a','M4':'#b33132','M5':'#16818b','M6':'#cf6227','T1':'#bd4b80','D1':'#4354a6'}
def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float) and not math.isfinite(x):return None
 return x
def js(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(clean(x),indent=2,ensure_ascii=False))
def figsave(fig,figs,name,svg=True):
 fig.savefig(figs/(name+'.png'),dpi=220,bbox_inches='tight')
 if svg:fig.savefig(figs/(name+'.svg'),bbox_inches='tight')
 plt.close(fig)

def fracture_area(p,g,end,events):
 X=g['X'];b=g['bonds'];alive=end['alive'].astype(bool);a,c=X[b[:,0]],X[b[:,1]]
 cross=((a[:,2]<-1e-12)!=(c[:,2]<-1e-12));broken=~alive
 where=np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(b)),where=abs(c[:,2]-a[:,2])>1e-15)
 at=a+where[:,None]*(c-a);dx=js(p/'test_plan.json')['material']['dx_m'];xs=np.arange(-.0001,.0003-dx/2,dx);ys=np.arange(-.0003,.0003-dx/2,dx)
 lost=np.zeros((len(xs),len(ys)),bool);occup=np.zeros_like(lost)
 for i,x in enumerate(xs):
  for j,y in enumerate(ys):
   ids=cross&(at[:,0]>=x-1e-12)&(at[:,0]<x+dx-1e-12)&(at[:,1]>=y-1e-12)&(at[:,1]<y+dx-1e-12)
   occup[i,j]=ids.any();lost[i,j]=ids.any() and not alive[ids].any()
 front=[]
 for j in range(len(ys)):
  count=0
  for i in range(len(xs)):
   if not lost[i,j]:break
   count+=1
  front.append(count*dx)
 return dict(convention='One projected plane z=0, thickness600um x newly advanced length; never bond_count*dx^2',method='Reference crossing intersections rasterized in dx cells; contiguous rows with zero surviving crossing bonds from original notch tip x=-100um; conservative connectivity diagnostic, not physical surface reconstruction',resolution_m=dx,notch_m=.0002,total_new_breaks=int(broken.sum()),cross_plane_breaks=int((broken&cross).sum()),off_plane_break_fraction=float((broken&~cross).sum()/max(broken.sum(),1)),front_advance_per_thickness_strip_m=front,projected_connectivity_loss_area_m2=float(sum(front)*dx),fully_lost_cross_plane_bins=int(lost.sum()),area_identification='NOT_IDENTIFIABLE_AS_SINGLE_CRACK' if ((broken&~cross).sum()/max(broken.sum(),1)>.2 or np.ptp(front)>2*dx) else 'CONNECTIVITY_PROXY_ONLY_NOT_VALIDATED_SURFACE',Gc_J_m2=None)

def material(base,figs):
 rows=[];series={};details={}
 for name in ['M1','M2','M3','M3P','M4','M5','M6','T1','D1']:
  p=base/name;assert js(p/'EXECUTION.json')['status']=='COMPLETED';assert not list(p.glob('*ERROR*'))
  a=np.genfromtxt(p/'history.csv',delimiter=',',names=True);plan=js(p/'test_plan.json');diag=js(p/'DIAGNOSTIC_CHECKS.json');series[name]=a
  assert len(a)==plan['steps']+1 and np.array_equal(a['step'],np.arange(plan['steps']+1));assert np.isfinite(np.column_stack([a[k] for k in a.dtype.names])).all()
  assert np.max(abs(a['time_s']-a['step']*plan['dt_s']))<1e-14
  row=dict(case=name,node_count=plan['nodes'],dx_um=plan['material']['dx_m']*1e6,dt_us=plan['dt_s']*1e6,duration_ms=plan['duration_s']*1e3,wall_seconds=diag['wall_seconds'],E_Pa=None,G_Pa=None,nu=None,M_Pa=None,K_Pa=None,peak_stress_Pa=None,peak_boundary_strain=None,peak_time_s=None,peak_observed=False,first_break_s=None,first_break_type=None,broken_bulk=int(a['broken_bulk'][-1]),broken_interface=int(a['broken_interface'][-1]),energy_residual_work_fraction=float(max(abs(a['R_energy_J']))/max(abs(a['W_ext_J']))),signed_topology_energy_J=float(a['D_break_independent_J'][-1]),positive_topology_energy_J=float(np.maximum(np.diff(a['D_break_independent_J']),0).sum()),negative_topology_energy_J=float(np.minimum(np.diff(a['D_break_independent_J']),0).sum()),negative_topology_events=len(diag['negative_topology_releases']),quasistatic_scope='Fit-window only; full startup and fracture inertia retained',fit_max_KE_U=None,full_max_KE_U=float(a['KE_over_U'].max()),Gc_J_m2=None)
  row.update(end_max_speed_m_s=float(a['max_speed_m_s'][-1]),end_K_free_J=float(a['K_free_J'][-1]),end_U_bulk_J=float(a['U_bulk_J'][-1]),max_U_contact_J=float(a['U_contact_J'].max()))
  if name in ['M1','M2','M3','M3P']:
   frac=abs(a['boundary_strain'])/.008;fit=(a['time_s']<=plan['ramp_s'])&(frac>=.2)&(frac<=.7);key='gamma_xz' if name=='M2' else 'eps_z';slope=np.polyfit(a[key][fit],a['stress_Pa'][fit],1)[0]
   row[{'M1':'E_Pa','M2':'G_Pa','M3':'M_Pa','M3P':'M_Pa'}[name]]=float(slope)
   if name=='M1':row['nu']=float(-np.polyfit(a['eps_z'][fit],.5*(a['eps_x'][fit]+a['eps_y'][fit]),1)[0])
   row['fit_max_KE_U']=float(a['KE_over_U'][fit].max());row['fit_time_start_s']=float(a['time_s'][fit][0]);row['fit_time_end_s']=float(a['time_s'][fit][-1]);row['fit_transverse_to_axial_strain_max']=float(np.max(abs(a['eps_x'][fit])+abs(a['eps_y'][fit]))/max(abs(a['eps_z'][fit])))
   assert row['broken_bulk']==0
  else:
   peak=int(np.argmax(a['stress_Pa']));row.update(peak_stress_Pa=float(a['stress_Pa'][peak]),peak_boundary_strain=float(a['boundary_strain'][peak]),peak_time_s=float(a['time_s'][peak]),peak_observed=bool(peak<len(a)-1 and a['stress_Pa'][-1]<.9*a['stress_Pa'][peak]))
   ev=np.load(p/'break_events.npz')['events']
   if len(ev):
    first_step=int(ev[0,0]);row.update(first_break_s=float(ev[0,0]*plan['dt_s']),first_break_type=int(ev[0,2]),KE_over_U_immediately_before_first_break=float(a['KE_over_U'][first_step-1]),maximum_unbroken_KE_U_including_startup=float(a['KE_over_U'][:first_step].max()))
   if len(ev):
    gg=np.load(p/'geometry.npz');first_bond=gg['bonds'][ev[0,1]];mid=gg['X'][first_bond].mean(0);current=np.load(p/'first_break.npz')['x_m'][first_bond].mean(0)
    row.update(first_break_atom_ids=(first_bond+1).tolist(),first_break_reference_midpoint_um=(mid*1e6).tolist(),first_break_current_midpoint_um=(current*1e6).tolist(),first_break_near_grip_transition=bool(abs(mid[2])>=.00035-150.5e-6))
   if name in ['M5','T1','D1']:
    area=fracture_area(p,np.load(p/'geometry.npz'),np.load(p/'end.npz'),ev);save(p/'CRACK_AREA_ATTEMPT.json',area);row['off_plane_break_fraction']=area['off_plane_break_fraction'];row['fracture_energy_status']='SIGNED_TOPOLOGY_TRANSFER_IDENTIFIED; NONNEGATIVE_Gc_NOT_IDENTIFIABLE'
   if name=='M6':
    g=np.load(p/'geometry.npz');initial=g['bonds'][g['bond_types']==2,0];alive=np.load(p/'end.npz')['alive'];surv=g['bonds'][(g['bond_types']==2)&alive.astype(bool),0];lost=np.setdiff1d(np.unique(initial),np.unique(surv));row['basal_interface_nodes_fully_debonded']=len(lost);row['basal_reference_patch_area_m2']=float(len(lost)*(plan['material']['dx_m']**2));row['interface_area_scope']='Assigned basal nodal patch area losing every anchor bond; not a fracture surface area from bond counts';row['failure_mode']='BULK_FIRST_OR_MIXED' if row['first_break_type']==1 else 'INTERFACE_FIRST_WITH_BULK_DAMAGE' if row['broken_bulk'] else 'INTERFACE'
  details[name]=dict(plan=plan,diagnostics=diag,metrics=row);rows.append(row)
 r={v['case']:v for v in rows};r['M3']['K_naive_formula_Pa']=r['M3']['M_Pa']-4*r['M2']['G_Pa']/3
 r['M3']['K_identification']='NOT_IDENTIFIED: central transverse strain is20.55% of axial'
 assert r['M3P']['fit_transverse_to_axial_strain_max']<1e-6
 r['M3P']['K_Pa']=r['M3P']['M_Pa']-4*r['M2']['G_Pa']/3
 r['M3P']['K_using_input_G_Pa']=r['M3P']['M_Pa']-4*G/3
 mapping=dict(target=dict(E_Pa=E,G_Pa=G,nu=NU,K_Pa=K),measured=dict(E_Pa=r['M1']['E_Pa'],G_Pa=r['M2']['G_Pa'],nu=r['M1']['nu'],K_Pa=r['M3P']['K_Pa']),relative_errors=dict(E=(r['M1']['E_Pa']/E-1),G=(r['M2']['G_Pa']/G-1),K=(r['M3P']['K_Pa']/K-1)),nu_absolute_error=abs(r['M1']['nu']-NU),scope='Finite free-interior specimens, nominal grip reaction / initial area versus central fitted strain; surface/grip effects included; not uniform continuum calibration')
 mapping['matches']={k:abs(v)<=.15 for k,v in mapping['relative_errors'].items()};mapping['matches']['nu']=mapping['nu_absolute_error']<=.02
 comparison={}
 for name,limit in [('T1',.10),('D1',.15)]:
  aa=series['M5'];bb=series[name];interp=np.interp(aa['time_s'],bb['time_s'],bb['stress_Pa']);norm=np.linalg.norm(aa['stress_Pa']);curve=float(np.linalg.norm(interp-aa['stress_Pa'])/norm)
  metrics={key:abs(r[name][key]/r['M5'][key]-1) for key in ['peak_stress_Pa','signed_topology_energy_J','first_break_s']}
  comparison[name]=dict(relative_differences=metrics,time_aligned_stress_relative_L2=curve,interpolation='Linear analysis-only alignment on M5 physical times; displayed states never interpolated',criterion=limit,main_metrics_pass=all(v<=limit for v in metrics.values()) and curve<=limit,first_failure_type_matches=r[name]['first_break_type']==r['M5']['first_break_type'],same_duration=details[name]['plan']['duration_s']==details['M5']['plan']['duration_s'],same_amplitude_relative_error=abs(details[name]['plan']['amplitude_m']/details['M5']['plan']['amplitude_m']-1),mesh_scope='Fixed nominal dimensions/horizon/notch/loading amplitude; nodal surface/grip sampling and contact distance vary with dx; two levels only')
 save(base/'ELASTIC_MAPPING.json',mapping);save(base/'MATERIAL_COMPARISONS.json',comparison);save(base/'MATERIAL_DETAILS.json',details)
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (base/'MACRO_TEST_RESULTS.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
 for j,name in enumerate(['M1','M2','M3']):
  a=series[name];key='gamma_xz' if name=='M2' else 'eps_z';ax=axs[0,j];ax.plot(a[key]*100,a['stress_Pa'],color=COLORS[name],lw=1.5,label='Accepted load-unload');frac=abs(a['boundary_strain'])/.008;fit=(a['time_s']<=.010)&(frac>=.2)&(frac<=.7);ax.plot(a[key][fit]*100,a['stress_Pa'][fit],color='#e4a53d',lw=3,label='Fixed fit: 20–70% load');target=[E,G,K+4*G/3][j];xx=np.array([a[key].min(),a[key].max()]);ax.plot(xx*100,target*xx,'--',color='black',label='Continuum target');ax.set(xlabel=('Central shear strain (%)' if name=='M2' else 'Central axial strain (%)'),ylabel='Nominal stress (Pa)',title=f'{name}: '+['E = %.1f Pa'%r['M1']['E_Pa'],'G = %.1f Pa'%r['M2']['G_Pa'],'M = %.1f Pa'%r['M3']['M_Pa']][j]);ax.legend()
  ax=axs[1,j];ax.plot(a['time_s']*1e3,a['KE_over_U']*100,color=COLORS[name]);ax.axhline(5,color='black',ls='--');ax.set(xlabel='Physical time (ms)',ylabel='K / U (%)',title='Full startup and unloading retained');ax.set_yscale('symlog',linthresh=1)
 fig.suptitle('Fixed LPS: G=940 Pa, K=33.33 kPa; free interiors; 500 s⁻¹ numerical damping');figsave(fig,figs,'01_elastic_response')
 fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
 a=series['M1'];axs[0].plot(a['eps_z']*100,a['eps_x']*100,label='x transverse');axs[0].plot(a['eps_z']*100,a['eps_y']*100,'--',label='y transverse');axs[0].set(xlabel='Central axial strain (%)',ylabel='Transverse strain (%)',title=f'M1 ν = {r["M1"]["nu"]:.4f} (target {NU:.4f})');axs[0].legend()
 for name in ['M4','M5']:a=series[name];axs[1].plot(a['boundary_strain']*100,a['stress_Pa'],label=name,color=COLORS[name]);axs[1].scatter(r[name]['peak_boundary_strain']*100,r[name]['peak_stress_Pa'],s=24,color=COLORS[name])
 axs[1].set(xlabel='Grip engineering strain (%)',ylabel='Nominal tensile stress (Pa)',title='New free-interior specimens');axs[1].legend()
 for name in ['M5','T1','D1']:a=series[name];axs[2].plot(a['time_s']*1e3,a['stress_Pa'],label=name,color=COLORS[name],lw=1)
 axs[2].set(xlabel='Physical time (ms)',ylabel='Notched nominal stress (Pa)',title='T1: dt/2; D1: 40 µm / fixed horizon');axs[2].legend();figsave(fig,figs,'02_strength_and_controls')
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for ax,name in zip(axs.ravel(),['M1','M2','M5','M6']):
  a=series[name]
  for key,label in [('W_ext_J','Boundary work'),('U_bulk_J','Source-consistent bulk U'),('U_contact_J','Conservative contact U'),('K_free_J','Free-DOF K'),('D_num_J','Numerical damping'),('D_break_independent_J','Signed topology transfer'),('R_energy_J','Independent residual')]:ax.plot(a['time_s']*1e3,a[key]*1e9,label=label,lw=1)
  ax.plot(a['time_s']*1e3,a['native_bond_pe_times_V_J']*1e9,'--',color='gray',label='Native bond pe × V (incomplete)');ax.set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title=name+'; Gc unavailable' if name in ['M5','M6'] else name+' no fracture');ax.legend(fontsize=8,ncol=2)
 figsave(fig,figs,'03_energy_ledger')
 fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained');a=series['M6'];axs[0].plot(a['displacement_m']*1e6,a['stress_Pa'],color=COLORS['M6']);axs[0].set(xlabel='Top grip opening (µm)',ylabel='Nominal traction (Pa)',title='M6: 3 fixed anchor layers')
 for key,label in [('broken_bulk','Bulk'),('broken_interface','Interface')]:axs[1].plot(a['time_s']*1e3,a[key],label=label)
 axs[1].legend();axs[1].set(xlabel='Physical time (ms)',ylabel='Actual broken bonds',title=r['M6']['failure_mode'])
 a=series['M5'];axs[2].plot(a['time_s']*1e3,a['D_break_independent_J']*1e9,label='Net topology transfer');neg=np.diff(a['D_break_independent_J'])<0;axs[2].scatter(a['time_s'][1:][neg]*1e3,a['D_break_independent_J'][1:][neg]*1e9,s=7,color='crimson',label='Negative release event');axs[2].legend();axs[2].set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title='M5: cannot call all releases dissipation');figsave(fig,figs,'04_interface_and_energy_limits')
 # All real nodes; one-plane cutaway bond display to prevent an opaque line mass.
 for name in ['M4','M5','M6']:
  p=base/name;g=np.load(p/'geometry.npz');X=g['X'];b=g['bonds'];anc=g['anchors'];roi=g['roi'];types=g['bond_types'];fig=plt.figure(figsize=(17,6));fig.subplots_adjust(left=.02,right=.98,bottom=.10,top=.80,wspace=.12)
  snapshots={lab:dict(np.load(p/(lab+'.npz'))) for lab in ['initial','first_break','peak','end']}
  allx=np.concatenate([snap['x_m'] for snap in snapshots.values()])*1e3;lo=allx.min(0);hi=allx.max(0);pad=.04*(hi-lo);lo-=pad;hi+=pad
  ev=np.load(p/'break_events.npz')['events'];first_nodes=b[int(ev[0,1])] if len(ev) else np.array([],int)
  for j,label in enumerate(['initial','first_break','peak','end']):
   ax=fig.add_subplot(1,4,j+1,projection='3d');snap=snapshots[label];x=snap['x_m']*1e3;alive=snap['alive'].astype(bool)
   degree=np.bincount(b[alive].ravel(),minlength=len(x));isolated=degree==0;ax.scatter(*x[~anc&~isolated].T,s=6,c='#b43c37',alpha=.75,depthshade=False);ax.scatter(*x[anc].T,s=6,c='#3378aa',depthshade=False)
   if isolated.any():ax.scatter(*x[isolated].T,s=8,c='#7b858e',depthshade=False)
   cut=(abs(X[b[:,0],1])<=55e-6)&(abs(X[b[:,1],1])<=55e-6)&alive
   seg=x[b[cut]];cs=np.where(types[cut]==2,'#c99b35',np.where(types[cut]==3,'#638bb0','#b3a5a1'));ax.add_collection3d(Line3DCollection(seg,colors=cs,linewidths=.25,alpha=.35))
   if label=='first_break':ax.scatter(*x[first_nodes].T,s=45,c='#e1b42d',edgecolors='black',linewidths=.5,depthshade=False)
   ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=f'{name} {label.replace("_"," ")}\nt={float(snap["time_s"])*1e3:.3f} ms',xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]),xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.view_init(elev=16,azim=-55);ax.set_box_aspect(hi-lo);ax.tick_params(labelsize=9);ax.grid(True)
  fig.suptitle('Actual accepted states; common physical axis limits; surviving bonds in central y slab\nRed=clot, blue=anchors, gray=isolated, gold=first-break nodes; no shell',y=.98,fontsize=13);figsave(fig,figs,'05_states_'+name,False)
 # Boundary/measurement geometry from actual registered node masks.
 fig=plt.figure(figsize=(13,4.8),layout='constrained')
 for j,name in enumerate(['M1','M2','M3']):
  geo=np.load(base/name/'geometry.npz');x=geo['X']*1e3;mask=geo['prescribed_mask'].any(1);roi=geo['roi'];ax=fig.add_subplot(1,3,j+1,projection='3d')
  ax.scatter(*x[~mask].T,s=5,color='#b43c37',alpha=.6,label='Free interior')
  ax.scatter(*x[mask].T,s=5,color='#3e78a4',alpha=.6,label='Boundary: selected DOFs')
  ax.scatter(*x[roi].T,s=9,facecolors='none',edgecolors='#71409a',linewidths=.5,label='Central measurement ROI')
  ax.view_init(elev=16,azim=-55);ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=name+': '+['axial-only grips','x/z grip motion; free interior','lateral boundary normals fixed'][j]);ax.set_box_aspect([.6,.6,1]);ax.tick_params(labelsize=8)
  if j==1:ax.legend(loc='upper center',bbox_to_anchor=(.5,-.14),fontsize=9)
 fig.suptitle('600 × 600 × 1000 µm; dx=50 µm; horizon=150.5 µm; selected boundary DOFs only')
 figsave(fig,figs,'08_actual_geometry_and_ROI',False)
 fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
 for name,color in [('M3','#888888'),('M3P','#226ba1')]:
  a=series[name];axs[0].plot(a['eps_z']*100,a['stress_Pa'],color=color,label=name)
  axs[1].plot(a['time_s']*1e3,a['eps_x']*100,color=color,label=name+' x');axs[1].plot(a['time_s']*1e3,a['eps_y']*100,'--',color=color,label=name+' y')
  axs[2].plot(a['time_s']*1e3,a['R_energy_J']*1e9,color=color,label=name)
 xx=np.array([series['M3P']['eps_z'].min(),0.]);axs[0].plot(xx*100,(K+4*G/3)*xx,'--',color='black',label='Continuum M target')
 axs[0].set(xlabel='Central axial strain (%)',ylabel='Nominal stress (Pa)',title=f'M3P: K = {r["M3P"]["K_Pa"] / 1000:.3f} kPa');axs[0].legend()
 axs[1].set(xlabel='Physical time (ms)',ylabel='Central lateral strain (%)',title='M3 is not locally uniaxial strain');axs[1].legend()
 axs[2].set(xlabel='Physical time (ms)',ylabel='Independent residual (nJ)',title='M3P: K estimate after boundary check');axs[2].legend()
 fig.suptitle('Registered diagnostic9: no parameter tuning; only transverse boundary condition changed')
 figsave(fig,figs,'10_constrained_compression_diagnostic')
 return rows,mapping,comparison

if __name__=='__main__':
 base=Path(sys.argv[1]);figs=Path(sys.argv[2]);figs.mkdir(parents=True,exist_ok=True);material(base,figs)
