"""Review-only metadata and clarity additions after inspection; never invokes LAMMPS."""
from pathlib import Path
import sys,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from finish import plot_shear
from audit import save,js
from postprocess import figsave
base=Path(sys.argv[1]);old=Path(sys.argv[2]);figs=Path(sys.argv[3])
s=js(base/'SHEAR_BOUNDARY_AUDIT.json');s['periodic']['available_full_fields']=list(s['periodic']['states']);s['periodic']['missing_full_fields']='Full fields only at declared accepted snapshots; native boundary reactions sampled every20 steps.'
save(base/'SHEAR_BOUNDARY_AUDIT.json',s);save(base/'analysis/M2P/SHEAR_AUDIT.json',s['periodic'])
report_inputs=js(base/'REPORT_INPUTS.json');report_inputs['shear']=s;save(base/'REPORT_INPUTS.json',report_inputs)
h={n:np.genfromtxt(p/'history.csv',names=True,delimiter=',') for n,p in [('M2',old/'M2'),('M2P',base/'M2P'),('M5',old/'M5'),('M5B',base/'M5B')]}
plot_shear(base,figs,h,s)
# Full original scalar curves stay on server; an explicit compact accepted-row selection accompanies figures.
for name in ['M2','M5']:
 a=h[name];ix=np.unique(np.r_[np.arange(0,len(a),20),len(a)-1,np.argmax(a['stress_Pa']),np.flatnonzero(np.diff(a['broken_bulk'])>0)+1]);np.savez_compressed(base/'baseline'/name/'review_curves.npz',records=a[ix],selection='Every20 accepted steps plus all topology events, exact stress peak and endpoint; no interpolated state')
fig,axs=plt.subplots(1,2,figsize=(11,7),layout='constrained')
for ax,name,sub in zip(axs,['M5','M5B'],['baseline/M5','analysis/M5B']):
 z=np.load(base/sub/'first_break_risk.npz');X=z['X'];b=z['bonds'];alive=z['alive'].astype(bool);eligible=alive.copy();eligible[z['first_bond_indices']]=True;u=np.zeros(len(X));np.maximum.at(u,b[eligible].ravel(),np.repeat(z['utilization'][eligible],2))
 # Same actual central-slab bonds; node color is maximum over all y at each X,Z.
 pts=np.unique(X[:,[0,2]],axis=0);u2=np.array([u[np.all(X[:,[0,2]]==q,axis=1)].max() for q in pts]);cut=alive&(abs(X[b[:,0],1])<55e-6)&(abs(X[b[:,1],1])<55e-6)
 ax.add_collection(LineCollection((X[b[cut]][:,:,[0,2]])*1e3,colors='#aaa',linewidths=.4,alpha=.4,zorder=1));sc=ax.scatter(*(pts*1e3).T,c=u2,s=28,cmap='magma',vmin=0,vmax=1.05,zorder=2)
 for zz in [-.5,.35]:ax.add_patch(Rectangle((-.3,zz),.6,.15,color='#236c9d',alpha=.09,zorder=0))
 for zz in [-.35,.35]:ax.axhline(zz,c='#236c9d',ls='--',lw=1.4)
 ax.plot([-.3,-.1],[0,0],c='green',lw=3,label='Initial pre-cut plane',zorder=4);mid=X[b[z['first_bond_indices']]].mean(1)*1e3;ax.scatter(mid[:,0],mid[:,2],marker='*',s=180,c='cyan',edgecolors='black',zorder=5,label='First native breaks')
 ax.annotate('Notch tip',(-.1,0),(-.27,.12),arrowprops=dict(arrowstyle='->'),fontsize=11)
 ax.set(xlim=(-.33,.33),ylim=(-.53,.53),xlabel='Reference X (mm)',ylabel='Reference Z (mm)',title=f'{name}: first accepted break at {float(z["time_s"])*1e3:.4f} ms');ax.set_aspect('equal');ax.legend(loc='lower left',fontsize=9)
fig.colorbar(sc,ax=axs,shrink=.7,label='Max incident live/new-break utilization, projected over y');fig.suptitle('B3 | Orthographic audit of real nodes and central-slab surviving bonds\nBlue shading=prescribed axial collars; no rendered continuous shell',fontsize=12)
figsave(fig,figs,'B3_notch_and_grip_projection')
# Reference intersections, not broken-bond count as area. Keep original half-open bin convention.
p=base/'M5B';g=np.load(p/'geometry.npz');X=g['X'];b=g['bonds'];a,c=X[b[:,0]],X[b[:,1]];cross=(a[:,2]<-1e-12)!=(c[:,2]<-1e-12);at=a+np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(b)),where=abs(c[:,2]-a[:,2])>1e-15)[:,None]*(c-a);dx=50e-6
fig,axs=plt.subplots(1,3,figsize=(12,5),layout='constrained');arrays={}
for ax,label,file in zip(axs,['Before forced collar threshold','Nominal reaction peak','End of prescribed loading'],['frame_012000.npz','peak.npz','end.npz']):
 z=np.load(p/file);alive=z['alive'].astype(bool);total=np.zeros((8,12),int);surv=np.zeros_like(total)
 for i in range(8):
  for j in range(12):
   xx=-.0001+i*dx;yy=-.0003+j*dx;ids=cross&(at[:,0]>=xx-1e-12)&(at[:,0]<xx+dx-1e-12)&(at[:,1]>=yy-1e-12)&(at[:,1]<yy+dx-1e-12);total[i,j]=ids.sum();surv[i,j]=(ids&alive).sum()
 ratio=np.divide(surv,total,out=np.zeros_like(surv,dtype=float),where=total>0);im=ax.imshow(ratio.T,origin='lower',extent=[-.1,.3,-.3,.3],vmin=0,vmax=1,cmap='viridis',interpolation='nearest');ax.set(xlabel='Reference X (mm)',ylabel='Reference Y (mm)',title=f'M5B {float(z["time_s"])*1e3:.3f} ms\n{label}');arrays[file]=dict(total=total,surviving=surv,time_s=float(z['time_s']))
fig.colorbar(im,ax=axs,shrink=.75,label='Surviving / initial crossing bonds in each 50µm bin');fig.suptitle('D2 | Original z=0 ligament: crossing-connectivity map, NOT physical crack area\nZero crossing connectivity alone does not establish a unique crack surface; Gc unavailable',fontsize=12)
figsave(fig,figs,'D2_crossing_connectivity_limits');np.savez_compressed(base/'analysis/M5B/AREA_GRID_SOURCE.npz',**{file.replace('.npz','')+'_'+key:val for file,d in arrays.items() for key,val in d.items()})
# Additional dynamic/energy and exact event onset diagnostics.
out={}
for name in ['M5','M5B']:
 a=h[name];npeak=int(np.argmax(a['stress_Pa']));pre=a['time_s']<6.8423750062e-3
 out[name]=dict(max_KE_U_full=float(a['KE_over_U'].max()),max_KE_U_before_collar_limit=float(a['KE_over_U'][pre].max()),energy_at_peak={k:float(a[k][npeak]) for k in ['W_ext_J','U_bulk_J','U_contact_J','K_free_J','D_num_J','D_break_independent_J','R_energy_J']},energy_end={k:float(a[k][-1]) for k in ['W_ext_J','U_bulk_J','U_contact_J','K_free_J','D_num_J','D_break_independent_J','R_energy_J']})
for name in ['M5','M5B']:
 a=h[name];first=np.flatnonzero(a['broken_bulk']>0)[0];ix=np.arange(first,len(a));k=int(ix[np.argmax(a['KE_over_U'][ix])]);i=int(np.argmax(a['KE_over_U']));out[name]['full_max_time_s']=float(a['time_s'][i]);out[name]['max_after_first_break']=dict(time_s=float(a['time_s'][k]),KE_U=float(a['KE_over_U'][k]),K_J=float(a['K_free_J'][k]),U_bulk_J=float(a['U_bulk_J'][k]))
save(base/'DYNAMICS_AND_ENERGY_DETAIL.json',out)
print('Review additions completed; no solver run.')
