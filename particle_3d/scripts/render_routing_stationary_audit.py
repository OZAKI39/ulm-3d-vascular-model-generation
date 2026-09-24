#!/usr/bin/env python3
"""Eight academic figures, generated exclusively from saved audit evidence."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from particle_3d.routing_stationary_audit import read,dump,ROLES,sha
from particle_3d.particle8_replay import REPO
R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data';F=R/'figures';F.mkdir(exist_ok=True)
COLORS=['#3479b5','#e69f00','#169b82','#8a8a8a'];LABELS=['O1','O2','O3','No exit']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})
manifest=[]
def save(fig,name,sources):
 fig.savefig(F/(name+'.png'),dpi=230,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
 manifest.append(dict(figure=name,sources=sources,png_sha256=sha(F/(name+'.png')),pdf_sha256=sha(F/(name+'.pdf'))))
def axes3(ax,points):
 lo=points.min(0);hi=points.max(0);center=(lo+hi)/2;half=(hi-lo).max()*.53
 ax.set(xlim=(center[0]-half,center[0]+half),ylim=(center[1]-half,center[1]+half),zlim=(center[2]-half,center[2]+half),xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)');ax.set_box_aspect([1,1,1]);ax.view_init(24,-62)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.pane.fill=False;a.set_major_locator(plt.MaxNLocator(4))

def main():
 s=read(D/'audit_summary.json');raw=read(D/'point_completed.json')['results'];st=read(D/'stationary_audit.json');geo=read(D/'plot_geometry.json');wall=np.array(geo['wall_triangles_m'])*1e6;tri=np.array(geo['inlet_triangles_m'])*1e6
 anchors=np.array([r['anchor_m'] for r in raw])*1e6;center=tri.reshape(-1,3).mean(0);_,_,vt=np.linalg.svd(tri.reshape(-1,3)-center,full_matrices=False);basis=vt[:2].T;xy=(anchors-center)@basis;polys=(tri-center)@basis
 dump(D/'inlet_projection.json',dict(origin_um=center,basis=basis,axis_role='ORTHONORMAL_INLET_PLANE_COORDINATES_FOR_DISPLAY_ONLY'))
 fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
 for ax in axs:
  ax.add_collection(PolyCollection(polys,facecolor='#eef2f5',edgecolor='#bbc3cb',linewidth=.3));ax.autoscale();ax.set_aspect('equal');ax.set(xlabel='Inlet coordinate 1 (µm)',ylabel='Inlet coordinate 2 (µm)')
 for k,role in enumerate(ROLES):
  mask=np.array([r['point_outlet']==role for r in raw]);axs[0].scatter(*xy[mask].T,s=13,c=COLORS[k],label=f'{LABELS[k]}: {mask.sum()}',alpha=.8,linewidths=0)
 accepted=np.array([r['accepted'] for r in raw]);axs[1].scatter(*xy[~accepted].T,s=13,marker='x',c='#b04f50',linewidth=.7,label=f'Rejected: {(~accepted).sum()}');axs[1].scatter(*xy[accepted].T,s=12,c='#246e9b',alpha=.8,linewidth=0,label=f'Accepted: {accepted.sum()}')
 axs[0].set_title('(a) Point-tracer outlet basins');axs[1].set_title('(b) Finite-size admission');axs[0].legend(fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,frameon=False);axs[1].legend(fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,frameon=False)
 save(fig,'01_inlet_basin_and_admission',['point_completed.json','plot_geometry.json','inlet_projection.json'])
 fig,axs=plt.subplots(1,2,figsize=(11.5,5),layout='constrained')
 for ax,key,title in zip(axs,['point_to_p65_transition','p65_to_p9a1_transition'],['(a) Point tracer → P6.5','(b) P6.5 → P9-A.1']):
  m=np.array(s[key]['counts']);fra=np.array(s[key]['fraction_of_source_basin']);ax.imshow(np.sqrt(m),cmap='Blues',vmin=0,vmax=np.sqrt(max(1,m.max())))
  for i in range(4):
   for j in range(4):ax.text(j,i,f'{m[i,j]}\n({fra[i,j]*100:.1f}%)',ha='center',va='center',fontsize=10,color='white' if m[i,j]>m.max()*.25 else '#17232d')
  ax.set(xticks=range(4),xticklabels=LABELS,yticks=range(4),yticklabels=LABELS,xlabel='Destination stage',ylabel='Source stage',title=title)
 fig.supxlabel('Cell: particle count (percentage of the source row); all 500 particles retained',fontsize=10)
 save(fig,'02_routing_transition',['paired_routing.json','audit_summary.json'])
 loss=s['outlet01_loss_chain'];fig,ax=plt.subplots(figsize=(10.8,4.2));ax.axis('off')
 values=[loss['raw_O1'],loss['accepted_O1'],loss['P65_keeps_O1'],loss['P9_keeps_O1']];names=['Raw O1 point basin','Admitted O1 seeds','P6.5 remains at O1','P9-A.1 remains at O1']
 for i,(value,name) in enumerate(zip(values,names)):
  x=.12+i*.255;ax.text(x,.63,f'{value}',fontsize=29,weight='bold',ha='center',transform=ax.transAxes,color='#275f85');ax.text(x,.43,name,ha='center',transform=ax.transAxes,fontsize=10)
  if i<3:ax.annotate('',xy=(x+.18,.64),xytext=(x+.075,.64),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.5,color='#6d7782'))
 ax.text(.245,.2,f"{loss['admission_rejected']} rejected\nat admission",ha='center',transform=ax.transAxes)
 ax.text(.5,.2,f"{loss['accepted_P65_reroutes']} rerouted; {loss['accepted_P65_stops']} stopped\nby P6.5",ha='center',transform=ax.transAxes)
 ax.text(.755,.2,f"{loss['P65_keeps_then_P9_reroutes']} additionally rerouted\nwith P9-A.1",ha='center',transform=ax.transAxes)
 ax.set_title('Outlet 01 loss chain | counts from the original O1 point basin',pad=20)
 fig.text(.5,.02,'Admission rejection and “not admitted” refer to the same particles and are counted once.',ha='center',fontsize=9)
 save(fig,'03_outlet01_loss_chain',['outlet01_loss_chain.json','paired_routing.json'])
 fig=plt.figure(figsize=(11,6.3));ax=fig.add_subplot(121,projection='3d',computed_zorder=False);ax.add_collection3d(Poly3DCollection(wall,facecolor='#718da2',edgecolor='none',alpha=.38,rasterized=True,zorder=1));axes3(ax,wall.reshape(-1,3));palette=plt.get_cmap('tab10')
 for h in s['stationary_hotspots']:
  items=[r for r in st if r['hotspot_id']==h['hotspot_id']];p=np.array([x['final_position_m'] for x in items])*1e6
  ax.scatter(*p.T,s=[55*(x['radius_m']*1e6)**2 for x in items],c=[palette((h['hotspot_id']-1)%10)],edgecolors='white',linewidth=.6,depthshade=False,zorder=3,label=f"H{h['hotspot_id']} (n={len(items)})")
 ax.legend(loc='upper left',fontsize=9);ax.set_title('(a) Stationary locations in the true vessel',fontsize=11)
 ax=fig.add_subplot(122);yy=np.arange(len(st));ax.scatter([x['radius_m']*1e6 for x in st],yy,c=[palette((x['hotspot_id']-1)%10) for x in st],s=60)
 ax.set(yticks=yy,yticklabels=[f"ID {x['particle_id']} | H{x['hotspot_id']}" for x in st],xlabel='Radius (µm)',title='(b) All 14 stationary particles');ax.invert_yaxis();ax.grid(axis='x',alpha=.18)
 fig.tight_layout();save(fig,'04_stationary_hotspots',['stationary_audit.json','plot_geometry.json','audit_summary.json'])
 sizes=list(csv.DictReader((D/'particle_sizes.csv').open()));finished=np.array([float(x['radius_um']) for x in sizes if x['stationary']=='False']);stopped=np.array([x['radius_m']*1e6 for x in st]);fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
 ax=axs[0];ax.boxplot([finished,stopped],tick_labels=['Completed (486)','Stationary (14)'],showfliers=False,widths=.5,medianprops=dict(color='#151f26',linewidth=2));ax.scatter(1+.13*np.sin(np.arange(len(finished))*2.399963229728653),np.sort(finished),s=7,c='#3479b5',alpha=.3,linewidth=0);ax.scatter(2+np.linspace(-.1,.1,len(stopped)),np.sort(stopped),s=26,c='#c44e52',linewidth=0);ax.set(ylabel='Radius (µm)',title='(a) Accepted-particle size distributions');ax.grid(axis='y',alpha=.18)
 ax=axs[1];ax.scatter([r['particle_id'] for r in st],[r['size_percentile'] for r in st],c='#c44e52',s=45);ax.set(xlabel='Stationary particle ID',ylabel='Percentile in accepted 500 (%)',title='(b) Position in the full size distribution',ylim=(95,100.3));ax.grid(alpha=.18)
 save(fig,'05_radius_completed_vs_stationary',['particle_sizes.csv','stationary_audit.json','audit_summary.json'])
 selected=[]
 for desc in ['ROBUST_RIGID_SIZE_EXCLUSION','HIGHLY_RADIUS_SENSITIVE','UNRESOLVED']:
  group=[r for r in st if r['radius_sensitivity']==desc]
  if group:selected.append(group[0])
 if not selected:selected=[st[0]]
 if not any(x['radius_sensitivity']=='ROBUST_RIGID_SIZE_EXCLUSION' for x in selected):
  group=[r for r in st if r['radius_sensitivity']=='INTERMEDIATE']
  if group and len(selected)<3:selected.insert(0,group[0])
 fig=plt.figure(figsize=(5.3*len(selected),5.9));examples=[]
 for k,r in enumerate(selected):
  ax=fig.add_subplot(1,len(selected),k+1,projection='3d',computed_zorder=False);c=np.array(r['final_position_m'])*1e6;a=r['radius_m']*1e6;relative=wall-c;mask=np.min(np.linalg.norm(relative,axis=2),axis=1)<3*a;local=relative[mask]
  ax.add_collection3d(Poly3DCollection(local,facecolor='#d4dfe7',edgecolor='#b4c0ca',linewidth=.15,alpha=.16,rasterized=True,zorder=1))
  u,v=np.mgrid[0:2*np.pi:40j,0:np.pi:25j];ax.plot_surface(a*np.cos(u)*np.sin(v),a*np.sin(u)*np.sin(v),a*np.cos(v),color='#e4a438',alpha=.35,linewidth=0,rasterized=True,zorder=2)
  kept=r['contact_redundancy']['kept_constraint_ids'];contact=[x for x in r['contacts'] if x['canonical_id'] in kept]
  for contact_row in contact:
   p=np.array(contact_row['wall_point_m'])*1e6-c;n=np.array(contact_row['normal']);ax.scatter(*p,c='#b3363d',s=23,zorder=12);ax.quiver(*p,*(n*a*.8),color='#b3363d',linewidth=2,arrow_length_ratio=.25,zorder=11)
  flow=np.array(r['free_fem_velocity_m_s']);flow/=np.linalg.norm(flow);ax.quiver(0,0,0,*(flow*a*2),color='#236fac',linewidth=2.4,arrow_length_ratio=.2,zorder=13)
  axes3(ax,np.array([[-2*a]*3,[2*a]*3]));ax.set(xlabel='Δx (µm)',ylabel='Δy (µm)',zlabel='Δz (µm)');ax.set_title(f"ID {r['particle_id']} | radius {a:.3f} µm\n{r['radius_sensitivity'].replace('_',' ').title()}",fontsize=10)
  examples.append(dict(particle_id=r['particle_id'],local_triangle_ids=np.flatnonzero(mask).tolist(),kept_contact_ids=kept,caption='Red: independent inward handoff normals (2 nm). Blue: local FEM velocity.'))
 fig.text(.5,.04,'Rigid spheres and actual wall triangles. Red arrows: retained handoff normals (2 nm gap). Blue arrow: FEM flow.\nVirtual radius tests describe model sensitivity, not measured bubble deformation.',ha='center',fontsize=10);fig.subplots_adjust(bottom=.17,wspace=.12)
 dump(D/'geometry_example_selection.json',examples);save(fig,'06_stationary_geometry_examples',['stationary_audit.json','plot_geometry.json','geometry_example_selection.json'])
 fig,ax=plt.subplots(figsize=(10,6.3));ratios=[1.,.99,.95,.90];statuscolors={'STATIONARY':'#bc5858','PASSED_LOCAL_HOTSPOT':'#20866a'}
 for i,r in enumerate(st):
  by={t['radius_ratio']:t for t in r['trials']}
  for j,ratio in enumerate(ratios):
   t=by.get(ratio);status=t['status'] if t else 'NOT_NEEDED';color=statuscolors.get(status,'#e8ebee' if t is None else '#dcab41')
   ax.add_patch(plt.Rectangle((j-.47,i-.4),.94,.8,facecolor=color,edgecolor='white'));ax.text(j,i,'Stopped' if status=='STATIONARY' else 'Passed' if status=='PASSED_LOCAL_HOTSPOT' else '—' if t is None else 'Unresolved',ha='center',va='center',color='white' if status in statuscolors else '#303030',fontsize=9)
 ax.set(xlim=(-.6,3.6),ylim=(len(st)-.5,-.8),xticks=range(4),xticklabels=['100%','99%','95%','90%'],yticks=range(len(st)),yticklabels=[f"ID {r['particle_id']} | {r['radius_m']*1e6:.3f} µm" for r in st],xlabel='Virtual radius / original radius',title='Local passage sensitivity | diagnostic only')
 ax.spines[['left','bottom']].set_visible(False);fig.text(.5,.01,'— = not needed for the final bracket; earlier surplus trials remain archived.\nVirtual radius sensitivity is not physical bubble shrinkage.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1))
 save(fig,'07_virtual_radius_sensitivity',['stationary_audit.json','virtual_radius_trials.csv'])
 fig,ax=plt.subplots(figsize=(11.6,5),layout='constrained');keys=['raw_candidate_point_split','accepted500_point_split','p65_500_split','p9a1_500_split'];counts=[None]+[s[k] for k in keys];fracs=[[s['flow_split']['outlet_fractions'].get(role,0) for role in ROLES]]+[[s[k][role]/sum(s[k].values()) for role in ROLES] for k in keys]
 for i,fr in enumerate(fracs):
  left=0.
  for j,frac in enumerate(fr):
   ax.barh(i,frac*100,left=left,color=COLORS[j],height=.64,label=LABELS[j] if i==0 else None)
   if frac>.055:ax.text(left+frac*50,i,f'{frac*100:.1f}%' if counts[i] is None else f'{counts[i][ROLES[j]]} ({frac*100:.1f}%)',ha='center',va='center',fontsize=10,color='white' if j!=1 else '#252525')
   left+=frac*100
  values=' / '.join(f'{v*100:.2f}%' for v in fr) if counts[i] is None else ' / '.join(str(counts[i][role]) for role in ROLES)
  ax.text(102,i,values,va='center',fontsize=10)
 ax.text(102,-.7,'O1 / O2 / O3 / no exit',fontsize=10,weight='bold');ax.set(yticks=range(5),yticklabels=['Fluid volume flow','Raw candidates | point','Accepted 500 | point','P6.5 | finite size','P9-A.1 | saved'],xlabel='Share of all candidates / particles (%)',xlim=(0,137),xticks=[0,25,50,75,100],title='Where does outlet routing change?');ax.invert_yaxis();ax.legend(loc='lower center',bbox_to_anchor=(.46,-.29),ncol=4,frameon=False);ax.spines[['left','right','top']].set_visible(False)
 save(fig,'08_routing_stage_summary',['audit_summary.json','stage_splits.csv'])
 dump(D/'figure_manifest.json',manifest);print('rendered',len(manifest),'figures')
if __name__=='__main__':main()
