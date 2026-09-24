#!/usr/bin/env python3
"""Six English academic figures from saved numerical tables (no simulation)."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from analyze_particle9a1 import ROOT,R,D,read,dump
F=R/'figures';F.mkdir(exist_ok=True);manifest=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42})
OLD='#bc4b51';NEW='#176b9b';REF='#368479'
def csvrows(name):return list(csv.DictReader((D/name).open()))
def save(fig,name,sources):
 for ext in ['png','pdf']:fig.savefig(F/(name+'.'+ext),dpi=220,bbox_inches='tight')
 plt.close(fig);manifest.append(dict(png='figures/'+name+'.png',pdf='figures/'+name+'.pdf',sources=sources))
def table():
 rows=read(D/'same12.json');ids=sorted({r['particle_id'] for r in rows});models=['P65','P9A','P9A1','P65_NEW']
 labels={'COMPLETED':'Completed','INLET_ESCAPE':'Inlet escape','HANDOFF_OR_SAFETY_STOP':'Safety stop','TIME_LIMIT':'Time limit','OUTSIDE_OTHER':'Other outside','OTHER':'Other'}
 colors={'COMPLETED':'#deeee8','INLET_ESCAPE':'#f5dfe0','HANDOFF_OR_SAFETY_STOP':'#fff0d6','TIME_LIMIT':'#e9dff3','OUTSIDE_OTHER':'#f5dfe0','OTHER':'#eeeeee'}
 fig,ax=plt.subplots(figsize=(12,7));ax.axis('off');body=[];cells=[]
 for pid in ids:
  rr=[next(r for r in rows if r['particle_id']==pid and r['model']==m) for m in models]
  body.append([str(pid)]+['Stationary stop' if 'STATIONARY' in str(r['failure_detail']) else labels[r['outcome']] for r in rr]);cells.append(['#edf1f5']+[colors[r['outcome']] for r in rr])
 tab=ax.table(cellText=body,cellColours=cells,colLabels=['ID','P6.5 original','P9-A original','P9-A.1','P6.5 + new handoff'],cellLoc='center',bbox=[0,0,1,.88],colWidths=[.08,.22,.22,.22,.26]);tab.auto_set_font_size(False);tab.set_fontsize(11)
 for (i,j),cell in tab.get_celld().items():
  cell.set_edgecolor('white')
  if i==0:cell.set_facecolor('#29475c');cell.set_text_props(color='white',weight='bold')
 counts=[sum(r['completed'] for r in rows if r['model']==m) for m in models]
 ax.set_title('Same 12 particles | Identical 2.0 mm/s FEM and birth states',fontsize=16,pad=28)
 ax.text(.5,.94,'Completed:  '+ '  |  '.join(f'{m}: {n}/12' for m,n in zip(['P6.5','P9-A','P9-A.1','P6.5 new'],counts)),ha='center',transform=ax.transAxes)
 save(fig,'01_same12_before_after',['data/same12.csv','data/same12.json'])
def inlet():
 rows=read(D/'inlet_3_15.json');fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
 for i,r in enumerate(rows):
  ax,bx=axes[i];values=np.array([r['FEM_inward_m_s'],r['old_inward_m_s'],r['new_inward_m_s']])*1e3
  ax.bar(['FEM','Old P9-A','P9-A.1'],values,color=['#666666',OLD,NEW],width=.55)
  for x,y in enumerate(values):ax.annotate(f'{y:.3g}',(x,y),xytext=(0,8 if y>=0 else -18),textcoords='offset points',ha='center')
  ax.axhline(0,color='#444444',lw=.7);ax.set_ylim(-.45,max(values)*1.25);ax.set_ylabel('Inward velocity (mm/s)');ax.set_title(f'ID {r["particle_id"]} | First solve')
  vv=np.array([r['bulk_tangential_m_s'],r['old_target_tangential_m_s'],r['new_solved_tangential_m_s']])*1e3
  bx.scatter(range(3),vv,s=65,c=['#666666',OLD,NEW]);bx.set_yscale('log');bx.set_xticks(range(3),['FEM bulk','Old target','New solved']);bx.set_ylim(min(vv)/5,max(vv)*6)
  for x,y in enumerate(vv):bx.annotate(f'{y:.3g}',(x,y),xytext=(0,9),textcoords='offset points',ha='center')
  bx.set_ylabel('Tangential speed (mm/s)');bx.set_title('Actual flow retained at the open rim');bx.grid(axis='y',alpha=.2)
 fig.suptitle('Inlet compatibility | Exact cap-triangle inward normal',fontsize=16)
 save(fig,'02_inlet_3_15',['data/inlet_3_15.csv','data/inlet_3_15.json'])
def canonical():
 rows=csvrows('canonical.csv');x=[float(r['gap_ratio']) for r in rows];fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
 for ax,quantity,label,scale in [(axes[0],'Vt_m_s','Tangential velocity (mm/s)',1e3),(axes[1],'omega_s_inv',r'Angular velocity (s$^{-1}$)',1)]:
  for prefix,name,color,ls,lw in [('reference_2d','2D reference',REF,'-',3),('old','Old P9-A',OLD,'--',2),('new','P9-A.1',NEW,':',2)]:
   ax.plot(x,[float(r[prefix+'_'+quantity])*scale for r in rows],label=name,color=color,ls=ls,lw=lw)
  ax.set_xscale('log');ax.set_xlabel('Wall gap / particle radius');ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend(fontsize=9)
 fig.suptitle('Canonical planar shear | Same mobility and reciprocity projection',fontsize=15)
 save(fig,'03_canonical_planar',['data/canonical.csv','data/canonical_summary.json'])
def handoff():
 rows=csvrows('handoff.csv');fig,axes=plt.subplots(3,2,figsize=(12,9),layout='constrained')
 for i,pid in enumerate([12,14,17]):
  for model,label,color in [('P9A','Old P9-A',OLD),('P9A1','P9-A.1',NEW)]:
   rs=[r for r in rows if int(r['particle_id'])==pid and r['model']==model];t=np.array([float(r['time_s']) for r in rs])*1000;g=np.array([float(r['g_nf_m']) for r in rs])*1e9;rat=np.array([float(r['dt_over_nominal']) for r in rs])
   axes[i,0].plot(t,g,lw=1.3,label=label,color=color);axes[i,1].plot(t[1:],rat[1:],lw=1.1,label=label,color=color)
  axes[i,0].set_yscale('symlog',linthresh=1e-7);axes[i,0].axhline(0,color='#888888',lw=.7);axes[i,1].set_yscale('log');axes[i,1].axhline(1,color='#888888',ls=':',lw=.7)
  axes[i,0].set_ylabel('Gap above handoff (nm)');axes[i,1].set_ylabel('Accepted dt / nominal dt')
  for ax in axes[i]:ax.set_title(f'ID {pid}');ax.set_xlabel('Elapsed physical time (ms)');ax.grid(alpha=.18);ax.legend(fontsize=8)
 fig.suptitle('Handoff histories | All accepted intervals; original safety guards retained',fontsize=15)
 save(fig,'04_handoff_before_after',['data/handoff.csv','data/same12.csv'])
def contact():
 rows=csvrows('contact_id7.csv');fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained');x=np.arange(2)
 axes[0].bar(x-.15,[float(r['retained_count']) for r in rows],width=.28,label='Retained constraints',color=[OLD,NEW]);axes[0].bar(x+.15,[float(r['numerical_rank']) for r in rows],width=.28,label='Solver numerical rank',color=['#d99499','#83b3ce']);axes[0].set_ylim(0,2.5);axes[0].set_yticks([0,1,2]);axes[0].set_ylabel('Count');axes[0].legend(fontsize=8)
 axes[1].bar(x,[float(r['condition']) for r in rows],color=[OLD,NEW],width=.55);axes[1].set_yscale('log');axes[1].set_ylabel('Contact condition estimate')
 axes[2].bar(x,[float(r['minimum_normal_speed_m_s'])*1e6 for r in rows],color=[OLD,NEW],width=.55);axes[2].axhline(0,color='#555555',lw=.7);axes[2].set_ylabel('Minimum normal velocity (µm/s)')
 for ax in axes:ax.set_xticks(x,['Original','Fixed']);ax.grid(axis='y',alpha=.2)
 fig.suptitle('ID 7 | Replay of the exact saved duplicate-contact system',fontsize=15)
 save(fig,'05_contact_redundancy_id7',['data/contact_id7.csv','data/id7_fixed_fixture.json','reference/id7_duplicate_contact.json'])
def trajectories():
 from particle_3d.particle81_simulation import environment
 env=environment();rows=csvrows('real_paths.csv');summary=read(D/'same12.json');gates=read(D/'gates.json');fig=plt.figure(figsize=(11,9));ax=fig.add_subplot(111,projection='3d',computed_zorder=False)
 mesh=env.wall.triangles*1e6;ax.add_collection3d(Poly3DCollection(mesh,facecolor='#778da0',edgecolor='#8397a8',linewidth=.035,alpha=.18,zorder=2))
 colors=plt.colormaps['tab20'](np.array([0,2,4,6,8,10,12,14,16,18,1,9])/19)
 for pid,color in zip(sorted({int(r['particle_id']) for r in rows}),colors):
  rr=[r for r in rows if int(r['particle_id'])==pid];xyz=np.array([[float(r[k])*1e6 for k in ['x_m','y_m','z_m']] for r in rr]);s=next(r for r in summary if r['model']=='P9A1' and r['particle_id']==pid)
  ax.plot(*xyz.T,color=color,lw=1.8,label=f'ID {pid}'+('' if s['completed'] else ' (incomplete)'),zorder=3)
  ax.scatter(*xyz[-1],color=color,s=24,marker='o' if s['completed'] else 'x',zorder=4)
 bounds=env.wall.triangles.reshape(-1,3)*1e6
 for setter,k in [(ax.set_xlim,0),(ax.set_ylim,1),(ax.set_zlim,2)]:setter(bounds[:,k].min(),bounds[:,k].max())
 ax.set_box_aspect(np.ptp(bounds,axis=0));ax.view_init(elev=22,azim=130)
 ax.set_xlabel('X (µm)',labelpad=12);ax.set_ylabel('Y (µm)',labelpad=12);ax.set_zlabel('Z (µm)',labelpad=12)
 ax.set_title('Real FEM trajectories | All 12 repaired-model paths',fontsize=16,pad=20);ax.legend(loc='upper left',bbox_to_anchor=(1.12,.93),fontsize=9,frameon=False)
 fig.text(.5,.025,'Same-12 key gates: '+('PASS' if gates['same12_pass'] else 'FAIL — incomplete paths retained'),ha='center',fontsize=11)
 dump(D/'trajectory_geometry.json',dict(wall_provenance=env.wall.provenance,rendered_triangle_stride=1,coordinate_display_scale=1e6,rendered_triangles_m=env.wall.triangles,all_paths_included=True))
 save(fig,'06_real_same12_trajectories',['data/real_paths.csv','data/trajectory_geometry.json','data/same12.json'])
if __name__=='__main__':
 for fn in [table,inlet,canonical,handoff,contact,trajectories]:fn()
 dump(D/'figure_manifest.json',manifest)
