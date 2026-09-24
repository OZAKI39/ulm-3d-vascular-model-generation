"""High-resolution figures from saved paired data; original vascular triangles."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv
R=Path(__file__).resolve().parents[1]
read=lambda p:json.loads(Path(p).read_text())
S=read(R/'data/analysis_summary.json');T=read(R/'data/paired_outcome_transition.json')
M={label:{r['bubble_id']:r for r in read(R/'data'/f'{label.lower()}_paired_metrics.json')} for label in ['OLD','NEW']}
COL={'O1':'#e69f00','O2':'#0072b2','O3':'#009e73','STATIONARY':'#cc3311','OTHER':'#666666'}
plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','font.size':11,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
flow=read(R/'data/old_new_flow_contract.json');wall=pv.read(flow['boundaries']['WALL']['paths'][0]);tri=np.asarray(wall.points)[wall.faces.reshape(-1,4)[:,1:]]*1e6
pts=np.asarray(wall.points)*1e6;lo=pts.min(0);hi=pts.max(0);center=(hi+lo)/2;span=hi-lo

def save(fig,name):
    if len(sys.argv)>1 and name not in sys.argv[1:]:
        plt.close(fig);return
    fig.savefig(R/'figures'/f'{name}.png',dpi=320,bbox_inches='tight')
    fig.savefig(R/'figures'/f'{name}.pdf',dpi=320,bbox_inches='tight');plt.close(fig)
    print(name,flush=True)
def geometry(ax,labels=True):
    ax.add_collection3d(Poly3DCollection(tri,facecolors=(.60,.64,.69,.07),edgecolors='none',rasterized=True))
    ax.set_xlim(lo[0]-2,hi[0]+2);ax.set_ylim(lo[1]-2,hi[1]+2);ax.set_zlim(lo[2]-2,hi[2]+2)
    ax.set_box_aspect(span);ax.view_init(elev=24,azim=-62);ax.set_proj_type('ortho');ax.grid(False)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.pane.set_edgecolor('white')
    if labels:ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_zlabel('z (µm)')
    else:ax.set_axis_off()
def path(label,pid):return np.load(R/'outputs'/label/'trajectories'/f'mb_{pid:06d}.npz')['samples'][:,1:4]*1e6
def draw(ax,label,pid,color,width=1.5,linestyle='-'):
    x=path(label,pid);ax.plot(*x.T,color=color,lw=width,ls=linestyle,zorder=20)
    ax.scatter(*x[0],c=color,s=9,zorder=21)
    ax.scatter(*x[-1],c=color,s=22,marker='x',zorder=21)

fig=plt.figure(figsize=(11,10));ax=fig.add_subplot(projection='3d');geometry(ax)
for t in T:draw(ax,'NEW',t['bubble_id'],COL.get(t['new_outcome'],COL['OTHER']),1.3)
for pid in S['stationary_new_ids']:
    x=path('NEW',pid)[-1];ax.text(*x,f'  ID {pid}: stationary',color=COL['STATIONARY'],fontsize=10,zorder=30)
ax.set_title('Network-derived flow: 30 paired microbubble trajectories',pad=22)
handles=[Line2D([0],[0],color=c,lw=3,label=(f'{k}: {S["new_outlet_counts"].get(k,0)}' if k.startswith('O') and k!='OTHER' else f'Stationary: {S["new_stationary"]}' if k=='STATIONARY' else 'Other: 0')) for k,c in COL.items()]
ax.legend(handles=handles,loc='upper left',frameon=False);fig.text(.5,.04,'Actual vascular WALL triangles • fixed historical cohort • N=30 paired validation',ha='center')
save(fig,'01_new_flow_30_tracks')

selected=[];reasons={}
def take(rows,reason):
    for t in rows:
        i=t['bubble_id']
        if i not in selected and len(selected)<8:selected.append(i);reasons[str(i)]=reason
for condition,reason in [
    (lambda t:t['old_status']=='STATIONARY' and t['new_status']=='COMPLETED','stationary -> completed'),
    (lambda t:t['old_status']=='COMPLETED' and t['new_status']=='STATIONARY','completed -> stationary'),
    (lambda t:t['old_outcome']=='O2' and t['new_outcome']=='O3','O2 -> O3'),
    (lambda t:t['old_outcome']=='O3' and t['new_outcome']=='O2','O3 -> O2'),
    (lambda t:t['old_status']==t['new_status']=='STATIONARY','persistent stationary')]:
    candidates=sorted([t for t in T if condition(t)],key=lambda t:t['bubble_id']);take(candidates[:1],reason)
for k in sorted({(t['old_outcome'],t['new_outcome']) for t in T}):
    candidates=sorted([t for t in T if (t['old_outcome'],t['new_outcome'])==k],key=lambda t:t['bubble_id']);take(candidates[:1],'one lowest-ID representative of each observed transition')
changed=sorted([t for t in T if t['outlet_changed'] and t['bubble_id'] not in selected],key=lambda t:(t['diameter_um'],t['bubble_id']))
if changed:
    take([changed[i] for i in np.unique(np.linspace(0,len(changed)-1,min(8-len(selected),len(changed))).round().astype(int))],'diameter-spaced changed-route examples, deterministic ID tie break')
take(sorted(T,key=lambda t:t['bubble_id']),'remaining slots by stable ID')
selection=dict(selected_ids=selected,reasons=reasons,rule='Prioritize stationary transitions, O2->O3, O3->O2, persistent stationary; then one representative per observed transition; fill changed-route diameter-spaced cases and stable IDs. No trajectory appearance used.',
    all_observed_transitions=S['transition_matrix'],absent_requested_categories=[k for k in ['STATIONARY -> O1','STATIONARY -> O2','STATIONARY -> O3','O3 -> O2'] if k not in S['transition_matrix']])
(R/'data/representative_selection.json').write_text(json.dumps(selection,indent=2)+'\n')
fig=plt.figure(figsize=(18,11));lookup={t['bubble_id']:t for t in T}
for k,pid in enumerate(selected):
    ax=fig.add_subplot(2,4,k+1,projection='3d');geometry(ax,False);draw(ax,'OLD',pid,'#0072b2',2,'--');draw(ax,'NEW',pid,'#d55e00',2)
    t=lookup[pid];ax.set_title(f'ID {pid} | D={t["diameter_um"]:.3f} µm\nOLD {t["old_outcome"]} → NEW {t["new_outcome"]}',fontsize=12,pad=1)
fig.suptitle('Same bubble, same initial state, different background flow',fontsize=18,y=.98)
fig.legend(handles=[Line2D([0],[0],color='#0072b2',lw=2,ls='--',label='OLD'),Line2D([0],[0],color='#d55e00',lw=2,label='NEW')],loc='lower center',ncol=2,frameon=False)
fig.subplots_adjust(wspace=.02,hspace=.13,top=.89,bottom=.08)
save(fig,'02_paired_old_new_trajectories')

def matrix_plot(matrix,xnames,ynames,xlabel,ylabel,title,name):
    a=np.zeros((len(ynames),len(xnames)),int)
    for k,n in matrix.items():
        y,x=k.split(' -> ');a[ynames.index(y),xnames.index(x)]+=n
    fig,ax=plt.subplots(figsize=(8,6.5));im=ax.imshow(a,cmap='Blues',vmin=0,vmax=max(1,a.max()))
    for i in range(len(ynames)):
        for j in range(len(xnames)):ax.text(j,i,str(a[i,j]),ha='center',va='center',fontsize=16,color='white' if a[i,j]>.55*a.max() else 'black')
    ax.set_xticks(range(len(xnames)),xnames,rotation=20);ax.set_yticks(range(len(ynames)),ynames);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel);ax.set_title(title,pad=15)
    fig.colorbar(im,ax=ax,label='Bubble count',shrink=.75);fig.tight_layout();save(fig,name)
names=['O1','O2','O3','STATIONARY','OTHER']
matrix_plot(S['transition_matrix'],names,names,'NEW outcome','OLD outcome','Paired outcome transitions | N=30','03_paired_outcome_transition')
v=read(R/'data/paired_transit_time.json');x=np.array([r['old_s'] for r in v])*1e3;y=np.array([r['new_s'] for r in v])*1e3
fig,ax=plt.subplots(figsize=(8,7));ax.scatter(x,y,c='#0072b2',s=45);bounds=[0,max(x.max(),y.max())*1.08];ax.plot(bounds,bounds,'--',color='.5',lw=1)
outliers=sorted(v,key=lambda r:(-abs(np.log(r['ratio'])),r['bubble_id']))[:4]
cluster=sorted([r for r in outliers if r['old_s']*1e3<bounds[1]/2],key=lambda r:r['new_s'])
for r in outliers:
    xy=(r['old_s']*1e3,r['new_s']*1e3)
    if r in cluster:
        target=(bounds[1]*.23,70+22*cluster.index(r))
        ax.annotate('ID '+str(r['bubble_id']),xy,xytext=target,textcoords='data',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='.4',lw=.7),va='center')
    else:ax.annotate('ID '+str(r['bubble_id']),xy,xytext=(-5,10),ha='right',textcoords='offset points',fontsize=10)
ax.set(xlim=bounds,ylim=bounds,xlabel='OLD transit time (ms)',ylabel='NEW transit time (ms)',title=f'Paired transit time | both completed, n={len(v)}');ax.set_aspect('equal');fig.tight_layout();save(fig,'04_transit_time_paired')
(R/'data/transit_outlier_selection.json').write_text(json.dumps(dict(rule='Top four absolute log(NEW/OLD transit time), stable ID tie break',ids=[r['bubble_id'] for r in outliers]),indent=2))
fig,axes=plt.subplots(1,2,figsize=(13,6))
for t in T:
    o,n=M['OLD'][t['bubble_id']],M['NEW'][t['bubble_id']];stationary='STATIONARY' in [t['old_status'],t['new_status']]
    for ax,key,scale in [(axes[0],'minimum_h_over_a',1),(axes[1],'nearwall_exposure_h_over_a_le_0p1_s',1000)]:
        ax.scatter(o[key]*scale,n[key]*scale,c='#cc3311' if stationary else '#0072b2',marker='*' if stationary else 'o',s=120 if stationary else 32)
        if stationary:ax.annotate(str(t['bubble_id']),(o[key]*scale,n[key]*scale),xytext=(6,4),textcoords='offset points')
for ax in axes:
    a,b=ax.get_xlim();c,d=ax.get_ylim();lim=[min(a,c),max(b,d)];ax.plot(lim,lim,'--',c='.5',lw=1);ax.set_xlim(lim);ax.set_ylim(lim)
axes[0].set(xlabel='OLD minimum h/a',ylabel='NEW minimum h/a',title='Minimum wall clearance | all 30')
axes[1].set(xlabel='OLD exposure (ms)',ylabel='NEW exposure (ms)',title='Saved-step exposure at h/a ≤ 0.1')
fig.suptitle('Near-wall exposure • red star: stationary in either flow');fig.tight_layout();save(fig,'05_minimum_gap_and_nearwall_exposure')
pm=S['point_vs_MB_transition'];yn=sorted({k.split(' -> ')[0] for k in pm});xn=names
matrix_plot(pm,xn,yn,'NEW finite-size MB outcome','NEW point-tracer outcome','Same initial centers | N=30','06_point_tracer_vs_microbubble')
fig,axes=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.6,1]});xx=np.arange(3);width=.34
axes[0].bar(xx-width/2,[S['old_fluid_split'][f'O{i}']*100 for i in range(1,4)],width,label='OLD',color='#0072b2')
axes[0].bar(xx+width/2,[S['new_fluid_split'][f'O{i}']*100 for i in range(1,4)],width,label='NEW',color='#d55e00')
axes[0].set(xticks=xx,xticklabels=['O1','O2','O3'],ylabel='Integrated fluid Qout / Qin (%)',ylim=(0,100),title='Background FEM fluid routing');axes[0].legend(frameon=False)
axes[1].axis('off');axes[1].text(.5,.75,str(S['outlet_changed_count']),ha='center',fontsize=55,color='#d55e00')
axes[1].text(.5,.53,f'paired bubbles changed outlet\n{len(v)} completed in both flows',ha='center',fontsize=13)
axes[1].text(.5,.18,'N=30 fixed paired cohort\nNot population split statistics',ha='center',fontsize=12)
fig.tight_layout();save(fig,'07_fluid_split_and_paired_route_changes')
