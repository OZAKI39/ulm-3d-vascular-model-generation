from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import sys
sys.path.insert(0,str(R/'analysis/softening_feasibility_v1'))
from model import SofteningLPS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.labelsize':14,'axes.titlesize':15,'legend.fontsize':11,'figure.facecolor':'white','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
B='#187F94';RED='#C34049';ORANGE='#C16A29';PURPLE='#7756A2';GRAY='#7F95A4';GREEN='#3B8262'
def save(fig,name):
 for ext in ['png','svg']:fig.savefig(F/(name+'.'+ext),dpi=190,bbox_inches='tight')
 plt.close(fig)
def grid(ax,x='Physical time (µs)'):
 ax.set_xlabel(x);ax.grid(True,color='#607080',alpha=.35,lw=.7)
def main():
 f=json.loads((O/'PLOT_AND_REPORT_FACTS.json').read_text());z=np.load(O/'full/trajectory.npz');s=z['series'];t={c:s[:,i] for i,c in enumerate(z['columns'])};tt=t['time_s']*1e6
 lt=json.loads((O/'local/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());ft=json.loads((O/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());accepted=[q for q in lt['trials'] if q['decision']=='ACCEPT'];rejected=[q for q in lt['trials'] if q['decision']=='REJECT'];at=np.array([q['end_s']*1e6 for q in accepted]);ah=np.array([q['h_s']*1e6 for q in accepted]);ad=np.array([q['max_dd'] for q in accepted])
 fig,ax=plt.subplots(2,2,figsize=(16,10),layout='constrained');fig.suptitle('Damage-event resolution | local PASS, complete 700 µs protocol NOT PASSED',fontsize=19)
 ax[0,0].step(tt[1:],t['dt_s'][1:]*1e6,where='pre',color=B);ax[0,0].axvline(640.5,color=ORANGE,ls='--');ax[0,0].set(title='Actual accepted step size',ylabel='Accepted h (µs)',xlim=(0,700),ylim=(.10,.27))
 ax[0,1].step(at,ah,where='pre',color=B,marker='o',label='Accepted physical substep')
 for k,q in enumerate(rejected):ax[0,1].scatter(q['end_s']*1e6,q['h_s']*1e6,c=RED,marker='x',s=90,label='Rejected 0.25 µs candidate' if k==0 else None)
 ax[0,1].set(title='Native rollback → two real 0.125 µs substeps',ylabel='h (µs)',ylim=(.10,.28));ax[0,1].legend(loc='upper right')
 ax[1,0].plot(at,ad,'o-',color=B,label='Accepted max Δd')
 ax[1,0].scatter([q['end_s']*1e6 for q in rejected],[q['max_dd'] for q in rejected],c=RED,marker='x',s=85,label='Rejected max Δd');ax[1,0].axhline(.10,color=RED,ls='--',label='Unchanged admission: 0.10');ax[1,0].set(title='Damage increments at real candidate coordinates',ylabel='max Δd');ax[1,0].legend()
 ax[1,1].plot(tt,t['max_dd'],color=B);ax[1,1].axhline(.10,color=RED,ls='--');ax[1,1].set(title='All 2755 accepted substeps pass Δd ≤ 0.10',ylabel='max Δd',xlim=(0,700));ax[1,1].axvspan(tt[-1],700,color=RED,alpha=.13);ax[1,1].text(.03,.86,'Last accepted: 687.75 µs\nNext native trial loses a point at 688 µs',transform=ax[1,1].transAxes,fontsize=12)
 for a in ax.flat:grid(a)
 save(fig,'Fig1_event_substeps')
 # Native endpoint and edge-force series; old candidate is never an accepted point.
 zz=np.load(O/'local/bond1813_evolution.npz');v=zz['values'];tv=v[:,0]*1e6;old=np.load(O/'local/transactions/rejected_trial_00001.npz');before=np.load(O/'local/transactions/before_trial_00001.npz');nodes=[23,244]
 fig,ax=plt.subplots(3,3,figsize=(20,14),layout='constrained');fig.suptitle('Original bond 1813 | atoms 24 and 245 | real native event history',fontsize=20)
 ax[0,0].plot(tv,v[:,1],color=B,label='Current stretch s');ax[0,0].plot(tv,v[:,2],color=PURPLE,ls='--',label='History κ');ax[0,0].axhline(.15,c=GRAY,ls=':');ax[0,0].set(title='Stretch and maximum history',ylabel='s, κ');ax[0,0].legend()
 ax[0,1].plot(tv,v[:,3],'o-',color=B,label='Accepted adaptive d');ax[0,1].plot([640.25,640.5],[before['d'][1812],old['d'][1812]],'x--',color=RED,label='Old rejected candidate');ax[0,1].set(title='Irreversible damage',ylabel='d');ax[0,1].legend()
 ax[0,2].plot(tv,v[:,-1]*1e9,color=B,label='Native conjugate edge force');ax[0,2].scatter([640.5],[float(old['bond1813_tension_N'])*1e9],color=RED,marker='x',s=70,label='Rejected candidate');ax[0,2].set(title='Endpoint-star edge force',ylabel='Bond force (nN)');ax[0,2].legend()
 for k in range(3):
  a=ax[1,k]
  for n,color in [(0,B),(1,ORANGE)]:a.plot(tv,v[:,4+3*n+k]*1e6,color=color,label=f'Atom {nodes[n]+1}');a.scatter([640.5],[old['x'][nodes[n],k]*1e6],c=color,marker='x',s=60)
  a.set(title=f'Native {"xyz"[k]} coordinate',ylabel=f'{"xyz"[k]} (µm)');a.legend()
 for n,color in [(0,B),(1,ORANGE)]:
  ax[2,0].plot(tv,np.linalg.norm(v[:,10+3*n:13+3*n],axis=1),color=color,label=f'Atom {nodes[n]+1}');ax[2,0].scatter([640.5],[np.linalg.norm(old['velocity'][nodes[n]])],c=color,marker='x',s=60)
  ax[2,1].plot(tv,np.linalg.norm(v[:,16+3*n:19+3*n],axis=1)*1e6,color=color,label=f'Atom {nodes[n]+1}');ax[2,1].scatter([640.5],[np.linalg.norm(old['material_force'][nodes[n]])*1e6],c=color,marker='x',s=60)
 ax[2,0].set(title='Native endpoint speed',ylabel='Speed (m/s)');ax[2,0].legend();ax[2,1].set(title='Native total material-force norm',ylabel='Material force (µN)');ax[2,1].legend();ax[2,2].step(at,ah,where='pre',color=B,marker='o');ax[2,2].set(title='Actual physical integration step',ylabel='h (µs)')
 for a in ax.flat:grid(a);a.set_xlim(640.2,643.05)
 save(fig,'Fig2_bond1813_history')
 fig,ax=plt.subplots(3,2,figsize=(17,13),layout='constrained');fig.suptitle('New t = 0 native path | last accepted 687.75 µs | 700 µs NOT COMPLETED',fontsize=19)
 u=np.array([30*(np.cos(np.pi*(0 if x<350 else min((x-350)/300,1)))-np.cos(np.pi*min(x/300,1))) for x in tt]);ax[0,0].plot(tt,u,color=B,label='Right grip');ax[0,0].plot(tt,-u,color=ORANGE,label='Left grip');ax[0,0].set(title='Unchanged half-cosine prescribed loading',ylabel='Grip displacement (µm)');ax[0,0].legend()
 ax[0,1].plot(tt,t['max_free_disp_m']*1e6,color=B);ax[0,1].set(title='Maximum actual free-node displacement',ylabel='Displacement (µm)')
 ax[1,0].plot(tt,t['max_speed_m_s'],color=B);ax[1,0].set(title='Maximum actual node speed',ylabel='Speed (m/s)')
 ax[1,1].plot(tt,t['new_soft'],color=ORANGE,label='New softening IDs');ax[1,1].plot(tt,t['new_fail'],color=B,label='New complete failures');ax[1,1].set(title='Natural events per accepted step',ylabel='Original-bond count');ax[1,1].legend()
 ax[2,0].plot(tt,t['total_failed'],color=RED);ax[2,0].set(title='4918 failed / 7129 live original bonds at stop',ylabel='Completely failed bonds')
 ax[2,1].plot(tt,t['contact_pairs'],color=PURPLE);ax[2,1].set(title='Unchanged native 1.0 dx contact',ylabel='Active contact pairs')
 for a in ax.flat:
  grid(a);a.set_xlim(0,700);a.axvspan(tt[-1],700,color=RED,alpha=.13)
  for x in [300,350,650]:a.axvline(x,c=GRAY,alpha=.55,ls=':')
 save(fig,'Fig3_actual_dynamic_protocol')
 fig,ax=plt.subplots(2,3,figsize=(20,11),layout='constrained');fig.suptitle('Complete variational energy ledger | accepted interval only: 0–687.75 µs',fontsize=19)
 for name,color in [('W_J',B),('U_J',ORANGE),('K_J',PURPLE),('D_J',RED)]:ax[0,0].plot(tt,t[name]*1e9,label=name[:-2],c=color)
 ax[0,0].set(title='W = ΔU + ΔK + D + ΔUc − Jc + R',ylabel='Energy (nJ)');ax[0,0].legend(ncol=2)
 for name,color in [('Uc_J',PURPLE),('Jc_J',GREEN)]:ax[0,1].plot(tt,t[name]*1e12,label=name[:-2],c=color)
 ax[0,1].set(title='Contact potential and topology injection',ylabel='Energy (pJ)');ax[0,1].legend()
 ax[0,2].plot(tt,t['R_J']*1e12,color=RED,label='R');ax[0,2].plot(tt,-t['switch_work_J']*1e12,color=B,ls='--',label='− damage-force switch work');ax[0,2].set(title='Signed numerical balance',ylabel='Energy (pJ)');ax[0,2].legend()
 ax[1,0].plot(tt,t['relative_R']*100,c=B);ax[1,0].axhline(2,c=RED,ls='--',label='Frozen 2% limit');ax[1,0].set(title='Maximum relative balance: 0.37243%',ylabel='|R| / running energy scale (%)',ylim=(0,2.2));ax[1,0].legend()
 comp=f['local']['comparison_same_640p25_to_640p5_window'];labels=['ΔR','VV defect','Switch work','Contact defect'];names=['R_J','VV_defect_J','switch_work_J','contact_defect_J'];xx=np.arange(4)
 ax[1,1].bar(xx-.18,[comp[k]['old_single_step']*1e15 for k in names],width=.36,color=RED,label='Rejected old 0.25 µs');ax[1,1].bar(xx+.18,[comp[k]['new_two_steps']*1e15 for k in names],width=.36,color=B,label='Two real 0.125 µs');ax[1,1].set_xticks(xx,labels,rotation=15);ax[1,1].set(title='Same physical window: 640.25 → 640.5 µs',ylabel='Increment (fJ)');ax[1,1].legend(fontsize=10)
 for name,color in [('VV_defect_J',B),('spatial_defect_J',ORANGE),('contact_defect_J',PURPLE)]:ax[1,2].plot(tt,t[name]*1e12,c=color,label=name.replace('_J','').replace('_',' '))
 ax[1,2].set(title='Other cumulative time-layer defects',ylabel='Energy (pJ)');ax[1,2].legend(fontsize=10)
 for a in ax.flat:grid(a)
 ax[1,1].set_xlabel('Actual event-window comparison')
 save(fig,'Fig4_energy_ledger')
 # Fixed camera/axes; exact points and original reference bond identities.
 refs=np.load(O/'provenance/reference.npz');b=refs['bonds'];choices=[(0,'Initial'),(340,'First softening'),(370,'First complete failure'),(2563,'Resolved event'),(2604,'Unloading complete'),(2755,'Last accepted')]
 states=[np.load(O/'full/keyframes'/f'step_{step:06d}.npz') for step,_ in choices];lim=np.ceil(max(abs(q['x']).max() for q in states)*1e6/50)*50+25
 fig=plt.figure(figsize=(20,13));fig.suptitle('Actual 341-point specimen | irreversible damage | fixed axes and camera',fontsize=21,y=.97);cm=LinearSegmentedColormap.from_list('actual_damage',[B,'#E9A969',RED]);norm=Normalize(0,1)
 near=refs['r0']<=1.01*DX
 for k,(q,(_,title)) in enumerate(zip(states,choices)):
  a=fig.add_subplot(2,3,k+1,projection='3d');x=q['x']*1e6;d=q['d'];nd=np.zeros(341);np.maximum.at(nd,b[:,0],d);np.maximum.at(nd,b[:,1],d)
  for mask,ls,alpha in [(near&(d<1),'-',.4),(near&(d==1),'--',.18)]:
   if np.any(mask):a.add_collection3d(Line3DCollection(x[b[mask]],colors=cm(d[mask]),linewidths=.65,linestyles=ls,alpha=alpha))
  a.scatter(x[:,0],x[:,1],x[:,2],s=20,c=nd,cmap=cm,norm=norm,depthshade=False,edgecolors='#2A3540',linewidths=.2)
  if k==5:a.scatter(*x[72],s=110,facecolors='none',edgecolors='#00A5C0',linewidths=2);a.text(*x[72],' 73',color=B,fontsize=13)
  a.set(xlim=(-lim,lim),ylim=(-lim,lim),zlim=(-lim,lim),xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title=f'{title} | {float(q["time"])*1e6:.3f} µs');a.view_init(elev=17,azim=-58);a.set_box_aspect((1,1,1));a.tick_params(labelsize=11,pad=1)
  for axis in [a.xaxis,a.yaxis,a.zaxis]:axis._axinfo['grid'].update(color='#738290',linewidth=.7);axis.pane.fill=True;axis.pane.set_facecolor('#F2F5F7')
 fig.subplots_adjust(left=.02,right=.90,bottom=.06,top=.91,wspace=.08,hspace=.10)
 cax=fig.add_axes([.93,.20,.014,.58]);cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cm),cax=cax);cb.set_label('Node: max incident d / bond: actual d',fontsize=14)
 fig.text(.04,.015,'All 341 real points. Short reference edges (r₀ ≤ 1.01 dx) shown for readability. Dashed failed edges are identity traces, not live forces.\n687.75 µs is the last accepted state; no 700 µs geometry is fabricated. Circled atom 73 crosses the fixed box on the next drift.',fontsize=12)
 save(fig,'Fig5_actual_material_points')
 print('FIVE_FIGURE_GROUPS_WRITTEN')
if __name__=='__main__':main()
