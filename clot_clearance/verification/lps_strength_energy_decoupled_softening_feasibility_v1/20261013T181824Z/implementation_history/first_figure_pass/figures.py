from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import time,resource
start=time.perf_counter();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'legend.fontsize':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.25,'grid.color':'#607080'})
colors=['#555D66','#187F94','#C16A29'];names=['instant_reference','soft_sc018','soft_sc024'];labels=['Instant break reference','Softening: sc = 0.18','Softening: sc = 0.24'];data={n:np.load(O/(n+'.npz')) for n in names+['planar_partial_sc024']};audit=json.loads((O/'STRENGTH_ENERGY_DECOUPLING_AUDIT.json').read_text());budget=json.loads((O/'ENERGY_TARGET_FEASIBILITY.json').read_text());figs=[]
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);figs.append(dict(name=name,caption=caption,png_sha256=sha(F/(name+'.png')),svg_sha256=sha(F/(name+'.svg'))))
fig,ax=plt.subplots(1,2,figsize=(13.4,5.2),layout='constrained')
for n,label,col in zip(names,labels,colors):
 z=data[n];a=z['series'];pick=(z['phase']<=2)|(z['phase']==6)
 ax[0].plot(a[pick,0],a[pick,2]*1e3,label=label,color=col,lw=2.2)
 ax[1].plot(a[pick,0],a[pick,5]*1e9,color=col,lw=2.2)
for aa in ax:aa.axvline(.15,color='#333',ls=':',lw=1.3);aa.set_xlabel('Compatible isotropic dilation s (dimensionless)')
ax[0].set(ylabel='Generalized reaction Q (mN)',title='A  Actual 341-node lattice reactions');ax[0].legend(loc='upper right');ax[0].annotate('onset = 0.15',(.15,38.36),(.075,43),arrowprops=dict(arrowstyle='->',color='#444'),fontsize=10)
for sc,col in [(.18,colors[1]),(.24,colors[2])]:ax[0].plot(sc,0,'o',color=col);ax[0].annotate(f'sc = {sc}',(sc,0),(sc-.027,5),fontsize=10,color=col)
ax[1].set(ylabel='Cumulative damage dissipation (nJ)',title='B  Signed global energy ledger');ax[1].text(.02,.95,'Stored energy and endpoint coupling included\nNo bond-spring area summation',transform=ax[1].transAxes,va='top',fontsize=10)
fig.suptitle('Fixed K, G and onset; one variational softening mechanism',fontsize=17)
fig.supxlabel('Offline prescribed motion; Q = −Σ F · (X / dx), conjugate to δ = s dx. Not native fracture or material Gc.',fontsize=10)
save(fig,'Fig1_instant_and_softening','Actual offline lattice reaction and global damage dissipation under compatible uniform expansion. Envelope excludes unloading branches, which appear in Fig3; instantaneous jump preserved.')
# Every displayed connection is an original real edge. Gray intact lines are optically sparse.
z=data['soft_sc024'];p=data['planar_partial_sc024'];s=z['series'][:,0];ph=z['phase'];ip=np.flatnonzero(np.isclose(s,.165)&(ph==2))[0];iu=np.flatnonzero((s==0)&(ph==3))[0];ir=np.flatnonzero(np.isclose(s,.165)&(ph==5))[0];pa=np.flatnonzero((p['phase']==1))[-1]
frames=[(z,0,'A  Original reference'),(z,ip,'B  Partial damage, s = 0.165'),(z,iu,'C  Unloaded: damage retained'),(z,ir,'D  Reloaded: same history'),(p,len(p['series'])-1,'E  Partial plane: 49 failed edges'),(z,len(s)-1,'F  Fully damaged uniform protocol')]
fig=plt.figure(figsize=(15.8,10.6));axes=[]
for k,(zz,it,title) in enumerate(frames):
 aa=fig.add_subplot(2,3,k+1,projection='3d');axes.append(aa);x=zz['x'][it]*1e6;d=zz['d'][it];b=zz['bonds'];intact=d<1e-10;part=(d>=1e-10)&(d<1-1e-10);full=d>=1-1e-10
 for sel,col,alpha,lw in [(intact,'#7F95A4',.028,.35),(part,'#E39438',.09,.50),(full,'#C34049',.065,.50)]:
  if sel.any():aa.add_collection3d(Line3DCollection(x[b[sel]],colors=col,alpha=alpha,linewidths=lw,linestyles='--' if col=='#C34049' else '-'))
 nd=np.bincount(b.ravel(),weights=np.repeat(d,2),minlength=len(x))/np.bincount(b.ravel(),minlength=len(x))
 aa.scatter(*x.T,c=nd,cmap='viridis',vmin=0,vmax=1,s=15,edgecolors='white',linewidth=.16,depthshade=False)
 aa.set_title(title,pad=8,fontsize=13);aa.set_xlabel('X (µm)',labelpad=2);aa.set_ylabel('Y (µm)',labelpad=2);aa.set_zlabel('Z (µm)',labelpad=2);aa.view_init(elev=17,azim=-58);aa.set_box_aspect((1,1,1));aa.set(xlim=(-205,200),ylim=(-200,200),zlim=(-200,200));aa.set_xticks([-150,0,150]);aa.set_yticks([-150,0,150]);aa.set_zticks([-150,0,150]);aa.tick_params(labelsize=9,pad=0)
 for axis in [aa.xaxis,aa.yaxis,aa.zaxis]:axis.pane.set_facecolor((.965,.974,.981,1));axis._axinfo['grid']['color']=(.52,.59,.64,.6)
 handles=[Line2D([0],[0],color='#7F95A4',lw=2,label='Intact original edge'),Line2D([0],[0],color='#E39438',lw=2,label='Partially damaged edge'),Line2D([0],[0],color='#C34049',ls='--',lw=2,label='Failed edge: reference connection only')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.035),ncol=3,frameon=False)
fig.suptitle('Actual reference IDs and bond histories: 341 nodes / 12,047 edges',fontsize=18,y=.98)
fig.text(.5,.017,'All states are computed offline. Dashed failed edges identify original pairs and carry zero force; panel E is not complete plane separation.',ha='center',fontsize=10)
fig.subplots_adjust(left=.025,right=.965,bottom=.12,top=.91,wspace=.12,hspace=.18)
sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='viridis');cax=fig.add_axes([.972,.3,.012,.42]);fig.colorbar(sm,cax=cax,label='Mean incident damage')
save(fig,'Fig2_reference_lattice_and_damage','Original full-precision coordinates and edge identities, partial damage, unloaded/reloaded irreversible states, actual partial plane, and fully damaged uniform protocol. Failed connections drawn dashed as identity markers, not active forces.')
fig,ax=plt.subplots(2,3,figsize=(15,8.5),layout='constrained');a=z['series'];step=np.arange(len(a));q=a[:,2]*1e3
ax[0,0].plot(step,a[:,1]*1e6,color=colors[2]);ax[0,0].set(title='A  Prescribed displacement',xlabel='Loading index (no physical time)',ylabel='δ = s dx (µm)')
ax[0,1].plot(a[:,1]*1e6,q,color=colors[2]);ax[0,1].set(title='B  Actual reaction loop',xlabel='δ = s dx (µm)',ylabel='Q (mN)')
for col,label,idx in [('#333','External work',4),('#287CA6','Recoverable energy',3),('#C16A29','Damage dissipation',5)]:ax[0,2].plot(step,a[:,idx]*1e9,label=label,color=col,lw=1.8)
ax[0,2].legend();ax[0,2].set(title='C  Global energy ledger',xlabel='Loading index',ylabel='Energy (nJ)')
ax[1,0].plot(step,a[:,6]*1e15,color='#735A8F');ax[1,0].set(title='D  Signed work-balance residual',xlabel='Loading index',ylabel='W − U − D (fJ)')
ax[1,1].plot(step,a[:,9],label='Mean damage d',color=colors[2]);ax[1,1].plot(step,a[:,11],label='Maximum history κ',color='#287CA6');ax[1,1].legend();ax[1,1].set(title='E  Irreversible state',xlabel='Loading index',ylabel='Dimensionless')
for pno,col,label in [(2,'#C16A29','Damage loading'),(3,'#287CA6','Unload'),(4,'#805790','Compression'),(5,'#3B8D61','Reload')]:
 ii=ph==pno;ax[1,2].plot(a[ii,1]*1e6,q[ii],color=col,label=label,lw=2,ls='--' if pno==5 else '-')
ax[1,2].legend();ax[1,2].set(title='F  Frozen damage on return path',xlabel='δ = s dx (µm)',ylabel='Q (mN)')
fig.suptitle('Loading → unloading → compression → reloading: sc = 0.24',fontsize=18)
fig.supxlabel('Complete energy balance; no inertia or damping. Splitting error vs analytical dissipation = 0.2615%, separate from work closure.',fontsize=10)
save(fig,'Fig3_work_energy_irreversibility','Full controlled displacement history and all-DOF external work; signed residual shown without clipping. Reversible unloading/reloading at fixed damage, no healing or unexplained energy source.')
fig,ax=plt.subplots(2,2,figsize=(12.8,8.4),layout='constrained');rr=audit['records'];xx=np.arange(3);short=['Instant','sc 0.18','sc 0.24'];stiff=np.array([v['initial_generalized_stiffness_N_m'] for v in rr]);peaks=np.array([v['peak_generalized_reaction_N'] for v in rr])
ax[0,0].bar(xx-.16,stiff/stiff[0],.3,color=colors[1],label='Initial stiffness / reference');ax[0,0].bar(xx+.16,peaks/peaks[0],.3,color=colors[2],label='Protocol peak / reference');ax[0,0].set_xticks(xx,short);ax[0,0].set(ylim=(0,1.25),ylabel='Normalized value',title='A  Same intact response and protocol peak');ax[0,0].legend(loc='upper center',fontsize=9)
ax[0,1].plot(xx,[.15]*3,'o-',color='#333',label='Onset s0');ax[0,1].plot(xx,[.15,.18,.24],'s-',color=colors[2],label='Complete-damage parameter');ax[0,1].set_xticks(xx,short);ax[0,1].set(ylabel='Stretch',title='B  One independent softening interval');ax[0,1].legend()
diss=np.array([v['total_dissipation_J'] for v in rr])*1e9;anal=np.array([v['analytic_continuous_dissipation_J'] for v in rr])*1e9
ax[1,0].bar(xx,diss,color=colors,width=.6,label='Actual split update');ax[1,0].plot(xx,anal,'k_',ms=28,mew=2,label='Continuous analytical value')
for i,v in enumerate(diss):ax[1,0].text(i,v+5,f'{v:.2f}',ha='center',fontsize=11)
ax[1,0].set_xticks(xx,short);ax[1,0].set(ylim=(0,285),ylabel='Total dissipation (nJ)',title='C  Energy control demonstrated');ax[1,0].legend(fontsize=9)
ax[1,1].bar([0,1],[budget['work_upper_bound_J_m2'],6],color=['#187F94','#C1C8CC'],width=.5);ax[1,1].set_yscale('log');ax[1,1].set_xticks([0,1],['Tested opening\nwork upper bound','6 J/m²\nreference only']);ax[1,1].set(ylabel='J/m² (log axis)',title='D  High-energy target not admitted');ax[1,1].text(0,.04,f"{budget['work_upper_bound_J_m2']:.5f}",ha='center');ax[1,1].text(1,8,'Not achieved',ha='center');ax[1,1].set_ylim(.006,30)
fig.suptitle('Fixed s0 = 0.15: numerical energy control, limited physical scope',fontsize=17)
fig.supxlabel('Uniform dilation is not a tensile strength test. The opening budget is protocol-specific; no material Gc or 6–8 J/m² calibration.',fontsize=10)
save(fig,'Fig4_energy_control_and_scope','Same initial generalized stiffness and onset, distinct observed irreversible energies. Protocol peak equality does not certify macroscopic strength. Energy budget shown with clearly unachieved reference.')
inputs={str(p.relative_to(O)):sha(p) for p in O.glob('*.npz')};dump(O/'FIGURE_DATA_PROVENANCE.json',dict(script_sha256=sha(Path(__file__)),figures=figs,inputs=inputs,illustrative_fake_data=False,all_scientific_curves_from_saved_states=True))
links=''.join(f'<section><h2>{v["name"]}</h2><p>{v["caption"]}</p><a href="{v["name"]}.svg">SVG</a><br><img src="{v["name"]}.png"></section>' for v in figs)
(F/'OPEN_RESULTS.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>LPS 渐进软化离线可行性</title><style>body{max-width:1320px;margin:30px auto;font:17px/1.65 system-ui;color:#21333f;padding:0 22px}img{width:100%;height:auto}section{margin:36px 0;border-top:1px solid #bdcbd3;padding-top:12px}.note{background:#edf4f7;padding:18px;border-left:5px solid #187f94}</style><h1>LPS 渐进软化：真实离线结果</h1><p class="note">341 节点、12,047 原始键；一个变分柔度机制、两个软化终点。未损伤恢复、能量梯度、不可逆性与耗能可调通过。宏观强度、材料 Gc、血栓实验均未认证。未修改生产 LAMMPS。</p><p><a href="../../../verification/'+C['task']+'/'+C['run']+'/REVIEW_LPS_STRENGTH_ENERGY_DECOUPLED_SOFTENING_V1.md">中文审查与35项回答</a> · <a href="../../../verification/'+C['task']+'/'+C['run']+'/FINAL_STATUS.json">分项状态</a></p>'+links+'</html>')
(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/'+C['task']+'/'+C['run']+'/OPEN_RESULTS.html"><a href="../../../results/'+C['task']+'/'+C['run']+'/OPEN_RESULTS.html">打开四组科研图</a>')
dump(O/'FIGURE_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,groups=4,PNG=4,SVG=4))
print('FIGURES_COMPLETE',flush=True)
