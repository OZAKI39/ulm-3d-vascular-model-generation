from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import time,resource,html
T=time.perf_counter();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.titlesize':18,'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize':13,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#8594a0','grid.alpha':.65,'grid.linewidth':.7,'savefig.facecolor':'white','svg.fonttype':'none'})
blue='#2b7099';red='#b84c50';orange='#d88927';gray='#9aa7b1';purple='#9067a8';dark='#263746'
rt=np.load(O/'reference_topology.npz');d=np.load(O/'final_candidate.npz');sp=np.load(O/'spatial_response.npz');X=d['X'];x=d['x'];b=d['bonds'];alive=d['alive'].astype(bool);s=d['stretch'];u=x-X
a=json.loads((O/'STATIC_EQUILIBRIUM_AUDIT.json').read_text());n=json.loads((O/'NATIVE_FORCE_GRADIENT_CROSSCHECK.json').read_text());regions=json.loads((O/'NOTCH_VS_LOADING_END_STRETCH.json').read_text());resp=json.loads((O/'LOW_LOAD_COMPLIANCE_AND_ENERGY.json').read_text());solver=json.loads((O/'INDEPENDENT_EQUILIBRIUM_CANDIDATE.json').read_text());rec=json.loads((O/'RECOVERY_SOLVER_EXIT.json').read_text());accepted=a['accepted'];label='ACCEPTED 2 Pa equilibrium' if accepted else 'STATIC_STATE_NOT_ACCEPTED';files=[];samples={}
def save(fig,name):
 for ext in ['png','svg']:fig.savefig(F/(name+'.'+ext),dpi=240,bbox_inches='tight');files.append(name+'.'+ext)
 plt.close(fig)
def spatial(ax,wide=False):
 ax.set(xlabel='X (mm)',ylabel='Z (mm)',xlim=(-2.56,2.56),ylim=(-.9,.9) if wide else (-.6,.6));ax.set_aspect('equal');ax.axvline(-1,color=orange,lw=1.5,ls='--')
sheet=np.isclose(X[:,1],-25e-6,atol=1e-12);ids=np.flatnonzero(alive&sheet[b[:,0]]&sheet[b[:,1]]);samples['XZ_sheet_atom_ids']=np.flatnonzero(sheet)+1;samples['XZ_sheet_alive_bond_ids']=ids+1
fig=plt.figure(figsize=(18,9));fig.subplots_adjust(left=.03,right=.98,bottom=.23,top=.84,wspace=.36);ax=fig.add_subplot(121,projection='3d');ay=fig.add_subplot(122)
ax.scatter(*(X*1e3).T,s=3,c=gray,alpha=.2,depthshade=False,rasterized=True)
for mask,col in [(rt['high'],red),(rt['low'],blue)]:ax.scatter(*(X[mask]*1e3).T,s=9,c=col,alpha=.8,depthshade=False,rasterized=True)
subset=ids[::7];ax.add_collection3d(Line3DCollection(X[b[subset]]*1e3,colors=blue,alpha=.3,linewidths=.5));samples['Fig1_displayed_3D_alive_bond_ids']=subset+1
outline=np.array([[-2.5,-.2,0],[-1,-.2,0],[-1,.2,0],[-2.5,.2,0],[-2.5,-.2,0]]);ax.plot(*outline.T,c=red,lw=2.5)
for xx in [-2,-.8,.5,1.8]:
 ax.quiver(xx,0,.5,0,0,.22,color=red,arrow_length_ratio=.3,lw=2);ax.quiver(xx,0,-.5,0,0,-.22,color=blue,arrow_length_ratio=.3,lw=2)
ax.set(xlabel='X (mm)',ylabel='Y (mm), periodic',zlabel='Z (mm)',xlim=(-2.6,2.6),ylim=(-.24,.24),zlim=(-.8,.8),xticks=[-2.5,0,2.5],yticks=[-.2,.2],zticks=[-.5,0,.5]);ax.view_init(elev=17,azim=-64);ax.set_proj_type('ortho');ax.set_box_aspect((5,1.1,1.6));ax.xaxis.labelpad=20;ax.yaxis.labelpad=14;ax.zaxis.labelpad=15
for aa in [ax.xaxis,ax.yaxis,ax.zaxis]:aa.pane.set_facecolor((.96,.97,.98,.6));aa._axinfo['grid']['color']=(.45,.52,.58,.7)
ax.set_title('(a) All 16,000 original material points',pad=20)
ay.scatter(X[sheet,0]*1e3,X[sheet,2]*1e3,s=9,c=gray);ay.add_collection(LineCollection(X[b[ids]][:,:,[0,2]]*1e3,colors=blue,linewidths=.45,alpha=.3));ay.plot([-2.5,-1],[0,0],c=red,lw=3);ay.add_patch(Circle((-1,0),.1505,fill=False,ec=orange,lw=2));ay.axhspan(.3245,.50,color=red,alpha=.08);ay.axhspan(-.50,-.3245,color=blue,alpha=.08)
for k,offset in [(4710,(-.5,-.48)),(4871,(.3,.5))]:
 pos=X[k-1,[0,2]]*1e3;ay.annotate(f'atom {k}',xy=pos,xytext=pos+offset,arrowprops={'arrowstyle':'-','color':dark},fontsize=13,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
for xx in [.25,1.25,2.25]:
 ay.annotate('',xy=(xx,.79),xytext=(xx,.51),arrowprops={'arrowstyle':'-|>','color':red,'lw':2.2});ay.annotate('',xy=(xx,-.79),xytext=(xx,-.51),arrowprops={'arrowstyle':'-|>','color':blue,'lw':2.2})
ay.text(.75,.82,'2 Pa (+Z)',color=red,fontsize=13);ay.text(.75,-.88,'2 Pa (−Z)',color=blue,fontsize=13)
spatial(ay,True);ay.set_title('(b) Actual Y = −25 µm slice; alive bonds')
fig.suptitle('Frozen Mode-I geometry and uniform 2 Pa outward traction',fontsize=23,y=.97)
fig.text(.06,.12,'5 × 0.4 × 1 mm; Y periodic. Initial pre-crack: 1.50 mm; front X = −1.00 mm.\n800 loaded nodes per face × 5 nN = 4 µN per face. Arrow lengths are schematic.\nFull fixed reference m; 909,392 reference bonds, 16,416 prescribed initial cuts.',fontsize=16)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=red,label='Top load / initial cut'),Line2D([],[],color=blue,label='Bottom load / alive bonds'),Line2D([],[],color=orange,label='Notch-front region')],loc='lower center',bbox_to_anchor=(.5,.012),ncol=3,frameon=False)
save(fig,'Fig1_specimen_and_2Pa_loading')

q=np.genfromtxt(O/'STATIC_SOLVER_TRACE.csv',delimiter=',',names=True,dtype=None,encoding='utf-8');it=q['accepted_state'];fig,axs=plt.subplots(2,3,figsize=(19,10),layout='constrained')
panels=[('max_force_N',1e-13,'Max free-node residual (N)'),('relative_force',1e-5,'L2 residual / external L2'),('Pi_J',None,'Total potential Π (J)'),('energy_span_last3_relative',1e-8,'Last 3 accepted states: ΔU / U'),('position_span_last3_m',1e-10,'Last 3 accepted states: Δx (m)'),('global_torque_relative',1e-10,'Current residual torque / (Fface W)')]
for ax,(key,gate,txt) in zip(axs.flat,panels):
 vals=q[key].astype(float);ax.plot(it,vals,c=blue,lw=1.6,label='Main solve: accepted states')
 if gate is not None:ax.set_yscale('log');ax.axhline(gate,c=red,ls='--',lw=1.5,label='Frozen gate')
 ax.axvline(it[-1],c=orange,ls=':',lw=1.7);ax.set(xlabel='Accepted state index (not physical time)',ylabel=txt);ax.set_title(txt,fontsize=16)
axs[0,0].legend(fontsize=12,loc='lower left');fig.suptitle(f'{label} — finite main solve + one bounded recovery',fontsize=23)
fig.text(.03,-.065,f'Main: wall-budget stop at 240 s. Recovery: {rec["HVP_calls"]} HVP calls; CG limit 600; no accepted correction.\nOrange marker = main endpoint / recovery starting point. Recovery has no additional accepted-state curve.',fontsize=16,color=red)
save(fig,'Fig2_solver_convergence_and_recovery')

mag=np.linalg.norm(u,axis=1)*1e6;im=int(np.argmax(mag));fig,axs=plt.subplots(2,1,figsize=(17,10),layout='constrained');fig.set_constrained_layout_pads(h_pad=.22,hspace=.12)
for ax,coords,title in [(axs[0],X,'(a) Original A0 reference, Y = −25 µm'),(axs[1],x,'(b) Saved 2 Pa candidate, deformation ×1')]:
 ax.add_collection(LineCollection(coords[b[ids]][:,:,[0,2]]*1e3,colors=blue,linewidths=.4,alpha=.3));sc=ax.scatter(coords[sheet,0]*1e3,coords[sheet,2]*1e3,c=mag[sheet] if coords is x else np.zeros(sheet.sum()),s=14,cmap='viridis',vmin=0,vmax=mag.max(),rasterized=True);spatial(ax);ax.set_title(title)
fig.colorbar(sc,ax=axs,label='Displacement magnitude (µm)',shrink=.7,pad=.02)
axs[1].annotate(f'Max |u| atom {im+1}: {mag[im]:.4f} µm\n(actual 3D max projected onto XZ)',xy=x[im,[0,2]]*1e3,xytext=(-.4,.30),arrowprops={'arrowstyle':'->','color':dark},fontsize=14,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
fig.suptitle(f'{label}\nCandidate geometry only; no certified elastic response',fontsize=22,color=red if not accepted else blue);samples['Fig3_all_displacement_atom_ids']=np.arange(1,len(X)+1)
save(fig,'Fig3_real_candidate_displacement_x1')

rn=list(regions['regions']);short=['Notch front','Notch mouth','Loading ends','Other X surface','Interior','Y seam*'];fig,axs=plt.subplots(2,2,figsize=(18,11),layout='constrained');ax=axs[0,0]
values=[regions['regions'][k]['max_s_over_s00'] for k in rn];ax.bar(np.arange(6),values,color=[orange,gray,red,gray,blue,purple]);ax.axhline(.03/.15,c=red,ls='--',label='Small-strain ceiling / s00');ax.set(xticks=np.arange(6),xticklabels=short,ylabel='Maximum positive stretch / s00',title='(a) Unaccepted-candidate regional diagnostic');ax.tick_params(axis='x',rotation=20);ax.legend(fontsize=12)
ax=axs[0,1]
for k,col in [(0,orange),(2,red),(4,blue)]:
 mask=alive&(sp['bond_regions']==k);ax.hist(s[mask]/S00,bins=100,density=True,histtype='step',lw=2,color=col,label=short[k])
ax.set(xlabel='Bond stretch / s00',ylabel='Probability density',title='(b) Actual distributions of surviving bonds');ax.legend()
ax=axs[1,0];ax.scatter(X[sheet,0]*1e3,X[sheet,2]*1e3,s=4,color=gray,alpha=.5)
for name,col in [('notch_front',orange),('loading_ends',red)]:
 r=regions['regions'][name];j=r['reference_bond_id']-1;v=X[b[j]][:,[0,2]]*1e3;ax.plot(v[:,0],v[:,1],c=col,lw=4,marker='o',ms=5);offset=(.45,-.68) if name=='notch_front' else (1.1,.27);pos=np.array(r['reference_midpoint_m'])[[0,2]]*1e3;ax.annotate(f'Bond {j+1}; atoms {b[j,0]+1}–{b[j,1]+1}',xy=pos,xytext=pos+offset,arrowprops={'arrowstyle':'-','color':col},fontsize=12,color=col,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
spatial(ax,True);ax.set_title('(c) True IDs; projected reference positions')
ax=axs[1,1];ax.axis('off');ax.set_title('(d) Certified response withheld',loc='left')
texts=[f'Notch / load-end candidate ratio: {regions["notch_to_loading_end_max_stretch_ratio"]:.6f}',f'Saved U_LPS = {solver["last"]["U_LPS_J"]*1e12:.6f} pJ',f'Saved endpoint F·u = {solver["last"]["Fdotu_J"]*1e12:.6f} pJ','F·u is not integrated loading work.','Accepted opening / secant compliance: UNAVAILABLE','GC_EFFECTIVE: NOT MEASURED','* Y seam overlaps the other five regions.']
for y,txt in zip(np.linspace(.9,.06,len(texts)),texts):ax.text(0,y,txt,fontsize=15,color=red if 'UNAVAILABLE' in txt or 'NOT MEASURED' in txt else dark)
fig.suptitle('STATIC_STATE_NOT_ACCEPTED — no certified compliance or loading curve',fontsize=23,color=red);save(fig,'Fig4_candidate_stretch_and_withheld_compliance')

fig,ax=plt.subplots(figsize=(16,10),layout='constrained');ax.axis('off');table=[['Small-strain admission',f'PASS: s_max = {a["max_positive_stretch"]:.8f} ≤ 0.03'],['Absolute free-node force',f'FAIL: {a["native_max_force_N"]:.6g} N > 1e-13 N'],['Relative free-node force',f'FAIL: {a["native_relative_force"]:.6g} > 1e-5'],['Energy / displacement platform','PASS last 3; last-10 energy audit does not pass'],['Current torque balance','FAIL frozen gate; see periodic-domain scope'],['Native force vs independent gradient',f'PASS: relative error {n["native_vs_independent_force_relative"]:.3g}'],['Low-load static baseline','NOT ESTABLISHED'],['Notch vs loading-end ranking','Candidate trend only; not a credible equilibrium result'],['Shear crossing times: H / W','1.031 ms / 5.157 ms; 24 ms is NOT certified quasi-static'],['New breaks / formal fracture runs','0 / 0'],['G0_THEORY / G0_DISCRETE','0.0717007083 / 0.0696962373 J/m² (nominal only)'],['GC_EFFECTIVE','NOT MEASURED; not zero']]
tbl=ax.table(cellText=table,colLabels=['Independent scientific item','Actual finding'],cellLoc='left',colLoc='left',loc='center',colWidths=[.40,.60]);tbl.auto_set_font_size(False);tbl.set_fontsize(15);tbl.scale(1,2.75)
for (r,c),cell in tbl.get_celld().items():
 cell.set_edgecolor('#d4dde4')
 if r==0:cell.set_facecolor(dark);cell.set_text_props(color='white',weight='bold')
 elif r in [2,3,5,7,12]:cell.set_facecolor('#faeeee')
 else:cell.set_facecolor('#f3f7f9')
ax.set_title('No accepted static baseline; material and production solver unchanged',fontsize=23,pad=24);save(fig,'Fig5_scientific_admission_and_limits')
np.savez_compressed(O/'FIGURE_DISPLAY_SUBSETS.npz',**samples)
dump(O/'FIGURE_DATA_PROVENANCE.json',dict(figures=files,data_sources={p.name:sha(p) for p in [O/'reference_topology.npz',O/'final_candidate.npz',O/'STATIC_SOLVER_TRACE.csv',O/'spatial_response.npz',O/'STATIC_EQUILIBRIUM_AUDIT.json']},exact_display_IDs='FIGURE_DISPLAY_SUBSETS.npz',scientific_status=label,PNG_DPI=240,no_synthetic_fracture_or_shell=True,all_displacement_scales=1))
caps=['原始16000点、存活键和2 Pa载荷','主求解及一次恢复的停止证据','未接受端点的真实倍率位移','未接受端点伸长诊断；柔度不认证','逐项准入与未测量范围'];cards=''.join(f'<section><h2>{cap}</h2><a href="{name}.png">PNG</a> · <a href="{name}.svg">SVG</a><img src="{name}.png" alt="{cap}"></section>' for name,cap in zip([Path(z).stem for z in files if z.endswith('.png')],caps))
page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>2 Pa Mode-I 静态准入审查</title><style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:30px auto;padding:0 24px;color:#263746;background:#f7f9fb}section{background:white;padding:24px;margin:24px 0;border:1px solid #d6dde4}img{width:100%}p{line-height:1.8}a{color:#2b7099}</style><h1>2 Pa Mode-I：静态基线尚未建立</h1><p>原始试样与真实节点外力核对通过；最大伸长0.015379满足小变形门槛，但绝对/相对力与冻结力矩门槛未通过。主求解240秒预算停止，一次恢复到CG 600次上限。没有追加求解、扫描或动态断裂。<b>GC_EFFECTIVE仍未测量。</b></p>'+cards+'</html>'
(F/'OPEN_RESULTS.html').write_text(page);rel='../../../results/'+C['task']+'/'+C['run']+'/';(O/'OPEN_RESULTS.html').write_text(page.replace('href="Fig','href="'+rel+'Fig').replace('src="Fig','src="'+rel+'Fig'))
dump(O/'FIGURE_RESOURCES.json',dict(wall_s=time.perf_counter()-T,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,new_solves=0));print('FIGURES_READY',len(files),flush=True)
