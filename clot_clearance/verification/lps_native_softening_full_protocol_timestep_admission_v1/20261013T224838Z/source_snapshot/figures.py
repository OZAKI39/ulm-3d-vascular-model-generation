"""Six publication-style groups drawn exclusively from saved actual states."""
from common import *
import time,resource
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
from mpl_toolkits.mplot3d.art3d import Line3DCollection
start=time.perf_counter()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':16,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.grid':True,'grid.color':'#607080','grid.alpha':.35,'grid.linewidth':.7,'figure.facecolor':'white','savefig.facecolor':'white'})
blue='#187F94';orange='#C16A29';red='#C34049';purple='#7756A2';gray='#7F95A4';green='#3B8262'
z=np.load(O/'trajectory.npz');a=z['series'];c={n:a[:,i] for i,n in enumerate(z['columns'])};t=c['time_s']*1e6;end=t[-1];x=z['x'];vel=z['velocity'];b=z['bonds'];X=z['reference_X'];ref=np.load(O/'provenance/reference.npz');ev=json.loads((O/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json').read_text());stop=json.loads((O/'NUMERICAL_SAFETY_STOP_RECORD.json').read_text());audit=json.loads((O/'FULL_PROTOCOL_NATIVE_TIMESTEP_AUDIT.json').read_text());top=json.loads((O/'FINAL_TOPOLOGY_AND_COMPONENTS.json').read_text());ct=json.loads((O/'DYNAMIC_CONTACT_CLOSURE_AUDIT.json').read_text());comparison=json.loads((O/'HISTORICAL_DT_COMPARISON.json').read_text());firstsoft=ev['first_softening']['time_s']*1e6;firstfail=ev['first_complete_failure']['time_s']*1e6
outputs=[]
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=190,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);outputs.append(dict(name=name,caption=caption))
def timeaxis(ax,planned=True):
 ax.set_xlabel('Physical time (µs)');ax.set_xlim(0,700 if planned else end)
 for q in [300,350]:ax.axvline(q,color=gray,lw=.8,ls=':')
 ax.axvline(end,color=red,lw=1.2,ls='--')
 if planned:ax.axvspan(end,700,color=gray,alpha=.12)
def heading(fig,text,sub):
 fig.suptitle(text,fontsize=21,fontweight='bold');fig.text(.5,-.015,sub,ha='center',fontsize=12)

fig,axs=plt.subplots(2,2,figsize=(16.5,10),layout='constrained');axs=axs.ravel();right=z['right_ids'][0]-1;left=z['left_ids'][0]-1
axs[0].plot(t,(x[:,right,0]-X[right,0])*1e6,color=blue,label='Right grip');axs[0].plot(t,(x[:,left,0]-X[left,0])*1e6,color=orange,label='Left grip');axs[0].set(ylabel='Prescribed x displacement (µm)',title='Actual analytic boundary motion');axs[0].legend()
axs[1].plot(t,c['max_free_disp_m']*1e6,color=blue,lw=2);axs[1].set(ylabel='Maximum free-node displacement (µm)',title='245 free nodes: native integration')
axs[2].plot(t,c['max_free_speed_m_s'],color=purple,lw=2,label='Free nodes');axs[2].plot(t,np.linalg.norm(vel[:,right],axis=1),color=orange,label='Right grip');axs[2].set(ylabel='Speed (m/s)',title='Boundary rest does not imply equilibrium');axs[2].legend()
axs[3].plot(t,c['K_J']*1e9,color=purple,label='All-node kinetic energy');axs[3].plot(t,c['Kfree_J']*1e9,color=blue,label='Free-node kinetic energy');axs[3].set(ylabel='Kinetic energy (nJ)',title='Inertia persists during unloading');axs[3].legend()
for ax in axs:
 timeaxis(ax)
 for q,color in [(firstsoft,orange),(firstfail,purple)]:ax.axvline(q,color=color,ls=':',lw=1)
heading(fig,'1  |  Actual loading and free-node dynamics',f'First softening {firstsoft:g} µs · First failure {firstfail:g} µs · STOP {end:g} µs · Shaded interval NOT_RUN')
save(fig,'Fig1_loading_and_dynamics','Only actual native states; no continuation through the unexecuted final interval.')

fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained');axs=axs.ravel()
axs[0].plot(t,c['max_dd'],color=blue,lw=1);axs[0].axhline(.1,color=red,ls='--',label='Frozen limit 0.10');axs[0].scatter([end],[c['max_dd'][-1]],s=60,color=red,zorder=5);axs[0].set(ylabel='Maximum step Δd',title=f'FAILED: {c["max_dd"][-1]:.6f} at {end:g} µs');axs[0].legend()
axs[1].plot(t,c['new_fail'],color=orange,lw=1);axs[1].set(ylabel='New complete failures / step',title='Unique original bond events')
axs[2].plot(t,c['total_failed'],color=red,lw=2);axs[2].set(ylabel='Cumulative failed bonds',title=f'{int(c["total_failed"][-1]):,} failed; {int(c["live_bonds"][-1]):,} alive')
axs[3].plot(t,c['step_disp_dx'],color=blue,label='Node step / dx');axs[3].plot(t,c['relative_endpoint_step_dx'],color=purple,label='Endpoint relative step / dx');axs[3].axhline(.05,color=blue,ls='--',lw=1,label='Node gate 0.05');axs[3].axhline(.1,color=purple,ls='--',lw=1,label='Relative gate 0.10');axs[3].set(ylabel='Step displacement / dx',title='Displacement gates remained satisfied');axs[3].legend(fontsize=10)
bid=stop['violating_bonds'][0]['bond_id'];k=bid-1;stretch=np.linalg.norm(x[:,b[k,1]]-x[:,b[k,0]],axis=1)/ref['r0'][k]-1;kappa=np.maximum.accumulate(np.maximum(stretch,0));damage=np.zeros(len(t));part=(kappa>S0)&(kappa<SC);damage[part]=SC*(kappa[part]-S0)/(kappa[part]*(SC-S0));damage[kappa>=SC]=1
axs[4].plot(t,damage,color=red,marker='.',ms=6);axs[4].scatter(t[-2:],damage[-2:],color=[blue,red],s=60,zorder=4);axs[4].set(xlim=(end-2,end+.12),xlabel='Physical time (µs)',ylabel=f'Original bond {bid}: d',title='Late softening event near onset')
axs[5].plot(t,stretch,color=blue,label='Current stretch');axs[5].plot(t,kappa,color=orange,ls='--',label='Irreversible κ');axs[5].axhline(.15,color=red,ls=':',label='s₀ = 0.15');axs[5].set(xlim=(end-2,end+.12),ylim=(.125,.157),xlabel='Physical time (µs)',ylabel=f'Original bond {bid}: stretch / κ',title='Native history matches independent model');axs[5].legend(fontsize=10)
for ax in axs[:4]:timeaxis(ax)
for ax in axs[4:]:ax.xaxis.set_major_locator(MultipleLocator(.5));ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
heading(fig,'2  |  Damage increment and frozen timestep gates',f'Bond {bid} (nodes 24, 245): d = 0.000281 → 0.113916 in one 0.25 µs step. Violating point is retained.')
save(fig,'Fig2_damage_timestep_gate','Safety stop is a late new-softening event; no curve clipping or threshold relaxation.')

keyidx=json.loads((O/'KEYFRAME_INDEX.json').read_text());keypaths={p['step']:p['path'] for p in keyidx};firstcontact=ct['events'][0]['step'];selected=[(0,'Initial'),(ev['first_softening']['step'],'First softening'),(ev['first_complete_failure']['step'],'First complete failure'),(firstcontact,'First closure contact'),(1400,'Unloading starts'),(len(t)-1,'Safety-stop state')]
fig=plt.figure(figsize=(19,12));cmap=LinearSegmentedColormap.from_list('damage',[blue,'#EDC795',red]);norm=Normalize(0,1);bound=np.ceil(np.max(abs(x))*1e6/50)*50+25;axs=[]
for n,(step,title) in enumerate(selected):
 ax=fig.add_subplot(2,3,n+1,projection='3d');axs.append(ax);q=np.load(O/keypaths[step]);xx=q['x']*1e6;d=q['d'];nodal=np.zeros(len(X));np.maximum.at(nodal,b[:,0],d);np.maximum.at(nodal,b[:,1],d);short=ref['r0']<=1.01*DX
 for mask,color,width,alpha,ls in [(short&(d==0),gray,.6,.25,'-'),(short&(d>0)&(d<1),orange,1.0,.6,'-'),(short&(d==1),red,.8,.28,'--')]:
  if np.any(mask):ax.add_collection3d(Line3DCollection(xx[b[mask]],colors=color,linewidths=width,alpha=alpha,linestyles=ls))
 ax.scatter(*xx.T,c=nodal,cmap=cmap,norm=norm,s=21,edgecolors='white',linewidths=.25,depthshade=False)
 if step==firstcontact:
  pair=np.array(ct['events'][0]['pairs'][0])-1;ax.plot(*xx[pair].T,color='#063CCB',lw=3);ax.scatter(*xx[pair].T,s=85,facecolors='none',edgecolors='#063CCB',lw=2)
 ax.set(xlim=(-bound,bound),ylim=(-bound,bound),zlim=(-bound,bound),xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)');ax.set_box_aspect((1,1,1));ax.view_init(elev=17,azim=-58);ax.set_title(f'{title} | {t[step]:g} µs\n{int(np.sum(d==1)):,} fully failed bonds',fontsize=14,pad=6)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_pane_color((.97,.98,.99,1));axis._axinfo['grid'].update(color='#8E9CA8',linewidth=.7)
 ax.tick_params(labelsize=10,pad=1)
fig.subplots_adjust(left=.01,right=.90,top=.92,bottom=.10,wspace=.04,hspace=.10);cbax=fig.add_axes([.93,.21,.013,.57]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax,label='Node color: maximum incident bond damage')
fig.legend(handles=[Line2D([0],[0],color=gray,label='Intact'),Line2D([0],[0],color=orange,label='Partially softened'),Line2D([0],[0],color=red,ls='--',label='Failed reference-edge trace'),Line2D([0],[0],color='#063CCB',lw=3,label='First actual contact pair')],loc='lower center',bbox_to_anchor=(.49,.038),ncol=4,frameon=False,fontsize=11)
fig.suptitle('3  |  Actual 341-node states: fixed physical axes and camera',fontsize=21,fontweight='bold',y=.99);fig.text(.49,.012,'Edges displayed only for original r₀ ≤ 1.01 dx; all 12,047 bonds are audited. Dashed traces are not live bonds. No reconstructed shell.',ha='center',fontsize=11)
save(fig,'Fig3_actual_point_cloud','True x and native d at six observed events. Fixed camera and equal physical scale. No unseen unloading/end state.')

fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained');axs=axs.ravel()
for name,label,color in [('W_J','Boundary work W',blue),('U_J','Variational U',green),('K_J','Physical K',purple),('D_J','Irreversible D',orange)]:axs[0].plot(t,c[name]*1e9,label=label,color=color)
axs[0].set(ylabel='Energy (nJ)',title='W = ΔU + ΔK + D + ΔUc − Jc + R');axs[0].legend(fontsize=11)
axs[1].plot(t,c['Uc_J']*1e12,color=blue,label='Contact Uc');axs[1].plot(t,c['Jc_J']*1e12,color=orange,ls='--',label='Topology Jc');axs[1].set(ylabel='Energy (pJ)',title='Actual conservative contact');axs[1].legend()
axs[2].plot(t,c['R_J']*1e12,color=red,label='Signed residual R');axs[2].plot(t,-c['switch_work_J']*1e12,color=purple,ls='--',label='− Damage-force switch work');axs[2].set(ylabel='Numerical balance (pJ)',title='Residual retained, not added to D');axs[2].legend(fontsize=10)
for name,label,color in [('VV_defect_J','VV quadrature',blue),('spatial_defect_J','Fixed-d spatial quadrature',orange),('contact_defect_J','Contact quadrature',purple)]:axs[3].plot(t,c[name]*1e12,color=color,label=label)
axs[3].set(ylabel='Cumulative defect (pJ)',title='Smaller independent time-layer terms');axs[3].legend(fontsize=10)
axs[4].plot(t,c['relative_R']*100,color=blue,label='Actual |R| / frozen scale');axs[4].axhline(2,color=red,ls='--',label='Frozen 2% gate');axs[4].set(ylabel='Relative energy balance (%)',title=f'Observed maximum {100*c["relative_R"].max():.4f}%');axs[4].legend(fontsize=10)
axs[5].axis('off');tab=[['W',f'{c["W_J"][-1]*1e9:.9f} nJ'],['U',f'{c["U_J"][-1]*1e9:.9f} nJ'],['K',f'{c["K_J"][-1]*1e9:.9f} nJ'],['D',f'{c["D_J"][-1]*1e9:.9f} nJ'],['Uc',f'{c["Uc_J"][-1]*1e9:.9f} nJ'],['Jc','0 J (measured)'],['R',f'{c["R_J"][-1]*1e12:.9f} pJ']];table=axs[5].table(cellText=tab,colLabels=['Stop-state term','Actual value'],loc='center',cellLoc='left',colWidths=[.35,.65]);table.auto_set_font_size(False);table.set_fontsize(12);table.scale(1,1.8);axs[5].set_title(f'{end:g} µs: stop-state ledger')
for ax in axs[:5]:timeaxis(ax)
heading(fig,'4  |  Complete variational energy and time-layer defects','Boundary path work includes exact grip ΔK. Signed topology = −D is display-only; damping is absent. No full700µs certification.')
save(fig,'Fig4_energy_and_time_layer','Every recorded step checked against frozen complete variational energy. Physical dissipation and numerical force-switch work are separate.')

old=np.load(OLD/'main_accepted_trajectory.npz');oldseries=old['series'];oc={n:oldseries[:,i] for i,n in enumerate(old['columns'])};ot=oc['time_s']*1e6;lim=ot[-1];mask=t<=lim+1e-6;overlap=np.load(O/'historical_overlap_comparison.npz')['series'];fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained');axs=axs.ravel()
for ax,on,nn,label,title,scale in [(axs[0],'total_failed','total_failed','Fully failed bonds','Natural failure count',1),(axs[1],'vmax_m_s','max_speed_m_s','Maximum speed (m/s)','Actual dynamic velocity',1),(axs[2],'max_dd','max_dd','Maximum step Δd','Different-sized physical steps',1),(axs[3],'residual_J','R_J','Signed residual (pJ)','Complete energy balance',1e12)]:
 ax.plot(ot,oc[on]*scale,color=orange,label='Saved dt = 0.50 µs');ax.plot(t[mask],c[nn][mask]*scale,color=blue,ls='--',label='New dt = 0.25 µs');ax.set(xlim=(0,lim),xlabel='Physical time (µs)',ylabel=label,title=title);ax.legend(fontsize=10)
axs[2].axhline(.1,color=red,lw=1,ls=':');axs[4].plot(overlap[:,0]*1e6,overlap[:,3]*1e9,color=purple,label='Largest coordinate component');axs[4].plot(overlap[:,0]*1e6,overlap[:,4]*1e9,color=blue,ls='--',label='Largest nodal vector difference');axs[4].set(xlabel='Identical physical time (µs)',ylabel='Position difference (nm)',title='No temporal interpolation');axs[4].legend(fontsize=10)
axs[5].axis('off');axs[5].text(0,.95,'Actual overlap: 0–122 µs\n\nFirst softening: 85 µs in both\nFirst complete failure: 92.5 µs in both\nSame seven first-failure IDs observed\n\nAt 122 µs: 520 failed in both\nFailed-ID sets identical at this endpoint\n\nMaximum component difference:\n'+f'x = {comparison["maximums"]["x_max_component_m"]*1e9:.6f} nm\nv = {comparison["maximums"]["v_max_component_m_s"]:.9f} m/s'+'\n\nNo old simulation re-executed.',va='top',fontsize=13)
heading(fig,'5  |  Saved 0.50 µs versus fresh 0.25 µs: actual overlapping times','TIMESTEP ADMISSION ≠ GLOBAL TIME CONVERGENCE  ·  Old short restart window is not used as a full-path reference')
save(fig,'Fig5_historical_timestep_comparison','Only matching physical times are compared; neither event sequences nor nonlinear trajectories are forced to coincide.')

fig,axs=plt.subplots(1,2,figsize=(17.5,10),gridspec_kw={'width_ratios':[1.35,1]},layout='constrained');axs[0].axis('off');status=[('Actual native free-node dynamics','PASS: observed window',green),('Full 700 µs / 2800-step protocol','FAILED_ADMISSION',red),('Maximum Δd ≤ 0.10','FAILED: 0.113635',red),('Complete energy ledger ≤ 2%','PASS: max 0.3724%',green),('Irreversible κ / d and momentum','PASS: observed window',green),('Unloading coverage','PARTIAL: 350–640.5 µs',orange),('Dynamic closure contact','OBSERVED: numerical checks pass',blue),('Contact penetration accuracy','NOT CERTIFIED: max 24.7% dx',orange),('Complete physical separation','NOT OBSERVED: 1 component',orange),('Natural area / material Gc','UNAVAILABLE',gray),('Global time convergence','NOT CERTIFIED',orange),('Next P4 specimen admission','NOT GRANTED; NOT_RUN',red)]
for i,(label,val,color) in enumerate(status):
 yy=.95-i*.077;axs[0].text(.0,yy,label,fontsize=12.5,fontweight='bold',va='top');axs[0].text(.52,yy,val,color=color,fontsize=12.5,va='top')
axs[0].set_title('Separate scientific admission decisions',pad=18)
axs[1].axis('off');axs[1].text(.02,.94,'OBSERVED SAFETY STOP',fontsize=18,fontweight='bold',color=red,va='top');axs[1].text(.02,.865,'Step 2562 · 640.5 µs\nBond 1813 · free nodes 24 and 245\nκ: 0.150007 → 0.152903\nd: 0.000281 → 0.113916\n\nNative and independent histories agree.\nNo direct contact force on these endpoints\nat the stop event.\n\nNext single priority:\nDamage-event-resolved time integration\nwith consistent force and energy accounting.\n\nNo automatic dt = 0.125 µs run.\nNo material retuning or added damping.\n\nSame verified C++ library; no recompilation.\nProduction and archived results preserved.',fontsize=14,va='top',linespacing=1.55)
heading(fig,'6  |  Scientific evidence and admission to further development','Program stopped correctly at the frozen gate. This is not a successfully admitted complete loading–unloading baseline.')
save(fig,'Fig6_scientific_admission','Implementation consistency, protocol failure, observed contact, convergence and material calibration remain separate decisions.')
dump(O/'FIGURE_MANIFEST.json',dict(figures=outputs,render_wall_s=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,source='Native actual trajectory and immutable saved historical overlap only',no_synthetic_shell=True,no_video=True))
report='REVIEW_LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1.md'
html='<!doctype html><html lang="zh"><meta charset="utf-8"><title>0.25 µs 全时程准入审查</title><style>body{font:18px system-ui;max-width:1550px;margin:36px auto;padding:0 20px;background:#f6f8fa;color:#22313a}img{width:100%;background:white}section{margin:40px 0}a{color:#187f94}.stop{color:#ad2d39}</style><h1>341 节点原生软化：0.25 µs 完整协议准入</h1><p class="stop"><b>FAILED_ADMISSION：</b>第2562步、640.5 µs，原始键1813的Δd=0.113635超过0.10，立即停止。最后合格状态640.25 µs；700 µs未完成。</p><p>完整力/能量/历史与动量在已观察窗口内一致；动态接触已发生，最大重叠24.7% dx的物理精度仍未认证。材料Gc、新生裂面面积及宏观强度不可用。</p>'
html+=f'<p><a href="../../../verification/{TASK}/{RUN}/{report}">完整中文报告（36问）</a> · <a href="../../../verification/{TASK}/{RUN}/FINAL_STATUS.json">独立验收状态</a> · <a href="../../../verification/{TASK}/{RUN}/NUMERICAL_SAFETY_STOP_RECORD.json">原始越界现场</a></p>'
for i,q in enumerate(outputs,1):html+=f'<section><h2>图{i}</h2><p>{q["caption"]}</p><a href="{q["name"]}.png">PNG</a> · <a href="{q["name"]}.svg">SVG</a><img src="{q["name"]}.png" loading="lazy"></section>'
(F/'OPEN_RESULTS.html').write_text(html+'</html>');(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/{TASK}/{RUN}/OPEN_RESULTS.html"><a href="../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">打开六组科研图</a>')
print('SIX_FIGURE_GROUPS_RENDERED',F,flush=True)
