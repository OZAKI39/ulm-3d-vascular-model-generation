"""Actual-result reduction and publication-style figures; no synthetic fracture geometry."""
from common import *
from traction_library import Couette,GAP,MU,RAMP
from traction_to_v3 import FrozenV3Provider,high_fidelity_gate
import csv,time,resource,html,platform,importlib.metadata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,cm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
def read(p):return json.loads(Path(p).read_text())
def history(p):return np.genfromtxt(p,delimiter=',',names=True)
def writecsv(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def style():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':16,'axes.labelsize':14,'legend.fontsize':11,'xtick.labelsize':12,'ytick.labelsize':12,'figure.facecolor':'white','axes.facecolor':'#fafafa','grid.color':'#8b929a','grid.alpha':.55,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220,'svg.fonttype':'none'})
def axis3(ax,X,unit='mm'):
 ax.set(xlabel=f'X ({unit})',ylabel=f'Y ({unit})',zlabel=f'Z ({unit})');ax.set_box_aspect(np.maximum(np.ptp(X,axis=0),.05));ax.view_init(19,-62)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.set_major_locator(MaxNLocator(4));a.labelpad=12;a.pane.set_facecolor((.95,.96,.97,1));a._axinfo['grid']['color']=(.45,.48,.52,.65)
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);return dict(name=name,caption=caption)
def figures(cases,rows,best,mevo,mapping):
 style();gallery=[];lib=load_npz(O/'PHYSICAL_TRACTION_LIBRARY.npz');mapped=load_npz(O/'V3_MAPPED_FORCES.npz');peak=int(np.argmax(lib['lower_wall_shear_Pa']));tri=lib['triangles_xyz_m']*1e3;cent=lib['face_centroid_m']*1e3;n=lib['n_clot'];traction=lib['total_traction_increment_Pa'][peak];mag=np.linalg.norm(traction,axis=1)
 fig=plt.figure(figsize=(15,7.5));fig.subplots_adjust(top=.78,wspace=.25);ax=fig.add_subplot(121,projection='3d');norm=colors.Normalize(0,.16);poly=Poly3DCollection(tri,facecolors=cm.viridis(norm(mag)),edgecolor='#73777c',linewidth=.08,alpha=.8);ax.add_collection3d(poly);ax.auto_scale_xyz(*cent.T)
 ix=np.arange(0,len(cent),55);vec=traction[ix]/.16*.12;ax.quiver(*cent[ix].T,*vec.T,color='#0a4264',arrow_length_ratio=.25,linewidth=.9);axis3(ax,cent);ax.set_title('V3 exposed CFD surface\nPHYSICS_APPROXIMATION',fontsize=14);ax.text2D(.02,.03,'0.16 Pa → 0.12 mm arrow',transform=ax.transAxes,fontsize=11)
 ax=fig.add_subplot(122,projection='3d');X=mapped['reference_coordinates_m']*1e3;fn=mapped['force_N'];normf=np.linalg.norm(fn,axis=1);ax.scatter(*X.T,c='#87939e',s=9,alpha=.55);sel=np.flatnonzero(normf>0)[::8];v=fn[sel]/max(normf.max(),1e-30)*.12;ax.quiver(*X[sel].T,*v.T,color='#a32638',linewidth=.8,arrow_length_ratio=.3);axis3(ax,X);ax.set_title('Mapped original material nodes\nFROZEN_SURFACE_ONE_WAY',fontsize=14);ax.text2D(.02,.03,f'{normf.max()*1e9:.3f} nN → 0.12 mm arrow',transform=ax.transAxes,fontsize=11);fig.suptitle('Frozen reference template: U = 4 mm/s, H = 25 μm',fontsize=17,y=.96);gallery.append(save(fig,'Fig01_V3_faces_and_nodes','Actual exposed CFD triangles and original 1927 material nodes. Traction is a local plane approximation; no new V3 fluid or solid solve. Signed affine weights are used.'))
 fig,axes=plt.subplots(1,2,figsize=(14,5));c=Couette(.004,.0008);y=np.linspace(0,GAP,80);xx=np.linspace(-100e-6,100e-6,90);uu=np.tile(c.velocity(y,.0004)[:,None],(1,len(xx)));im=axes[0].pcolormesh(xx*1e6,y*1e6,uu*1e3,shading='auto',cmap='viridis');axes[0].streamplot(xx*1e6,y*1e6,uu,np.zeros_like(uu),color='white',density=.7,linewidth=.8);fig.colorbar(im,ax=axes[0],label='Flow speed (mm/s)');axes[0].set(xlabel='Tangential position (μm)',ylabel='Fluid gap coordinate (μm)',title='Computed plane Stokes flow at 0.4 ms\nPHYSICS_APPROXIMATION');tt=np.linspace(0,c.end,600)
 axes[1].plot(tt*1e3,c.shear(tt),color='#a32638',lw=2.2,label='Lower-wall shear');axes[1].axhline(MU*.004/GAP,ls='--',color='#6a737d',label='Steady μU/H');axes[1].axvline(.8,ls=':',color='#444');axes[1].set(xlabel='Time (ms)',ylabel='Shear traction (Pa)',title='Resolved startup and decay\np increment = 0 Pa');axes[1].grid();axes[1].legend();fig.tight_layout();gallery.append(save(fig,'Fig02_Stokes_flow_and_startup','Infinite parallel planes with a stationary lower wall and moving upper wall. This is not a computed bubble flow; the finite specimen edge and V3 curvature are not resolved.'))
 fig=plt.figure(figsize=(14,6))
 for j,(val,title) in enumerate([(mag,'Tangential shear magnitude (Pa)'),(np.zeros(len(mag)),'Pressure increment = 0 Pa')]):
  ax=fig.add_subplot(1,2,j+1,projection='3d');poly=Poly3DCollection(tri,facecolors=cm.viridis(norm(val)),edgecolor='#6f7880',linewidth=.08);ax.add_collection3d(poly);ax.auto_scale_xyz(*cent.T);axis3(ax,cent);ax.set_title(title+'\nPHYSICS_APPROXIMATION')
 fig.subplots_adjust(left=.02,right=.84,wspace=.32);cbax=fig.add_axes([.92,.24,.015,.5]);fig.colorbar(cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cbax,label='Increment (Pa)');gallery.append(save(fig,'Fig03_surface_shear_pressure','Zero pressure is the explicit plane-flow result, not missing pressure replaced with a fitted field. Historical steady background pressure is not added again.'))
 fig,ax=plt.subplots(1,2,figsize=(14,5));nums=np.array([r['trial'] for r in rows]);u=np.array([r['drive_speed_m_s'] for r in rows])*1e3;dur=np.array([r['hold_until_s'] for r in rows])*1e3;obj=np.array([r['objective'] for r in rows])
 sc=ax[0].scatter(u,dur,c=obj*1e6,cmap='plasma',s=150,edgecolors='#333');fig.colorbar(sc,ax=ax[0],label='Frozen proxy ×10⁶')
 for i,x,y in zip(nums,u,dur):ax[0].annotate(str(i),(x,y),xytext=(7,7),textcoords='offset points')
 ax[0].set(xlabel='Upper-plane drive (mm/s)',ylabel='Start of switch-off (ms)',title='Actual optimizer candidates\nLOW_FIDELITY_PROXY');ax[0].grid();ax[1].bar(nums,obj*1e6,color=['#447b9c']*3+['#b2474a']);ax[1].set(xlabel='Optuna trial number',ylabel='Frozen proxy ×10⁶',title='3 random + 1 BoTorch evaluation\nNo damage: strain precursor only');ax[1].grid(axis='y');fig.tight_layout();gallery.append(save(fig,'Fig04_optimizer_feedback','Numbers are persisted Optuna trial IDs. Trial 0 failed during initialization and has no score. Samples are insufficient to validate GP accuracy or an optimum.'))
 ref=load_npz(SPEC/'reference.npz');state=load_npz(best['path']/'keys/final.npz');X=state['x']*1e6;nb=int(ref['clot_count']);b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];plane=np.abs(ref['X'][:,1])<=DX*.51;draw=(bt==1)&plane[b[:,0]]&plane[b[:,1]];live=draw&(state['d']<1);broken=draw&cut
 fig=plt.figure(figsize=(15,6));ax=fig.add_subplot(121,projection='3d');ax.scatter(*X[:nb].T,s=7,c='#972b45',alpha=.4,label='Attached material');ax.scatter(*X[nb:].T,s=5,c='#646c75',alpha=.3,label='Fixed anchor');ax.add_collection3d(Line3DCollection(X[b[live]],colors='#627a8b',linewidths=.35,alpha=.45));ax.add_collection3d(Line3DCollection(X[b[broken]],colors='#e09922',linewidths=.6,linestyles='dashed',alpha=.6));axis3(ax,X,'μm');ax.set_title('Actual native final point/bond state\nACTUAL_SOLVE / small specimen');ax.legend(loc='upper left');ax=fig.add_subplot(122);h=history(best['path']/'history.csv');ax.plot(h['time_s']*1e3,h['max_tensile_kappa'],color='#27647d',label='Maximum new tensile history');ax.axhline(.15,color='#a32638',ls='--',label='Bulk softening onset');ax.set_yscale('symlog',linthresh=1e-6);ax.set(xlabel='Time (ms)',ylabel='Tensile κ (dimensionless)',title='No newly softened or failed bonds\nNo actual free component');ax.grid();ax.legend();fig.tight_layout();gallery.append(save(fig,'Fig05_actual_small_specimen','Original point coordinates and IDs; no displacement magnification, fabricated shell or crack surface. For readability, displayed bonds lie in the central y slice. Orange dashed bonds are the unchanged original precrack, not new damage.'))
 fig,ax=plt.subplots(1,2,figsize=(14,5));h=history(best['path']/'history.csv')
 for field,label,color in [('W_J','External work','#a32638'),('U_J','Stored energy','#246580'),('K_J','Kinetic energy','#e49b2a'),('D_J','Damage dissipation','#363b40')]:ax[0].plot(h['time_s']*1e3,h[field]*1e12,label=label,color=color,lw=1.8)
 ax[0].set(xlabel='Time (ms)',ylabel='Energy (pJ)',title='Actual native energy ledger\nACTUAL_SOLVE');ax[0].legend();ax[0].grid();ax[1].plot(h['time_s']*1e3,h['relative_R']*100,color='#6c437b');ax[1].set(xlabel='Time (ms)',ylabel='Relative balance residual (%)',title='W − (ΔU + ΔK + D + ΔUc − Jc)\nResidual is not damage dissipation');ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[1].grid();fig.tight_layout();gallery.append(save(fig,'Fig06_actual_energy','Trapezoidal actual-force work on actual native displacement; original contact activation term retained. No new fracture means D = 0.'))
 evidence=load_npz(O/'MEVO_ACTUAL_EVIDENCE.npz');mt=evidence['trajectory_time_s'];mx=evidence['positions_m'];ft=evidence['force_sample_time_s'];p=abs(evidence['force_sample_pressure_Pa']);fig,ax=plt.subplots(1,3,figsize=(17,5));ax[0].fill_between(ft*1e3,p.min(1),p.max(1),color='#417f9e',alpha=.4,label='32-particle range');ax[0].set(xlabel='Time (ms)',ylabel='Acoustic |P| (Pa)',title='Actual F10S pressure samples\nACTUAL_SOLVE');ax[0].grid();ax[1].plot(mt*1e3,(mx.mean(1)-mx[0].mean(0))[:,0]*1e6,color='#a32638');ax[1].set(xlabel='Time (ms)',ylabel='Cluster COM X shift (μm)',title='Actual center translation\nFlow configured STRICT_ZERO');ax[1].grid();ax[2].axis('off');ax[2].text(.04,.92,'NOT_MEASURED',fontsize=19,weight='bold',transform=ax[2].transAxes);ax[2].text(.04,.79,'Bubble → mean flow → clot traction\n\nNo validated transfer law.\n\nAPECSS interpolation exists;\nstreaming shear is not inferred\nfrom particle translation.\n\nOverall: NOT_YET_EVALUABLE',fontsize=14,va='top',linespacing=1.5,transform=ax[2].transAxes);fig.tight_layout();gallery.append(save(fig,'Fig07_MEVO_available_and_missing','Read actual current F10S records. Its homogeneous-water proxy differs from the historical k-Wave geometry. Radial table support does not establish realizable clot shear.'))
 fig,ax=plt.subplots(1,2,figsize=(14,5));ax[0].bar(nums,[r['wall_s'] for r in rows],color='#3b7796',label='Wall');ax[0].plot(nums,[r['cpu_s'] for r in rows],'o-',color='#a32638',label='CPU');ax[0].set(xlabel='Native trial number',ylabel='Measured seconds',title='Native evaluation costs\nACTUAL_SOLVE');ax[0].legend();ax[0].grid(axis='y');ax[1].bar(nums,[r['peak_RSS_bytes']/1024**3 for r in rows],color='#6b7d8c');ax[1].set(xlabel='Native trial number',ylabel='Peak RSS (GiB)',title='Measured process peak memory\nNo new full V3 run');ax[1].grid(axis='y');fig.tight_layout();gallery.append(save(fig,'Fig08_measured_resources','CPU and wall are measured per native subprocess. Different geometry and physical duration preclude a validated speedup claim versus historical V3.'))
 fig,ax=plt.subplots(1,2,figsize=(14,5));vals=[mapping[k] for k in ['force_relative_error','moment_relative_error','power_relative_error']];ax[0].bar(['Force','Moment','Power'],vals,color='#377b85');ax[0].set_yscale('log');ax[0].axhline(1e-10,ls='--',color='#a32638',label='Acceptance threshold');ax[0].set(ylabel='Relative discrepancy',title='Conservation / transpose pairing\nACTUAL_SOLVE: algebraic map');ax[0].legend();ax[1].axis('off');ax[1].text(.02,.94,f"Signed extrapolation weights\n\nMinimum: {mapping['min_weight']:.3f}\nMaximum: {mapping['max_weight']:.3f}\nMax column absolute sum: {mapping['max_column_absolute_sum']:.3f}\n\nConservation passed.\nLocal load accuracy is not certified.\nCurrent-surface FSI is not implemented.",fontsize=15,va='top',linespacing=1.4);fig.tight_layout();gallery.append(save(fig,'Fig09_mapping_conservation_and_scope','Affine weights preserve reference force and moment; the same transpose preserves virtual mechanical power. Negative weights can amplify local nodal forces, so these checks do not certify local traction accuracy.'))
 return gallery
def main():
 start=time.perf_counter();cpu=time.process_time();F.mkdir(parents=True,exist_ok=True)
 cases=[]
 for p in sorted((O/'trials').glob('*/SUMMARY.json')):
  z=read(p)
  if z['status']=='COMPLETED':cases.append(dict(path=p.parent,summary=z))
 assert len(cases)==4,'Only finalize the four actual completed native evaluations'
 rows=[]
 for c in cases:
  z=c['summary'];h=history(c['path']/'history.csv');params=z['params'];rows.append(dict(trial=int(c['path'].name.split('_')[-1]),**params,fluid_Re=params['drive_speed_m_s']*25,actual_peak_shear_Pa=float(h['shear_Pa'].max()),mechanical_end_s=z['time_s'],objective=z['objective'],max_tensile_kappa=float(h['max_tensile_kappa'].max()),new_softened=int(h['new_softened'].max()),new_bulk_failed=int(h['bulk_failed'].max()),new_interface_failed=int(h['interface_failed'].max()),historical_largest_free=z['history_max_free_nodes'],final_largest_free=int(h['largest_detached_nodes'][-1]),N_ge14=int(h['N_ge14'][-1]),N_ge34=int(h['N_ge34'][-1]),N_ge114=int(h['N_ge114'][-1]),detached_mass_kg=float(h['detached_nodes'][-1]*RHO*V),debris_1_2_mass_kg=float(h['debris_1_2_mass_kg'][-1]),W_J=float(h['W_J'][-1]),D_J=float(h['D_J'][-1]),max_relative_R=z['max_relative_energy_residual'],accepted_steps=z['accepted_steps'],wall_s=z['wall_s'],cpu_s=z['cpu_s'],peak_RSS_bytes=z['peak_RSS_bytes'],result_class='REAL_NO_FRACTURE' if not h['new_softened'].max() else 'REAL_DAMAGE',required_acoustic_pressure_Pa='NOT_YET_EVALUABLE',fracture_threshold='NOT_IDENTIFIED'))
 writecsv(O/'LOAD_REQUIREMENT_TABLE.csv',rows);best=max(cases,key=lambda c:c['summary']['objective']);bestrow=max(rows,key=lambda r:r['objective']);mapping=read(O/'FORCE_MAPPING_AUDIT.json');mevo=read(O/'MEVO_ASSET_AUDIT.json');frozen=read(O/'DESIGN_FREEZE.json')
 for path,digest in frozen['reused_source_sha256'].items():assert sha(R/path)==digest,path
 assert sha(LIB)==frozen['library_sha256']
 gate=high_fidelity_gate([c['summary'] for c in cases]);provider=FrozenV3Provider(best['summary']['params']);ids=np.arange(1,3489)[::-1];force=provider.force(.0004,ids);assert not force[ids>1927].any();base=np.ones((1927,3))*1e-9;old=dict(load=base,external_force=base*3);new=provider.replace_load(old,.0004);assert np.max(abs(new['external_force']-new['load']-2*base))<1e-22
 np.savez_compressed(O/'NEXT_V3_LOAD_PORT.npz',original_atom_ids=ids,reference_force_N=force,time_s=.0004);dump(O/'V3_PROVIDER_PORT_CHECK.json',dict(status='OFFLINE_PORT_PASS_NOT_NATIVE_V3_RUN',nonsequential_original_ids=True,anchors_zero=True,preserved_background_difference_N=float(np.max(abs(new['external_force']-new['load']-2*base))),candidate=best['summary']['params']))
 numerical=dict(status='FOUR_COMPLETED_NATIVE_EVALUATIONS',rows=rows,trial_summaries=[dict(path=str(c['path'].relative_to(R)),summary=c['summary']) for c in cases],failed_initialization_trial=0,native_unchanged=True,reused_input_hashes_match=True,material_Gc='UNVALIDATED',interface='UNCALIBRATED',macro_chunk_size_convergence='UNVALIDATED',full_V3=gate)
 dump(O/'NUMERICAL_RESULTS.json',numerical)
 status=dict(task=TASK,run_id=RUN,TRACTION_MODEL_STATUS='PASS_LOCAL_UNSTEADY_STOKES_APPROXIMATION',FORCE_MAPPING_STATUS='PASS_REFERENCE_FORCE_MOMENT_POWER_SIGNED_WEIGHTS',OPTIMIZER_STATUS='REAL_CLOSED_LOOP_INITIAL_DESIGN_BOTORCH_EXPLORATORY',REAL_NATIVE_EVALUATION_STATUS='FOUR_COMPLETED_ONE_RETAINED_INIT_FAILURE',CHUNK_FORMATION_STATUS='NOT_OBSERVED',MEVO_FEASIBILITY_STATUS='NOT_YET_EVALUABLE',FULL_V3_RUN_STATUS='NOT_RUN_EVIDENCE_INSUFFICIENT',SCIENTIFIC_CONCLUSION='Within the bounded low-Re plane-flow surrogate, four real native evaluations show no new damage or detached components. A persisted physical-load optimization loop works, but neither chunk-producing demand nor MEVO streaming realization has been identified.',material_Gc='UNVALIDATED',interface_parameters='UNCALIBRATED',macro_chunk_size_convergence='UNVALIDATED',curved_V3_flow='UNVALIDATED_LOCAL_PATCH_APPROXIMATION',best_proxy_trial=bestrow['trial'],new_full_V3_cases=0)
 dump(O/'FINAL_STATUS.json',status)
 gallery=figures(cases,rows,best,mevo,mapping);dump(O/'FIGURE_INDEX.json',gallery)
 old=read(R/'verification/v3_chunk_fragmentation_load_direction_control_v1/20261016T105642Z/RUN_RESOURCES.json');stages=[read(p) for p in O.glob('OPTIMIZER_STAGE_*.json')];resources=dict(native_wall_sum_s=sum(r['wall_s'] for r in rows),native_cpu_sum_s=sum(r['cpu_s'] for r in rows),native_max_RSS_bytes=max(r['peak_RSS_bytes'] for r in rows),accepted_steps=sum(r['accepted_steps'] for r in rows),optimizer_stages=stages,historical_full_V3_wall_s=old['main_wall_s'],historical_full_V3_CPU_s=old['main_CPU_s'],historical_full_V3_duration_s=.00233,new_full_V3_count=0,validated_speedup=None,saved_counterfactual_wall_s=None,cost_scope='Actual small subprocess cost versus historical different-size/duration V3 cost. No equal-work timing or hypothetical saved seconds claimed.',install_and_early_preflight_cost='NOT_MEASURED; package log and timestamps retained',analysis_wall_s=time.perf_counter()-start,analysis_CPU_s=time.process_time()-cpu,analysis_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
 dump(O/'RUN_RESOURCES.json',resources)
 report(cases,rows,bestrow,mevo,mapping,resources,status)
 cards=''.join(f'<section><h2>{html.escape(g["name"])}</h2><a href="{g["name"]}.png"><img src="{g["name"]}.png"></a><p>{html.escape(g["caption"])}</p><a href="{g["name"]}.svg">SVG</a></section>' for g in gallery)
 doc='<!doctype html><html><meta charset="utf-8"><title>Physics constrained clot inverse design</title><style>body{font:18px system-ui;max-width:1450px;margin:30px auto;padding:0 24px;background:#f5f7f9;color:#202932}section{background:white;padding:22px;margin:24px 0;border:1px solid #ccd3da}img{width:100%}a{color:#24617b}p{line-height:1.6}</style><h1>物理约束血栓载荷反向设计</h1><p>4 次有效原生小试样计算；没有新增损伤或真实游离块体。Stokes 局部平面近似；尚无微泡 → 声流 → 血栓牵引的已验证转换。没有运行新完整 V3。</p><p><a href="../../../verification/'+TASK+'/'+RUN+'/REVIEW_PHYSICS_CONSTRAINED_LOAD_INVERSE_DESIGN_V1.md">完整科研报告</a> · <a href="../../../verification/'+TASK+'/'+RUN+'/LOAD_REQUIREMENT_TABLE.csv">实际试验表</a></p>'+cards+'</html>'
 (F/'OPEN_RESULTS.html').write_text(doc);(F/'index.html').write_text(doc);(O/'OPEN_RESULTS.html').write_text('<meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html"><a href="../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html">打开科研图表</a>')
 print(json.dumps(status,indent=2),flush=True)
def report(cases,rows,best,mevo,mapping,res,status):
 table='\n'.join('| '+ ' | '.join(str(v) for v in [r['trial'],f"{r['drive_speed_m_s']*1e3:.4f}",f"{r['hold_until_s']*1e3:.4f}",f"{r['actual_peak_shear_Pa']:.5f}",f"{r['objective']:.6g}",r['new_softened'],r['historical_largest_free'],f"{r['wall_s']:.2f}"])+' |' for r in rows)
 text=f'''# PHYSICS_CONSTRAINED_CLOT_LOAD_INVERSE_DESIGN_V1

本轮实现并运行了局部流体牵引 → 原始节点映射 → Optuna 提交 → 原生 LAMMPS 求解 → 实际拓扑/能量读回 → BoTorch 下一候选的最小闭环。完成四个有效小试样；保留一次 run 0 初始化接口拒绝。**没有新软化、新断键或真实脱离分量；没有运行新完整 V3。** 这不是已找到块状碎裂载荷的结果。

## 实际观测

| Trial | 驱动速度 mm/s | 开始关断 ms | 实际峰值剪切 Pa | 冻结代理目标 | 新软化 | 历史最大自由分量 | 原生墙钟 s |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

代理目标最高的 trial {best['trial']} 仍只有弹性伸长前兆。所有 ≥14、≥34、≥114 节点分量、脱离质量、1–2 节点碎屑质量和不可逆损伤耗能均为零；不存在可报告的自由块体内部活键比例或二次裂解。不能将这些零值转化为声流普遍不能碎裂血栓的结论。

## 十五项科学问题

1. **是否有流体牵引模板？** 有：直接解一维非稳态 Stokes 扩散方程的 Fourier 模态 ODE，形成可查询时间场和真实 V3 面片 NPZ。下壁静止、上壁平移，侧向无限，初始流体静止，零压力梯度。物理来源为明确的流体边界值问题；微泡驱动尚未求解。
2. **与 Gaussian 的差别？** 没有人工指定应力幅度/空间 Gaussian 宽度或任意混合方向；参数为上壁速度和关断时刻，剪切由 μ∂u/∂y 得到，启动由黏性扩散决定。上壁运动是理想化平均流驱动，不能称作已实现的微泡 Streaming。V3 上投影为独立切平面，尚未满足全曲面流动连续性。
3. **压力、剪切、速度和时间依据？** μ=0.001 Pa·s、ρ=1000 kg/m³，H=25 μm，速度范围0.5–4 mm/s，Re_H≤0.1。稳态剪切0.02–0.16 Pa；p增量严格0 Pa。驱动0时刻开始，0.1 ms半余弦上升，0.4–0.8 ms开始0.1 ms关断，随后继续计算0.3 ms黏性尾流。H²/ν=0.625 ms，最慢模态时间0.063326 ms，5倍为0.316629 ms。不是Cho实验建立时程，也没有第一步施加稳态峰值。独立有限差分与模态解的峰值归一化应力差{read(O/'FLUID_TRACTION_CHECK.json')['finite_difference_vs_series_max_scaled_difference']:.6g}。
4. **真实V3映射？** 使用含栓CFD的3611个真实暴露面片、原始网格三角形节点ID及1927原始材料ID。它与上一轮4918个几何法向面片的离散来源不同，不进行身份混同。复用历史稀疏支持；共面/病态支持仅补入邻近原始材料点。原V3库、材料和历史数据未改变。
5. **力、矩、功率？** 每面片 Σw=1、ΣwX=c，力为W(A t)，面片速度为Wᵀv。相对误差依次为{mapping['force_relative_error']:.3g}、{mapping['moment_relative_error']:.3g}、{mapping['power_relative_error']:.3g}。同参考点、参考面积/法向；标记FROZEN_SURFACE_ONE_WAY。权重约{mapping['negative_weight_fraction']*100:.2f}%为负，范围[{mapping['min_weight']:.3f},{mapping['max_weight']:.3f}]，列绝对和最大{mapping['max_column_absolute_sum']:.3f}。守恒不等于局部载荷精度认证，变形后真实空间力矩/真实壁面功率尚未保证。
6. **Optuna提交与恢复？** 实际SQLite ask/tell，Python独立子进程运行原生LAMMPS。完成两条有效记录后进程退出，再加载DB及pickle采样器状态继续；每个阶段从磁盘重新加载并核对状态和值。失败初始化未作为零分样本。支持中断后重新挂接RUNNING试验；已有SUMMARY直接读回。未完成的单个原生轨迹重启从原参考态执行，不宣称本轮实现了半途机械checkpoint续算。
7. **BoTorch样本够吗？** 不够。三条有效随机初始观测后实际拟合SingleTaskGP，运行带Re及能量约束的qLogExpectedImprovement，提出并完成第四条真实原生计算；四条观测后另存NEXT_CANDIDATE.json但不运行。只验证可运行反馈接口，不声称GP预测可靠、已全局优化或高质量贝叶斯优化。
8. **真实小试样？** 四组，原始1792材料点+336固定锚点、原始预裂纹与endpoint-star、G/K/s0/sc及库完全复用。加载改为平顶面Stokes剪切力，解除旧位移夹持，材料点全部积分；底部原三层锚点固定。这是局部有限试样机制筛选，不能预测V3碎块尺寸。没有为结果新增预切或调整材料。
9. **裂纹集中/凝聚候选？** 未观察到新增损伤，因此未得到凝聚裂纹证据。冻结上半试样作为诊断区域，剩余z=0韧带为候选边界；代理为边界平均损伤减内部平均损伤，加0.01倍边界平均clip(κ/s0,0,1)，取全历史最大值，预裂键排除。只剩伸长前兆时明确如此标注。没有事后改变目标。
10. **独立多节点块体？** 没有。每个新断键事件都会基于原始活键做锚点可达性分析；本轮无新断键，末态原生gather_bonds与参考存活端点逐一核对通过。原有预裂键不计作新增断裂。
11. **有效载荷需求合理吗？** 已试速度/剪切在规定局部低Re范围内，但没有识别出能产生块体的有效需求或断裂阈值。不无限线性放大，不反推一个未经运行的临界剪切。表中是实际试验需求，不是已达成碎裂的需求。
12. **MEVO能直接实现吗？** 尚不能评价。实际读取历史F1 k-Wave复基底及当前F10S的复压力、APECSS有效叶单元、32粒子位置/有限差分速度和一/二次Bjerknes力。当前F10S为HOMOGENEOUS_WATER_CONTROL_PROXY、STRICT_ZERO背景流，覆盖3 ms。当前径向数据库外包范围R=1.779–1.781 μm、f=2.499–2.501 MHz、P=0.1–40000 Pa，还须落入有效叶。缺少机械剪切需求→声压的逆映射，不能凭粒子平移把牵引标为已实现。历史k-Wave与当前F10S装置不同。接口分别返回SUPPORTED_BY_AVAILABLE_MODEL（仅有效径向插值）、CONDITIONALLY_PLAUSIBLE（局部流体/时间范围）、OUTSIDE_VALIDATED_RANGE、NOT_YET_EVALUABLE；总体为最后一项。
13. **主要不确定性？** 微泡/群集至近壁平均速度梯度的转换、曲面流动及边缘误差、局部有符号节点映射的载荷精度、材料Gc/界面参数、尺寸收敛。当前全体弱弹性结果无法区分更强但仍物理允许的局部流动是否有效。
14. **下一轮具体机制？** 优先以现有APECSS谐波和群集状态为输入，建立单一局部近壁Stokes平均流/受迫流动问题，包含实际壁面与血栓几何，再生成空间非均匀剪切。先验证源项、时间建立及表面牵引，不将中心速度代替流速；若需要惯性射流，应另立Navier–Stokes模型，不扩大本轮低Re标记范围。
15. **成本节约多少？** 本轮四个有效小试样实测墙钟合计{res['native_wall_sum_s']:.3f} s、CPU合计{res['native_cpu_sum_s']:.3f} s、最大RSS {res['native_max_RSS_bytes']/1024**3:.3f} GiB，共{res['accepted_steps']}步。上一轮单个2.33 ms完整V3实测{res['historical_full_V3_wall_s']:.3f} s；几何、物理时长和检查开销不同，不报告等工作量加速倍数，也不把未执行的V3成本当作实测节约。本轮新增完整V3次数0。

## 状态与复现

```json
{json.dumps(status,ensure_ascii=False,indent=2)}
```

见 `REPRODUCE.md`、`OPTIMIZATION_STUDY.db`、`OPTIMIZATION_TRIALS.csv`、`NUMERICAL_RESULTS.json`、`TRACTION_MODEL_SCOPE.md`、`MEVO_FEASIBILITY.json` 和图表索引。设计冻结后只有初始化回调名称适配修正，记录于PORT_ADAPTER_FIX.json；物理域、目标、材料和原生库不变。完整V3载荷提供器已离线验证原始ID、锚点零力及背景力保留，但没有宣称做过本轮V3动力学验收。

## 文献与适用边界

Oh 2021本地全文已读：第4.1节以指定点对伸长作为裂纹路径代理，第4.2节以有限差分核对伴随灵敏度。本轮借鉴预定义区域/伸长前兆，未移植二维形状伴随到三维不可逆原生系统。[论文](https://doi.org/10.1016/j.engfracmech.2021.107837)

Cho 2015仅摘要/出版记录核实：PIV测得群集诱导射流与对称涡结构；这些测量不能使本项目粒子速度自动成为流体速度。全文建立时程和全部实验边界未核实，不采用其中峰值作首步加载。[摘要](https://pubmed.ncbi.nlm.nih.gov/26025507/)

Zhang 2022仅出版社摘要核实：FEM/PIV比较入射角对流速、表面剪切及溶栓的影响。未将含药物/微泡因素的溶栓率映射成本材料断裂阈值。[出版社](https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.15874)

Karmakar 2025本地全文已读：第1节及界面建模讨论本体凝聚破裂与界面脱粘竞争。保持既有LPS变分软化，不移植其NOSB/Ogden及界面标定参数。[论文](https://doi.org/10.1002/cnm.70118)

实现参考官方[Optuna持久化](https://optuna.readthedocs.io/en/stable/tutorial/20_recipes/001_rdb.html)与[BoTorchSampler](https://optuna-integration.readthedocs.io/en/stable/reference/generated/optuna_integration.BoTorchSampler.html)，实际使用Optuna 5.0.0的set_constraint接口。依赖版本以本轮安装清单为准。
'''
 (O/'REVIEW_PHYSICS_CONSTRAINED_LOAD_INVERSE_DESIGN_V1.md').write_text(text)
if __name__=='__main__':main()
