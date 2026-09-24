#!/usr/bin/env python3
"""Render 12 figures from immutable saved CSV/JSON evidence, never rerun physics."""
from pathlib import Path
import sys,json,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent;sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
REPORT=PACKAGE/'reports/particle6_5';DATA=REPORT/'data'
NAMES=['00_particle6_5_scope_and_literature','01_nearfield_activation_weight','02_h_lower_vs_particle_radius',
'03_sphere_wall_regularized_resistance','04_sphere_sphere_regularized_resistance','05_sphere_wall_continuum_handoff',
'06_sphere_sphere_continuum_handoff','07_real_fem_subnanometer_replay','08_handoff_scale_sensitivity',
'09_handoff_physical_time_refinement','10_particle6_5_lammps_bridge_parity','11_particle6_5_timestep_comparison']
SOURCES=[['00_scope'],['01_activation'],['02_lower'],['03_wall_static'],['04_pair_static'],
['05_wall_dt1','05_wall_old_p5'],['06_pair_dt1','06_pair_old_p5'],['07_historical_p5','07_real_matched_old_p5','07_real_v1_dt1'],
['08_sensitivity','08_real_floor1.5_dt1','07_real_v1_dt1','08_real_floor3_dt1'],['09_large_step_refinement'],['10_bridge'],['11_timestep']]
NOTES=[
('区分文献先例、项目选择和尚未冻结的模型。','0.05a 与 10⁻³a 有悬浮液先例；2 nm、0.01 与平滑函数按批准的项目合同使用。','文献没有直接测得 SonoVue–小鼠血管 2 nm 接触距离；图中尺度顺序以 1 µm 球为例。'),
('看权重是否有跳变，两个端点斜率是否归零。','权重从 1 平滑降到 0，导数在两个端点为 0，区间内不回升。','无数值异常；smoothstep 是数值耦合选择。'),
('看分子尺度与粒径尺度在哪相交。','2 µm 半径处相交；历史约 0.588 与 1.266 µm 两个例子的下限都为 2 nm。','没有将下限硬编码为总是 2 nm，也没有改采样半径。'),
('看原始 gap、权重、阻力和未施加约束时的法向速度。','V1 外区贡献为零，完整区与 P5 leading term 一致，分母在下限封顶。','下限以下仅为静态诊断扫描；虚线 P5 leading 是公式参照，点线 P5 default 保留原启用范围。'),
('看等半径、1:2 和 1:4 配对的权重、阻力与相对法向速度。','参考长度均取较小半径，阻力仍使用原 R_eff；交换配对顺序不改变结果。','图示为未约束的静态响应；下限以下不属于可接受的动态状态。'),
('看靠墙球到交棒时刻后的间隙、速度和切向位移。','V1 保留约 2 nm 原始间隙并约束法向速度；切向位移继续增长，P5 对照继续进入更小间隙。','无位置回推；P5 对照沿用旧 0.01 硬切换，初期差异也包含旧启用边界效应。交棒时刻需结合图 11。'),
('看两球是否对称靠近，交棒后是否继续处理切向运动。','两球中心关于中线对称，间隙停在下限，公共切向运动继续。','图中旧 P5 初始 gap 超过它的默认配对启用范围；该历史默认行为保留。'),
('区分旧双 MB 历史、同初态 P5 重放和 V1 重放。','旧 0.02354 nm 历史保留；同一 3.233 nm 保存状态重放中，P5 最小约 0.0815 nm，V1 最小约 2 nm。','旧双 MB 起点已低于下限而被拒绝；本图的有效起点未移动，交棒次数对步长敏感。'),
('分别看 1.5、2、3 nm 原始轨迹及差异，留意未观察到的事件。','三个尺度独立显示，记录最小 gap、交棒时间、位移与终态差异；未发生的交棒不填零。','事件及终态速度有步长和尺度依赖；没有用人为百分比宣布科学敏感性通过。'),
('看跨越下限的大试探是否细分，剩余时间是否处理。','越界试探被拒绝，接受子区间完整覆盖请求窗口，原始 gap 保留，切向运动继续。','时间细分保证约束安全，不能替代整体时间精度收敛验证。'),
('看多余 LAMMPS 候选如何筛选，以及 R、b、U、约束和位置差。','候选包含额外配对，精确近场筛选后与独立路径相同；所有矩阵及状态差为零，实际二进制重启后继续一致。','只验证单 MPI rank；邻居 cutoff 与 skin 均为验证设置。'),
('看 dt、dt/2、dt/4 的事件、最小间隙、终态位置和速度。','各组完整覆盖时间并保持下限；真实场景交棒次数由 6 变为 0、0，合成事件时间也随步长变化。','尚未证明时间精度或事件拓扑收敛；本图不选择生产步长。')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'figure.dpi':100,'savefig.dpi':140})
C=['#176a96','#ca6a23','#29834e'];
def read(name):return json.loads((DATA/(name+'.json')).read_text())
def xy(rows,key,scale=1):return np.array([r[key] for r in rows])*scale
def time(rows):return xy(rows,'time_s',1e3)
def style(ax,xlabel,ylabel,title):ax.set(xlabel=xlabel,ylabel=ylabel,title=title);ax.grid(alpha=.2)
def handoff_line(ax,d):
 t=d['summary']['handoff_time_s']
 if t is not None:ax.axvline(t*1e3,color='#ad3333',ls=':',label=f'first handoff {t*1e3:.4g} ms')
def trajectory_panels(case):
 index='05' if case=='wall' else '06';d=read(f'{index}_{case}_dt1');old=read(f'{index}_{case}_old_p5');r=d['interactions']
 fig,axs=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
 for data,label,color,ls in [(old,'P5 default',C[1],'--'),(d,'V1',C[0],'-')]:axs[0,0].plot(time(data['interactions']),xy(data['interactions'],'h_geom_m',1e9),ls,color=color,label=label)
 axs[0,0].axhline(2,color='gray',ls=':',label='h_lower = 2 nm');axs[0,0].set_yscale('log');style(axs[0,0],'time (ms)','raw h_geom (nm)','Original geometry, no clipping');axs[0,0].legend(fontsize=8)
 axs[0,1].plot(time(r),xy(r,'g_nf_m',1e9));axs[0,1].axhline(0,color='gray',ls=':');style(axs[0,1],'time (ms)','g_NF (nm)','Distance remaining to handoff')
 for key,label,ls in [('free_vn_m_s','free',':'),('lubricated_vn_m_s','lubricated','--'),('constrained_vn_m_s','constrained','-')]:axs[0,2].plot(time(r[1:]),xy(r[1:],key,1e6),ls,label=label)
 axs[0,2].set_yscale('symlog',linthresh=.01);style(axs[0,2],'time (ms)','normal speed (um/s)','Interval-held speeds');axs[0,2].legend(fontsize=8)
 p=d['positions'];ids=sorted({r['particle_id'] for r in p})
 for i in ids:
  a=[q for q in p if q['particle_id']==i];axis=2 if case=='wall' else 0
  axs[1,0].plot(time(a),[(q['center_m'][axis]-a[0]['center_m'][axis])*1e9 for q in a],label=f'ID {i}')
 style(axs[1,0],'time (ms)','normal center displacement (nm)','Exact center motion from initial position');axs[1,0].legend(fontsize=8)
 q=[q for q in p if q['particle_id']==ids[0]];axis=0 if case=='wall' else 1
 axs[1,1].plot(time(q),[(v['center_m'][axis]-q[0]['center_m'][axis])*1e9 for v in q]);style(axs[1,1],'time (ms)','tangential displacement (nm)','Motion continues after handoff')
 axs[1,2].step(time(r),[int(q['interaction_state']=='CONTINUUM_HANDOFF_CONTACT') for q in r],where='post');axs[1,2].set_yticks([0,1],['continuum','handoff']);style(axs[1,2],'time (ms)','accepted state','Continuum boundary, not molecular contact')
 for ax in axs.flat:handoff_line(ax,d)
 fig.suptitle(f"{case}: first handoff = {d['summary']['handoff_time_s']*1e3:.6f} ms | raw geometry, no position projection",fontsize=15)
 return fig

def make(index):
 if index==0:
  fig,axs=plt.subplots(1,3,figsize=(15,8),layout='constrained')
  text=[('LITERATURE PRECEDENT','A: suspension simulation\nouter = 0.05 a\nlower order = 0.001 a\n\nB / C: confined water\nnanoscale boundary evidence\n\nD / E: lipid / shell context\n4 nm shell is NOT a radius offset'),('PROJECT MODEL CHOICE','Example: sphere radius 1 um\n\n0.05 a: near-field off\n        | C1 transition |\n0.01 a: full leading term\n        | continuum |\nh_lower = max(2 nm, 0.001 a)\n        | handoff constraint |\nCapped lubrication remains'),('NOT FROZEN','Glycocalyx wall model\nRBC lubrication\nProduction timestep\nProduction neighbor cutoff / skin\n\nNo molecular solid-contact model\nNo full production hydrodynamics\nNo Particle-7')]
  for ax,(title,body),c in zip(axs,text,C):
   ax.axis('off');ax.text(.02,.94,title,fontsize=17,weight='bold',color=c,va='top');ax.text(.02,.83,body,fontsize=13,va='top',linespacing=1.8)
  fig.suptitle('Sphere normal near-field V1: evidence roles and model boundary',fontsize=18)
 elif index==1:
  d=read('01_activation');fig,axs=plt.subplots(2,1,figsize=(11,8),layout='constrained')
  for ax,key,title in [(axs[0],'w','C1 activation weight'),(axs[1],'derivative','Derivative vanishes at both endpoints')]:
   ax.plot(xy(d,'chi'),xy(d,key));ax.axvline(.01,color=C[1],ls=':',label='FULL = 0.01');ax.axvline(.05,color=C[2],ls=':',label='OFF = 0.05');style(ax,'chi = h_geom / a_ref',key,title);ax.legend()
 elif index==2:
  d=read('02_lower');fig,ax=plt.subplots(figsize=(11,7),layout='constrained')
  for key,label,ls in [('h_molecular_component_m','2 nm molecular floor',':'),('h_suspension_component_m','0.001 a_ref','--'),('h_lower_m','h_lower = maximum','-')]:ax.plot(xy(d,'a_ref_m',1e6),xy(d,key,1e9),ls,label=label)
  for a,label in [(.588015722765549,'0.588 um'),(1.2658269815352818,'1.266 um')]:ax.plot(a,2,'o');ax.annotate(label,(a,2),xytext=(a,3.5 if a<1 else 5.5),arrowprops={'arrowstyle':'->'})
  ax.axvline(2,color='gray',ls=':',label='crossover: 2 um');ax.set_xscale('log');ax.set_yscale('log');style(ax,'reference radius (um)','scale (nm)','Current historical SonoVue radii: molecular component dominates');ax.legend(loc='upper left')
 elif index in [3,4]:
  rows=read('03_wall_static' if index==3 else '04_pair_static');fig,axs=plt.subplots(1,3,figsize=(16,6.5),layout='constrained')
  for ratio in ([None] if index==3 else [1,2,4]):
   d=[r for r in rows if r['radius_ratio']==ratio];label='sphere-wall' if ratio is None else f'radius ratio 1:{ratio}'
   axs[0].plot(xy(d,'chi'),xy(d,'w'),label=label)
   axs[1].loglog(xy(d,'h_geom_m',1e9),xy(d,'coefficient_kg_s'),label=label+' V1')
   axs[2].semilogx(xy(d,'h_geom_m',1e9),xy(d,'v1_unconstrained_vn_ratio'),label=label+' V1')
   if index==3:
    axs[1].loglog(xy(d,'h_geom_m',1e9),xy(d,'old_p5_leading_uncapped_zeta_kg_s'),'--',label='P5 leading formula (uncapped)')
    axs[1].loglog(xy(d,'h_geom_m',1e9),xy(d,'old_p5_default_zeta_kg_s'),':',label='P5 default eligibility')
    axs[2].semilogx(xy(d,'h_geom_m',1e9),xy(d,'old_p5_default_vn_ratio'),':',label='P5 default')
  for ax in axs[1:]:ax.axvline(2,color='gray',ls=':',label='h_lower');ax.axvspan(.01,2,color='gray',alpha=.08)
  axs[0].axvline(.01,color='gray',ls=':');axs[0].axvline(.05,color='gray',ls=':')
  style(axs[0],'chi = raw h / a_ref','w','Same smooth activation');style(axs[1],'raw h_geom (nm)','normal resistance (kg/s)','Capped coefficient; raw geometry retained');style(axs[2],'raw h_geom (nm)','V_n / U_free,n','Unconstrained normal response')
  for ax in axs:ax.legend(fontsize=8)
  fig.suptitle('Below 2 nm shaded: static diagnostic ONLY; not accepted V1 dynamics',fontsize=14)
 elif index in [5,6]:fig=trajectory_panels('wall' if index==5 else 'pair')
 elif index==7:
  history=read('07_historical_p5');old=read('07_real_matched_old_p5');new=read('07_real_v1_dt1');r=new['interactions'];fig,axs=plt.subplots(2,2,figsize=(15,9),layout='constrained')
  for i in [101,203]:
   d=[q for q in history['rows'] if q['particle_id']==i];axs[0,0].semilogy(time(d),xy(d,'h_geom_m',1e9),label=f'historical MB {i}')
  axs[0,0].axhline(2,color='gray',ls=':',label='2 nm');style(axs[0,0],'historical elapsed time (ms)','raw gap (nm)',f"Original P5 two-MB history preserved | min {history['minimum_h_geom_m']*1e9:.5f} nm");axs[0,0].legend(fontsize=8)
  for d,label,ls in [(old,'matched P5','--'),(new,'V1','-')]:axs[0,1].semilogy(time(d['interactions']),xy(d['interactions'],'h_geom_m',1e9),ls,label=label)
  axs[0,1].axhline(2,color='gray',ls=':',label='2 nm = h_lower');handoff_line(axs[0,1],new);style(axs[0,1],'replay elapsed time (ms)','raw gap (nm)','Same saved 3.233 nm state, same FEM, dt and window');axs[0,1].legend(fontsize=8)
  for key,label,ls in [('free_vn_m_s','free',':'),('lubricated_vn_m_s','lubricated','--'),('constrained_vn_m_s','constrained','-')]:axs[1,0].plot(time(r[1:]),xy(r[1:],key,1e6),ls,label=label)
  axs[1,0].set_yscale('symlog',linthresh=.01);style(axs[1,0],'accepted interval end (ms)','normal speed (um/s)','Interval-held speeds; exact endpoint normal');axs[1,0].legend()
  axs[1,1].plot(time(r),xy(r,'h_geom_m',1e9));axs[1,1].axhline(2,color='gray',ls=':',label='handoff');axs[1,1].set_ylim(1.8,3.4);style(axs[1,1],'replay elapsed time (ms)','raw gap (nm)','Close view: raw values, no smoothing');handoff_line(axs[1,1],new);axs[1,1].legend(fontsize=8)
 elif index==8:
  sensitivity=read('08_sensitivity');fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
  for col,(nm,name,c) in enumerate(zip([1.5,2,3],['08_real_floor1.5_dt1','07_real_v1_dt1','08_real_floor3_dt1'],C)):
   d=read(name);r=d['interactions'];ax=axs[0,col];ax.semilogy(time(r),xy(r,'h_geom_m',1e9),color=c);ax.axhline(nm,color=c,ls=':');style(ax,'time (ms)','raw gap (nm)',f'{nm:g} nm floor | SENSITIVITY_ONLY');ax.set_ylim(1,10);handoff_line(ax,d)
  for nm,c in zip([1.5,2,3],C):
   d=sorted([r for r in sensitivity['rows'] if r['h_molecular_floor_m']=={1.5:1.5e-9,2:2e-9,3:3e-9}[nm]],key=lambda r:r['dt_divisor']);x=[r['dt_divisor'] for r in d]
   axs[1,0].plot(x,xy(d,'minimum_h_geom_m',1e9),'o-',color=c,label=f'{nm:g} nm')
   for r in d:
    if r['handoff_time_s'] is not None:axs[1,1].scatter(r['dt_divisor'],r['handoff_time_s']*1e3,color=c)
    else:axs[1,1].text(r['dt_divisor'],12.7-.4*nm,f'{nm:g}: none',color=c,fontsize=8,ha='center')
   axs[1,2].plot(x,xy(d,'final_position_m_absolute_difference_from_2nm',1e9),'o-',color=c,label=f'{nm:g} nm')
  style(axs[1,0],'dt divisor','minimum raw gap (nm)','All three validation timesteps');axs[1,0].legend()
  style(axs[1,1],'dt divisor','first handoff (ms)','None = not observed in this window')
  style(axs[1,2],'dt divisor','final position difference (nm)','Absolute difference from same-dt 2 nm');axs[1,2].legend()
  for ax in axs[1]:ax.set_xticks([1,2,4],['dt','dt/2','dt/4'])
  fig.suptitle('Scientific sensitivity: report differences; user review pending',fontsize=16)
 elif index==9:
  d=read('09_large_step_refinement');r=d['interactions'];l=d['ledger'];fig,axs=plt.subplots(2,2,figsize=(14,8),layout='constrained')
  for a in l:axs[0,0].plot([a['t0_s']*1e3,a['t1_s']*1e3],[a['depth']]*2,linewidth=2,color=C[0])
  axs[0,0].plot([0,32],[-2,-2],color=C[1],linewidth=3,label='requested 32 ms');style(axs[0,0],'physical time (ms)','subdivision depth','Accepted child intervals');axs[0,0].legend()
  axs[0,1].semilogy([a['t1_s']*1e3 for a in l],[a['dt_s'] for a in l],'o-',ms=3);style(axs[0,1],'physical time (ms)','accepted interval (s)','Actual physical-time refinement')
  axs[1,0].plot(time(r),xy(r,'h_geom_m',1e9),'o-',ms=2);axs[1,0].axhline(2,color='gray',ls=':');style(axs[1,0],'physical time (ms)','raw gap (nm)','Accepted geometry')
  axs[1,1].plot([a['t1_s']*1e3 for a in l],np.cumsum([a['dt_s'] for a in l])-np.array([a['t1_s'] for a in l]));style(axs[1,1],'physical time (ms)','coverage error (s)','No time skipped; no endpoint projection')
  for ax in axs.flat:handoff_line(ax,d)
  fig.suptitle(f"First handoff = {d['summary']['handoff_time_s']*1e3:.6f} ms | validation dt = 32 ms",fontsize=15)
 elif index==10:
  d=read('10_bridge');fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained');r=d['matrix_rows'][0];ids=[17,203,901];ix={v:i for i,v in enumerate(ids)}
  for ax,key,title in zip(axs[0],['standalone_eligible_pairs','lammps_candidates','filtered_nearfield_pairs'],['Standalone eligible','LAMMPS broadphase','Filtered near-field']):
   m=np.zeros((3,3))
   for i,j in r[key]:m[ix[i],ix[j]]=m[ix[j],ix[i]]=1
   ax.imshow(m,vmin=0,vmax=1,cmap='Blues');ax.set_xticks(range(3),ids);ax.set_yticks(range(3),ids);ax.set_title(title)
  for ax,keys,title in [(axs[1,0],['R','b'],'Resistance matrix / rhs'),(axs[1,1],['U','J'],'Velocity / constraints')]:
   for k in keys:ax.plot([q['step'] for q in d['matrix_rows']],[q[k] for q in d['matrix_rows']],'o-',label=k)
   ax.set_ylim(-1,1);style(ax,'step','absolute maximum difference',title);ax.legend();ax.text(.5,.8,'exact zero',ha='center',transform=ax.transAxes)
  axs[1,2].plot([r['step'] for r in d['errors']],[r['position_m'] for r in d['errors']],'o-');axs[1,2].axvline(3,color=C[1],ls=':',label='binary restart');axs[1,2].set_ylim(-1,1);style(axs[1,2],'step','position difference (m)','State parity through actual restart');axs[1,2].legend()
 elif index==11:
  rows=read('11_timestep');fig,axs=plt.subplots(4,3,figsize=(15,12),layout='constrained')
  for col,case in enumerate(['wall','pair','real']):
   d=[r for r in rows if r['case']==case];x=[1,2,4];base=d[-1]
   for k,r in enumerate(d):
    if r['handoff_time_s'] is not None:axs[0,col].scatter(x[k],r['handoff_time_s']*1e3,color=C[col])
    else:axs[0,col].text(x[k],.2,'none',ha='center',transform=axs[0,col].get_xaxis_transform())
   axs[1,col].plot(x,xy(d,'minimum_g_nf_m',1e9),'o-',color=C[col]);axs[1,col].axhline(0,color='gray',ls=':')
   for row,key,scale in [(2,'final_position_m',1e9),(3,'final_velocity_m_s',1e6)]:axs[row,col].plot(x,[np.linalg.norm(np.asarray(r[key])-base[key])*scale for r in d],'o-',color=C[col])
   for row,(title,y) in enumerate([('first handoff','time (ms)'),('minimum g_NF','gap minus lower (nm)'),('final x vs dt/4','difference norm (nm)'),('final V vs dt/4','difference norm (um/s)')]):
    style(axs[row,col],'validation timestep',y,case+' | '+title);axs[row,col].set_xticks(x,['dt','dt/2','dt/4']);axs[row,col].set_xlim(.7,4.3)
  maxerr=max(r['time_coverage_error_s'] for r in rows)
  fig.suptitle(f'NOT PRODUCTION TIMESTEP SELECTION\nMaximum coverage error = {maxerr:.2g} s | event convergence NOT ESTABLISHED',fontsize=16)
 return fig


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path);parser.add_argument('--figures-only',action='store_true');args=parser.parse_args()
 target=args.output_dir or REPORT/'figures';target.mkdir(parents=True,exist_ok=True);manifest=[]
 for i,name in enumerate(NAMES):
  fig=make(i);file=target/(name+'.png');fig.savefig(file,metadata={'Software':'Particle-6.5 saved-evidence renderer'});plt.close(fig)
  inputs={f'data/{s}.json':sha256(DATA/(s+'.json')) for s in SOURCES[i]}
  manifest.append(dict(index=i,filename=file.name,sha256=sha256(file),data_sha256=inputs,review_notes=dict(zip(['应该看什么','实际看到什么','有没有异常'],NOTES[i]))))
 if not args.figures_only:
  write_json(REPORT/'FIGURE_MANIFEST.json',dict(renderer=str(Path(__file__).relative_to(REPO)),figures=manifest))
  print('12 figures and hash manifest written')
if __name__=='__main__':main()
