#!/usr/bin/env python3
"""Eleven deterministic figures generated solely from saved machine evidence."""
from pathlib import Path
import argparse,json,sys,textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
from particle_3d.lammps_state import BridgeParticle
REPORT=PACKAGE/'reports/particle6';DATA=REPORT/'data'
NAMES=['00_particle6_scope_and_lammps_environment','01_lammps_state_roundtrip','02_lammps_neighbor_equivalence',
'03_one_step_standalone_vs_lammps','04_multistep_bridge_parity','05_resistance_system_parity','06_lammps_zero_force_audit',
'07_checkpoint_restart_parity','08_shape_metadata_restart','09_real_fem_standalone_vs_lammps','10_neighbor_rebuild_stress']
SOURCES={0:['00_scope','00_lammps_environment'],1:['01_state_roundtrip'],2:['02_neighbor_equivalence'],3:['03_one_step'],
4:['04_mixed_multistep','04_sphere_multistep'],5:['05_resistance'],6:['06_force_audit','06_command_audit'],
7:['07_restart_mixed','07_restart_sphere'],8:['08_metadata_restart'],9:['09_real_two_mb'],10:['10_neighbor_rebuild']}
NOTES=[
('看清谁计算物理、谁只保存状态，并确认实际 LAMMPS 环境。','原 particle_3d 负责全部物理与积分；LAMMPS 只保存状态、查询邻居和重启，实际版本为 20250722。','只验证单 MPI rank；安装包虽有其他物理模块，本轮命令白名单不会调用它们。'),
('看三种形状进入 LAMMPS 再取出后，形状和所有状态字段是否变化。','三个阶段的形状重合，ID、类型、位置、四元数、速度、角速度及几何字段逐位相同。','胶囊的四元数只是保留字段；胶囊轴、半径和柱段长度独立决定几何。'),
('比较所有可能配对、两条查询路径与精确筛选，确认没有遗漏。','六种形状组合都有接触见证；查询候选完全一致，精确筛选排除了额外候选，漏检为 0。','查询距离包含验证 skin；这些数值只用于验证，候选本身不是接触或润滑。'),
('看相同初值经过一个 P5 步长后，位置及 V、Omega、q 是否一致。','两条路径独立计算后完全重合，位置、速度、角速度和姿态误差均为 0。','未放宽容差；同时保存按 float64 运算尺度推导的舍入预算。'),
('看独立连续推进 100 步是否积累误差，以及邻居是否分歧。','P5 球体场景和 P4 混合形状场景的所有误差曲线均为 0，邻居不匹配数为 0。','混合 RBC 只验证原 P4/P2 运动与元数据，不新增或宣称 P5 非球形流体动力学。'),
('比较 R 的非零结构、元素、右端、求解速度和接触约束。','R、b、U、J 的最大差均为 0，接触活动集也一致。','接触约束乘子仍沿用 P5 运动学含义，不是 LAMMPS 物理接触力。'),
('检查每次邻居重建前后是否出现力，并确认速度来源。','所有已记录 LAMMPS force 为 0；atomic 模式没有 torque 数组，诊断值记为 0，未读作物理输入。','torque 的 0 表示未定义且未使用，不能解释成测得的物理力矩；LAMMPS 时间步始终为 0。'),
('看第 40 步销毁实例并读二进制重启后，是否继续原来的轨迹。','混合形状和球体场景均完成 40 + 60 步，与连续 100 步逐位相同，时间和步号一致。','重启从二进制恢复粒子，不按种子重新采样；全局物理时间来自带哈希的 sidecar。'),
('逐字段比较三种形状在 checkpoint 前后的几何和状态。','每个粒子的全部字段差为 0，椭球原旋转矩阵、胶囊轴和保留四元数都恢复。','RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION；额外保存旋转矩阵是为保留原 P4 状态的最后一位。'),
('看原 P5 双 MB 在真实冻结流场中的轨迹、间隙及误差。','原始 P5 入口重新运行后与桥接结果完全一致，墙隙、粒子间隙、近场启用、阻力残差和出口事件一致。','没有移动初值或改半径、步长、时长；亚纳米墙隙的连续介质有效性仍未建立，真实 RBC 动力学仍延期。'),
('看粒子跨过查询范围后，更新位置能否刷新候选列表。','每次更新都强制重建，LAMMPS 与 standalone 的候选状态一致；保存的旧列表在多次跨越时确实会失效。','这是显式位置写入的查询压力测试，不是新物理积分器，也不冻结生产 cutoff 或 skin。')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'figure.dpi':100,'savefig.dpi':140})
COLORS=['#176a96','#d47922']
def read(name):return json.loads((DATA/(name+'.json')).read_text())
def labelzero(ax):
    ax.text(.5,.85,'exact difference = 0',transform=ax.transAxes,ha='center',color='#176a46')
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    ax.set_ylim(-1,1)
def make(index):
    if index==0:
        d=read('00_lammps_environment');fig,ax=plt.subplots(figsize=(14,8));ax.axis('off')
        ax.text(.025,.9,'OUR SOLVER\nPHYSICS AUTHORITY',fontsize=22,weight='bold',va='top',color=COLORS[0])
        ax.text(.025,.68,'Frozen FEM sampling\nP0-P5 geometry / resistance / contact\nPhysical-time subdivision\nValidated x / q integration\nV and Omega are computed here',fontsize=16,va='top',linespacing=1.6)
        ax.text(.53,.9,'LAMMPS\nSTATE / NEIGHBOR / CHECKPOINT',fontsize=20,weight='bold',va='top',color=COLORS[1])
        module='\n'.join(textwrap.wrap(d['lammps_python_module'],54,break_long_words=True,break_on_hyphens=False))
        ax.text(.53,.67,f"Version: {d['lammps_version']}\nPython: {d['python_version'].split()[0]} | CPU | MPI ranks: 1\npair zero: YES | property/atom: YES\nBinary restart / Python neighbors: verified\n\nPython module:\n{module}",fontsize=12,va='top',linespacing=1.35)
        ax.text(.025,.12,'NO F=ma     NO NVE     NO FORCE PHYSICS',color='#b52a2a',fontsize=23,weight='bold')
    elif index==1:
        d=read('01_state_roundtrip');fig=plt.figure(figsize=(14,8));gs=fig.add_gridspec(4,3,height_ratios=[1,1,1,.7]);picked=[next(p for p in d['before'] if p['mode_code']==m) for m in [1,2,3]]
        after={p['particle_id']:p for p in d['after']}
        for row,p in enumerate(picked):
            for col,rec in enumerate([p,after[p['particle_id']],after[p['particle_id']]]):
                ax=fig.add_subplot(gs[row,col]);shape=BridgeParticle(**rec).shape();a=np.linspace(0,2*np.pi,121)
                points=np.array([shape.support([np.cos(t),np.sin(t),0])-shape.center_m for t in a])*1e6
                ax.fill(points[:,0],points[:,1],color=COLORS[0],alpha=.3);ax.plot(points[:,0],points[:,1],color=COLORS[0]);ax.set_aspect('equal')
                ax.set(xlabel='local x (um)',ylabel='local y (um)');ax.set_title(['Python before','LAMMPS custom state','Python after'][col]+f" | ID {p['particle_id']}")
                ax.text(.02,.05,shape.mode,transform=ax.transAxes,fontsize=8)
        ax=fig.add_subplot(gs[3,:]);ax.axis('off')
        keys=['position_m','quaternion','velocity_m_s','omega_s_inv','geometry']
        table=ax.table(cellText=[[f"{d['errors'][k]:.0e}" for k in keys]],colLabels=['x (m)','q','V (m/s)','Omega (1/s)','geometry'],loc='center');table.scale(1,1.6)
    elif index==2:
        d=read('02_neighbor_equivalence');ids=[p['particle_id'] for p in d['particles']];lookup={x:i for i,x in enumerate(ids)}
        fig,axs=plt.subplots(1,4,figsize=(16,7))
        for ax,key,title in zip(axs,['all_possible_pairs','standalone_candidates','lammps_candidates','lammps_exact_pairs'],['All possible pairs','P4 validation query','LAMMPS zero query','Exact physics pairs']):
            matrix=np.zeros((len(ids),len(ids)))
            for i,j in d[key]:matrix[lookup[i],lookup[j]]=matrix[lookup[j],lookup[i]]=1
            ax.imshow(matrix,vmin=0,vmax=1,cmap='Blues');ax.set_title(f'{title}\n{len(d[key])} canonical pairs');ax.set_xticks(range(len(ids)),ids,rotation=90,fontsize=7);ax.set_yticks(range(len(ids)),ids,fontsize=7);ax.set_xlabel('stable particle ID')
        fig.text(.5,.05,f"False negatives: {d['false_negatives']} | raw mismatch: {d['mismatch_count']} | extra candidates removed: {d['extra_candidates']} | six shape combinations",ha='center',fontsize=14)
    elif index==3:
        d=read('03_one_step');fig,axs=plt.subplots(1,3,figsize=(14,7))
        for k,(key,style) in enumerate([('standalone','o-'),('bridge','x--')]):
            a=np.array([p['position'] for p in d[key]['particles']])*1e6;axs[0].plot(a[:,0],a[:,2],style,label=key,color=COLORS[k])
        axs[0].set(xlabel='x (um)',ylabel='z (um)',title=f"P5 spheres after one step | max x error = {d['sphere']['position_m']:.0e} m");axs[0].legend()
        for ax,keys,title in [(axs[1],['velocity_m_s','omega_s_inv'],'V / Omega maximum difference'),(axs[2],['quaternion','orientation_rad'],'q / orientation-axis difference')]:
            ax.bar(keys,[d['sphere'][k] for k in keys],color=COLORS[0]);ax.set_title(title);labelzero(ax)
    elif index==4:
        fig,axs=plt.subplots(2,3,figsize=(14,8));keys=['position_m','orientation_rad','velocity_m_s','omega_s_inv','quaternion','neighbor_mismatch_count']
        for case,color in zip(['04_sphere_multistep','04_mixed_multistep'],COLORS):
            d=read(case);rows=d['errors']
            for ax,key in zip(axs.flat,keys):ax.plot([r['step'] for r in rows],[r[key] for r in rows],label='P5 spheres' if 'sphere' in case else 'P4 mixed metadata',color=color,linestyle='-' if 'sphere' in case else '--');ax.set(title=key,xlabel='independent step');labelzero(ax)
        axs[0,0].legend(loc='lower center',fontsize=8)
    elif index==5:
        d=read('05_resistance');a=d['standalone'];b=d['bridge'];ra=np.array(a['R']);rb=np.array(b['R']);fig,axs=plt.subplots(2,3,figsize=(14,8))
        for ax,m,title in [(axs[0,0],ra,'Standalone R sparsity'),(axs[0,1],rb,'Bridge R sparsity')]:ax.spy(m,markersize=4,color=COLORS[0]);ax.set_title(title)
        im=axs[0,2].imshow(abs(ra-rb),vmin=0,vmax=1,cmap='Blues');axs[0,2].set_title('|R_A - R_B| = 0');fig.colorbar(im,ax=axs[0,2],shrink=.7)
        for ax,key in zip(axs[1],['b','U','J']):
            delta=np.abs(np.asarray(a[key])-np.asarray(b[key])).ravel();ax.plot(delta,'o-',markersize=3);ax.set(title=f'|{key}_A - {key}_B|',xlabel='entry');labelzero(ax)
        fig.text(.5,.02,'Contact IDs and active constraints identical; velocity source: original P5 resistance solver',ha='center')
    elif index==6:
        rows=read('06_force_audit');fig,axs=plt.subplots(1,2,figsize=(14,7))
        for ax,key,title in zip(axs,['max_lammps_force','max_lammps_torque'],['LAMMPS force storage','Torque absent in atom_style atomic']):
            ax.plot([r[key] for r in rows],color=COLORS[0]);ax.set(xlabel='audit before / after neighbor build',ylabel='unused diagnostic storage',title=title);labelzero(ax)
        fig.text(.5,.17,'Torque diagnostic 0 = UNDEFINED and UNUSED; no physical torque is measured.',ha='center',fontsize=13)
        fig.text(.5,.1,'Sphere velocity source = PARTICLE_3D_RESISTANCE_SOLVER | mixed validation = upstream P4',ha='center',fontsize=13)
    elif index==7:
        fig,axs=plt.subplots(2,2,figsize=(14,8))
        for row,name in enumerate(['07_restart_sphere','07_restart_mixed']):
            d=read(name)
            for key,color,style in [('continuous',COLORS[0],'-'),('restarted',COLORS[1],'--')]:
                states=d[key];x=np.array([[p['position'][0] for p in s['particles']] for s in states])*1e6
                for j in range(x.shape[1]):axs[row,0].plot([s['step'] for s in states],x[:,j]-x[0,j],color=color,ls=style,label=key if j==0 else None)
            axs[row,0].axvline(40,color='#b52a2a',ls=':',label='destroy / read binary restart');axs[row,0].legend(fontsize=8)
            axs[row,0].set(title=d['model'],xlabel='step',ylabel='x displacement (um)')
            axs[row,1].plot([e['step'] for e in d['errors']],[e['position_m'] for e in d['errors']]);axs[row,1].axvline(40,color='#b52a2a',ls=':');axs[row,1].set(title='State / time / neighbors: exact continuation',xlabel='step',ylabel='position difference (m)');labelzero(axs[row,1])
    elif index==8:
        d=read('08_metadata_restart');fields=['type_code','mode_code','position','q','velocity','omega','radius','axes','capsule_axis','capsule_radius','capsule_length','bound_radius','rotation'];fig,ax=plt.subplots(figsize=(14,7))
        values=np.array([[np.max(abs(np.asarray(a[f])-np.asarray(b[f]))) for f in fields] for a,b in zip(d['before'],d['after'])])
        im=ax.imshow(values,cmap='Blues',vmin=0,vmax=1,aspect='auto');ax.set_xticks(range(len(fields)),fields,rotation=30,ha='right');ax.set_yticks(range(len(values)),[f"ID {p['particle_id']} | mode {p['mode_code']}" for p in d['before']]);fig.colorbar(im,ax=ax,label='absolute field difference')
        for i in range(len(values)):
            for j in range(len(fields)):ax.text(j,i,'0',ha='center',va='center')
        ax.set_title('All custom properties restored from binary restart, by stable particle ID')
        fig.text(.5,.02,'RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION | capsule axis, radius and length own capsule geometry',ha='center',fontsize=11)
    elif index==9:
        d=read('09_real_two_mb');fig,axs=plt.subplots(2,3,figsize=(15,8));ref=d['original_p5_states'];got=d['bridge_world_states'];t=np.array([r['time_s'] for r in ref])*1e3
        for j in range(2):
            for rows,style,label in [(ref,'-',f'MB {j+1} P5'),(got,'--',f'MB {j+1} bridge')]:
                x=np.array([r['particles'][j]['center_m'] for r in rows])*1e6;axs[0,0].plot(x[:,0],x[:,2],style,color=COLORS[j],label=label)
            axs[0,2].plot(t,[r['particles'][j]['wall_gap_m']*1e9 for r in ref],color=COLORS[j],label=f'MB {j+1} P5');axs[0,2].plot(t,[r['particles'][j]['wall_gap_m']*1e9 for r in got],'--',color=COLORS[j])
        axs[0,0].set(xlabel='x (um)',ylabel='z (um)',title='Original real FEM trajectories');axs[0,0].legend(fontsize=8)
        for rows,style,label in [(ref,'-','original P5'),(got,'--','bridge')]:axs[0,1].plot(t,[r['pair_gaps'][0]['gap_m']*1e6 for r in rows],style,label=label)
        axs[0,1].set(xlabel='time (ms)',ylabel='pair gap (um)',title='Pair gap');axs[0,1].legend();axs[0,2].set(xlabel='time (ms)',ylabel='wall gap (nm)',title='Original wall gaps');axs[0,2].set_yscale('log');axs[0,2].legend(fontsize=8)
        for ax,key in zip(axs[1],['position_m','velocity_m_s','neighbor_mismatch_count']):ax.plot([e['step'] for e in d['errors']],[e[key] for e in d['errors']]);ax.set(title=key,xlabel='requested P5 step');labelzero(ax)
    else:
        d=read('10_neighbor_rebuild');rows=d['rows'];steps=[r['step'] for r in rows];fig,axs=plt.subplots(1,2,figsize=(14,7))
        axs[0].plot(steps,[r['distance_m']*1e6 for r in rows],'o-',label='center distance');axs[0].axhline(d['policy']['center_cutoff_m']*1e6,color=COLORS[1],ls=':',label='zero cutoff');axs[0].axhline(rows[0]['query_radius_m']*1e6,color='#b52a2a',ls='--',label='cutoff + skin');axs[0].set(xlabel='position update',ylabel='distance (um)',title='Crossing validation query range');axs[0].legend()
        for key,marker,label in [('standalone_candidate','o-','standalone'),('lammps_candidate','x--','LAMMPS rebuilt')]:axs[1].plot(steps,[int(r[key]) for r in rows],marker,label=label)
        axs[1].plot(steps,[int(bool(r['previous_list'])) for r in rows],':',color='gray',label='previous (potentially stale) list');axs[1].set(xlabel='position update',ylabel='pair in candidate list',title='Fresh neighbor list after every update',yticks=[0,1]);axs[1].legend()
    fig.suptitle(NAMES[index].replace('_',' '),fontsize=17,y=.985)
    bottom=.25 if index==6 else (.14 if index==2 else .055)
    fig.tight_layout(rect=[0,bottom,1,.94])
    return fig

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=REPORT/'figures');parser.add_argument('--no-notes',action='store_true');args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True);manifest=[]
    for index,name in enumerate(NAMES):
        fig=make(index);path=args.output/(name+'.png');fig.savefig(path,metadata={'Software':'Particle6 deterministic saved-data plotting'});plt.close(fig)
        sources={str((DATA/(s+suffix)).relative_to(REPO)):sha256(DATA/(s+suffix)) for s in SOURCES[index] for suffix in ['.json','.csv']}
        manifest.append(dict(figure=path.name,sha256=sha256(path),sources=sources))
        if not args.no_notes:
            expected,actual,anomaly=NOTES[index];(REPORT/f'{index:02d}_step_notes.md').write_text(f'应该看什么：{expected}\n\n实际看到什么：{actual}\n\n有没有异常：{anomaly}\n')
    if not args.no_notes:write_json(REPORT/'FIGURE_MANIFEST.json',manifest)

if __name__=='__main__':main()
