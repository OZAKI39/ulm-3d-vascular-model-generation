#!/usr/bin/env python3
"""Twelve reproducible audit figures drawn exclusively from saved P5 evidence."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle,Polygon
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json
REPORT=PACKAGE/'reports/particle5';DATA=REPORT/'data'
NAMES=['00_particle5_scope_and_resistance_contract','01_isolated_sphere_stokes_resistance',
'02_sphere_wall_lubrication','03_sphere_sphere_lubrication','04_rbc_mb_resistance_symmetry',
'05_sparse_resistance_matrix','06_resistance_dissipation_validation','07_resistance_metric_contact',
'08_lubrication_transit_delay','09_real_fem_near_wall_mb_lubrication',
'10_real_fem_two_mb_resistance_smoke','11_particle5_timestep_comparison']
SOURCES={0:['00_scope'],1:['01_stokes'],2:['02_wall'],3:['03_pair'],4:['04_algebraic'],5:['05_sparse'],6:['06_dissipation','05_sparse'],
7:['07_contact'],8:['08_delay','08_p4_wall_states','11_wall_dt1_states'],9:['09_real_near_wall'],
10:['10_initialization','10_real_summary','10_real_two_mb_states','10_pair_eligibility','09_real_near_wall'],
11:['11_timestep']+[f'11_{kind}_dt{divisor}_states' for kind in ['wall','pair'] for divisor in [1,2,4]]}
NOTES=[
('看清背景流、几何与新增阻力之间的数据关系。','6N 求解包含平移和自转；右端只由背景流的 self 项提供。','这里只定量验证球形法向润滑，生产截断和生产步长没有冻结。'),
('看平移和转动阻力是否分别随半径的一次方、三次方变化。','三个黏度、四个半径的计算点都落在解析线上，孤立粒子保持自由速度。','测试黏度扫描仅作公式验证，真实回放使用冻结来源核对的黏度。'),
('看间隙减小时阻力是否增大、法向速度是否连续降低。','阻力按间隙倒数增长，法向速度符合解析解，切向速度与自转保持原值。','0.1 等指定比值只验证渐近公式代数，不代表已验证任意间距的墙面流体动力学。'),
('看等大和不等大两球的相对法向速度是否随间隙缩小而下降。','两组都符合独立二元解析解；成对项作用相反，共同平移不受影响。','没有人为固定减速比例，也没有把远处配对当作近场。'),
('只看矩阵对称、交换互易和非负耗散。','使用测试单位下的正系数 1，交换球和原 P2 椭球的块顺序不改变矩阵。','这张图不是 RBC 润滑模型；测试系数不能进入真实物理接口。'),
('看 12 个球的非零块分布，以及候选筛选是否漏掉近场项。','复用 P4 候选查询后，稀疏矩阵和求解结果与穷举完全一致。','查询包围盒扩展只服务本轮近场验证，不改变球尺寸，也不是生产邻居层。'),
('看固定种子的大量速度下，各类耗散是否出现负值。','1024 个向量的 self、wall、pair 与总耗散均非负，self 缩放矩阵最小特征值为正。','舍入界与特征值原符号一起保存，没有对负特征值取绝对值。'),
('看相同流体速度起点下，两种接触修正是否都阻止继续入侵。','两球近墙和三球案例均满足约束；P5 的阻力加权修正代价更小。','乘子仍只是运动学约束量，不是物理接触力；微小负速度在原舍入界内。'),
('看同一初始间隙下，P4 恒速靠墙与 P5 逐渐减速的区别。','P4 到墙后停止法向运动；P5 随间隙降低而持续变慢，在验证时长内仍未碰墙。','没有为了产生接触而添加间隙下限，也不要求有限时间一定碰上。'),
('分别看原最小间隙样本和原轨迹中的靠墙样本，注意法向速度正负。','原最小间隙样本正在离墙；另一个原始靠墙样本的法向速度显著减小，切向变化仅为浮点误差。','没有移动中心或缩放球；这是局部法向修正，尚未包含切向、转动耦合、曲率和多墙效应。'),
('看原 P4 双球初值、尺寸和时长下的新轨迹，以及每一步的近场资格。','完整烟雾回放状态有限且无穿透；墙面近场启用，双球间隙比仍远大于 0.01，pair 润滑没有启用。','新轨迹产生更小墙隙是积分结果，不是改初值；验证窗口没有出口事件，也不代表完整悬浮液。'),
('看球墙和双球在步长减半后的间隙、速度和耗散差异。','两类案例完整覆盖物理时间，保持正间隙；末间隙差与独立隐式解析关系的误差随细化减小。','图中所有步长仅用于验证，不能据此选定生产步长。')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'figure.dpi':110,'savefig.dpi':145})


def read(name):return json.loads((DATA/(name+'.json')).read_text())


def curve(states,kind):
    t=np.array([s['time_s'] for s in states])
    h=np.array([s['particles'][0]['wall_gap_m'] if kind=='wall' else s['pair_gaps'][0]['gap_m'] for s in states])
    # Stored velocity is the value used over the preceding accepted interval.
    v=np.array([-s['particles'][0]['velocity_m_s'][2] if kind=='wall' else s['particles'][0]['velocity_m_s'][0]-s['particles'][1]['velocity_m_s'][0] for s in states[1:]])
    d=np.array([s['projection'].get('dissipation',{}).get('total',0.) for s in states[1:]])
    return t,h,v,d


def make(index):
    if index==0:
        fig,ax=plt.subplots(figsize=(13,8));ax.axis('off')
        boxes=[(.24,.82,'P0 frozen FEM\nU_free = [u, curl(u)/2]'),(.76,.82,'P3 / P4 original geometry\ngap, normals, touching constraints'),
               (.5,.53,'P5: R = R_self + R_wall + R_pair\nb = R_self U_free\nSolve R U = b, then simultaneous J U >= 0'),
               (.5,.22,'Sphere normal lubrication: quantitative V0\nNon-spherical lubrication: NOT FROZEN\nNo production cutoff / timestep\nNo LAMMPS / no CFD / no Particle-6')]
        for x,y,label in boxes:ax.text(x,y,label,ha='center',va='center',fontsize=13,bbox=dict(boxstyle='round,pad=.8',fc='#eaf1f6',ec='#225c83'))
        for x in [.24,.76]:ax.annotate('',xy=(.5,.65),xytext=(x,.73),arrowprops=dict(arrowstyle='->',lw=2,color='#225c83'))
        ax.annotate('',xy=(.5,.35),xytext=(.5,.42),arrowprops=dict(arrowstyle='->',lw=2,color='#225c83'))
        ax.set_title('Particle-5 | affine block resistance balance, full 6N',fontsize=17,pad=20)
    elif index==1:
        rows=read('01_stokes');fig,axes=plt.subplots(1,2,figsize=(13,5.5))
        for mu in sorted({r['mu_pa_s'] for r in rows}):
            rr=[r for r in rows if r['mu_pa_s']==mu];a=np.array([r['radius_m'] for r in rr]);dense=np.geomspace(min(a),max(a),100)
            for ax,key,power,factor in [(axes[0],'translation',1,6),(axes[1],'rotation',3,8)]:
                line,=ax.loglog(dense*1e6,factor*np.pi*mu*dense**power,label=f'analytic, mu={mu:g} Pa s')
                ax.loglog(a*1e6,[r[key] for r in rr],'o',mfc='white',color=line.get_color())
        for ax in axes:
            ax.set_xlabel('sphere radius (µm)');ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=8)
            ax.set_xticks([.25,.5,1,2,5],['0.25','0.5','1','2','5']);ax.xaxis.set_minor_formatter(NullFormatter())
        axes[0].set_ylabel('translation resistance (kg/s)');axes[1].set_ylabel('rotation resistance (kg m²/s)')
        axes[0].set_title('Stokes translation: 6 pi mu a');axes[1].set_title('Stokes rotation: 8 pi mu a³')
        fig.suptitle('Isolated spheres retain all six free-velocity components',fontsize=15)
    elif index==2:
        rows=read('02_wall');x=np.array([r['gap_ratio'] for r in rows]);fig,axes=plt.subplots(1,3,figsize=(16,5.3))
        axes[0].loglog(x,[r['coefficient'] for r in rows],'o-');axes[0].set_ylabel('normal resistance (kg/s)');axes[0].set_title('Leading resistance ~ 1/h')
        axes[1].loglog(x,[r['normal_attenuation'] for r in rows],'o',label='6N solve');axes[1].loglog(x,x/(1+x),'-',label='analytic h/(a+h)');axes[1].set_ylabel('Vn / U_free,n');axes[1].set_title('Normal motion slows');axes[1].legend()
        axes[2].semilogx(x,[r['velocity'][0]/r['free'][0] for r in rows],'o-',label='tangential translation')
        axes[2].semilogx(x,[r['velocity'][3]/r['free'][3] for r in rows],'x--',label='rotation');axes[2].set_ylim(.98,1.02);axes[2].set_ylabel('solved / free');axes[2].set_title('Unchanged modes');axes[2].legend(fontsize=8)
        for ax in axes:ax.set_xlabel('h/a [VALIDATION_GAPS_ONLY]');ax.grid(alpha=.2)
    elif index==3:
        rows=read('03_pair');fig,axes=plt.subplots(1,2,figsize=(13,5.5))
        for ai in sorted({r['radius_i_m'] for r in rows}):
            rr=[r for r in rows if r['radius_i_m']==ai];x=[r['gap_ratio'] for r in rr];label='equal radii' if ai==1e-6 else 'unequal radii (1:4)'
            axes[0].loglog(x,[r['coefficient'] for r in rr],'o-',label=label)
            axes[1].loglog(x,[r['relative_attenuation'] for r in rr],'o',label=label)
        x=np.logspace(-4,-1,80);axes[1].loglog(x,x/(1+x),'-',color='black',label='independent analytic relative solution')
        for ax in axes:ax.set_xlabel('h/R_eff [VALIDATION_GAPS_ONLY]');ax.grid(alpha=.2);ax.legend(fontsize=8)
        axes[0].set_ylabel('pair normal resistance (kg/s)');axes[1].set_ylabel('relative solved / free normal speed')
        fig.suptitle('Pair normal lubrication | common translation is a null mode',fontsize=15)
    elif index==4:
        row=read('04_algebraic');fig,axes=plt.subplots(1,2,figsize=(13,6.2));r=np.array(row['matrix'])
        im=axes[0].imshow(r,cmap='RdBu_r',vmin=-1,vmax=1);fig.colorbar(im,ax=axes[0],label='internally normalized TEST units')
        axes[0].axvline(5.5,color='black');axes[0].axhline(5.5,color='black');axes[0].set_title('RBC / MB algebraic 6 + 6 block');axes[0].set_xlabel('generalized component');axes[0].set_ylabel('generalized component')
        axes[1].axis('off');axes[1].text(.04,.8,f"zeta_test = 1 (validation only)\n\nSymmetry error: {row['symmetry_error']:.1e}\nSwap error: {row['swap_error']:.1e}\nMinimum eigenvalue: {row['minimum_eigenvalue']:.2e}\nRoundoff bound: {row['eigenvalue_roundoff_bound']:.2e}\n\nPhysical API rejects this block.\nNo RBC lubrication magnitude inferred.",va='top',fontsize=12)
        fig.suptitle('ALGEBRAIC ASSEMBLY TEST ONLY\nNOT PHYSICAL RBC LUBRICATION MODEL',fontsize=18,color='#b02424',weight='bold')
    elif index==5:
        row=read('05_sparse');fig,axes=plt.subplots(1,2,figsize=(13,6));r=np.array(row['matrix'])
        axes[0].spy(r,markersize=2,color='#1e6794');axes[0].set_title('12 spheres: 72 × 72 sparse resistance');axes[0].set_xlabel('generalized component');axes[0].set_ylabel('generalized component')
        axes[1].axis('off');axes[1].text(.03,.88,f"P4 AABB candidate reuse\n\nCandidate pairs: {len(row['candidate_pairs'])}\nAll pairs: {row['all_pairs_count']}\n\nmax |R_sparse - R_all| = {row['matrix_difference']:.1e}\nmax |U_sparse - U_all| = {row['velocity_difference']:.1e}\n\nSelf-scaled condition: {row['solver']['condition_estimate']:.5g}\nResidual: {row['solver']['relative_backward_residual']:.2e}\n\nQuery dilation: validation only\nParticle dimensions unchanged",va='top',fontsize=12)
        fig.suptitle('Sparse near-field assembly agrees exactly with all pairs',fontsize=15)
    elif index==6:
        rows=read('06_dissipation');matrix=read('05_sparse');fig,axes=plt.subplots(1,2,figsize=(13,5.8))
        for key in ['self','wall','pair','total']:
            values=np.array([r[key] for r in rows]);axes[0].hist(np.log10(values[values>0]),bins=30,histtype='step',lw=1.6,label=key)
        axes[0].set_xlabel('log10 dissipation (W)');axes[0].set_ylabel('fixed-seed sample count');axes[0].legend();axes[0].set_title('1024 random generalized velocities')
        keys=['self','wall','pair','total'];axes[1].bar(keys,[min(r[k] for r in rows) for k in keys],color=['#3477aa','#e8a43a','#3e9a79','#8c58a5']);axes[1].set_yscale('log');axes[1].set_ylabel('minimum sampled dissipation (W)');axes[1].set_title('All components are nonnegative')
        fig.suptitle(f"Minimum self-scaled eigenvalue = {matrix['solver']['minimum_scaled_eigenvalue']:.6g} (signed, no absolute value)",fontsize=14)
    elif index==7:
        rows=read('07_contact');fig,axes=plt.subplots(2,2,figsize=(13,8))
        for k,row in enumerate(rows):
            x=np.arange(len(row['p5_normal_speeds']));axes[0,k].bar(x-.16,np.array(row['p4_normal_speeds'])*1e6,width=.32,label='P4 geometry metric');axes[0,k].bar(x+.16,np.array(row['p5_normal_speeds'])*1e6,width=.32,label='P5 resistance metric');axes[0,k].axhline(0,color='black',lw=.8)
            axes[0,k].set_xticks(x,[f'contact {i+1}' for i in x]);axes[0,k].set_ylabel('corrected normal speed (µm/s)');axes[0,k].set_title(row['case'].replace('_',' '));axes[0,k].legend(fontsize=8)
            scale=max(np.max(np.abs(np.array(row['hydro']).reshape(-1,6)[:,:3]))*1e6,np.max(row['p5_normal_speeds'])*1e6)
            axes[0,k].set_ylim(-.1*scale,1.3*scale)
            axes[0,k].text(.03,.83,f"min P5 = {min(row['p5_normal_speeds']):.2e} m/s\nroundoff budget = {row['p5_audit']['contact_kkt']['velocity_budget_m_s']:.2e} m/s",transform=axes[0,k].transAxes,fontsize=9)
            axes[1,k].bar(['P4 metric','P5 metric'],[row['p4_resistance_objective'],row['p5_resistance_objective']],color=['#327dab','#efaa49']);axes[1,k].set_ylabel('0.5 deltaU.T R deltaU (W)');axes[1,k].set_title('Same unconstrained hydro velocity as reference')
        fig.suptitle('Simultaneous hard contact | multiplier is KINEMATIC_CONSTRAINT_MULTIPLIER',fontsize=14)
    elif index==8:
        fig,axes=plt.subplots(1,2,figsize=(13,5.5))
        for label,source in [('P4 hard contact','08_p4_wall_states'),('P5 self + wall resistance','11_wall_dt1_states')]:
            t,h,v,_=curve(read(source),'wall');axes[0].plot(t,h*1e9,'o-',ms=2,label=label)
            if source=='08_p4_wall_states':axes[1].step(h[:-1]*1e9,v*1e6,where='post',label=label)
            else:axes[1].plot(h[:-1]*1e9,v*1e6,'o-',ms=3,label=label)
        axes[0].set_xlabel('physical time (s)');axes[0].set_ylabel('surface gap (nm)');axes[1].set_xlabel('gap at velocity evaluation (nm)');axes[1].set_ylabel('normal approach speed (µm/s)')
        for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle('Transit delay comes from resistance | no fixed velocity multiplier',fontsize=15)
    elif index==9:
        rows=read('09_real_near_wall');fig,axes=plt.subplots(2,2,figsize=(14,9))
        for k,row in enumerate(rows):
            ax=axes[0,k];values=np.array([row['free_normal_velocity_m_s'],row['lubricated_normal_velocity_m_s']])*1e6
            ax.bar(['original free','P5 normal correction'],values,color=['#b5b9be','#2784ad']);ax.set_yscale('symlog',linthresh=.001);ax.axhline(0,color='black',lw=.7);ax.set_ylabel('signed inward-normal speed (µm/s)')
            if k==0:ax.set_ylim(0,abs(values).max()*5)
            else:ax.set_ylim(-abs(values).max()*5,0)
            ax.set_title('Original minimum gap: SEPARATING' if k==0 else 'Original approaching state: APPROACHING')
            for j,val in enumerate(values):ax.annotate(f'{val:.6g}',(j,val),xytext=(0,8 if val>0 else -15),textcoords='offset points',ha='center',fontsize=9)
            ax.text(.53 if k==0 else .03,.94,f"h = {row['gap_m']:.6e} m\nh/a = {row['gap_ratio']:.6g}\nVn/U_free,n = {row['attenuation_ratio']:.6g}",transform=ax.transAxes,va='top',fontsize=10)
        row=rows[1];ax=axes[1,0];center=np.array(row['center_m']);n=np.array(row['normal']);e=np.eye(3)[np.argmin(abs(n))];e=(e-e@n*n);e/=np.linalg.norm(e);basis=np.array([e,n])
        for tri in row['wall_patch_m']:
            points=(np.array(tri)-center)@basis.T*1e6;ax.add_patch(Polygon(points,facecolor='#aab6c1',edgecolor='#788592',alpha=.15,lw=.5))
        a=row['radius_m']*1e6;ax.add_patch(Circle((0,0),a,ec='#247da7',fc='#c4e1ee',alpha=.7));ax.plot(0,0,'o',color='#247da7');ax.set_xlim(-2*a,2*a);ax.set_ylim(-2*a,2*a);ax.set_aspect('equal');ax.set_xlabel('local tangent projection (µm)');ax.set_ylabel('local inward normal (µm)');ax.set_title('Unmoved original sphere / original WALL triangles')
        axes[1,1].axis('off');axes[1,1].text(.03,.9,f"Approaching sample: P4 row {row['source_row_index']}, ID {row['original_particle']['particle_id']}\nOriginal time: {row['source_time_s']:.9g} s\n\nTangential change: {row['tangential_difference_m_s']:.3e} m/s\nAngular change: {row['angular_difference_s_inv']:.3e} 1/s\n\nOriginal center, radius and positive gap retained.\nLocal leading normal correction only.\nTangential / rotation coupling / curvature /\nmultiple-wall hydrodynamics are not covered.",va='top',fontsize=11)
        fig.suptitle('Real frozen-state replay | do not confuse separation with approach',fontsize=15)
    elif index==10:
        states=read('10_real_two_mb_states');summary=read('10_real_summary');fig,axes=plt.subplots(2,2,figsize=(14,9));t=np.array([s['time_s'] for s in states])*1e3
        colors=['#2582b0','#d77630']
        for k in range(2):
            points=np.array([s['particles'][k]['center_m'] for s in states])*1e6
            axes[0,0].plot(points[:,0],points[:,2],color=colors[k],label=f"MB {states[0]['particles'][k]['particle_id']}")
            axes[0,0].plot(*points[0,[0,2]],'o',color=colors[k]);axes[0,0].plot(*points[-1,[0,2]],'x',color=colors[k])
            axes[1,0].semilogy(t,[s['particles'][k]['wall_gap_m']*1e9 for s in states],color=colors[k],label=f'MB {k+1}')
        axes[0,0].set_aspect('equal',adjustable='datalim');axes[0,0].set_xlabel('x (µm)');axes[0,0].set_ylabel('z (µm)');axes[0,0].set_title('Same P4 initial states; original radii');axes[0,0].legend()
        ratio=[s['pair_gaps'][0]['gap_m']/summary['relevant_pair_radius_m'] for s in states];axes[0,1].semilogy(t,ratio,label='actual gap / R_eff');axes[0,1].axhline(.01,color='#b32c2c',ls='--',label='validation eligibility only');axes[0,1].set_ylabel('pair gap / R_eff');axes[0,1].set_xlabel('physical time (ms)');axes[0,1].set_title('PAIR LUBRICATION NOT ACTIVE');axes[0,1].legend(fontsize=8)
        axes[1,0].set_xlabel('physical time (ms)');axes[1,0].set_ylabel('actual WALL gap (nm)');axes[1,0].legend();axes[1,0].set_title('Accepted positive gaps, no artificial floor')
        axes[1,1].axis('off');axes[1,1].text(.03,.92,f"Accepted states: {summary['rows']}\nPhysical time: {summary['final_time_s']:.9g} s\nTime coverage error: {summary['time_coverage_error_s']:.1e} s\nMinimum pair gap: {summary['minimum_pair_gap_m']:.6e} m\nMinimum pair ratio: {summary['minimum_pair_gap_ratio']:.6g}\nActive pair blocks: {summary['pair_active_block_count']}\nActive wall blocks: {summary['wall_active_block_count']}\nMaximum scaled condition: {summary['maximum_condition']:.6g}\nMaximum residual: {summary['maximum_residual']:.2e}\nNo outlet event in this window.\nRBC hydrodynamics: DEFERRED",va='top',fontsize=11)
        fig.suptitle('Real two-MB resistance smoke | one-way frozen FEM, finite validation window',fontsize=15)
    else:
        fig,axes=plt.subplots(3,2,figsize=(14,11))
        for col,kind in enumerate(['wall','pair']):
            for divisor in [1,2,4]:
                t,h,v,d=curve(read(f'11_{kind}_dt{divisor}_states'),kind);label='dt' if divisor==1 else f'dt/{divisor}'
                axes[0,col].semilogy(t,h*1e9,label=label);axes[1,col].semilogy(t[1:],v*1e6,label=label);axes[2,col].plot(t[1:],d,label=label)
            axes[0,col].set_title('Sphere–wall' if kind=='wall' else 'Sphere–sphere')
            for row in range(3):axes[row,col].set_xlabel('physical time (s)');axes[row,col].legend(fontsize=8);axes[row,col].grid(alpha=.2)
            axes[0,col].set_ylabel('gap (nm)');axes[1,col].set_ylabel('normal approach speed (µm/s)');axes[2,col].set_ylabel('total dissipation (W)')
        fig.suptitle('NOT PRODUCTION TIMESTEP SELECTION\nPositive gaps throughout; no P5 hard-contact event within these horizons',fontsize=15)
    fig.tight_layout(rect=(0,0,1,.94 if index not in [0,4,11] else .90))
    return fig


def generate(output, *, notes=True):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);manifest=[]
    for i,name in enumerate(NAMES):
        fig=make(i);path=output/(name+'.png');fig.savefig(path,metadata={'Software':'Particle-5 reproducible audit'});plt.close(fig)
        sources=[DATA/(name+ext) for name in SOURCES[i] for ext in ['.json','.csv']]
        if notes:
            a,b,c=NOTES[i]
            (REPORT/f'{i:02d}_step_notes.md').write_text(f"### 图 {i:02d}：{name}\n\n![图 {i:02d}](figures/{name}.png)\n\n应该看什么：{a}\n\n实际看到什么：{b}\n\n有没有异常：{c}\n\n数据："+'、'.join(f'[{p.name}](data/{p.name})' for p in sources)+'。\n')
        manifest.append(dict(figure=path.name,sha256=sha256(path),sources={str(p.relative_to(REPO)):sha256(p) for p in sources}))
    if notes:write_json(REPORT/'FIGURE_MANIFEST.json',manifest)
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=REPORT/'figures');p.add_argument('--no-notes',action='store_true');args=p.parse_args();generate(args.output,notes=not args.no_notes)
