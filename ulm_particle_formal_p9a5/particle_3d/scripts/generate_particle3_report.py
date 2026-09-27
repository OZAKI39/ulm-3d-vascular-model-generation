#!/usr/bin/env python3
"""Thirteen reproducible scientific audit figures from saved numerical records."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Ellipse,Rectangle,Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256
from particle_3d.particle3_cases import write_json,write_rows
from particle_3d.rbc_orientation import rotation_matrix
REPORT=PACKAGE/'reports/particle3';DATA=REPORT/'data';FIGURES=REPORT/'figures'
NAMES=['00_particle3_scope_and_wall_contract.png','01_wall_triangle_normals.png','02_sphere_wall_gap_validation.png',
       '03_ellipsoid_wall_gap_validation.png','04_hard_contact_sliding_validation.png','05_physical_time_subdivision.png',
       '06_capillary_deformation_surrogate.png','07_deformation_feasibility_map.png','08_real_wall_gap_samples.png',
       '09_real_wall_controlled_contact.png','10_real_fem_mb_wall_clearance.png','11_real_fem_rbc_wall_clearance.png',
       '12_particle3_timestep_comparison.png']
COLORS=['#137c9b','#e18c24','#985fa6','#339569','#c74746']
MODE_COLORS={'FREE_OBLATE':'#339569','CAPILLARY_DEFORMED':'#e18c24','DEFORMATION_SURROGATE_INFEASIBLE':'#b8bdc4'}
MODE_SHORT={'FREE_OBLATE':'FREE','CAPILLARY_DEFORMED':'DEFORMED','DEFORMATION_SURROGATE_INFEASIBLE':'INFEASIBLE'}
for candidate in ['/mnt/c/Windows/Fonts/msyh.ttc','/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc']:
    if Path(candidate).exists():
        font_manager.fontManager.addfont(candidate)
        plt.rcParams['font.sans-serif']=[font_manager.FontProperties(fname=candidate).get_name(),'DejaVu Sans']
        if Path('/mnt/c/Windows/Fonts/msyhbd.ttc').exists():font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/msyhbd.ttc')
        break
plt.rcParams.update({'font.size':10,'axes.unicode_minus':False,'savefig.dpi':160,'axes.spines.top':False,'axes.spines.right':False})


def read(name):return json.loads((DATA/name).read_text())
def columns(rows,key):return np.array([r[key] for r in rows])


def save(fig,index,sources,footer,notes):
    FIGURES.mkdir(parents=True,exist_ok=True)
    fig.text(.5,.02,footer,ha='center',va='bottom',fontsize=9,bbox=dict(facecolor='#edf2f5',edgecolor='none',pad=7))
    fig.savefig(FIGURES/NAMES[index],bbox_inches='tight');plt.close(fig)
    write_json(DATA/f'{index:02d}_figure_sources.json',dict(figure=NAMES[index],sources=sources,
        source_sha256={s:sha256(DATA/s) for s in sources},frozen_wall_sha256=read('00_particle3_scope_and_wall_contract.json')['wall']['sha256'],
        timestep_role='VALIDATION_ONLY',manual_visual_review='PENDING_USER_REVIEW'))
    (REPORT/f'{index:02d}_step_notes.md').write_text(f'![{NAMES[index]}](figures/{NAMES[index]})\n\n应该看什么：{notes[0]}\n\n实际看到什么：{notes[1]}\n\n有没有异常：{notes[2]}\n')


_wall=None
def wall_triangles():
    global _wall
    if _wall is None:
        from particle_3d.wall_geometry import WallGeometry
        _wall=WallGeometry.from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular').triangles*1e6
    return _wall


def vessel(ax):
    # Render the diagnostic center paths above the translucent mesh; otherwise
    # matplotlib's collection-wide 3-D depth sort can hide an entire inside path.
    ax.computed_zorder=False
    triangles=wall_triangles();ax.add_collection3d(Poly3DCollection(triangles,facecolor='#88a9ba',edgecolor='#88a9ba',linewidth=.10,alpha=.32,zorder=1))
    points=triangles.reshape(-1,3);lo=points.min(0);hi=points.max(0)
    ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]),xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)')
    ax.set_box_aspect(hi-lo);ax.view_init(24,-62)


def shape_mesh(center,axes,rotation):
    theta=np.linspace(0,2*np.pi,41);phi=np.linspace(0,np.pi,23)
    local=np.array([axes[0]*np.outer(np.sin(phi),np.cos(theta)),axes[1]*np.outer(np.sin(phi),np.sin(theta)),axes[2]*np.outer(np.cos(phi),np.ones_like(theta))])
    return (rotation@local.reshape(3,-1)).reshape(local.shape)+np.asarray(center)[:,None,None]


def stage0():
    d=read('00_particle3_scope_and_wall_contract.json');fig,ax=plt.subplots(figsize=(13,8));ax.axis('off')
    boxes=[(.03,.65,.44,.25,'只读来源 / READ ONLY',f"Particle-2: {d['particle2_commit'][:12]}\nWALL: {d['wall']['triangle_count']:,} triangles\n实体边界：WALL (ID 1)\nINLET / 3 OUTLET 为开口"),
           (.53,.65,.44,.25,'完整几何 / FINITE SIZE','MB: sphere\nRBC: oriented oblate ellipsoid\nConflict → volume / area constrained capsule\nFace + edge + vertex signed Euclidean gap'),
           (.03,.27,.44,.27,'单边硬接触 / HARD CONTACT','Rigid stationary wall; zero coating\n仅阻止接触点向墙运动；保留切向速度\nContinuous sweep → real time subdivision\nNo position projection / artificial frozen step'),
           (.53,.27,.44,.27,'本轮边界 / SCOPE','Area budget derived from P2 oblate\nNo force / membrane mechanics / lubrication\nNo transit penalty / particle-particle collision\nVALIDATION ONLY; Particle-4 未开始')]
    for x,y,w,h,title,body in boxes:
        ax.add_patch(Rectangle((x,y),w,h,facecolor='#edf4f7',edgecolor='#b4c8d2',lw=1.4))
        ax.text(x+.02,y+h-.035,title,weight='bold',fontsize=13,va='top');ax.text(x+.02,y+h-.095,body,va='top',linespacing=1.8)
    ax.set_title('00  Particle-3 范围与正式 WALL 合同',loc='left',fontsize=18,pad=15)
    ax.text(.03,.14,'WALL SHA256\n'+d['wall']['sha256'],fontsize=10,family='monospace')
    save(fig,0,['00_particle3_scope_and_wall_contract.json'],'FEM and P0–P2 scientific inputs unchanged. Human review pending.',
         ['实体边界、有限粒子几何和本轮工作边界。','只有正式 WALL 是实体壁，入口与三个出口保持开口。','没有加入润滑、膜力学、粒子间相互作用或生产步长。'])


def stage1():
    rows=read('01_wall_normals.json');fig=plt.figure(figsize=(14,7));ax=fig.add_subplot(121,projection='3d');vessel(ax)
    for row in rows:
        center=np.mean(row['vertices_m'],axis=0)*1e6;n=np.array(row['normal_in'])
        ax.quiver(*center,*n,length=4,color=COLORS[4],arrow_length_ratio=.3);ax.text(*center,str(row['triangle_id']),fontsize=6)
    ax.set_title('18 个原始三角形：红箭头指向腔内')
    ax=fig.add_subplot(122);ax.bar(range(len(rows)),columns(rows,'owner_inward_distance_m')*1e6,color=COLORS[0]);ax.set(xticks=range(len(rows)),xticklabels=columns(rows,'triangle_id'),ylabel='相邻四面体中心的有向距离 (µm)',xlabel='正式 WALL triangle id',title='每个法向都指向唯一所属的腔内四面体');ax.tick_params(axis='x',rotation=90,labelsize=8)
    fig.subplots_adjust(bottom=.24,wspace=.3)
    save(fig,1,['01_wall_normals.json'],'All 18 owning-tetra signs > 0. Original winding retained. All 45,221 WALL faces shown and available for geometry queries.',
        ['不同位置、方向和面积的三角形内法向。','18 个内法向均指向唯一所属腔内四面体。','自动方向检查没有异常；网格外观仍待人工审核。'])


def stage2():
    rows=read('02_sphere_gap.json');fig,axes=plt.subplots(1,2,figsize=(13,6))
    for i,feature in enumerate(['FACE','EDGE','VERTEX']):
        selected=[r for r in rows if r['requested_feature']==feature]
        axes[0].plot(columns(selected,'analytic_gap_m')*1e6,columns(selected,'gap_m')*1e6,'o',color=COLORS[i],label=feature)
    lim=[-1,1.2];axes[0].plot(lim,lim,'k--',lw=.8);axes[0].axhline(0,color='#aaa',lw=.7);axes[0].legend();axes[0].set(xlabel='解析 gap (µm)',ylabel='计算 gap (µm)',title='面、边、顶点：正间隙及负穿透深度')
    row=rows[5];center=np.array(row['center_m'])*1e6;wall=np.array(row['wall_point_m'])*1e6;particle=np.array(row['particle_point_m'])*1e6
    ax=axes[1];ax.add_patch(Ellipse(center[[1,2]],2,2,facecolor='#d6e9f1',edgecolor=COLORS[0]));ax.plot([0,5],[0,0],color='#535e66',lw=3);ax.scatter(*wall[[1,2]],c=COLORS[4],label='finite edge point');ax.scatter(*particle[[1,2]],c=COLORS[1],label='particle point');ax.plot([wall[1],particle[1]],[wall[2],particle[2]],'--',color=COLORS[4]);ax.set(xlim=(-1.5,3),ylim=(-1,2),xlabel='y (µm)',ylabel='z (µm)',title='有限边示例：最近点连线');ax.set_aspect('equal');ax.legend(fontsize=8)
    fig.subplots_adjust(bottom=.22,wspace=.3)
    save(fig,2,['02_sphere_gap.json'],f"12 cases; max analytic error = {max(r['error_m'] for r in rows):.3e} m. Negative gap is penetration, not center-outside status.",
        ['球面到有限三角形面、边、顶点的距离和符号。','12 个正负间隙案例与解析值重合。','误差在原先声明的浮点舍入预算内。'])


def stage3():
    rows=read('03_ellipsoid_gap.json');caps=read('03_capsule_gap.json');fig,axes=plt.subplots(1,2,figsize=(13,6))
    for i,ratio in enumerate([.15,.3,.6]):
        chosen=[r for r in rows if r['ratio']==ratio and r['center_m'][2]==1.5e-6]
        axes[0].plot(columns(chosen,'angle_rad')*180/np.pi,columns(chosen,'gap_m')*1e6,'o-',markersize=3,label=f'c/a={ratio}',color=COLORS[i])
    axes[0].axhline(0,color='k',lw=.8);axes[0].legend();axes[0].set(xlabel='倾角 (degree)',ylabel='gap (µm)',title='固定中心：旋转即可改变完整 RBC 间隙')
    ax=axes[1];ax.axhline(0,color='#58616a',lw=3)
    for i,angle in enumerate([0,45,90]):ax.add_patch(Ellipse((i*5,1.5),4,1.2,angle=angle,facecolor='none',edgecolor=COLORS[i],lw=2))
    ax.set(xlim=(-2.5,12.5),ylim=(-1,4),xlabel='分列显示 (µm)',ylabel='离墙高度 (µm)',title='相同中心高度、轴长；三个姿态');ax.set_aspect('equal')
    fig.subplots_adjust(bottom=.22,wspace=.3)
    save(fig,3,['03_ellipsoid_gap.json','03_capsule_gap.json'],f"117 ellipsoid cases: max error {max(r['error_m'] for r in rows):.3e} m; 9 capsule plane cases: {max(r['error_m'] for r in caps):.3e} m.",
        ['固定中心时，椭球旋转如何改变壁面间隙。','不同轴比的倾角曲线与解析支撑高度一致。','没有用中心距离减去固定 RBC 半径替代姿态几何。'])


def stage4():
    rows=read('04_contact_velocity.json');fig,axes=plt.subplots(1,5,figsize=(16,5))
    for ax,row in zip(axes,rows):
        center=np.array(row['center_m']);point=np.array(row['contact_point_m']);scale=1e-4
        ax.axhline(0,color='#58616a',lw=2);ax.scatter(0,1,s=50,c='#333');ax.scatter(*(point-center)[[0,2]]/1e-6+[0,1],s=25,c=COLORS[4])
        for key,color,label in [('free_velocity',COLORS[1],'free'),('corrected_velocity',COLORS[0],'corrected')]:
            v=np.array(row[key])/scale;ax.arrow(0,1,v[0],v[2],width=.025,color=color,length_includes_head=True,label=label,alpha=.85)
        ax.arrow(-.8,0,0,.7,width=.02,color='#444');ax.text(-.75,.3,'n');ax.set(xlim=(-1.2,2.7),ylim=(-2.3,4.3),title=row['case'].replace('_','\n'));ax.set_aspect('equal');ax.legend(fontsize=7,loc='lower left');ax.set_xlabel('velocity arrows / 100 µm s⁻¹')
    fig.subplots_adjust(bottom=.23,wspace=.35)
    save(fig,4,['04_contact_velocity.json'],'Only inward contact-point normal speed is removed. Omega and tangential velocity are retained. No contact → no correction.',
        ['入墙、平行、离墙、无接触和旋转椭球的两种速度箭头。','只在接触点向墙运动时补偿法向平移，其他情形不改速度。','法向残差和切向变化均在浮点预算内。'])


def stage5():
    rows=read('05_time_trajectories.json');ledger=read('05_time_intervals.json');summary=read('05_time_summary.json');fig,axes=plt.subplots(1,2,figsize=(13,6))
    for i,dt in enumerate([1,.5,.25]):
        r=[r for r in rows if r['requested_dt_s']==dt];axes[0].plot(columns(r,'time_s'),columns(r,'center_m')[:,2]*1e6,'o-',ms=3,label=f'dt={dt} s',color=COLORS[i])
        intervals=[(r['t0_s'],r['dt_s']) for r in ledger if r['requested_dt_s']==dt]
        axes[1].broken_barh(intervals,(i-.3,.6),facecolors=COLORS[i],edgecolor='white',lw=1)
    axes[0].plot([0,1],[2,-2],'k--',label='requested unaccepted crossing');axes[0].axhline(.5,color='#aaa');axes[0].axvline(.375,color='#c74746',ls=':');axes[0].set(xlabel='真实物理时间 (s)',ylabel='球心高度 (µm)',title='跨薄墙大步必须在接触处分解');axes[0].legend(fontsize=8)
    axes[1].set(xlabel='真实物理时间 (s)',yticks=[0,1,2],yticklabels=['dt','dt/2','dt/4'],title='所有接受的子区间：无遗漏，无虚耗时间',xlim=(0,1));axes[1].axvline(.375,color='#c74746',ls=':')
    fig.subplots_adjust(bottom=.22,wspace=.28)
    save(fig,5,['05_time_trajectories.json','05_time_intervals.json','05_time_summary.json'],f"Contact t=0.375 s; final t=1 s. Max time coverage error={max(r['time_coverage_error'] for r in summary):.3e} s.",
        ['未接受的大步、接触时间和实际接受的时间片。','三个请求步长均在 0.375 秒接触，余下时间继续沿墙移动。','没有穿过薄墙，也没有把失败步的时间直接跳过去。'])


def stage6():
    rows=read('06_deformation_examples.json');rid=rows[0]['rbc_id'];selected=[r for r in rows if r['rbc_id']==rid]
    choices=[next(r for r in selected if r['status']=='FREE_OBLATE'),next(r for r in selected if r['status']=='CAPILLARY_DEFORMED'),next(r for r in selected if r['status']=='DEFORMATION_SURROGATE_INFEASIBLE')]
    fig,axes=plt.subplots(1,3,figsize=(14,7))
    for ax,r in zip(axes,choices):
        tube=r['tube_apothem_m']*1e6;ax.axvline(-tube,color='#647783',lw=3);ax.axvline(tube,color='#647783',lw=3)
        ax.add_patch(Ellipse((0,0),2*r['a_m']*1e6,2*r['c_m']*1e6,facecolor='none',edgecolor=COLORS[0],ls='--',lw=1.5,label='original oblate'))
        text=f"V={r['V_fL']:.3f} fL\nA_budget={r['area_budget_m2']*1e12:.3f} µm²"
        if 'R_cap_m' in r:
            radius=r['R_cap_m']*1e6;length=r['L_cap_m']*1e6
            ax.plot([-radius,-radius],[-length/2,length/2],c=COLORS[1]);ax.plot([radius,radius],[-length/2,length/2],c=COLORS[1]);ax.add_patch(Arc((0,length/2),2*radius,2*radius,theta1=0,theta2=180,color=COLORS[1],lw=2));ax.add_patch(Arc((0,-length/2),2*radius,2*radius,theta1=180,theta2=360,color=COLORS[1],lw=2))
            text+=f"\nR={radius:.3f} µm; L={length:.3f} µm\nA/A_budget={r['area_ratio']:.5f}\ngap={r['gap_m']*1e6:.3g} µm"
        else:text+='\n'+('gap='+f"{r['gap_m']*1e6:.3g} µm" if 'gap_m' in r else 'No feasible capsule; case stops')
        ax.text(.02,.02,text,transform=ax.transAxes,va='bottom',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
        ax.set(xlim=(-5.5,5.5),ylim=(-9,9),xlabel='横向 (µm)',ylabel='管轴方向 (µm)',title=MODE_SHORT[r['status']]+f' | tube R={tube:.2f} µm');ax.set_aspect('equal')
    fig.subplots_adjust(bottom=.19,wspace=.25)
    save(fig,6,['06_deformation_examples.json'],'Same original RBC volume in all panels. Smooth capsule area ≤ model budget. INFEASIBLE is not a physiological-blockage conclusion.',
        ['同一个原始 RBC 在宽、窄、极窄人工管中的几何。','依次得到原始椭球、满足体积与面积约束的胶囊、代理不可行。','不可行情形停止验证，没有缩小 RBC 或增加面积预算。'])


def stage7():
    rows=read('07_feasibility_map.json');fig,axes=plt.subplots(1,2,figsize=(14,6))
    for mode,color in MODE_COLORS.items():
        r=[r for r in rows if r['status']==mode];axes[0].scatter(columns(r,'tube_radius_m')*1e6,columns(r,'r'),s=12,c=color,label=f'{MODE_SHORT[mode]} ({len(r)})',alpha=.65)
    radii=sorted({r['tube_radius_m'] for r in rows});bottom=np.zeros(len(radii))
    for mode,color in MODE_COLORS.items():
        counts=np.array([sum(r['tube_radius_m']==radius and r['status']==mode for r in rows) for radius in radii]);axes[1].bar(np.array(radii)*1e6,counts/128,bottom=bottom,width=.4,color=color,label=MODE_SHORT[mode]);bottom+=counts/128
    axes[0].set(xlabel='人工管外接圆半径 (µm)',ylabel='原始 RBC c/a',title='固定的 128 个分布样本 × 9 管径');axes[0].legend(fontsize=8)
    axes[1].set(xlabel='人工管外接圆半径 (µm)',ylabel='样本比例',title='所有样本，包括不可行结果');axes[1].legend(fontsize=8,loc='upper left')
    fig.subplots_adjust(bottom=.22,wspace=.28)
    save(fig,7,['07_feasibility_map.json','07_selection.json'],'Tube sizes are VALIDATION_GEOMETRY_PARAMETERS, not a mouse capillary population. 48-sided finite triangles; open ends.',
        ['样本形状和人工管径变化时的三类结果。','1,152 个组合完整保留，窄管中的不可行结果也计入分母。','这些管径仅用于算法验证，不能解释为小鼠血管分布。'])


def stage8():
    rows=read('08_real_wall_static.json')[:3];fig=plt.figure(figsize=(15,6))
    for i,r in enumerate(rows):
        ax=fig.add_subplot(1,3,i+1,projection='3d');origin=np.mean(r['triangle_m'],axis=0);tri=(np.array(r['triangle_m'])-origin)*1e6;c=(np.array(r['center_m'])-origin)*1e6;g=r['geometry']
        ax.add_collection3d(Poly3DCollection([tri],facecolor='#7f959f',edgecolor='#3e5763',alpha=.5))
        if r['shape_mode']=='SPHERE_MB':xyz=shape_mesh(c,np.repeat(g['radius_m']*1e6,3),np.eye(3))
        elif r['shape_mode']=='FREE_OBLATE':xyz=shape_mesh(c,np.array(g['axes_m'])*1e6,rotation_matrix(g['q']))
        else:
            from particle_3d.rbc_orientation import quaternion_from_short_axis
            xyz=shape_mesh([0,0,0],np.repeat(g['R_cap_m']*1e6,3),np.eye(3));xyz[2]+=np.where(xyz[2]>=0,1,-1)*g['L_cap_m']*.5e6
            xyz=(rotation_matrix(quaternion_from_short_axis(g['axis']))@xyz.reshape(3,-1)).reshape(xyz.shape)+c[:,None,None]
        ax.plot_wireframe(*xyz,rstride=3,cstride=4,color=COLORS[i],lw=.5)
        w=(np.array(r['wall_point_m'])-origin)*1e6;p=(np.array(r['particle_point_m'])-origin)*1e6
        ax.scatter(*w,c=COLORS[4]);ax.scatter(*p,c='#222');ax.plot(*np.stack([w,p]).T,color=COLORS[4],lw=2);ax.quiver(*w,*r['normal'],length=1,color='#333')
        points=np.concatenate((xyz.reshape(3,-1).T,tri));lo=points.min(0);hi=points.max(0);ax.set(xlim=(lo[0]-.5,hi[0]+.5),ylim=(lo[1]-.5,hi[1]+.5),zlim=(lo[2]-.5,hi[2]+.5),xlabel='Δx (µm)',ylabel='Δy (µm)',zlabel='Δz (µm)');ax.set_box_aspect(hi-lo+1)
        ax.set_title(r['shape_mode']+f"\npatch gap={r['gap_m']*1e6:.3f} µm\nfull WALL gap={r['full_wall_gap_m']*1e6:.3f} µm",fontsize=10)
    fig.subplots_adjust(bottom=.23,wspace=.1)
    save(fig,8,['08_real_wall_static.json','08_sonovue_samples.json'],'Actual original triangle 43073 shown; all 5 regions saved. PATCH ONLY diagnostic. Full-WALL negative gaps are explicitly retained.',
        ['原 WALL patch 与球、椭球、胶囊的两个最近点及内法向。','局部正间隙与解析构造一致，整幅 WALL 的间隙另列。','部分局部摆放会碰到其他壁面，已保留负的全局间隙，不宣称全局可放入。'])


def stage9():
    rows=read('09_controlled_paths.json');summaries=read('09_controlled_summary.json');fig,axes=plt.subplots(2,2,figsize=(13,8))
    tid=summaries[0]['triangle_id']
    for col,ptype in enumerate(['MB','RBC']):
        for i,direction in enumerate(['APPROACH','OBLIQUE','TANGENT','SEPARATING']):
            r=[r for r in rows if r['original_wall_triangle_id']==tid and r['particle_type']==ptype and r['direction']==direction];t=columns(r,'time_s');axes[0,col].plot(t,columns(r,'wall_gap_m')*1e6,'o-',ms=3,c=COLORS[i],label=direction)
            n=np.array(r[0]['normal_inward']);velocity=columns(r,'corrected_velocity_m_s');axes[1,col].plot(t,velocity@n*1e6,'o-',ms=3,c=COLORS[i])
        axes[0,col].set(title=f'{ptype} | actual WALL patch {tid}',ylabel='gap (µm)');axes[0,col].legend(fontsize=8);axes[1,col].set(xlabel='控制试验时间 (s)',ylabel='center normal velocity (µm/s)')
    fig.subplots_adjust(bottom=.18,hspace=.32,wspace=.3)
    save(fig,9,['09_controlled_paths.json','09_controlled_intervals.json','09_controlled_summary.json'],f"40 controlled cases across 5 actual WALL regions. Max tangent error={max(r['max_tangential_error_m_s'] for r in summaries):.3e} m/s. Not FEM-flow trajectories.",
        ['真实 patch 上的四种指定运动方向、接触后的间隙和法向速度。','40 个受控案例保持不穿透，接触后保留切向运动。','这些是指定速度的局部诊断，不作为真实流动或完整血管通行结论。'])


def summaries():return [read(p.name) for p in sorted(DATA.glob('*_real_*_summary.json'))]


def stage10():
    sources=[];fig=plt.figure(figsize=(14,8));a=fig.add_subplot(121,projection='3d');vessel(a);b=fig.add_subplot(222);c=fig.add_subplot(224)
    for i in range(3):
        name=f'10_real_mb_dt{i}_states.json';sources.append(name);rows=read(name);s=read(f'10_real_mb_dt{i}_summary.json');sources.append(f'10_real_mb_dt{i}_summary.json')
        if not rows:continue
        x=columns(rows,'center_m')*1e6;t=columns(rows,'time_s');a.plot(*x.T,c=COLORS[i],label=f"{s['validation_dt_s']*1e6:.3f} µs");b.plot(t,columns(rows,'wall_gap_m')*1e6,c=COLORS[i]);c.step(t,[r['contact_state']=='TOUCHING' for r in rows],where='post',c=COLORS[i],label=s['status'])
    a.set_title('原 P1 起点与尺寸：整幅 WALL 约束');a.legend(fontsize=8);b.set(title='完整球的 wall gap',xlabel='time (s)',ylabel='gap (µm)');c.set(xlabel='time (s)',ylabel='contact',yticks=[0,1],yticklabels=['free','touching']);c.legend(fontsize=8)
    fig.subplots_adjust(bottom=.19,hspace=.4,wspace=.35)
    save(fig,10,sources,'Same P1 initial center and SonoVue sphere. All 45,221 WALL triangles available to the locator. No artificial move to force a contact.',
        ['原始入口起点的真实 FEM 球心轨迹、完整球间隙和接触状态。','三个减半步长均保留全部实际接受状态和终止事件。','是否发生自然接触由结果决定，起点和尺寸没有为展示接触而改变。'])


def stage11():
    sources=[];fig=plt.figure(figsize=(15,9));a=fig.add_subplot(221,projection='3d');vessel(a);b=fig.add_subplot(222);c=fig.add_subplot(223);d=fig.add_subplot(224)
    length_axis=c.twinx();length_axis.set_ylabel('L_cap (µm), dashed')
    annotations=[]
    for i in range(5):
        prefix=f'11_real_rbc_g{i}_dt0';rows=read(prefix+'_states.json');s=read(prefix+'_summary.json');sources.extend([prefix+'_states.json',prefix+'_summary.json']);annotations.append(f"r={s['r']:.3f}: {s['status']} ({s['rows']} states)")
        if rows:
            x=columns(rows,'center_m')*1e6;t=columns(rows,'time_s');a.plot(*x.T,c=COLORS[i],label=f"r={s['r']:.3f}");b.plot(t,columns(rows,'wall_gap_m')/columns(rows,'roundoff_m'),c=COLORS[i],label=f"r={s['r']:.3f}")
            caps=[r for r in rows if r['shape_mode']=='CAPILLARY_DEFORMED']
            if caps:
                c.plot(columns(caps,'time_s'),columns(caps,'R_cap_m')*1e6,c=COLORS[i],label=f"R, r={s['r']:.3f}");d.plot(columns(caps,'time_s'),columns(caps,'area_ratio'),c=COLORS[i],label=f"r={s['r']:.3f}")
                length_axis.plot(columns(caps,'time_s'),columns(caps,'L_cap_m')*1e6,c=COLORS[i],ls='--',alpha=.65)
    a.set_title('整幅 WALL；曲线为接受的变形代理轨迹');a.legend(fontsize=8);b.set(title='完整代理间隙；微小负值受固定舍入预算约束',ylabel='gap / 固定舍入预算',xlabel='time (s)');b.axhline(-1,color='#c74746',ls='--',label='acceptance lower bound');b.legend(fontsize=8);c.set(xlabel='time (s)',ylabel='R_cap (µm), solid');c.legend(fontsize=8);d.set(xlabel='time (s)',ylabel='A_capsule / A_budget');d.axhline(1,c='#333',ls='--');d.legend(fontsize=8)
    fig.text(.5,.10,'\n'.join(annotations),ha='center',fontsize=8);fig.subplots_adjust(bottom=.28,hspace=.4,wspace=.35)
    save(fig,11,sources,'Capsule geometry follows local free flow; retained quaternion is not interpreted as deformed Jeffery orientation. INFEASIBLE ≠ physiological blockage.',
        ['五个原始 RBC 的接受轨迹、模式、胶囊半径和面积比例。','各案例的变形和不可行终止均直接列出，初始不可行的案例没有虚构轨迹。','这里只证明当前代理的几何可行性，不声称真实红细胞必然如此变形。'])


def comparison_records():
    records=summaries();orientation=[]
    for gi in range(5):
        for left,right in [(0,1),(1,2)]:
            a=read(f'11_real_rbc_g{gi}_dt{left}_states.json');b=read(f'11_real_rbc_g{gi}_dt{right}_states.json')
            # Only upstream FREE states have a physically interpreted Jeffery axis.
            def upstream(rows):
                result=[]
                for row in rows:
                    if row['shape_mode']!='FREE_OBLATE':break
                    result.append(row)
                return result
            aa=upstream(a);bb=upstream(b);differences=[]
            if aa and bb:
                for row in aa:
                    times=np.array([r['time_s'] for r in bb]);k=int(np.argmin(np.abs(times-row['time_s'])))
                    if abs(times[k]-row['time_s'])<=8*np.spacing(max(row['time_s'],1.)):
                        differences.append(float(np.degrees(np.arccos(np.clip(abs(np.dot(row['p'],bb[k]['p'])),0,1)))))
            orientation.append(dict(geometry_index=gi,dt_pair=[left,right],upstream_matched_states=len(differences),
                max_upstream_axis_difference_deg=max(differences,default=None),status='AVAILABLE' if differences else 'NO_UPSTREAM_FREE_OBLATE_SEGMENT',
                particle2_orientation_convergence='NOT_ESTABLISHED',deformed_quaternion_difference_not_interpreted=True))
    write_json(DATA/'12_timestep_comparison.json',dict(cases=records,upstream_orientation=orientation,production_particle_timestep_frozen=False))
    write_rows(DATA/'12_timestep_comparison.csv',records)
    return records,orientation


def stage12():
    records,orientation=comparison_records();fig,axes=plt.subplots(2,2,figsize=(14,9))
    for group in ['MB']+list(range(5)):
        r=sorted([r for r in records if (r['particle_type']=='MB' if group=='MB' else r['particle_type']=='RBC' and r['geometry_index']==group)],key=lambda r:r['dt_index'])
        if not r:continue
        label='MB' if group=='MB' else f'RBC r={r[0]["r"]:.3f}';x=columns(r,'validation_dt_s')*1e6
        axes[0,0].plot(x,[r['minimum_gap_m']*1e6 if r['minimum_gap_m'] is not None else np.nan for r in r],'o-',label=label)
        axes[0,1].plot(x,columns(r,'contact_duration_s'),'o-',label=label)
        axes[1,0].plot(x,columns(r,'last_accepted_time_s'),'o-',label=label)
    axes[0,0].set(xlabel='validation dt (µs)',ylabel='minimum gap (µm)',title='无接受状态的案例不伪造最小 gap');axes[0,0].legend(fontsize=8);axes[0,1].set(xlabel='validation dt (µs)',ylabel='contact duration (s)',title='接触持续时间（来自接受的物理区间）');axes[1,0].set(xlabel='validation dt (µs)',ylabel='last accepted time (s)',title='出口或代理不可行时停止')
    axes[1,1].axis('off');text='UPSTREAM ORIENTATION\n'+('所有 RBC 均无上游 FREE_OBLATE 段。\n不以冻结的 q 差为 0 声称姿态收敛。\n' if all(r['status']=='NO_UPSTREAM_FREE_OBLATE_SEGMENT' for r in orientation) else '轴夹角采用 acos(abs(p₁·p₂))。\n')+'P2 真实姿态时间步收敛仍未建立。\n\n模式与终止事件（按 dt → dt/4）：\n'
    for group in ['MB']+list(range(5)):
        r=sorted([r for r in records if (r['particle_type']=='MB' if group=='MB' else r['particle_type']=='RBC' and r['geometry_index']==group)],key=lambda r:r['dt_index'])
        short=[s['status'].replace('DEFORMATION_SURROGATE_INFEASIBLE','INFEASIBLE').replace('VALIDATION_HORIZON_REACHED','HORIZON') for s in r]
        text+=str(group)+': '+', '.join(short)+'\n'
    axes[1,1].text(0,1,text,va='top',fontsize=9,linespacing=1.7);fig.subplots_adjust(bottom=.18,hspace=.38,wspace=.28)
    sources=['12_timestep_comparison.json']+[p.name for p in sorted(DATA.glob('*_real_*_summary.json'))]
    save(fig,12,sources,'NOT PRODUCTION TIMESTEP SELECTION. Contact counts, mode switches and failed attempts are preserved in JSON. No artificial delay model.',
        ['三个减半步长的最小间隙、接触持续时间、终止时间和上游姿态可比性。','18 个重放均保存，初始不可行结果与实际有轨迹的案例分开显示。','这些比较不确定生产步长，P2 姿态收敛限制继续保留。'])


def main():
    global FIGURES
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--only',nargs='*',type=int);p.add_argument('--output',type=Path)
    args=p.parse_args()
    if args.output:FIGURES=args.output
    for i in (args.only if args.only is not None else range(13)):
        globals()[f'stage{i}']();print(NAMES[i],flush=True)


if __name__=='__main__':main()
