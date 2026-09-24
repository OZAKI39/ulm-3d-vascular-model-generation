#!/usr/bin/env python3
"""Sixteen reproducible audit figures from saved CSV/JSON; no trajectory changes."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule,unit
from particle_3d.particle3_cases import write_json
from particle_3d.audit import sha256
REPORT=PACKAGE/'reports/particle4';DATA=REPORT/'data'
NAMES=['00_particle4_scope_and_contact_contract','01_pair_gap_geometry_matrix','02_mb_mb_head_on_contact','03_mb_mb_glancing_contact',
'04_rbc_mb_contact_geometry','05_rbc_rbc_orientation_contact','06_off_center_angular_contact','07_capsule_pair_contact',
'08_pair_physical_time_subdivision','09_three_particle_simultaneous_contact','10_contact_order_invariance','11_wall_plus_particle_multicontact',
'12_all_pairs_vs_broadphase','13_mixed_four_particle_contact_trajectory','14_real_fem_two_mb_contact_smoke','15_particle4_timestep_comparison']
SOURCES={0:['00_scope'],1:['01_geometry_matrix','01_independent_references'],2:['02_head_on_projection','02_head_on_states'],3:['03_glancing_projection','03_glancing_states'],
4:['04_rbc_mb_geometry'],5:['05_rbc_orientation_geometry'],6:['06_offcenter_projection'],7:['07_capsule_projection'],
8:['08_no_tunnel_summary']+[f'08_no_tunnel_dt{i}_{s}' for i in range(3) for s in ['states','intervals']],9:['09_simultaneous_projections'],10:['10_order_invariance'],11:['11_wall_pair_projection'],
12:['12_broadphase_comparison','12_broadphase_summary'],13:['13_mixed_dt2_states','13_mixed_dt2_intervals'],14:['14_real_two_mb_states','14_real_two_mb_summary','14_initialization'],15:['15_timestep_comparison']}
NOTES=[
('看清 P4 继承哪些输入、只新增什么接触规则。','本阶段复用 P0–P3，新增粒子间有限尺寸硬接触和同时求解。','真实 RBC 通行仍未建立；没有润滑、黏附、新变形或生产步长。'),
('看六类配对在分离、接触、穿透时是否返回正确符号。','六类均有固定方向和旋转几何检查，输入交换后间隙不变、法向反向。','负间隙只用于几何诊断，不作为轨迹初始状态。'),
('看两球正碰时两边速度修正是否对称。','接触后相对法向速度消失，两边中心修正相反。','没有反弹或自旋修正；箭头仅作显示。'),
('看擦边接触时切向运动是否保留。','只去掉继续入侵的法向分量，切向速度保留，轨迹随后分开。','分开来自原有切向运动，不是人为反弹。'),
('看五个分位 RBC 与球的真实接触点和法向。','使用原 P2 样本及三种姿态，不用一个平均 RBC 代替分布。','轮廓是标明平面的几何投影，接触判断仍为三维。'),
('看同样中心位置下转动 RBC 怎样改变间隙。','两个不同 P2 样本随相对姿态变化出现不同间隙。','这里是静态几何检查，可能显示穿透，不是被接受的轨迹。'),
('看偏心接触时中心速度和角速度如何共同修正。','自由椭球的角速度修正来自同一个几何加权求解，接触法向速度满足约束。','乘子不是接触力，图中角速度变化不代表碰撞冲量。'),
('看胶囊与三类粒子接触时哪些量会变化。','胶囊只接受平移速度修正，轴、尺寸和保存的自旋保持不变。','胶囊不通过碰撞进一步压扁或拉长。'),
('看端点分离的大步是否仍检测到中途相撞。','三个请求步长均在 0.375 秒接触，并完整处理到 1 秒。','所有细分区间都有真实起止时间，失败尝试不消耗时间。'),
('看直线和三角形三粒子能否同时满足全部约束。','所有接触一起求解，图中列出修正前后法向速度。','没有按粒子输入顺序逐对修正。'),
('看保持粒子 ID 不变时，排列输入是否改变结果。','三粒子的全部排列和四粒子的全部排列给出相同速度、角速度及接触集合。','没有按位置重新编号。'),
('看粒子被 WALL 与另一粒子夹住时是否两边都满足约束。','WALL 和粒子接触进入同一个求解，球和 RBC 混合案例均通过。','固定 WALL 的作用不要求所有粒子的平移修正总和为零。'),
('看候选筛选有没有漏掉真正接触的配对。','同一批混合形状的候选结果与穷举精确接触结果一致。','候选可以多给，但不能漏给；耗时只作烟雾检查。'),
('看四粒子轨迹中的三类接触事件及最小间隙。','两球和两个原 P2 椭球产生球–球、球–RBC、RBC–RBC 接触，所有接受状态满足舍入界。','此例使用中心线接触及零自由角速度，不能据此宣称真实姿态收敛。'),
('看真实 WALL 内两只原尺寸 MB 的轨迹、间隙和出口记录。','冻结流场驱动两只 MB，初始壁面及粒子间间隙均为正；自然接触结果直接标在图中。','这是 256 个验证步的烟雾检查，不是完整悬浮液；未到出口会明确记录。'),
('看步长减半后接触时间、间隙、末位置、速度和短轴差。','三个步长均只用于人工混合场景验证，比较的是短轴而非 quaternion 分量。','NOT PRODUCTION TIMESTEP SELECTION；真实 RBC 通行限制仍保留。')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'figure.dpi':120,'savefig.dpi':140})


def read(name):return json.loads((DATA/(name+'.json')).read_text())


def restore(record):
    if record['type']=='Sphere':return Sphere(record['center_m'],record['radius_m'])
    if record['type']=='Ellipsoid':return Ellipsoid(record['center_m'],record['axes_m'],record['rotation'])
    return Capsule(record['center_m'],record['axis_world'],record['radius_m'],record['cylindrical_length_m'])


def outline(ax,shape,basis,origin,color,label=None):
    points=[]
    for angle in np.linspace(0,2*np.pi,361):
        n=np.cos(angle)*basis[0]+np.sin(angle)*basis[1]
        points.append(basis@(shape.support(n)-origin)*1e6)
    points=np.array(points);ax.fill(points[:,0],points[:,1],color=color,alpha=.12);ax.plot(points[:,0],points[:,1],color=color,lw=1.5,label=label)
    center=basis@(shape.center_m-origin)*1e6;ax.plot(*center,'o',color=color,ms=4)


def contact_picture(ax,record,title='',wall=False):
    shapes={int(i):restore(s) for i,s in record['shapes'].items()};ids=sorted(shapes);origin=shapes[ids[0]].center_m
    basis=np.array([[1,0,0],[0,0,1.]]) if wall else np.array([[1.,0,0],[0,1.,0]])
    gaps=record['gaps']
    if len(ids)==2 and not wall:
        n=-np.array(gaps[0]['normal_j_to_i']);relative=np.array(gaps[0]['point_i_m'])-origin
        tangent=relative-(relative@n)*n
        if np.linalg.norm(tangent)<1e-12*max(s.bounding_radius_m for s in shapes.values()):
            axis=np.eye(3)[np.argmin(np.abs(n))];tangent=axis-(axis@n)*n
        basis=np.array([n,unit(tangent)])
    colors=['#1776b6','#dc7633','#28926d','#9256aa'];radius=max(s.bounding_radius_m for s in shapes.values())
    vscale=max(np.linalg.norm(v) for v in record['free_velocities'].values());arrow_scale=.65*radius/max(vscale,np.finfo(float).tiny)
    for k,i in enumerate(ids):
        s=shapes[i];outline(ax,s,basis,origin,colors[k],f'ID {i}: {s.mode}')
        start=basis@(s.center_m-origin)*1e6
        for key,color,style in [('free_velocities','#777777','--'),('velocities','#b92236','-')]:
            velocity=np.array(record[key].get(str(i),record[key].get(i)))
            end=start+basis@velocity*arrow_scale*1e6
            ax.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',color=color,lw=1.6,linestyle=style))
    for gap in gaps:
        if len(ids)>2 and gap['state']!='TOUCHING':continue
        pi=basis@(np.array(gap['point_i_m'])-origin)*1e6;pj=basis@(np.array(gap['point_j_m'])-origin)*1e6
        ax.plot(*pi,'o',mfc='white',mec='black',ms=5);ax.plot(*pj,'x',color='black',ms=5)
        ax.plot([pi[0],pj[0]],[pi[1],pj[1]],':',color='black',lw=.8)
        normal=basis@np.array(gap['normal_j_to_i'])*radius*.45*1e6
        ax.annotate('n',xy=pi+normal,xytext=pi,arrowprops=dict(arrowstyle='->',color='#278c52'),color='#278c52')
    if wall:
        height=(basis@(np.array([0,0,0])-origin))[1]*1e6;ax.axhline(height,color='black',lw=3,label='fixed WALL')
    ax.set_aspect('equal',adjustable='datalim');ax.set_title(title);ax.set_xlabel('projection coordinate 1 (µm)');ax.set_ylabel('projection coordinate 2 (µm)')
    ax.margins(.22);ax.text(.01,.01,'DISPLAY SCALE ONLY: gray free / red corrected arrows',transform=ax.transAxes,fontsize=7)
    ax.legend(fontsize=7,loc='upper right')
    return basis,origin


def static_picture(ax,row,title):
    record=dict(shapes={1:row['shape_i'],2:row['shape_j']},gaps=[row['gap']],free_velocities={1:[0,0,0],2:[0,0,0]},velocities={1:[0,0,0],2:[0,0,0]})
    contact_picture(ax,record,title)
    ax.text(.02,.94,f"gap = {row['gap']['gap_m']*1e6:.4g} µm",transform=ax.transAxes,fontsize=8)


def plot(number):
    if number==0:
        fig,ax=plt.subplots(figsize=(13,6));ax.axis('off');scope=read('00_scope')
        ax.text(.04,.86,'P0 frozen field → P1 MB → P2 RBC → P3 WALL',fontsize=22,color='#205f80')
        ax.text(.04,.65,'P4: finite-size pair gap + simultaneous hard constraints',fontsize=20,weight='bold')
        ax.text(.04,.41,'Frictionless · symmetric translation · FREE_OBLATE angular correction\nGeometry metric only; multiplier is NOT physical force',fontsize=16,linespacing=1.8)
        ax.text(.04,.19,'NO lubrication / adhesion / new deformation / LAMMPS / production timestep\nCarry forward: REAL RBC PASSAGE NOT ESTABLISHED',fontsize=14,color='#a33b30',linespacing=1.6)
        ax.text(.04,.015,'P3 tested commit: '+scope['particle3_dependency_commit'],fontsize=10)
    elif number==1:
        rows=read('01_geometry_matrix');labels=list(dict.fromkeys(r['pair_type'] for r in rows));matrix=np.array([[r['gap']['gap_m']*1e6 for r in rows if r['pair_type']==label] for label in labels])
        fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'width_ratios':[1.5,1]});im=axes[0].imshow(matrix,cmap='coolwarm',vmin=-1,vmax=1,aspect='auto')
        axes[0].set_yticks(range(6),labels);axes[0].set_xticks(range(9),['sep','touch','overlap']*3,rotation=45);axes[0].set_title('Three support directions × signed gap (µm)');fig.colorbar(im,ax=axes[0])
        refs=read('01_independent_references');axes[1].plot([r['error_m'] for r in refs],'.');axes[1].set_title('Independent analytic / convex-reference errors');axes[1].set_ylabel('error / bracket width (m)');axes[1].set_xlabel('reference case');axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    elif number in [2,3,6]:
        name={2:'02_head_on',3:'03_glancing',6:'06_offcenter'}[number];record=read(name+'_projection')
        fig,axes=plt.subplots(1,2,figsize=(14,6));basis,origin=contact_picture(axes[0],record,{2:'Head-on: symmetric correction',3:'Glancing: tangent retained',6:'Off-center FREE_OBLATE + MB'}[number])
        if number==6:
            ids=sorted(record['shapes']);values=record['angular_corrections'][ids[0]]
            axes[1].bar(['delta Omega x','delta Omega y','delta Omega z'],values,color='#1776b6');axes[1].set_ylabel('angular correction (1/s)');axes[1].set_title('Signed angular correction: FREE_OBLATE only')
            axes[1].text(.03,.97,'r × n / ell² determines correction direction\nKinematic multiplier; NO force / impulse',transform=axes[1].transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
            gap=record['gaps'][0];i=str(gap['particle_i_id']);center=np.array(record['shapes'][i]['center_m']);point=np.array(gap['point_i_m']);r=point-center
            start=basis@(center-origin)*1e6;contact=basis@(point-origin)*1e6
            axes[0].plot([start[0],contact[0]],[start[1],contact[1]],'--',color='#8b4d9d',lw=1.4);axes[0].text(*((start+contact)/2),'r',color='#8b4d9d',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
            vc_free=np.array(record['free_velocities'][i])+np.cross(record['free_omegas'][i],r)
            vc_new=np.array(record['velocities'][i])+np.cross(record['omegas'][i],r)
            scale=.6*restore(record['shapes'][i]).bounding_radius_m/max(np.linalg.norm(vc_free),np.linalg.norm(vc_new))
            for velocity,color,label in [(vc_free,'#555555','v_c free'),(vc_new,'#b92236','v_c corrected')]:
                end=contact+basis@velocity*scale*1e6;axes[0].annotate('',xy=end,xytext=contact,arrowprops=dict(arrowstyle='->',color=color))
                axes[0].annotate(label,xy=end,xytext=(3,9 if label=='v_c corrected' else -16),textcoords='offset points',fontsize=8,color=color)
        else:
            rows=read(name+'_states')
            for i in range(2):axes[1].plot([r['time_s'] for r in rows],[r['particles'][i]['center_m'][0]*1e6 for r in rows],label=f"ID {rows[0]['particles'][i]['particle_id']}")
            axes[1].set_xlabel('physical time (s)');axes[1].set_ylabel('center x (µm)');axes[1].set_title('Accepted trajectories');axes[1].legend()
        free_label=', '.join(f'{v:.3e}' for v in record['audit']['free_normal_speeds_m_s']);new_label=', '.join(f'{v:.3e}' for v in record['audit']['corrected_normal_speeds_m_s'])
        fig.text(.51,.025,f'normal speed (m/s): free {free_label} → corrected {new_label}',fontsize=8)
    elif number==4:
        rows=read('04_rbc_mb_geometry');fig,axes=plt.subplots(2,3,figsize=(15,9))
        for k in range(5):static_picture(axes.flat[k],rows[3*k+2],f"r quantile {rows[3*k+2]['quantile']:.0%}; original ID {rows[3*k+2]['rbc_id']}")
        axes.flat[5].axis('off');axes.flat[5].text(.03,.7,'5 original P2 population quantiles\n3 orientations each\nTrue 3D support contact\nNo average-RBC replacement',fontsize=14,linespacing=1.8)
    elif number==5:
        rows=read('05_rbc_orientation_geometry');fig,axes=plt.subplots(2,2,figsize=(14,10))
        for ax,index in zip(axes.flat,[0,6,12]):static_picture(ax,rows[index],f"relative tilt {rows[index]['angle_rad']*180/np.pi:.0f}°")
        axes.flat[3].plot([r['angle_rad']*180/np.pi for r in rows],[r['gap']['gap_m']*1e6 for r in rows],'o-');axes.flat[3].set_xlabel('tilt (degrees)');axes.flat[3].set_ylabel('pair gap (µm)');axes.flat[3].axhline(0,color='black',lw=.8);axes.flat[3].set_title('Centers fixed; geometry only')
    elif number==7:
        records=read('07_capsule_projection');fig,axes=plt.subplots(1,3,figsize=(17,6))
        for ax,row in zip(axes,records):contact_picture(ax,row,'CAPSULE + '+row['other_shape'])
        fig.text(.3,.01,'Capsule axis, dimensions and retained spin unchanged by contact',fontsize=12)
    elif number==8:
        fig,axes=plt.subplots(1,3,figsize=(16,6));cases=read('08_no_tunnel_summary')
        for i,c in enumerate(cases):
            rows=read(f'08_no_tunnel_dt{i}_states');ledger=read(f'08_no_tunnel_dt{i}_intervals')
            axes[0].plot([r['time_s'] for r in rows],[r['pair_gaps'][0]['gap_m']*1e6 for r in rows],'.-',label=f"dt={c['validation_dt_s']}")
            axes[1].hlines([i]*len(ledger),[r['t0_s'] for r in ledger],[r['t1_s'] for r in ledger],lw=4);axes[1].plot([r['t1_s'] for r in ledger],[i]*len(ledger),'|',color='white',ms=8)
        axes[0].legend();axes[0].set(xlabel='physical time (s)',ylabel='pair gap (µm)',title='Contact at 0.375 s');axes[1].set(xlabel='physical time (s)',yticks=range(3),yticklabels=['dt','dt/2','dt/4'],title='Accepted physical intervals')
        t=np.linspace(0,1,101);axes[2].plot(t,np.abs(4-8*t)-1,'--',label='unchecked free trial');axes[2].plot(t,np.maximum(3-8*t,0),label='hard contact');axes[2].axhline(0,color='black');axes[2].set(xlabel='time (s)',ylabel='gap (µm)',title='Separated endpoints can hide collision');axes[2].legend()
    elif number==9:
        records=read('09_simultaneous_projections')[:2];fig,axes=plt.subplots(2,2,figsize=(14,10))
        for row,record in enumerate(records):
            contact_picture(axes[row,0],record,record['case']+' simultaneous constraints');a=record['audit'];x=np.arange(len(a['multipliers']))
            axes[row,1].bar(x-.15,np.array(a['free_normal_speeds_m_s'])*1e6,.3,label='free');axes[row,1].bar(x+.15,np.array(a['corrected_normal_speeds_m_s'])*1e6,.3,label='corrected');axes[row,1].set(xlabel='canonical contact',ylabel='normal speed (µm/s)');axes[row,1].legend()
    elif number==10:
        rows=read('10_order_invariance');fig,axes=plt.subplots(1,2,figsize=(14,6));axes[0].plot([r['velocity_error_m_s'] for r in rows],'o');axes[1].plot([r['angular_error_s_inv'] for r in rows],'o')
        axes[0].set(ylabel='velocity difference (m/s)',title='Stable IDs; all 6 + 6 + 24 permutations');axes[1].set(ylabel='angular difference (1/s)',title='Contact sets identical')
        for ax in axes:ax.set_xlabel('deterministic input permutation');ax.grid(alpha=.2)
    elif number==11:
        records=read('11_wall_pair_projection');fig,axes=plt.subplots(1,2,figsize=(15,7))
        for ax,record in zip(axes,records):contact_picture(ax,record,record['case'],wall=True)
        fig.text(.25,.025,'All WALL and pair rows enter ONE kinematic projection',fontsize=14)
    elif number==12:
        rows=read('12_broadphase_comparison');s=read('12_broadphase_summary');fig,axes=plt.subplots(1,2,figsize=(14,6))
        axes[0].bar(['all pairs','candidates','true contacts'],[s['all_pairs'],s['candidates'],s['true_contact_pairs']],color=['#aaa','#1776b6','#28926d']);axes[0].set_title(f"AABB sweep; false negatives = {s['false_negatives']}")
        selected=[r for r in rows if r['candidate']];axes[1].scatter([r['all_pairs_gap_m']*1e6 for r in selected],[r['narrowphase_gap_m']*1e6 for r in selected],s=18);axes[1].set(xlabel='all-pairs gap (µm)',ylabel='candidate exact gap (µm)',title='Same exact narrowphase geometry')
    elif number==13:
        rows=read('13_mixed_dt2_states');fig,axes=plt.subplots(1,3,figsize=(17,6));times=[r['time_s'] for r in rows]
        for k,p in enumerate(rows[0]['particles']):axes[0].plot(times,[r['particles'][k]['center_m'][0]*1e6 for r in rows],label=f"{p['particle_id']} {p['shape_mode']}")
        axes[0].set(xlabel='time (s)',ylabel='center x (µm)',title='Two MB + two original P2 RBC');axes[0].legend(fontsize=8)
        for k,g in enumerate(rows[0]['pair_gaps']):
            if k in [0,3,5]:axes[1].plot(times,[r['pair_gaps'][k]['gap_m']*1e6 for r in rows],label=str(g['canonical_pair_id']))
        axes[1].set(xlabel='time (s)',ylabel='pair gap (µm)',title='Three contact types');axes[1].legend(fontsize=8)
        axes[2].step(times,[sum(g['state']=='TOUCHING' for g in r['pair_gaps']) for r in rows],where='post');axes[2].set(xlabel='time (s)',ylabel='touching pair count',title='Every accepted state checked')
    elif number==14:
        rows=read('14_real_two_mb_states');result=read('14_real_two_mb_summary');fig=plt.figure(figsize=(15,10));ax=fig.add_subplot(221,projection='3d')
        from particle_3d.wall_geometry import WallGeometry
        wall=WallGeometry.from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular');paths=[np.array([r['particles'][k]['center_m'] for r in rows]) for k in range(2)]
        allpoints=np.concatenate(paths);margin=3*max(p['radius_m'] for p in rows[0]['particles']);lo=allpoints.min(axis=0)-margin;hi=allpoints.max(axis=0)+margin
        mask=np.all(wall.triangles.max(axis=1)>=lo,axis=1)&np.all(wall.triangles.min(axis=1)<=hi,axis=1);origin=allpoints.mean(axis=0)
        ax.add_collection3d(Poly3DCollection((wall.triangles[mask]-origin)*1e6,facecolor='#aac9d2',edgecolor='#779aa5',linewidth=.15,alpha=.10))
        for k,path in enumerate(paths):ax.plot(*((path-origin)*1e6).T,lw=2.5,label=f"MB {rows[0]['particles'][k]['particle_id']}")
        for setlim,l,h,o in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],lo,hi,origin):setlim((l-o)*1e6,(h-o)*1e6)
        ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title='Actual WALL local view; display crop only');ax.legend()
        times=np.array([r['time_s'] for r in rows])*1e3;ax=fig.add_subplot(222);ax.plot(times,[r['pair_gaps'][0]['gap_m']*1e6 for r in rows]);ax.set(xlabel='time (ms)',ylabel='pair gap (µm)',title=result['pair_contact_status'].replace('_',' '))
        ax=fig.add_subplot(223)
        for k in range(2):ax.plot(times,[r['particles'][k]['wall_gap_m']*1e6 for r in rows],label=f'MB {k+1}')
        ax.set(xlabel='time (ms)',ylabel='finite-size WALL gap (µm)',title='Full WALL queried at every accepted state');ax.legend()
        ax=fig.add_subplot(224);ax.step(times,[int(r['pair_gaps'][0]['state']=='TOUCHING') for r in rows],where='post');ax.set(xlabel='time (ms)',ylabel='pair touching (0 / 1)',title='Outlet events: '+('NONE WITHIN SMOKE HORIZON' if not result['outlet_events'] else str(result['outlet_events'])))
    else:
        data=read('15_timestep_comparison');cases=data['cases'];fig,axes=plt.subplots(2,3,figsize=(16,9));x=['dt','dt/2','dt/4']
        for ax,key,label in zip(axes.flat[:3],['first_contact_time_s','minimum_pair_gap_m','contact_duration_s'],['first contact (s)','minimum gap (m)','contact duration (s)']):ax.plot(x,[r[key] for r in cases],'o-');ax.set_ylabel(label)
        for ax,key,label in zip(axes.flat[3:],['final_position_max_difference_m','final_velocity_max_difference_m_s','max_short_axis_difference_rad'],['final center difference (m)','final velocity difference (m/s)','short-axis difference (rad)']):ax.bar(['dt vs dt/2','dt/2 vs dt/4'],[r[key] for r in data['comparisons']]);ax.set_ylabel(label)
        fig.suptitle('NOT PRODUCTION TIMESTEP SELECTION',color='#a33b30',fontsize=16)
    fig.tight_layout(rect=(0,.055,1,.96));return fig


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--only',type=int,nargs='+');parser.add_argument('--output',type=Path);args=parser.parse_args()
    output=args.output or REPORT/'figures';output.mkdir(parents=True,exist_ok=True)
    for i in args.only if args.only is not None else range(16):
        fig=plot(i);path=output/(NAMES[i]+'.png');fig.savefig(path,metadata={'Software':'Particle4 deterministic matplotlib evidence'});plt.close(fig)
        if args.output is None:
            source={name+'.json':sha256(DATA/(name+'.json')) for name in SOURCES[i]}
            write_json(DATA/f'{i:02d}_figure_sources.json',dict(figure=path.name,source_sha256=source,manual_visual_review='PENDING_USER_REVIEW',arrows='DISPLAY_SCALE_ONLY',geometry_queries='3D_EXACT_SUPPORT_GEOMETRY'))
            notes=NOTES[i];(REPORT/f'{i:02d}_step_notes.md').write_text(f'![{NAMES[i]}](figures/{path.name})\n\n应该看什么：{notes[0]}\n\n实际看到什么：{notes[1]}\n\n有没有异常：{notes[2]}\n')
        print(path.name,flush=True)


if __name__=='__main__':main()
