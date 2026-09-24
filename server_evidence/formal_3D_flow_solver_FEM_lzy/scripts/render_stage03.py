#!/usr/bin/env python3
"""WSL review artifacts for the measured singular Stage 3 attempt; no fabricated fields."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pyvista as pv
from fem3d.audit import write_json,sha256,timestamp
from fem3d.vascular import load_config
R=ROOT/'reports/stage03';B=ROOT/'outputs/stage03/reference'
read=lambda p:json.loads(p.read_text())
config,_=load_config(ROOT);fail=read(B/'metadata/failure.json');assert fail['status']=='FAIL'
diag=read(B/'qc/singularity_diagnosis.json');a=read(ROOT/'outputs/stage03/preflight/assembly_verified.json')
r=read(B/'metadata/resource_accounting.json');raw=np.load(ROOT/config['mesh']['source']/'mesh/volume_mesh.npz')
ports=read(ROOT/config['mesh']['source']/'planar_port_contract_v2.json')['ports']
font=Path('/mnt/c/Windows/Fonts/msyh.ttc');assert font.is_file()
font_manager.fontManager.addfont(str(font));family=font_manager.FontProperties(fname=font).get_name()
plt.rcParams.update({'font.family':family,'axes.unicode_minus':False,'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
ink='#26384b';muted='#607184';red='#aa354c';teal='#168b80';orange='#d28a38';figures={}
points=raw['points_m']*1e6;tri=raw['boundary_triangles'];tags=raw['facet_tags']
surf=pv.PolyData(points,np.c_[np.full(len(tri),3),tri].ravel())
wall=surf.extract_cells(tags==1)
tetra=pv.UnstructuredGrid(np.c_[np.full(len(raw['tetra']),4),raw['tetra']].ravel(),np.full(len(raw['tetra']),pv.CellType.TETRA,dtype=np.uint8),points)

def screenshot(kind):
    p=pv.Plotter(off_screen=True,window_size=(1550,1050),border=False);p.set_background('white')
    if kind=='geometry':
        p.add_mesh(wall,color='#b4c6cd',smooth_shading=False)
        anchors=[];labels=[]
        for name,port in ports.items():
            color=orange if name=='inlet' else teal
            p.add_mesh(surf.extract_cells(tags==port['entity_id']),color=color,lighting=False)
            center=np.array(port['plane_origin_m'])*1e6;normal=np.array(port['outward_normal'])
            if name=='inlet':p.add_mesh(pv.Arrow(start=center+normal*10,direction=-normal,scale=9),color=color)
            else:p.add_mesh(pv.Arrow(start=center,direction=normal,scale=9),color=color)
            anchors.append(center+normal*14);labels.append(name.upper()+' | '+('fixed Q' if name=='inlet' else 'zero traction'))
        p.add_point_labels(np.array(anchors),labels,font_size=19,text_color=ink,shape_color='white',shape_opacity=.95,show_points=False,always_visible=True)
    else:
        p.add_mesh(wall,color='#b4c6cd',opacity=.17)
        centers=np.array([row['centroid_m'] for row in read(ROOT/'inputs/stage03/residual_source_evidence.json')['records']])*1e6
        ids=[]
        from scipy.spatial import cKDTree
        dist,ids=cKDTree(points[raw['tetra']].mean(axis=1)).query(centers);assert max(dist)<1e-7
        p.add_mesh(tetra.extract_cells(ids),color=red,show_edges=True,line_width=2)
        p.add_points(centers,color=red,point_size=13,render_points_as_spheres=True)
        labels=['cells 1 + 2: two null pressure modes','cell 3: no flow samples']
        anchor=np.array([centers[:2].mean(axis=0),centers[2]])
        p.add_point_labels(anchor,labels,font_size=19,text_color=ink,shape_color='white',shape_opacity=.95,show_points=False,always_visible=True)
    p.view_isometric();p.enable_parallel_projection();p.reset_camera();p.camera.zoom(.88)
    p.add_axes(xlabel='x',ylabel='y',zlabel='z',line_width=2)
    image=p.screenshot(return_img=True);p.close();return image

def save(fig,name,role,kind):
    fig.savefig(R/name,dpi=150,bbox_inches='tight');plt.close(fig)
    figures[name]={'sha256':sha256(R/name),'description':role,'kind':kind,'valid_solution_field':False}

def footer(fig,text):
    fig.text(.055,.045,text,color=muted,fontsize=11)
    fig.text(.055,.013,'参考数值条件 · 非实验泵流量  |  STAGE 3: FAIL — SINGULAR_PRESSURE_SUPPORT',color=red,fontsize=10)

fig=plt.figure(figsize=(13,9));ax=fig.add_axes([.02,.16,.96,.75]);ax.imshow(screenshot('geometry'));ax.axis('off')
fig.suptitle('真实血管：从一个入口泵入，三个出口自然流出',fontsize=19,y=.97)
footer(fig,'灰色壁面：速度为零。橙色入口：仅约束总流量。绿色出口：σn = 0。\n箭头用于标识端口与边界方向，不代表已求得的速度；出口箭头不表示非零牵引。')
save(fig,'real_geometry_and_bc.png','Actual unchanged Stage 1.7 geometry, all four labeled ports; BC annotations only','MEASURED_GEOMETRY')

fig=plt.figure(figsize=(13,9));ax=fig.add_axes([.02,.18,.96,.73]);ax.imshow(screenshot('residual'));ax.axis('off')
fig.suptitle('三个低质量单元在哪里？',fontsize=20,y=.97)
footer(fig,'3 个物理质心已在当前 DOLFINx 网格中重新定位，位移均为 0 m。\n相邻的单元 1、2 支撑了两个已证实的压力零模式；第三个位置仅完成定位。')
save(fig,'residual_cell_locations.png','Actual three residual tetra locations; markers enlarged for visibility, geometry unmodified','MEASURED_GEOMETRY')

fig,axes=plt.subplots(1,3,figsize=(14,7),gridspec_kw={'width_ratios':[1.1,1,1]})
axes[0].bar(['仅组装','正式尝试\n分解失败'],[a['assembly_wall_time_s'],r['elapsed_time_s']],color=[teal,red])
for i,v in enumerate([a['assembly_wall_time_s'],r['elapsed_time_s']]):axes[0].text(i,v+.35,f'{v:.2f} s',ha='center')
axes[0].set(ylabel='耗时 / 秒',ylim=(0,r['elapsed_time_s']*1.25),title='真实记录的时间')
peak=r['valid_peak_rss_kib']/1024**2
axes[1].bar(['组装\n最大 rank','正式尝试\n最大子进程'],[a['max_rank_peak_rss_kib']/1024**2,peak],color=[teal,red])
axes[1].set(ylabel='峰值 RSS / GiB',ylim=(0,peak*1.3),title='不是四个 rank 的同时总和')
axes[1].text(1,peak+.07,f'{peak:.3f} GiB',ha='center')
axes[2].axis('off');axes[2].text(.05,.96,f"{a['total_dofs']:,}\n总未知量",fontsize=26,color=ink,va='top')
axes[2].text(.05,.63,f"速度 {a['velocity_dofs']:,}\n压力 {a['pressure_dofs']:,}\n全局约束 {a['real_dofs']}\nnnz {a['nnz']:,}",fontsize=13,linespacing=1.7,va='top')
axes[2].text(.05,.10,'4 MPI ranks · 单线程\nGPU 未参与\n没有 OOM / 内存终止',fontsize=13,color=muted)
fig.suptitle('真实计算用了多少时间和内存？',fontsize=20);fig.subplots_adjust(top=.80,bottom=.29,wspace=.35)
footer(fig,f"{r['elapsed_time_s']:.2f} s 是失败尝试的整体耗时，不是成功求解时间。RSS 来自 Linux 子进程资源记录。\n进程组采样漏掉 MPICH worker，已排除该无效总量；详见资源计量说明。")
save(fig,'solver_resource_usage.png','Measured assembly and failed-attempt time; valid individual-child peak RSS with aggregation limitation','MEASURED_RESOURCES')

records=read(ROOT/'inputs/stage03/residual_source_evidence.json')['records']
fig,ax=plt.subplots(figsize=(14,8));ax.axis('off')
fig.suptitle('剩下的三个低质量网格附近，流场有没有异常？',fontsize=20,y=.94)
ax.text(.5,.86,'当前无法做流场对比：矩阵奇异，未得到有效解',ha='center',color=red,fontsize=19,transform=ax.transAxes)
rows=[]
for i,row in enumerate(records):
    locked=any(np.linalg.norm(np.array(row['centroid_m'])-r['centroid_m'])<1e-15 for r in diag['topological_evidence']['wall_velocity_fully_constrained_cells'])
    rows.append([f'单元 {i+1}',f"{row['min_sicn']:.6f}",'全部被壁面固定' if locked else '存在自由节点','未获得','未获得','未获得'])
table=ax.table(cellText=rows,colLabels=['物理定位','minSICN','P2 速度节点','速度','压力','应变率'],loc='center',cellLoc='center',colWidths=[.13,.13,.27,.14,.14,.14]);table.auto_set_font_size(False);table.set_fontsize(13);table.scale(1,2.6)
ax.text(.5,.19,'两个压力基函数仅支撑在前两个单元上；实际矩阵验证：两条零行，且 A·e1 = A·e2 = 0。',ha='center',fontsize=13,transform=ax.transAxes,color=ink)
footer(fig,'没有生成 residual / neighbor median 比值，也没有进行平滑或添加压力固定点。\n本图呈现失败证据，不能作为局部流场有限性或离散误差通过的证明。')
save(fig,'residual_cell_flow_check.png','Topology and matrix failure evidence for residual cells; requested flow/neighbor monitor unavailable','DIAGNOSIS_FLOW_UNAVAILABLE')

cards={
'velocity_global.png':('血管中的速度有多大？','P2 速度场未生成，无法绘制真实三维速度分布。'),
'velocity_slices.png':('内部切面上的速度是什么样？','没有有效速度系数，不能生成切面或用假定剖面代替。'),
'pressure_global.png':('血管中的压力如何分布？','P1 压力系统存在两个局部零模式，没有可验收的压力场。'),
'pressure_sections.png':('各个截面的平均压力是多少？','入口、三个出口及内部截面的压力积分均未计算。'),
'flux_balance.png':('流量是否守恒？',f"目标参考流量 Q = {config['physics']['inlet_volume_flow_m3_s']:.12e} m³/s\n实际入口、三个出口、总出口与闭合误差：均未评估。"),
'outlet_flow_split.png':('三个出口各分到多少流量？','没有实际 FEM 流量积分；Q1、Q2、Q3 及比例均未知。'),
'velocity_gradient_slice.png':('速度变化最剧烈的区域在哪里？','未生成速度梯度样本。没有孤立峰值判断，也没有计算 WSS。')}
for name,(title,message) in cards.items():
    fig,ax=plt.subplots(figsize=(12,7));ax.axis('off')
    fig.suptitle(title,fontsize=23,y=.90)
    ax.text(.5,.68,'未生成流场图',color=red,fontsize=30,ha='center',transform=ax.transAxes)
    ax.text(.5,.47,message,ha='center',va='center',fontsize=15,linespacing=1.8,transform=ax.transAxes)
    ax.text(.5,.20,'原因：MUMPS 数值分解失败；实际矩阵已证实有两个独立压力零模式。',ha='center',fontsize=13,color=muted,transform=ax.transAxes)
    footer(fig,'这是缺失结果状态页，不是流场、零值数据或通过验收的证据。')
    save(fig,name,message,'UNAVAILABLE_STATUS_PAGE')
write_json(R/'visualization_manifest.json',{'timestamp':timestamp(),'generated_on':'WSL','stage_status':'FAIL',
    'figures':figures,'actual_field_plots':0,'unavailable_status_pages':len(cards),
    'source_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in [ROOT/config['mesh']['source']/'mesh/volume_mesh.npz',B/'qc/singularity_diagnosis.json',B/'metadata/resources.json',B/'metadata/resource_accounting.json',ROOT/'configs/stage03_reference_vascular.yaml']},
    'font_source':str(font),'reference_numerical_condition':True,'experimental':False})
print('11 WSL artifacts: 3 geometry/resource figures, 1 diagnostic table, 7 explicitly unavailable flow pages.')
