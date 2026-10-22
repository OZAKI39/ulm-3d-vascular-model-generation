#!/usr/bin/env python3
"""Render measured native fields, including explicit failed-state labels."""
import json
import sys
from scipy.spatial import cKDTree
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'mesh_generate/src'))
from vascular_validation.visuals import pv, np, plt, FONT, save_scene, save_plot, COLORS
from vascular_validation.provenance import sha256, write_json
from vascular_validation.execution import read_resource_log
R = ROOT/'reports/mesh_and_flow'
def read(name): return json.loads((R/(name+'.json')).read_text())
q = read('flow_qc')
grid = pv.read(ROOT/q['path'])
grid['Speed'] = np.linalg.norm(grid['Velocity'], axis=1)
wall = pv.read(ROOT/'outputs/mesh_and_flow/solver_mesh/mesh-complete.exterior.vtp')
failed = q['status'] != 'PASS'
note = ('未验收的瞬态状态；' if failed else '') + '第 {} 步，t={:.4e} s；数值参考工况'.format(q['step'], q['time_s'])
scene = pv.Plotter(off_screen=True, window_size=(1440,1080))
scene.add_mesh(wall, color='#a6b8c4', opacity=.08)
wall_only = pv.read(ROOT/'outputs/mesh_and_flow/solver_mesh/mesh-surfaces/WALL.vtp')
wall_distance, wall_index = cKDTree(np.asarray(grid.points,float)).query(np.asarray(wall_only.points,float))
assert wall_distance.max() < 1e-14
visible = np.ones(grid.n_points,dtype=bool)
visible[wall_index] = False
cloud = pv.PolyData(grid.points[visible])
cloud['Speed'] = grid['Speed'][visible]
scene.add_mesh(cloud, scalars='Speed', cmap='viridis', point_size=3,
               scalar_bar_args={'title':'Speed (m/s)'}, render_points_as_spheres=False)
scene.add_axes()
capabilities = save_scene(scene, 'velocity_global.png', '血管中的速度分布', note+'；隐藏零速壁面节点，显示内部及端口真实节点速度')

points = np.asarray(grid.points, float)
center = points.mean(axis=0)
_,_,vectors = np.linalg.svd(points-center, full_matrices=False)
axis = vectors[0]
if axis[np.argmax(abs(axis))] < 0: axis = -axis
projection = (points-center)@axis
limits = [float(projection.min()),float(projection.max())]
locations = np.linspace(*limits,14)[1:-1]
section_rows=[]
for distance in locations:
    cut = grid.slice(normal=axis, origin=center+distance*axis).triangulate()
    if cut.n_cells == 0: continue
    tri = cut.faces.reshape(-1,4)[:,1:]
    xyz = np.asarray(cut.points,float)[tri]
    area = .5*np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)
    mean = float(np.dot(area,np.asarray(cut['Pressure'])[tri].mean(axis=1))/area.sum())
    section_rows.append({'offset_m':float(distance),'area_m2':float(area.sum()),'mean_pressure_pa':mean})
cuts = [grid.slice(normal=axis,origin=center+(limits[0]+fraction*(limits[1]-limits[0]))*axis) for fraction in (.2,.5,.8)]
slice_max = max(float(np.max(cut['Speed'])) for cut in cuts)
scene = pv.Plotter(shape=(1,3),off_screen=True,window_size=(1800,800))
scene.set_background('white',all_renderers=True)
for index, (label, fraction) in enumerate(zip(('S1','S2','S3'),(.2,.5,.8))):
    scene.subplot(0,index)
    distance = limits[0]+fraction*(limits[1]-limits[0])
    cut = cuts[index]
    scene.add_mesh(cut,scalars='Speed',cmap='viridis',clim=(0,slice_max),
                   scalar_bar_args={'title':'Speed (m/s)', 'fmt':'%.1e', 'n_labels':3})
    scene.add_text(label+' (axis fraction '+str(fraction)+')',position='upper_left',font_size=12,color='#334155')
    focus=np.asarray(cut.center)
    span=max(float(np.ptp(cut.points,axis=0).max()),1e-6)
    scene.camera_position=(focus+axis*3*span,focus,vectors[1])
    scene.enable_parallel_projection()
    scene.reset_camera()
# Keep each camera perpendicular to its section rather than resetting to the global view.
frame=scene.screenshot(return_img=True)
scene.close()
fig,ax=plt.subplots(figsize=(15,7))
ax.imshow(frame);ax.axis('off')
fig.subplots_adjust(top=.9,bottom=.08,left=0,right=1)
save_plot(fig,'velocity_slices.png','血管内部几个截面的速度',note+'；沿几何主轴 20%、50%、80% 位置，分别放大，色标统一采用三个截面的最大值')
scene = pv.Plotter(off_screen=True,window_size=(1440,1080))
scene.add_mesh(grid.extract_surface(algorithm='dataset_surface'),scalars='Pressure',cmap='coolwarm',
               scalar_bar_args={'title':'Pressure (Pa)'})
scene.add_axes()
save_scene(scene,'pressure_global.png','血管中的压力分布',note+'；出口为自然牵引参考')
fig,ax=plt.subplots(figsize=(10,6))
ax.plot([a['offset_m'] for a in section_rows],[a['mean_pressure_pa'] for a in section_rows],'o-',color='#3b6f9b')
ax.set_xlabel('沿几何主轴的平面偏移 / m',fontproperties=FONT)
ax.set_ylabel('截面面积加权平均压力 / Pa',fontproperties=FONT)
ax.grid(alpha=.2)
fig.tight_layout(rect=(0,.09,1,.92))
save_plot(fig,'pressure_sections.png','沿血管不同位置的平均压力',note+'；平面可跨多个分支，不代表单根血管中心线压降')

fig,ax=plt.subplots(figsize=(10,6))
values=[q['Q_in_m3_s'],q['Q_out_total_m3_s']]
ax.bar(['Inlet','Total outlet'],np.array(values)/q['Q_target_m3_s'],color=['#146fcb','#e98624'],width=.55)
ax.axhline(1,color='#45566c',ls='--',label='Target inflow')
for i,value in enumerate(values): ax.text(i,value/q['Q_target_m3_s']+.025,f'{value:.7e} m³/s',ha='center')
ax.set_ylabel('Flow / Q_target');ax.set_ylim(0,max(values)/q['Q_target_m3_s']*1.25);ax.legend(frameon=False)
fig.tight_layout(rect=(0,.09,1,.92))
save_plot(fig,'flux_balance.png','流入和流出是否守恒？',note+f"；相对质量差={q['epsilon_mass']:.6g}，门限=1e-6")
fig,ax=plt.subplots(figsize=(10,6))
roles=['OUTLET_01','OUTLET_02','OUTLET_03']
fractions=[q['outlet_fractions'][n] for n in roles]
ax.bar(roles,np.array(fractions)*100,color=[COLORS[n] for n in roles],width=.6)
for i,(role,fraction) in enumerate(zip(roles,fractions)):
    ax.text(i,fraction*100+1,f'{fraction:.2%}\n{q["outlet_flows_m3_s"][role]:.5e} m³/s',ha='center')
ax.set_ylabel('出口流量占总流出量 / %',fontproperties=FONT);ax.set_ylim(0,max(fractions)*125)
fig.tight_layout(rect=(0,.09,1,.92))
save_plot(fig,'outlet_flow_split.png','三个出口分别流出多少？',note+'；从实际场积分得出，尚不能作为稳态分流结论')

steady=read('steady_state')
table=ROOT/'outputs/mesh_and_flow/vascular_flow/4-procs/B_NS_Velocity_flux.txt'
rows=[list(map(float,line.split())) for line in table.read_text().splitlines() if line.split() and line.split()[0].isdigit()]
a=np.array(rows)
fig,axes=plt.subplots(1,2,figsize=(13,6))
if steady['velocity_errors']:
    axes[0].semilogy(np.arange(len(steady['velocity_errors']))+1,steady['velocity_errors'],'o-')
else:
    axes[0].text(.5,.5,'仅保存一个实际速度场\n无法计算连续保存状态的速度变化',fontproperties=FONT,ha='center',va='center',transform=axes[0].transAxes)
axes[0].axhline(1e-5,color='#b54343',ls='--',label='Required <= 1e-5')
axes[0].set_xlabel('连续保存区间',fontproperties=FONT);axes[0].set_ylabel('Velocity relative change');axes[0].legend()
if len(a)>1:
    errors=np.max(np.abs(np.diff(a[:,2:],axis=0)),axis=1)/q['Q_target_m3_s']
    axes[1].semilogy(a[1:,0],errors,'o-',label='Native per-step flow change')
axes[1].axhline(1e-6,color='#b54343',ls='--',label='Required <= 1e-6')
axes[1].set_xlabel('时间步',fontproperties=FONT);axes[1].set_ylabel('Max boundary flow change / Q_target');axes[1].legend(fontsize=9)
for ax in axes: ax.grid(alpha=.2)
fig.tight_layout(rect=(0,.09,1,.92))
save_plot(fig,'steady_convergence.png','计算是否达到稳态？','未达到验收要求；右图使用逐步边界积分日志，不替代连续保存场的联合稳态判据')

native=read('native_solver');run=read('flow_execution')['runs'][-1]
measures=[read('native_dependencies')['VTK']['build'],native['build'],read('official_smoke'),run]
labels=['VTK dependency','Native solver build','Official smoke','Vascular solve']
resources = [{**read_resource_log(ROOT/m['resource_log']), 'label':label,
              'python_monotonic_elapsed_s':m['elapsed_s'], 'raw_log':m['resource_log']}
             for label,m in zip(labels,measures)]
write_json(R/'solver_resource_usage.json', {'measurements':resources,
           'memory_scope':'Maximum individual process RSS, not the sum of MPI ranks',
           'clock_note':'GNU wall-clock elapsed and Python monotonic elapsed are both preserved; no performance inference from their discrepancy'})
fig,axes=plt.subplots(1,2,figsize=(13,6))
for ax,key,factor,label in zip(axes,['wall_elapsed_s','peak_rss_kib'],[1,1/1024],['Wall-clock elapsed (s)','Peak process RSS (MiB)']):
    data=[m[key]*factor for m in resources]
    ax.barh(labels,data,color=['#bcc9d5','#7d9db6','#52a58d','#ce8255'])
    ax.invert_yaxis();ax.set_xlabel(label)
    for i,value in enumerate(data): ax.text(value,i,f' {value:.1f}',va='center')
    ax.set_xlim(0,max(data)*1.22)
fig.tight_layout(rect=(0,.10,1,.92))
save_plot(fig,'solver_resource_usage.png','这次计算用了多少时间和内存？','内存为 GNU time 记录的最大单进程常驻量，并非所有并行进程之和；真实血管计算因数值失败停止')
write_json(R/'visual_provenance.json',{'flow_result':q['path'],'flow_result_sha256':q['sha256'],
    'artifact_kind':q['artifact_kind'],'graphics_force_software':True,'global_velocity_wall_nodes_hidden':True,'velocity_slice_fractions':[.2,.5,.8],'velocity_slice_common_max_m_s':slice_max,
    'graphics_renderer':[line.strip() for line in capabilities.splitlines() if 'renderer' in line.lower()],
    'pressure_sections':{'method':'Area-weighted P1 pressure on geometric principal-axis planes; may cross multiple branches','axis':axis.tolist(),'origin_m':center.tolist(),'sections':section_rows},
    'native_flow_table_sha256':sha256(table),
    'figures':{p.name:sha256(p) for p in sorted(R.glob('*.png'))}})
print('Actual field figures complete; acceptance:',q['status'])
