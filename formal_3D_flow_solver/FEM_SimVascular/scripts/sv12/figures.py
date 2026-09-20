#!/usr/bin/env python3
"""SV1.2 diagnostics and, only when accepted, unchanged native field views."""
import os
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'
os.environ['MPLCONFIGDIR'] = str(ROOT/'outputs/sv1_2/plot_cache')
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sv_validation.sv12 import REPORT, load
from sv_validation.provenance import sha256, write_json

FONT = FontProperties(fname='/mnt/c/Windows/Fonts/msyh.ttc')
plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 11, 'figure.facecolor': 'white'})
flow, qc = load('flow_execution'), load('saved_state_qc')
history, resources = load('solver_history'), load('solver_resource_usage')
accepted = flow['accepted_solution_available'] and load('solution_reload')['status'] == 'PASS'
note = ('数值验收通过，待人工审核' if accepted else '未获得 accepted steady solution，仅作诊断')
figures = []

def save(fig, name, title, caption):
    fig.suptitle(title, fontproperties=FONT, fontsize=19, y=.97)
    fig.text(.5, .02, caption, fontproperties=FONT, fontsize=10, ha='center', va='bottom')
    fig.tight_layout(rect=(.015, .095, .985, .92))
    path = REPORT/name
    fig.savefig(path, dpi=180)
    plt.close(fig)
    figures.append({'path': str(path.relative_to(ROOT)), 'sha256': sha256(path),
                    'kind': 'accepted native field' if name in FORMAL else 'run diagnostic'})

FORMAL = {'velocity_global.png', 'velocity_slices.png', 'pressure_global.png',
          'pressure_sections.png', 'flux_balance.png', 'outlet_flow_split.png'}
states, intervals = qc['states'], qc['intervals']
steps = [s['step'] for s in states]
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
axes[0].plot(steps, [s['Q_in_m3_s']/s['Q_target_m3_s'] for s in states], 'o-', label='Qin / Qtarget')
axes[0].plot(steps, [s['Q_out_total_m3_s']/s['Q_target_m3_s'] for s in states], '.-', label='Qout total / Qtarget')
axes[0].set_ylabel('Normalized flow')
axes[0].legend()
axes[1].plot(steps, [s['epsilon_mass'] for s in states], 'o-', label='Mass error')
axes[1].plot(steps, [s['epsilon_Q'] for s in states], '.-', label='Inlet error')
axes[1].set_yscale('symlog', linthresh=1e-16)
axes[1].axhline(1e-6, ls='--', color='#ad3c44', label='Acceptance: 1e-6')
axes[1].set_ylabel('Relative error')
axes[1].legend()
for ax in axes:
    ax.set_xlabel('Saved timestep'); ax.grid(alpha=.2)
save(fig, 'mass_balance_over_time.png', '流入和流出什么时候开始真正守恒？',
     note+'；每个保存状态独立积分，step 10 是 SV1.1 的续算起点；早期瞬态误差不用于求解器排名')

fig, ax = plt.subplots(figsize=(12, 6))
linear = history['linear']
for iteration in sorted({r['nonlinear_iteration'] for r in linear}):
    subset = [r for r in linear if r['nonlinear_iteration'] == iteration]
    ax.plot([r['step'] for r in subset], [r['linear_iterations'] for r in subset], '.',
            ms=3, label=f'Nonlinear iteration {iteration}')
ax.axvspan(0, 10, color='#cbd5e1', alpha=.4, label='Inherited SV1.1')
ax.set_xlabel('Timestep'); ax.set_ylabel('PETSc KSP iterations'); ax.grid(alpha=.2)
ax.legend(ncol=2, fontsize=9)
save(fig, 'linear_iterations_over_time.png', '完整计算过程中线性求解是否一直稳定？',
     f"{note}；新增求解失败 {history['linear_failures']} 次；全部 KSP reason 和 true residual 保留在原始日志与逐次记录")

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
nonlinear = history['nonlinear']
axes[0].plot([r['step'] for r in nonlinear], [r['iterations'] for r in nonlinear], '.-')
axes[0].axhline(12, color='#ad3c44', ls='--', label='Maximum allowed: 12')
axes[0].set_ylabel('Outer nonlinear iterations'); axes[0].legend()
for key, label in [('initial_residual', 'Initial scaled residual'), ('final_residual', 'Final scaled residual')]:
    values = [r for r in nonlinear if r[key] is not None]
    axes[1].plot([r['step'] for r in values], [r[key] for r in values], '.-', label=label)
axes[1].set_yscale('symlog', linthresh=1e-26)
axes[1].set_ylabel('Scaled nonlinear residual norm'); axes[1].legend(fontsize=9)
ratios = [r for r in nonlinear if r['Ri_over_R0'] is not None and r['Ri_over_R1'] is not None]
axes[2].plot([r['step'] for r in ratios], [min(r['Ri_over_R0'], r['Ri_over_R1']) for r in ratios], '.-', label='min(Ri/R0, Ri/R1)')
axes[2].set_yscale('symlog', linthresh=1e-16)
axes[2].axhline(1e-10, color='#ad3c44', ls='--', label='Frozen tolerance: 1e-10')
axes[2].set_ylabel('Final relative nonlinear residual'); axes[2].legend(fontsize=9)
for ax in axes:
    ax.set_xlabel('Timestep'); ax.grid(alpha=.2)
save(fig, 'nonlinear_convergence.png', '外层流体方程是否一直正常收敛？',
     f"{note}；新增非线性失败 {history['nonlinear_failures']} 步；依据实际相对残差门限验收，不能只看迭代次数")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, key, threshold in zip(axes, ('E_u', 'E_Q'), (1e-5, 1e-6)):
    ax.plot([r['step'] for r in intervals], [r[key] for r in intervals], '.-')
    ax.set_yscale('symlog', linthresh=1e-16)
    ax.axhline(threshold, color='#ad3c44', ls='--', label=f'Threshold: {threshold:g}')
    if len(intervals) >= 5:
        ax.axvspan(intervals[-5]['previous_step'], intervals[-1]['step'], alpha=.12, color='#368d80', label='Last five saved intervals')
    ax.set_xlabel('Saved timestep'); ax.set_ylabel(key); ax.grid(alpha=.2); ax.legend(fontsize=9)
save(fig, 'steady_convergence.png', '流场什么时候达到稳定？',
     f"{note}；连续共同通过 {qc['consecutive_passing_intervals']} 个保存区间；E_u 使用体积 L2，E_Q 以 Qtarget 归一化")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
measurements = resources['measurements']
labels = [r['name'] for r in measurements]
positions = np.arange(len(labels))
axes[0].barh(positions-.15, [r['solver_wall_s'] for r in measurements], .3, label='GNU wall time')
axes[0].barh(positions+.15, [r['python_monotonic_s'] for r in measurements], .3, label='Python monotonic')
axes[0].set_yticks(positions, labels); axes[0].set_xlabel('Elapsed time (s)'); axes[0].legend()
rss = [r['peak_single_process_RSS_KiB']/1024 for r in measurements]
axes[1].barh(labels, rss, color='#419989'); axes[1].set_xlabel('Peak single-process RSS (MiB)')
for i, value in enumerate(rss): axes[1].text(value, i, f' {value:.1f}', va='center')
axes[1].set_xlim(0, max(rss)*1.2)
save(fig, 'solver_resource_usage.png', '这次完整计算用了多少时间和内存？',
     '仅统计 SV1.2 新增求解；内存为最大单进程 RSS，非 MPI 总和；原生继承计时另存，未计算加速比')

if accepted:
    import pyvista as pv
    pv.OFF_SCREEN = True
    final = load('accepted_solution')
    path = ROOT/final['path']
    assert sha256(path) == final['sha256']
    grid = pv.read(path)
    speed = np.linalg.norm(np.asarray(grid['Velocity']), axis=1)
    final_note = f"SV1.2 accepted step {final['step']}，t={final['time_s']:.6e} s；原生场，未平滑、裁剪或修改压力参考"
    wall = pv.read(ROOT/'outputs/sv1/SV_MESH/mesh-complete.exterior.vtp')
    renderer = []
    def scene_save(scene, name, title, caption, reset=True):
        scene.set_background('white', all_renderers=True)
        if reset:
            scene.view_isometric(); scene.reset_camera()
        frame = scene.screenshot(return_img=True)
        renderer.extend(line.strip() for line in scene.render_window.ReportCapabilities().splitlines() if 'renderer' in line.lower())
        scene.close()
        fig, ax = plt.subplots(figsize=(13, 8))
        ax.imshow(frame); ax.axis('off')
        save(fig, name, title, caption)
    scene = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    scene.add_mesh(wall, color='#9fb4c7', opacity=.08)
    cloud = pv.PolyData(grid.points.copy())
    cloud['Speed'] = speed
    scene.add_mesh(cloud, scalars='Speed', cmap='viridis', point_size=3, opacity=.10,
                   render_points_as_spheres=False, scalar_bar_args={'title': 'Speed (m/s)', 'fmt': '%.2e'})
    scene.add_axes()
    scene_save(scene, 'velocity_global.png', '稳态血管中的速度分布',
               final_note+'\n全部原生节点均参与显示，包括零速壁面；统一显示透明度 0.10 用于观察内部，色标为实际速度模长')
    section = load('pressure_sections')
    axis, center = np.array(section['axis']), np.array(section['origin_m'])
    lower, upper = section['projection_limits_m']
    cuts = []
    for fraction in (.2, .5, .8):
        cut = grid.slice(normal=axis, origin=center+(lower+fraction*(upper-lower))*axis)
        # First interpolate the P1 velocity vector, then compute its magnitude.
        cut['Speed'] = np.linalg.norm(np.asarray(cut['Velocity']), axis=1)
        cuts.append(cut)
    common_max = max(float(np.max(cut['Speed'])) for cut in cuts)
    scene = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1800, 700))
    for i, (cut, fraction) in enumerate(zip(cuts, (.2, .5, .8))):
        scene.subplot(0, i)
        scene.add_mesh(cut, scalars='Speed', cmap='viridis', clim=(0, common_max),
                       scalar_bar_args={'title': 'Speed (m/s)', 'fmt': '%.1e', 'n_labels': 3})
        scene.add_text(f'S{i+1}: axis fraction {fraction}', font_size=11, color='#334155')
        focus = np.array(cut.center)
        span = max(float(np.ptp(cut.points, axis=0).max()), 1e-6)
        scene.camera_position = (focus+axis*3*span, focus, section['view_up'])
        scene.enable_parallel_projection(); scene.reset_camera()
    scene_save(scene, 'velocity_slices.png', '血管内部三个截面的稳态速度',
               final_note+'\n几何主轴 20%、50%、80% 截面；速度向量线性插值后取模，三图使用同一色标', reset=False)
    scene = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    scene.add_mesh(grid.extract_surface(algorithm='dataset_surface'), scalars='Pressure', cmap='coolwarm',
                   scalar_bar_args={'title': 'Pressure (Pa)', 'fmt': '%.2e'})
    scene.add_axes()
    scene_save(scene, 'pressure_global.png', '稳态血管中的压力分布', final_note)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    rows = section['sections']
    axes[0].plot([r['offset_m']*1e6 for r in rows], [r['mean_pressure_pa'] for r in rows], 'o-')
    axes[0].set_xlabel('Principal-axis plane offset (um)'); axes[0].set_ylabel('Area-average pressure (Pa)')
    ports = ['INLET', 'OUTLET_01', 'OUTLET_02', 'OUTLET_03']
    values = [section['port_area_averages_pa'][r] for r in ports]
    axes[1].bar(ports, values, color=['#3a74ad', '#d99451', '#4caa8b', '#966db2'])
    axes[1].set_ylabel('Port area-average pressure (Pa)'); axes[1].tick_params(axis='x', rotation=20)
    for i, value in enumerate(values): axes[1].annotate(f'{value:.5g}', (i, value), ha='center', xytext=(0, 5), textcoords='offset points')
    for ax in axes: ax.grid(axis='y', alpha=.2)
    save(fig, 'pressure_sections.png', '截面与端口的面积平均压力',
         final_note+'\n主轴平面可能跨多个分支，不代表单根中心线压降；端口平均使用完整三角面面积权重')
    fig, ax = plt.subplots(figsize=(10, 6))
    values = [final['Q_in_m3_s'], final['Q_out_total_m3_s']]
    ax.bar(['Inlet', 'Total outlets'], np.array(values)/final['Q_target_m3_s'], color=['#3a74ad', '#d99451'])
    ax.axhline(1, ls='--', color='#334155', label='Qtarget')
    ax.set_ylim(0, 1.25*max(values)/final['Q_target_m3_s']); ax.set_ylabel('Flow / Qtarget'); ax.legend()
    for i, value in enumerate(values): ax.text(i, value/final['Q_target_m3_s']+.025, f'{value:.8e} m³/s', ha='center')
    save(fig, 'flux_balance.png', '最终稳态流入与流出是否守恒？',
         final_note+f"\nεQ={final['epsilon_Q']:.5e}，εmass={final['epsilon_mass']:.5e}；门限均为 1e-6")
    fig, ax = plt.subplots(figsize=(11, 6))
    ports = ['OUTLET_01', 'OUTLET_02', 'OUTLET_03']
    fractions = [final['outlet_fractions'][r] for r in ports]
    ax.bar(ports, 100*np.array(fractions), color=['#d99451', '#4caa8b', '#966db2'])
    ax.set_ylim(0, 130*max(fractions)); ax.set_ylabel('Fraction of total outflow (%)')
    for i, (port, fraction) in enumerate(zip(ports, fractions)):
        ax.text(i, 100*fraction+1, f"{fraction:.4%}\n{final['outlet_flows_m3_s'][port]:.7e} m³/s", ha='center')
    save(fig, 'outlet_flow_split.png', '三个出口最终各流出多少？', final_note+'\n分流比例来自 accepted 最终原生场的实际端面积分')
    assert sha256(path) == final['sha256']
else:
    renderer, common_max = [], None
    assert not any((REPORT/name).exists() for name in FORMAL), 'Unaccepted run must not carry formal field figures'

write_json(REPORT/'visuals.json', {
    'accepted_solution_available': accepted, 'figures': figures,
    'field_source': load('accepted_solution')['path'] if accepted else None,
    'field_source_sha256': load('accepted_solution')['sha256'] if accepted else None,
    'raw_fields_modified': False, 'velocity_global_all_native_nodes_displayed': accepted,
    'velocity_global_uniform_display_opacity': .10 if accepted else None,
    'pressure_offset_applied_pa': 0, 'slice_velocity_method': 'Interpolate vector before magnitude',
    'slice_common_max_m_s': common_max, 'graphics_renderer': sorted(set(renderer)),
    'software_rendering_requested': True,
    'diagnostic_source_sha256': {name: sha256(REPORT/(name+'.json')) for name in ('saved_state_qc', 'solver_history', 'solver_resource_usage')}})
print(f'SV1.2 figures complete: {len(figures)}; accepted field views: {accepted}', flush=True)
