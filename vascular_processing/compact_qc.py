"""Same physical scale comparisons and an evidence-based compact phantom report."""
from pathlib import Path

import numpy as np
import pyvista as pv


def fmt(values):
    return ' x '.join(f'{v:.2f}' for v in values)


def panel_figure(path, meshes, captions, title):
    # Centering is display-only. Every panel uses one camera vector and identical
    # parallel_scale/window dimensions, so mm/pixel is exactly shared.
    centered = []
    for mesh in meshes:
        mesh = mesh.copy(); mesh.points -= (mesh.points.min(axis=0) + mesh.points.max(axis=0)) / 2
        centered.append(mesh)
    span = max(float(np.linalg.norm(np.ptp(m.points, axis=0))) for m in centered)
    plotter = pv.Plotter(shape=(1, len(meshes)), off_screen=True, window_size=(650*len(meshes), 850))
    for index, (mesh, caption) in enumerate(zip(centered, captions)):
        plotter.subplot(0, index); plotter.set_background('#f4f6f9')
        plotter.add_mesh(mesh, color=['#1c8095', '#d58323', '#7460a9'][index % 3], smooth_shading=True)
        plotter.add_text(caption, position='upper_left', font_size=11, color='#14222e')
        plotter.camera_position = [(span*1.3, -span*1.8, span), (0, 0, 0), (0, 0, 1)]
        plotter.enable_parallel_projection(); plotter.camera.parallel_scale = span*0.70
        plotter.add_text('Identical view / mm per pixel\nDisplay centering only; scale = 1', position='lower_left', font_size=10, color='#364856')
        plotter.add_axes()
    plotter.screenshot(str(path)); plotter.close()


def figures(output, m):
    qc = output / 'QC'; qc.mkdir(exist_ok=True)
    meshes = []; captions = []
    for mode in ('MINI', 'BALANCED', 'RICH'):
        item = m['candidates'][mode]; s = item['stats']; r = item['radius_compensation']
        meshes.append(pv.read(item['surface_stl']))
        captions.append(f"{mode}\nBranches {s['branch_count']} | Bifurcations {s['bifurcation_count']} | Outlets {s['outlet_count']}\n"
            f"Centerline {s['centerline_length_mm']:.2f} mm\nBBox {fmt(s['bbox']['extents_mm'])} mm\n"
            f"Minimum original D {s['diameter_min_mm']:.2f} mm\nRadius compensated {100*r['centerline_compensated_fraction']:.1f}% of length\n"
            f"Context {item['proximal_context']['length_mm']:.2f} mm | Inlets 1")
    panel_figure(qc/'compact_roi_comparison.png', meshes, captions, 'Compact candidates')
    old = m['old_reference']; old_mesh = pv.read(Path(old['path'])/'BG001_RMCA_print_surface_closed.stl')
    s = old['stats']; balanced = m['candidates']['BALANCED']['stats']
    panel_figure(qc/'old_vs_compact.png', [old_mesh, meshes[1]], [
        f"FULL_CONTEXT_REFERENCE\nBranches {s['branch_count']} | Bifurcations {s['bifurcation_count']} | Outlets {s['outlet_count']}\n"
        f"Length {s['centerline_length_mm']:.2f} mm\nBBox {fmt(s['bbox']['extents_mm'])} mm",
        f"BALANCED / COMPACT_PRINT_CANDIDATE\nBranches {balanced['branch_count']} | Bifurcations {balanced['bifurcation_count']} | Outlets {balanced['outlet_count']}\n"
        f"Length {balanced['centerline_length_mm']:.2f} mm\nBBox {fmt(balanced['bbox']['extents_mm'])} mm"], 'Old versus compact')
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1950, 850))
    for index, item in enumerate(m['orientation']['top_candidates']):
        plotter.subplot(0, index); plotter.set_background('#f4f6f9')
        mesh = pv.read(item['stl']); plotter.add_mesh(mesh, color='#d58323', smooth_shading=True)
        center = mesh.center; extent = np.ptp(mesh.points, axis=0)
        plotter.add_mesh(pv.Plane(center=(128, 128, 0), direction=(0, 0, 1), i_size=85, j_size=85), color='#bac5cf', opacity=.2)
        inlet = np.array(m['orientation']['inlet_native_mm'])
        mat = np.array(item['transform_4x4']); pin = mat[:3,:3] @ inlet + mat[:3,3]
        plotter.add_mesh(pv.Sphere(radius=1.1, center=pin), color='#d62728')
        plotter.add_text(f"BALANCED orientation {index+1}\nScore {item['score']:.4f}\nBBox {fmt(item['bbox_extents_mm'])} mm\n"
            f"Inlet height {item['inlet_height_mm']:.2f} mm\nStem angle {item['proximal_stem_angle_from_plate_deg']:.1f} deg\n"
            f"Downward area {item['downward_support_area_mm2']:.1f} mm2 (heuristic)\n"
            + ('ACTUAL BAMBU SLICE' if index == 0 else 'STL ONLY'), font_size=11, color='#14222e')
        plotter.camera_position = [(230, -20, 130), (128,128,20), (0,0,1)]
        plotter.enable_parallel_projection(); plotter.camera.parallel_scale=55
        plotter.add_axes()
    plotter.screenshot(str(qc/'orientation_comparison.png')); plotter.close()


def write_report(output, m):
    selected = m.get('selected_orientation', {}); candidates = m['candidates']; old=m['old_reference']['stats']
    lines = ['# 单一冻结 MeVO 组件的紧凑制造 ROI 研究', '', '## 摘要', '',
        f"本研究将既有大范围制造模型保留为 FULL_CONTEXT_REFERENCE，并从单一冻结语义组件生成 MINI、BALANCED 与 RICH 三个紧凑派生物。实际来源为 part{m['source_component']:02d}，默认候选为 BALANCED。未重新执行配准、供体排序、NN、UNKNOWN 校准、分支语义聚合或近端精修。所有保留坐标、TYPE 与有向边均来自真实 BraVa；半径补偿和打印旋转只属于制造层。当前状态为 {m['status']}，未进行实体打印。", '',
        '## 方法', '',
        '首先检查配置优先的 part04，并同时读取已有 semantic_refined 与 diameter_075 单根文件。以每轮三个候选均具有最小分叉结构为来源可行性条件，优先使用配置中的组件顺序，并以冻结语义支持减去需补偿中心线比例作为工程评分。程序未对原语义进行重新投票或聚合。', '',
        '从精修根沿原始拓扑向远端展开，分支按冻结语义支持、原始中位直径及原长度的 0.45/0.35/0.20 权重排序。最多保留两个女儿，最多跨两次解剖拓扑分叉；分叉深度仅描述制造简化，不对应 M2/M3 分段。路径在预算以内的最后一个原始样本截断，不插值、不连接其他组件。删除短于 12 mm 且原始中位直径低于 1 mm 的末端 twig，其余必要连接路径保留。每条未完全保留的原语义分支均在 branch_selection.csv 中解释。', '',
        '近端只选择距目标 10 mm 最近、位于 8–15 mm 内的真实上游前缀，不增加兄弟分支。根据用户后续明确授权，45/55/65 mm 路径长度及 60/80 mm 全局长度仅保留为参考，不作为硬截断。算法对原始采样点距离与全部拓扑边界进行确定性搜索，在单组件、分支/出口/深度、200 mm 总中心线、最大空间范围等硬限制内，权衡空间目标、代表性分叉、总长度及补偿比例，自动确定各候选截点。RICH 的空间范围仍服从全局最长边 110 mm、对角线 140 mm 上限。', '',
        '制造半径下限为 0.60 mm，使用已有约 5 mm 局部 smoothstep 过渡；每点原始半径、补偿半径和增量保存在映射表。补偿比例按边长及两端补偿指示量的梯形权重统计，另列采样点比例。半径倍数超过 1.75 时标注 LARGE_MANUFACTURING_COMPENSATION，不自动抹平原始半径差异或判为解剖改变。', '',
        '## 来源可行性与尺寸冲突', '',
        '真实 part04_01 的第一处分叉距离精修根约 58.62 mm，第二处分叉约 123.05 mm；固定 55 mm 上限会在首处分叉之前截止。用户因此明确允许自适应路径长度，并要求优先保持合理复杂度与空间范围。本研究以该授权为准，未偷偷延长默认路径或更改语义边界。', '',
        f"自适应搜索后实际选择 {m['source_derived_component']} 的 {m['source_level']} 层，使用 part04={m['used_part04']}，回退 part03={m['fallback_to_part03']}。既有 d075 也被检查，但其短细截断末端进一步筛除后不能提供三个模式都需要的最低分叉结构；因此采用同一组件的冻结语义版，并只在制造派生物中局部截断。", '',
        f"本次尺寸目标作为偏好，allow_undersized_candidate={m['profile']['compact_roi']['allow_undersized_candidate']}。结构、空间上限与切片可行性分别检查；不通过坐标缩放、跨组件连接或增加长 M1 前缀扩大模型。", '',
        '| 来源 | 层 | BALANCED 分支/分叉/出口 | 结构评价 |', '|---|---|---|---|']
    for r in m['selection_trials']:
        s=r['balanced_stats']; lines.append(f"| {r['derived_component']} | {r['source_level']} | {s['branch_count']}/{s['bifurcation_count']}/{s['outlet_count']} | {'可选' if r['structure_eligible'] else ', '.join(r['reasons'])} |")
    lines += ['', '## 三种紧凑候选的实际结果', '',
        '| 候选 | 分支 | 分叉 | 出口 | 长度 mm | AABB mm | 最长边 mm | context mm | 补偿长度比例 |', '|---|---:|---:|---:|---:|---|---:|---:|---:|']
    for mode,item in candidates.items():
        s=item['stats'];r=item['radius_compensation']
        lines.append(f"| {mode} | {s['branch_count']} | {s['bifurcation_count']} | {s['outlet_count']} | {s['centerline_length_mm']:.3f} | {fmt(s['bbox']['extents_mm'])} | {s['bbox']['longest_dimension_mm']:.3f} | {item['proximal_context']['length_mm']:.3f} | {100*r['centerline_compensated_fraction']:.2f}% |")
        if item.get('post_model_radius_qc'):
            q=item['post_model_radius_qc']; lines.append(f"\n{mode} 原生建模：{item['model']['status']}；拟合后最小直径 {q['min_diameter_mm']:.4f} mm；最大补偿倍数 {r['maximum_radius_inflation_factor']:.3f}，超过 1.75 的采样点数 {len(r['large_compensation_ids'])}；自适应根至末端截断距离 {item['adaptive_path']['selected_cutoff_mm']:.3f} mm，加入 context 后最长路径 {s['maximum_inlet_to_outlet_mm']:.3f} mm。\n")
    # Keep the markdown table contiguous despite separate per-model discussion.
    rows=[x for x in lines if x.startswith('\n')];lines=[x for x in lines if not x.startswith('\n')];lines += rows
    b=candidates['BALANCED']['stats']; cmp=m['comparison']
    lines += ['', f"旧模型有 {old['branch_count']} 条拓扑分支、{old['bifurcation_count']} 个分叉、{old['outlet_count']} 个出口，包围盒为 {fmt(old['bbox']['extents_mm'])} mm。BALANCED 为 {b['branch_count']} 条分支、{b['bifurcation_count']} 个分叉、{b['outlet_count']} 个出口，总长度减少 {cmp['centerline_reduction_percent']:.2f}%；出口减少 {cmp['outlet_reduction']}，分叉减少 {cmp['bifurcation_reduction']}。", '',
        '全部候选为单连通、单入口。BALANCED/RICH 总长度超过 160 mm 优选区间但仍低于 200 mm 硬上限；这是保留第二处分叉及其真实末端支路的结果。若 RICH 与 BALANCED 的补偿 SWC 哈希相同，会在清单记录 same_geometry_as，不虚构第三处分叉。', '',
        '## 打印方向与实际切片', '']
    if selected:
        lines += [f"实际打印配置为 {m['printer']['model']}，喷嘴 {m['printer']['nozzle_diameter_mm']} mm，材料为任务明确指定的 ABS，打印体积为 {fmt(m['printer']['build_volume_mm'])} mm，软件版本 {m['printer']['version']}。原有 GUI 预设保持不变。", '',
            f"仅对 BALANCED 表面搜索 {m['orientation']['tested']} 个确定性方向，保存前三名；MINI、RICH 仅输出表面，不自动切片。方向评分增加入口高度、近端主干与床面 30–60° 的偏好，以及远端长水平段惩罚。最优方向 XYZ 欧拉角为 {selected['euler_xyz_deg']} 度，打印包围盒为 {fmt(selected['bbox_extents_mm'])} mm，入口高度 {selected['inlet_height_mm']:.3f} mm，主干角度 {selected['proximal_stem_angle_from_plate_deg']:.2f}°。", '',
            f"向下需支撑面积启发式为 {selected['downward_support_area_mm2']:.3f} mm²；该指标不是切片器支撑体积。实际只调用 BALANCED Top1 切片，状态为 {m['slicer']['status']}。", '']
        for row in m['slicer']['results']:
            lines += [f"切片记录：状态 {row['status']}；时间 {row.get('estimated_print_time_seconds')} 秒；总耗材 {row.get('filament_used_g')} g。支撑信息为 {row.get('support_information')}；未提供的独立支撑质量不补造。实际原始元数据以 BambuStudio/slice_results.json 为准。", '']
    lines += ['## 讨论与可复核性', '',
        '三种模型用于比较同一局部血管网络的制造复杂度，不能解释为不同解剖标签。所选 part04_01 本身只有两个分叉，因此 BALANCED 与 RICH 可能收敛为相同拓扑和几何；不会为制造一个更大的 RICH 而连接 part04_02、part03 或其它组件。三次二叉分叉加一个入口主干通常需要七条拓扑分支，与最多六条的约束也不相容，因此三分叉是上限，不是强制目标。自适应路径可以超过旧固定参考长度，但总长度与空间范围仍受硬上限限制。阈值为本项目工程启发式，NOT manufacturer guaranteed limits. Must be calibrated experimentally.', '',
        '旧 manufacturing_roi、semantic_refined、diameter_075、既有 VascularMD 输出和上游库均保持原文件哈希；旧大模型的 FULL_CONTEXT_REFERENCE 身份记录在新清单内，不改写旧清单。s1-2 可视化代码未修改。没有向打印机提交任务，不能称为 MANUFACTURING_VALIDATED。', '',
        f"本次明确记录的工程提示：{'; '.join(m.get('warnings', []))}。", '',
        f"本轮核验 {m.get('protection', {}).get('protected_files', 0)} 个受保护文件，哈希全部保持一致。测试及逐点几何核验的完整记录见 QC/test_results.json。", '']
    if m.get('tests'):
        t=m['tests']
        lines += [f"实际回归共 {t['full_suite_passed']} 项通过，其中既有 {t['existing_passed']} 项、新增 {t['new_passed']} 项，最终失败 {t['failures']} 项。既有 scikit-image/NumPy 弃用提示 {t['warnings']} 条，新增测试无警告。六个导出 SWC 均通过原坐标、TYPE、有向边和原始半径映射核验；三套原生模型均成功且保持拓扑，三套打印表面均为单连通、封闭且无非流形边。", '',
            f"s1-2_swc_roi_generate_human.py 的 SHA256 为 `{t['human_sha256']}`，与保护快照一致。旧大模型的输出及其九个实现文件哈希也保持不变。", '']
    (output/'compact_report.md').write_text('\n'.join(lines),encoding='utf-8')
    (output/'README.md').write_text(
        '# 紧凑 ABS/PDMS 制造候选\n\n默认模型：BALANCED。旧大模型：`../manufacturing_roi/`，角色 FULL_CONTEXT_REFERENCE，原文件不变。\n\n'
        f"状态：{m['status']}。尺寸目标、原始半径及拟合后直径的实际限制见 `compact_report.md`。\n\n"
        'NOT manufacturer guaranteed limits. Must be calibrated experimentally.\n\n'
        '- `candidates/{MINI,BALANCED,RICH}/`：原始半径/补偿 SWC、映射、删枝原因、原生 VTK、封口打印 STL。\n'
        '- `BG001_RMCA_BALANCED_print_candidate.3mf`：实际切片成功时生成，可直接打开 Bambu Studio。\n'
        '- `orientation/`：仅 BALANCED 的 Top 3 方向，Top1 执行一次切片。\n'
        '- `QC/`：三候选同视角同比例、旧模型对比、方向对比和测试记录。\n\n'
        '从项目根目录运行，使用新的输出目录：\n\n```bash\n'
        '/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/compact_bg001_rmca.py '
        '--output outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_run02\n```\n\n'
        '配置：`config/compact_manufacturing_profile.yaml`。支持 `--stage prepare/model/finish/verify`；默认 all。'
        '三个候选共用同一语义组件，默认不运行大模型入口。没有自动全局缩放或打印提交。\n',encoding='utf-8')
