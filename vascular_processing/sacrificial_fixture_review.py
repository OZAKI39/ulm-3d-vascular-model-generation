"""Chinese review figures and report for the independent fixture design stage."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageColor
import pyvista as pv
from scipy.spatial import cKDTree
import trimesh

from .sacrificial_fixture import NORMALS, ROOT, face_axis, enabled_faces, unit

COLORS = dict(vascular_core='#10a7b5', inlet_ports='#23ac62', outlet_ports='#f08c28',
              casting_box='#9099a5', holes='#d7333f')


def finite_json(value):
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite_json(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Path):
        return str(value)
    return value


def save_json(path, data):
    Path(path).write_text(json.dumps(finite_json(data), indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')


def polydata(mesh):
    return pv.PolyData(mesh.vertices, np.column_stack([np.full(len(mesh.faces), 3), mesh.faces]).ravel())


def core_face_labels(core, source, routes):
    """Color labels only; every exported face comes from the actual Boolean core."""
    centers = np.vstack([source.triangles_center] + [r.mesh.triangles_center for r in routes])
    labels = np.concatenate([np.zeros(len(source.faces), dtype=int)] +
        [np.full(len(r.mesh.faces), 1 if r.endpoint.role == 'inlet' else 2, dtype=int) for r in routes])
    return labels[cKDTree(centers).query(core.triangles_center)[1]]


def assembly_export(out, basename, source, routes, box_mesh, core):
    out.mkdir(parents=True, exist_ok=True)
    blocks = pv.MultiBlock()
    if core is not None:
        labels = core_face_labels(core, source, routes)
        for label, key in enumerate(['vascular_core', 'inlet_ports', 'outlet_ports']):
            indices = np.flatnonzero(labels == label)
            blocks[key] = polydata(core.submesh([indices], append=True, repair=False)) if len(indices) else pv.PolyData()
    else:
        blocks['vascular_core'] = polydata(source)
        for role, key in [('inlet', 'inlet_ports'), ('outlet', 'outlet_ports')]:
            meshes = [r.mesh for r in routes if r.endpoint.role == role]
            blocks[key] = polydata(trimesh.util.concatenate(meshes)) if meshes else pv.PolyData()
    blocks['casting_box'] = polydata(box_mesh) if box_mesh is not None else pv.PolyData()
    blocks.save(out / 'assembly_preview.vtm')
    result = dict(vtm=str(out / 'assembly_preview.vtm'), vtm_blocks=list(blocks.keys()),
                  vascular_layers='DISJOINT_PATCHES_OF_BOOLEAN_EXTERNAL_SURFACE' if core is not None else 'UNMERGED_PARTIAL_REVIEW',
                  internal_attachment_caps_displayed=core is None,
                  three_mf=None, three_mf_status='NOT_EXPORTED_INCOMPLETE_DESIGN')
    if core is not None and box_mesh is not None:
        scene = trimesh.Scene()
        # Exactly two physical bodies in 3MF; the diagnostic VTM has four layers.
        scene.add_geometry(core, geom_name='sacrificial_core_with_ports', node_name='sacrificial_core_with_ports')
        scene.add_geometry(box_mesh, geom_name='casting_box', node_name='casting_box')
        path = out / f'{basename}_assembly.3mf'
        try:
            scene.export(path)
            restored = trimesh.load_scene(path)
            if len(restored.geometry) != 2:
                raise ValueError('3MF round trip did not preserve two bodies')
            np.testing.assert_allclose(restored.bounds, scene.bounds, atol=1e-4)
            result.update(three_mf=str(path), three_mf_status='TWO_BODIES_ROUND_TRIP_VERIFIED')
        except (ImportError, ValueError, TypeError, AttributeError) as exc:
            if path.exists():
                path.unlink()
            result['three_mf_status'] = 'VTM_FALLBACK: ' + str(exc)
    return result


def chinese_font(cfg, size):
    candidates = [Path(cfg['visualization']['chinese_font']),
                  Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'),
                  Path('/mnt/c/Windows/Fonts/msyh.ttc')]
    for path in candidates:
        if path.is_file():
            return ImageFont.truetype(str(path), size=size)
    raise ValueError('A CJK font is required for readable Chinese review figures; set visualization.chinese_font')


def render_review(out, source, endpoints, routes, box, box_mesh, core, cfg, status, start_checks=None):
    out.mkdir(parents=True, exist_ok=True)
    titles = [
        ('01_source_vascular_core.png', '原始血管芯：保持已验收的打印姿态', 'source'),
        ('02_endpoints_and_tangents.png', '入口和出口位置：箭头表示向外几何方向', 'endpoints'),
        ('03_box_margin_preview.png', '血管与盒体之间的预留空间', 'margin'),
        ('04_port_face_assignment.png', '端口分配：仅使用 X、Y 四个侧壁', 'assignment'),
        ('05_port_routes_viewA.png', '端口延长路线：视角一', 'routeA'),
        ('05_port_routes_viewB.png', '端口延长路线：视角二', 'routeB'),
        ('06_core_with_ports.png', '血管芯与已通过检查的端口延长', 'core'),
        ('07_casting_box.png', '开放式浇注盒：红色标记为穿墙孔', 'box'),
        ('08_assembly_preview.png', '盒体与血管芯装配：两个独立部件', 'assembly'),
        ('09_top_view.png', '从上方检查端口与血管的间距', 'top'),
        ('10_side_view.png', '从侧面检查底部空间和端口高度', 'side'),
    ]
    conflicts = {key: check for key, check in (start_checks or {}).items() if not check['passed']}
    if conflicts:
        titles.append(('11_endpoint_clearance_conflict.png', '原始出口附近的间距：需单独核查', 'conflict'))
    relaxed = [r for r in routes if r.required_vessel_clearance < r.preferred_vessel_clearance]
    if relaxed:
        titles.append(('12_local_clearance_review.png', '间距放宽后的局部延长：保持血管互不相交', 'local'))
    source_poly = polydata(source)
    if core is not None:
        palette = np.array([ImageColor.getrgb(COLORS[k]) for k in ['vascular_core', 'inlet_ports', 'outlet_ports']])
        actual_core = polydata(core)
        actual_core.cell_data['body_color'] = palette[core_face_labels(core, source, routes)].astype(np.uint8)
    center = box.outer.mean(axis=0)
    scale = float(np.max(box.outer_size))
    font = chinese_font(cfg, cfg['visualization']['title_font_size'])
    small = chinese_font(cfg, cfg['visualization']['legend_font_size'])
    for name, title, mode in titles:
        plotter = pv.Plotter(off_screen=True, window_size=cfg['visualization']['window_size'])
        plotter.set_background('#fafbfc')
        use_union = core is not None and mode in ('core', 'routeA', 'routeB', 'assembly', 'top', 'side', 'local')
        if use_union:
            plotter.add_mesh(actual_core, scalars='body_color', rgb=True, smooth_shading=True)
        elif mode != 'box':
            plotter.add_mesh(source_poly, color=COLORS['vascular_core'], smooth_shading=True)
        if mode not in ('source', 'endpoints', 'core', 'conflict', 'local') and box_mesh is not None:
            plotter.add_mesh(polydata(box_mesh), color=COLORS['casting_box'], opacity=.23 if mode != 'box' else .60,
                             smooth_shading=False)
        if mode == 'margin':
            b = box.inner
            plotter.add_mesh(pv.Box(bounds=tuple(b[:, k].tolist()[j] for k in range(3) for j in (0, 1))).outline(),
                             color=COLORS['casting_box'], line_width=2)
        if mode in ('endpoints', 'assignment', 'routeA', 'routeB', 'core', 'assembly', 'top', 'side'):
            for e in endpoints:
                color = COLORS['inlet_ports' if e.role == 'inlet' else 'outlet_ports']
                if mode in ('endpoints', 'routeA', 'routeB'):
                    plotter.add_arrows(e.position[None, :], e.tangent[None, :], mag=cfg['visualization']['arrow_length_mm'], color=color)
                plotter.add_point_labels(e.position[None, :], [e.endpoint_id], font_size=19, text_color='#202733',
                    point_color=color, point_size=9, shape_color='white', shape_opacity=.75, always_visible=True)
        if mode in ('assignment', 'margin'):
            for face in enabled_faces(cfg):
                p = box.inner.mean(axis=0)
                axis, side = face_axis(face)
                p[axis] = box.inner[side, axis]
                plotter.add_point_labels(p[None, :], [face], font_size=22, text_color='#333333', show_points=False, always_visible=True)
        for r in routes:
            color = COLORS['inlet_ports' if r.endpoint.role == 'inlet' else 'outlet_ports']
            if mode == 'assignment':
                plotter.add_mesh(pv.Line(r.endpoint.position, r.target), color=color, line_width=4)
            if not use_union and (mode in ('routeA', 'routeB', 'assembly', 'top', 'side', 'local') or (mode == 'core' and core is None)):
                plotter.add_mesh(polydata(r.mesh), color=color, smooth_shading=True)
            if mode not in ('source', 'endpoints', 'margin', 'core'):
                ring = pv.Disc(center=r.target, inner=r.endpoint.radius + cfg['ports']['assembly_clearance_mm'],
                               outer=r.endpoint.radius + cfg['ports']['assembly_clearance_mm'] + .25,
                               normal=NORMALS[r.face], c_res=48)
                plotter.add_mesh(ring, color=COLORS['holes'])
                if mode in ('box', 'assignment'):
                    plotter.add_point_labels(r.target[None, :], [r.endpoint.endpoint_id], font_size=19,
                        text_color=COLORS['holes'], show_points=False, always_visible=True)
        direction = np.array([1, -1.3, 1.1])
        up = (0, 0, 1)
        if mode == 'routeB':
            direction = np.array([-1, -.8, .65])
        if mode == 'top':
            direction, up = np.array([0, 0, 1]), (0, 1, 0)
        if mode == 'side':
            direction = np.array([0, -1, 0])
        plotter.camera_position = [center + scale * 2 * direction, center, up]
        plotter.enable_parallel_projection()
        plotter.reset_camera()
        plotter.camera.zoom(.86)
        if mode == 'conflict':
            for endpoint_id, check in conflicts.items():
                e = next(e for e in endpoints if e.endpoint_id == endpoint_id)
                p, q = np.array(check['endpoint_disk_point_mm']), np.array(check['vessel_point_mm'])
                plotter.add_mesh(pv.Disc(center=e.position, inner=0, outer=e.radius, normal=e.tangent, c_res=96),
                                 color=COLORS['outlet_ports'], opacity=.9)
                plotter.add_mesh(pv.Line(p, q), color=COLORS['holes'], line_width=6)
                plotter.add_point_labels(np.array([e.position, (p + q) / 2]),
                    [endpoint_id, f"{check['distance_mm']:.3f} mm < {check['required_clearance_mm']:g} mm"],
                    font_size=22, point_color=COLORS['holes'], text_color='#a51622', always_visible=True)
                midpoint = (p + q) / 2
                plotter.camera_position = [midpoint + np.array([2, -6, 3]), midpoint, (0, 0, 1)]
                plotter.camera.parallel_scale = 5.5
        if mode == 'local':
            r = relaxed[0]
            p = r.endpoint.position
            plotter.add_point_labels(p[None, :], [f'{r.endpoint.endpoint_id}: min gap {r.vessel_clearance:.3f} mm'],
                font_size=22, point_color=COLORS['holes'], text_color='#a51622', always_visible=True)
            plotter.camera_position = [p + np.array([2, -6, 3]), p, (0, 0, 1)]
            plotter.camera.parallel_scale = 9
        plotter.add_axes(line_width=2)
        screenshot = plotter.screenshot(return_img=True)
        plotter.close()
        canvas = Image.fromarray(screenshot).convert('RGB')
        draw = ImageDraw.Draw(canvas)
        draw.rectangle([0, 0, canvas.width, 85], fill='#fafbfc')
        draw.text((30, 18), title, font=font, fill='#152b3b')
        draw.rectangle([0, canvas.height - 100, canvas.width, canvas.height], fill='#fafbfc')
        x = 30
        for text, key in [('原血管', 'vascular_core'), ('入口延长', 'inlet_ports'), ('出口延长', 'outlet_ports'), ('盒体', 'casting_box'), ('穿墙孔', 'holes')]:
            draw.rectangle([x, canvas.height - 82, x + 22, canvas.height - 60], fill=COLORS[key])
            draw.text((x + 32, canvas.height - 88), text, font=small, fill='#23394a')
            x += 230
        footer = f'单位：mm    {status}    已绘制延长：{len(routes)}/{len(endpoints)}；全部原血管保留'
        draw.text((30, canvas.height - 43), footer, font=small, fill='#465565')
        canvas.save(out / name)
    files = [str(out / name) for name, _, _ in titles]
    if core is not None and any(e.cap_ring is not None for e in endpoints):
        files.extend(render_attachment_comparison(out, source, endpoints, routes, core, cfg))
    return files


def render_attachment_comparison(out, source, endpoints, routes, core, cfg):
    """Same camera and colors on the actual old/new union; optional mesh edges."""
    previous = ROOT / cfg['input']['compact_root'] / cfg['output']['previous_directory']
    previous_qc = json.loads((previous / 'geometry_qc.json').read_text())
    old_core = trimesh.load_mesh(previous_qc['core']['stl'], process=True)
    old_layers = pv.read(previous_qc['assembly']['vtm'])
    centers, ids = [], []
    for index, key in enumerate(['vascular_core', 'inlet_ports', 'outlet_ports']):
        centers.append(old_layers[key].cell_centers().points)
        ids.append(np.full(old_layers[key].n_cells, index, dtype=int))
    palette = np.array([ImageColor.getrgb(COLORS[k]) for k in ['vascular_core', 'inlet_ports', 'outlet_ports']])
    old_poly = polydata(old_core)
    old_poly.cell_data['body_color'] = palette[np.concatenate(ids)[cKDTree(np.vstack(centers)).query(old_core.triangles_center)[1]]].astype(np.uint8)
    new_poly = polydata(core)
    new_poly.cell_data['body_color'] = palette[core_face_labels(core, source, routes)].astype(np.uint8)
    endpoint = next(e for e in endpoints if e.cap_ring is not None)
    normal = endpoint.tangent
    direction = np.cross(endpoint.attachment_qc['reference_swc_mean_tangent'], normal)
    if np.linalg.norm(direction) < 1e-6:
        direction = np.cross(endpoint.cap_ring[0] - endpoint.position, normal)
    direction = unit(direction)
    # Look from the side away from the known neighboring vessel, otherwise O2
    # can hide the O3 seam even though the camera is centered on the right cap.
    neighbor = previous_qc.get('endpoint_start_clearance', {}).get(endpoint.endpoint_id, {}).get('vessel_point_mm')
    if neighbor is not None and direction @ (np.asarray(neighbor) - endpoint.position) > 0:
        direction = -direction
    focus = endpoint.position + .5 * normal
    files = []
    for edges, name in [(False, '13_O3_attachment_comparison.png'), (True, '14_O3_attachment_mesh_comparison.png')]:
        plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 1050))
        for index, model in enumerate([old_poly, new_poly]):
            plotter.subplot(0, index)
            plotter.set_background('#fafbfc')
            plotter.add_mesh(model, scalars='body_color', rgb=True, smooth_shading=True,
                show_edges=edges, edge_color='#343b41', line_width=.6)
            ring = pv.lines_from_points(np.vstack([endpoint.cap_ring, endpoint.cap_ring[0]]))
            plotter.add_mesh(ring, color=COLORS['holes'], line_width=4)
            plotter.camera_position = [focus + 10 * direction, focus, normal]
            plotter.enable_parallel_projection()
            plotter.camera.parallel_scale = 3.5
        screenshot = plotter.screenshot(return_img=True)
        plotter.close()
        canvas = Image.fromarray(screenshot).convert('RGB')
        draw = ImageDraw.Draw(canvas)
        draw.rectangle([0, 0, 1800, 100], fill='#fafbfc')
        font = chinese_font(cfg, 30)
        small = chinese_font(cfg, 23)
        draw.text((25, 20), f"修复前：平均方向偏离端面 {endpoint.attachment_qc['before_axis_error_deg']:.2f}°", font=font, fill='#152b3b')
        draw.text((925, 20), '修复后：沿原端面轮廓和法向平滑接出', font=font, fill='#152b3b')
        draw.rectangle([0, 955, 1800, 1050], fill='#fafbfc')
        draw.text((30, 970), '青色：原血管   橙色：延长段   红线：原端面边界   左右使用相同视角和比例', font=small, fill='#23394a')
        draw.text((30, 1008), '显示实际合并后外表面；原模型、端口半径与全局打印姿态保持不变。', font=small, fill='#465565')
        path = out / name
        canvas.save(path)
        files.append(str(path))
    return files


def number(value):
    return '不适用' if value is None or not math.isfinite(value) else f'{value:.2f}'


def write_report(out, inputs, box, endpoints, routes, cfg, qc):
    n_in = sum(e.role == 'inlet' for e in endpoints)
    n_out = len(endpoints) - n_in
    dims = lambda v: ' × '.join(f'{x:.2f}' for x in v)
    lines = [
        '# BALANCED 血管芯与浇注盒设计审核', '',
        '## 1. 这次做了什么', '',
        '本轮围绕已经验收的血管芯设计了一个顶部开放的浇注盒。'
        '入口和出口的位置来自原有中心线及其编号对应表，延长段沿端点几何方向起步，再通向允许的盒体边界。'
        '原血管与延长段构成牺牲芯，盒体单独保存为另一个部件。'
        '图片用于人工检查位置、弯曲和装配关系，本轮没有切片或提交打印。'
        + (f'当前仅完成 {len(routes)}/{len(endpoints)} 个端口延长，其余端口因间距约束被保留为待调整。' if len(routes) < len(endpoints) else ''), '',
        '## 2. 输入是什么', '',
    ]
    for key, title in [('source_stl', '已验收血管 STL'), ('compact_swc', '紧凑血管中心线'),
                       ('fitted_swc', '生成最终表面的拟合中心线'), ('manifest', '来源清单'), ('transform', '固定打印坐标变换')]:
        path = inputs['paths'][key]
        lines.append(f'- {title}：`{path}`。')
    lines += ['', 'SWC 用来确认真实端点、父子关系、原分支编号和半径；STL 用来计算盒体尺寸、碰撞与最终实体。'
              '本轮不重新识别 MeVO，也不重新运行配准、NN、半径补偿、VascularMD 建模或朝向搜索。'
              '端口使用拟合中心线中与原封口记录一致的实际半径，早期补偿值和原始值一并保存在端点表中。', '',
              '![处理前：原始血管芯](QC/01_source_vascular_core.png)', '',
              '## 3. 盒子多大', '',
              f'- 内部：**{dims(box.inner_size)} mm**（X × Y × Z）。',
              f'- 外部：**{dims(box.outer_size)} mm**。',
              f"- X、Y 两侧分别留 {cfg['box']['margin_x_mm']:g}、{cfg['box']['margin_y_mm']:g} mm；原血管下方留 {cfg['box']['margin_bottom_mm']:g} mm，上方留 {cfg['box']['margin_top_mm']:g} mm。",
              f"- 侧壁厚 {cfg['box']['wall_thickness_mm']:g} mm，底板厚 {cfg['box']['bottom_thickness_mm']:g} mm，顶部无盖。", '',
              '本轮仅允许 ±X、±Y 四个侧壁；−Z 底板和 +Z 顶部均不放置端口。顶部保持开放。', '',
              '“顶部开放”指盒内空腔可从上方进入；盒壁和底板作为一个有厚度的实体，其三角网格仍应闭合，不能有破面。'
              '沿用原打印坐标，因此盒底可能位于 Z=0 以下；这不是重新选择打印姿态，后续切片时再分别定位两个部件。', '',
              '## 4. 有几个入口和出口', '', f'来源拓扑确认 **{n_in} 个入口、{n_out} 个出口**，未根据 STL 外观猜测身份。'
              f'本次仅有 **{len(routes)}/{len(endpoints)} 个端口**的延长通过检查；端点没有被删除。', '',
              'I 表示入口，O 表示出口。端点箭头是中心线向模型外侧的几何方向，不是实测血流方向。', '',
              '## 5. 每个端口在哪里', '',
              '| 端口 | 原 SWC 编号 / 分支 | 盒壁 | 路线 | 延长长度 mm | 半径 mm | 最小弯曲半径 mm |',
              '|---|---|---|---|---:|---:|---:|']
    route_map = {r.endpoint.endpoint_id: r for r in routes}
    for e in endpoints:
        r = route_map.get(e.endpoint_id)
        lines.append(f"| {e.endpoint_id}（{'入口' if e.role == 'inlet' else '出口'}） | {e.original_swc_id} / {e.branch_id} | " +
                     (f"{r.face} | {'直线' if r.kind == 'ROUTE_A_STRAIGHT' else '平滑曲线'} | {r.length:.2f} | {e.radius:.4f} | {number(r.minimum_bend_radius)} |"
                      if r else f'未找到可接受布局 | 待调整 | — | {e.radius:.4f} | — |'))
    lines += ['', '盒壁名称属于当前打印坐标，与左右脑解剖方向无关。每段长度从原端点算到最外端，不包含向原封口内侧的布尔连接重叠。', '',
              f"端口在盒外继续延伸 {cfg['ports']['outside_stub_mm']:g} mm。穿墙孔的余量为**径向 {cfg['ports']['assembly_clearance_mm']:g} mm**，即孔直径比端口直径大 {2*cfg['ports']['assembly_clearance_mm']:g} mm；孔不是软管接头。",
              f"同壁孔中心间距至少 {cfg['ports']['min_port_spacing_mm']:g} mm，且不少于两半径加两倍径向余量；孔中心距顶部、底部和侧角至少 {cfg['ports']['min_wall_edge_clearance_mm']:g} mm。", '',
              '![处理后：血管芯与通过检查的端口](QC/06_core_with_ports.png)', '',
              '## 6. 有没有风险', '']
    for warning in qc['warnings'] + qc['failures']:
        lines.append('- ' + warning)
    for endpoint_id, alignment in qc.get('attachment_alignment', {}).items():
        check = alignment.get('generated_mesh_section', {})
        lines += ['', '**本轮端口对准修复。** '
            f"{endpoint_id} 原先沿约 5 mm 中心线的平均方向起步，该方向与已验收端面的法向相差 {alignment['before_axis_error_deg']:.2f}°，造成斜接。"
            '本轮通过已有 SWC 和端口记录定位原封口，取其实际轮廓及向外法向作为连接基准；没有从 STL 猜测新的出口。'
            f"起始轮廓沿用原端面的 {alignment['cap_vertex_count']} 个点，半径仍为 {alignment['radius_mm']:.6f} mm。"
            f"向内连接改为长度 {alignment['collar_length_mm']:g} mm 的收窄锥段，其位于原血管外的数值体积为 {alignment['collar_outside_source_volume_mm3']:.3g} mm³，"
            '避免原来的完整半径倒插管从弯曲管壁侧面露出。']
        if check.get('evaluated'):
            lines.append(f"实际生成网格在端面外 {check['section_offset_mm']:g} mm 处进行截面核对：中心偏差 {check['center_offset_mm']:.3g} mm，"
                f"轮廓最大差异 {check['contour_hausdorff_mm']:.3g} mm，起始轴向误差 {check['axis_error_deg']:.6f}°。"
                '这项检查独立于闭合、连通和体积检查，专门防止斜接或错位。')
        lines += ['', '![O3 接口修复前后](QC/13_O3_attachment_comparison.png)', '',
                  '![O3 实际三角表面对照](QC/14_O3_attachment_mesh_comparison.png)', '',
                  '装配预览沿用现有配色和视角，数据改为合并后的实际外表面，不再显示两实体交叠处的内部封口。'
                  'VTM 中三个血管图层是同一外表面的互不重叠分组；两个物理部件仍分别导出。'
                  '详细数值见 [接口对准检查](attachment_alignment_qc.json)。']
    if routes:
        lines += [f"- 最大端点方向转角：{max(r.turn for r in routes):.2f}°。该值是起始切线与盒壁法线的总方向差，并非接口处的尖角。",
                  f"- 路线至非连接区血管的最小表面间距：{min(r.vessel_clearance for r in routes):.3f} mm。",
                  f"- 任意两个端口的最小表面间距：{min(r.nearest_other_port for r in routes):.3f} mm。",
                  f"- 同壁端口最小中心间距：{number(min(r.wall_spacing for r in routes))} mm。"]
    if any(not check['passed'] for check in qc.get('endpoint_start_clearance', {}).values()):
        lines += ['', '![起点间距冲突](QC/11_endpoint_clearance_conflict.png)', '',
                  '橙色圆盘是端口必须包含的起始截面，红线表示它到旁边血管的最短距离。'
                  '由于这段距离在延长起点就已经不足，改变目标盒壁、加大盒子或修改后续曲线都无法消除该冲突。'
                  '图中显示原始起点与原血管的关系；实际延长实体（含连接重叠）的最小间距另行逐口核验。']
    if qc.get('native_restoration'):
        restoration = qc['native_restoration']
        lines += ['', '**先检查能否补回原始血管。** 本轮读取完整原始 SWC，并核对源文件哈希、所有保留节点坐标及有向父子边。'
                  '补回长度应由真实下游路径和间距需要决定，不能凭空外推血管，也不能把旁边的血管接到当前出口。']
        for row in restoration['endpoints']:
            if row['restoration_needed']:
                lines.append(f"{row['endpoint_id']} 对应原始节点 {row['original_swc_id']}，原始子节点为 {row['raw_child_ids']}，"
                    f"可补回的真实下游长度为 **{row['available_arc_length_mm']:.3f} mm**。"
                    + ('该节点已是完整原始模型的末端，所以自适应搜索也没有可继续选取的原血管；本次没有伪造补回段。' if row['status'] == 'ORIGINAL_TERMINAL_NO_CONTINUATION' else '可用路径与停止原因保存在检查记录中。'))
        lines += ['', '完整依据见 [原始血管补回检查](native_restoration_audit.json)。']
    if qc.get('clearance_fallback', {}).get('used'):
        fallback = qc['clearance_fallback']
        lines += ['', '**最后手段：仅对受限端口放宽间距。** 在确认原始血管无法补回后，按用户允许的顺序，'
                  f"以 {cfg['clearance_fallback']['decrement_mm']:g} mm 步长由 {cfg['ports']['route_clearance_mm']:g} mm 向下搜索。"
                  f"首次完成全部端口布局的门槛为 **{fallback['selected_clearance_mm']:g} mm**，"
                  f"仅适用于 {', '.join(fallback['affected_endpoints'])} 到原血管的距离。"
                  f"其他端口以及端口两两之间仍要求 {cfg['ports']['route_clearance_mm']:g} mm，血管相交、壁边距不足或弯曲过急的候选仍被拒绝。"
                  '这不代表已经满足原来的全程 3 mm 条件，需人工确认局部 PDMS 壁厚是否可接受。']
        for name in fallback['affected_endpoints']:
            profile = qc['clearance_profiles'][name]
            end = profile['last_outward_sample_below_preferred_mm']
            lines.append(f"{name} 完整实体的最小间距为 {profile['minimum_complete_mesh_clearance_mm']:.4f} mm；"
                f"从原封口向外的截面抽样最小值为 {profile['minimum_outward_sampled_gap_mm']:.4f} mm。"
                + (f"最后一个低于 3 mm 的抽样截面位于起点后约 {end:.2f} mm。" if end is not None else '向外抽样截面均达到首选间距。')
                + '抽样仅用于定位局部狭窄位置，是否相交仍由完整三角表面检查决定。')
        lines += ['', '[逐级搜索记录](clearance_fallback_audit.json)；[沿程间距表](tables/route_clearance_profile.csv)。']
        lines += ['', '![间距例外位置的实际延长模型](QC/12_local_clearance_review.png)']
    if qc.get('core'):
        lines += ['', f"已输出的芯体包含 {qc['core']['connected_components']} 个连通实体，闭合检查为 {qc['core']['watertight']}，"
                  f"退化面为 {qc['core']['degenerate_faces']} 个；原血管体积减损的数值检验结果为 {qc['core']['source_subtracted_volume_mm3']:.3g} mm³，低于 {cfg['geometry']['volume_tolerance_mm3']:g} mm³ 数值容差。"
                  f"盒体已开 {qc['box'].get('port_hole_count', 0)} 个侧壁穿孔；完整设计要求全部 {len(endpoints)} 个端点均有可用通路。"]
    lines += ['', '**设计参数的含义。**', '',
              f"端点方向通过约 {cfg['routing']['tangent_window_mm']:g} mm 的中心线估计，按 {cfg['routing']['tangent_resample_step_mm']:g} mm 重采样后作直线回归，真实窗口和原采样点数逐口记录。"
              f"只有方向差不超过 {cfg['ports']['preferred_max_turn_angle_deg']:g}° 且间距满足要求时才考虑沿切线直连，否则使用 CadQuery 曲线扫掠。"
              f"弯曲半径至少 {cfg['ports']['minimum_bend_radius_mm']:g} mm，并以 {cfg['routing']['bend_radius_safety_factor']:g} 倍安全系数筛选；曲率按不超过 {cfg['routing']['curvature_step_mm']:g} mm 的弧长间隔查询 OCC。"
              f"首选血管及端口之间保持 {cfg['ports']['route_clearance_mm']:g} mm 表面间距；仅上面明确记录的最后手段例外可降低下限。无论是否放宽，均另留两倍 {cfg['geometry']['mesh_tolerance_mm']:g} mm 网格离散余量。", '',
              f"未启用端面对准修复的端口沿用向原封口内侧重叠 {cfg['routing']['attachment_overlap_mm']:g} mm 的方法；已修复端口使用上文记录的内埋连接段。碰撞检查仅豁免同一父分支端点附近 {cfg['routing']['attachment_exemption_mm']:g} mm 的连接区；其他分支及更远父分支仍参与检查。"
              f"表面归属依据已有 SWC，以 {cfg['geometry']['centerline_label_sample_step_mm']:g} mm 间隔关联；一个面仅在三个顶点都属于此连接区时才豁免。"
              '表面最短距离由 FCL 计算，不以最近网格顶点距离代替。', '',
              '允许的候选边界为 ' + '、'.join(enabled_faces(cfg)) + '，按距离、方向变化、碰撞和边缘条件评分；权重依次为 ' +
              '、'.join(str(v) for v in cfg['routing']['score_weights'].values()) + '。距离除以盒内对角线，角度除以 180°，不满足条件的候选被拒绝。'
              f"优先入口和出口相对，若曲线总转向超过偏好值 {cfg['routing']['preferred_max_curved_turn_angle_deg']:g}° 则尝试其他布局；出口最多占 {cfg['routing']['max_outlet_faces']} 个相邻面。"
              '所有候选（含被拒绝者）的分数、坐标和原因保存在 CSV。', '',
              '为避免把目标孔强制放在端点同一高度，程序还尝试沿端点方向顺势转弯。'
              f"试用的弯曲半径为 {cfg['routing']['natural_bend_radii_mm']} mm，各用 {cfg['routing']['natural_bend_interpolation_points']} 个圆弧导点交给 CadQuery 拟合，转弯后直达盒壁。"
              f"另一组候选在壁面横向偏移 {cfg['routing']['target_lateral_offsets_mm']} mm、竖向偏移 {cfg['routing']['target_vertical_offsets_mm']} mm；两点样条的端切线幅度取端点间距的 {cfg['routing']['spline_tangent_scales']} 倍。"
              f"每面保留至多 {cfg['routing']['options_per_face']} 个不同目标位置进入组合比较，曲线在穿墙前至少直行 {cfg['routing']['wall_approach_mm']:g} mm。", '',
              f"几何离散角度精度为 {cfg['geometry']['mesh_angular_tolerance_rad']:g} rad，来源坐标核对容差为 {cfg['geometry']['coordinate_tolerance_mm']:g} mm，面积不大于 {cfg['geometry']['minimum_triangle_area_mm2']:g} mm² 的三角形按退化面处理。"
              f"仅当合并结果出现数值退化面时，才使用 Manifold 官方简化功能，容差为 {cfg['geometry']['boolean_cleanup_tolerance_mm']:g} mm；不修改输入 STL，并重新检查体积保留和连通性。"
              '完整数值和图像显示参数保存在 [本次配置](resolved_design_config.yaml)。表内空的距离或弯曲半径表示不适用，不代表零。', '',
              'VMTK 检查：' + qc['vmtk']['decision'] + '。', '',
              '## 7. 当前能不能进入下一步', '', f"**{qc['status']}**", '',
              ('几何检查通过，可先人工查看图像和装配文件，再决定是否开展后续打印验证。' if qc['status'] == 'READY_FOR_HUMAN_REVIEW'
               else '当前仍有约束未满足，需调整本轮结构设计后再审核；失败原因保存在 failure_summary.md。'), '',
              '这不是打印批准。两件套没有增加可拆壁或卡扣，刚性芯能否穿入多个固定侧孔、穿墙余量如何密封，以及实际打印误差，仍需人工确认。', '',
              ('本次 STL/STEP 文件名含 PARTIAL_REVIEW_ONLY，仅供检查已完成部分，不能作为全部端口完成的设计。'
               '没有导出完整最终芯体或完整最终盒体，3MF 暂不生成。' if len(routes) != len(endpoints)
               else '所有端口均已包含在本次审核模型中。'), '',
              '审核入口：[装配图](QC/08_assembly_preview.png)；[顶视图](QC/09_top_view.png)；'
              '[交互装配](assembly/assembly_preview.vtm)；[端口表](tables/port_layout_summary.csv)；'
              '[完整检查](geometry_qc.json)；[源文件哈希](protected_source_hashes.json)。', '',
              'VTM 在 ParaView 中打开后可分别控制 vascular_core、inlet_ports、outlet_ports、casting_box 四层。3MF 若成功导出则只含两个物理部件。', '']
    (out / 'design_report.md').write_text('\n'.join(lines), encoding='utf-8')
