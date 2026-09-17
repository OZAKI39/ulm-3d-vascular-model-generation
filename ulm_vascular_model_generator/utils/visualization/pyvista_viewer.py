"""
PyVista 渲染器。

数据流是：
SWC 文件 -> `read_swc` 解析成树 -> `build_cylinder_mesh` 生成圆柱 mesh
-> PyVista 显示、截图或导出 mesh 文件。
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from .cylinder_mesh import CapStyle, MeshData, RadiusStyle, build_cylinder_mesh
from .publication_style import (
    add_publication_axes,
    add_publication_effects,
    add_smooth_junctions,
)
from .swc_loader import SWCTree, read_swc

ViewMode = Literal["xz", "xy", "yz", "iso"]


def _require_pyvista():
    """延迟导入 PyVista；没有安装时给出清晰的安装提示。"""
    # PyVista 只在真正可视化时需要；延迟导入可让纯生成流程不依赖它。
    try:
        import pyvista as pv
    # 如果环境里没有 pyvista，就给出更容易看懂的错误提示。
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "PyVista is required for this viewer. Install it with: "
            "d:\\anaconda3\\envs\\pmp\\python.exe -m pip install pyvista"
        ) from exc
    # 返回模块对象，调用方可以继续使用 pv.PolyData、pv.Plotter 等 API。
    return pv


def swc_to_polydata(
    swc_path: Path,
    *,
    radius_style: RadiusStyle = "cylindrical",
    sides: int = 18,
    radius_scale: float = 1.0,
    cap_style: CapStyle = "terminal",
    publication_style: bool = False,
):
    """读取 SWC 并转换成单个 PyVista PolyData。"""
    # 先确保 PyVista 可用；这样错误会在进入转换前就暴露。
    pv = _require_pyvista()
    # 读取 SWC 文本，得到节点字典和父子边列表。
    tree = read_swc(swc_path)
    # 把每条父子边转成圆柱面片。
    mesh = build_cylinder_mesh(
        tree,
        radius_style=radius_style,
        sides=sides,
        radius_scale=radius_scale,
        cap_style=cap_style,
    )
    # 如果没有任何顶点，说明 SWC 没有可绘制的有效边。
    if mesh.points.size == 0:
        raise ValueError(f"{swc_path} does not contain drawable SWC edges.")

    # PyVista PolyData 需要顶点数组和面片编码数组。
    polydata = pv.PolyData(mesh.points, mesh.faces)
    # 把半径写入 cell_data 后，可以在渲染时按半径上色。
    polydata.cell_data["radius_um"] = mesh.cell_radii_um
    # edge_id 让面片能追溯回原始 SWC 边，便于排查问题。
    polydata.cell_data["edge_id"] = mesh.cell_edge_ids
    # 生成平滑法线，改善圆柱在 PBR/Phong 光照下的高光连续性。
    try:
        polydata = polydata.compute_normals(
            cell_normals=False,
            point_normals=True,
            split_vertices=False,
            auto_orient_normals=True,
        )
    except Exception:
        # 老版本 VTK/PyVista 在法线生成失败时仍可继续显示原始 mesh。
        pass
    # 同时返回 PolyData、原始树和 mesh 缓冲区，兼顾显示和日志统计。
    if publication_style:
        polydata = add_smooth_junctions(
            pv,
            polydata,
            tree,
            radius_scale,
            sides,
        )
    return polydata, tree, mesh


def _set_camera(
    plotter,
    mesh,
    view: ViewMode,
    parallel_scale_factor=0.62,
) -> None:
    """根据配置选择 XZ/XY/YZ/等轴测视角。"""
    # bounds 是 mesh 的包围盒，用于估算相机需要离对象多远。
    bounds = mesh.bounds
    # center 是可视化对象中心点，相机会朝这里看。
    center = mesh.center
    # span 取三个方向中的最大尺寸，保证整个血管网络能被画面包含。
    span = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )
    # 相机距离略大于对象尺寸，避免模型贴到画面边缘。
    distance = span * 1.7

    # XZ 视图：沿负 y 方向看向模型，适合当前生成器的平面血管。
    if view == "xz":
        plotter.camera_position = [
            (center[0], center[1] - distance, center[2]),
            center,
            (0.0, 0.0, 1.0),
        ]
        # 平行投影消除透视缩放，便于观察血管拓扑。
        plotter.camera.parallel_projection = True
        # parallel_scale 控制画面显示范围。
        plotter.camera.parallel_scale = parallel_scale_factor * span
    # XY 视图：沿正 z 方向俯视。
    elif view == "xy":
        plotter.camera_position = [
            (center[0], center[1], center[2] + distance),
            center,
            (0.0, 1.0, 0.0),
        ]
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = parallel_scale_factor * span
    # YZ 视图：沿负 x 方向侧视。
    elif view == "yz":
        plotter.camera_position = [
            (center[0] - distance, center[1], center[2]),
            center,
            (0.0, 0.0, 1.0),
        ]
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = parallel_scale_factor * span
    # 其他情况使用 PyVista 默认等轴测视角，适合看三维形态。
    else:
        plotter.view_isometric()
        plotter.camera.parallel_projection = False


def _set_background(plotter, background: str, background_top: str | None) -> None:
    """设置纯色或渐变背景，兼容不同 PyVista 版本。"""
    if background_top:
        try:
            plotter.set_background(background, top=background_top)
            return
        except TypeError:
            pass
    plotter.set_background(background)


def _enable_anti_aliasing(plotter, mode: str) -> None:
    """启用抗锯齿；不支持时静默退化为默认渲染。"""
    mode = mode.lower()
    if mode == "none":
        return
    try:
        plotter.enable_anti_aliasing(mode)
    except Exception:
        pass


def _enable_depth_and_shadow(plotter, *, eye_dome_lighting: bool, shadows: bool) -> None:
    """启用提高空间层次感的 VTK 后处理。"""
    if eye_dome_lighting:
        try:
            plotter.enable_eye_dome_lighting()
        except Exception:
            pass
    if shadows:
        try:
            plotter.enable_shadows()
        except Exception:
            pass


def _configure_lighting(plotter, pv, mesh, enabled: bool) -> None:
    """配置三点布光，让圆柱血管有稳定的高光和轮廓。"""
    if not enabled:
        return

    bounds = mesh.bounds
    center = mesh.center
    span = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )

    try:
        plotter.remove_all_lights()
    except Exception:
        pass

    light_specs = [
        ((center[0] - 0.45 * span, center[1] - 1.25 * span, center[2] + 1.15 * span), "#fff4e5", 0.95),
        ((center[0] + 1.20 * span, center[1] - 0.75 * span, center[2] + 0.45 * span), "#c7dcff", 0.38),
        ((center[0] + 0.20 * span, center[1] + 1.05 * span, center[2] + 0.95 * span), "#ffffff", 0.50),
    ]
    for position, color, intensity in light_specs:
        try:
            light = pv.Light(
                position=position,
                focal_point=center,
                color=color,
                intensity=float(intensity),
            )
            plotter.add_light(light)
        except Exception:
            return


def _add_mesh_with_fallback(plotter, mesh, **kwargs) -> None:
    """
    添加 mesh，并在旧 PyVista 版本不支持 PBR/材质参数时自动降级。

    这样配置可以默认开启高级渲染，但不会因为某个环境的 PyVista 版本偏旧而
    影响基础可视化。
    """
    try:
        plotter.add_mesh(mesh, **kwargs)
        return
    except TypeError:
        pass

    pbr_keys = {"pbr", "metallic", "roughness"}
    reduced = {key: value for key, value in kwargs.items() if key not in pbr_keys}
    try:
        plotter.add_mesh(mesh, **reduced)
        return
    except TypeError:
        pass

    material_keys = {
        "ambient",
        "diffuse",
        "specular",
        "specular_power",
        "opacity",
        "lighting",
    }
    minimal = {
        key: value
        for key, value in reduced.items()
        if key not in material_keys
    }
    plotter.add_mesh(mesh, **minimal)


def render_swc(
    swc_path: Path,
    *,
    radius_style: RadiusStyle = "cylindrical",
    sides: int = 18,
    radius_scale: float = 1.0,
    cap_style: CapStyle = "terminal",
    color: str = "#0b37b7",
    color_map: str = "turbo",
    color_by_radius: bool = False,
    background: str = "#d7e4ec",
    background_top: str | None = None,
    show_edges: bool = False,
    show_axes: bool = True,
    show_title: bool = True,
    scalar_bar: bool = True,
    publication_style: bool = False,
    font_size: int = 18,
    pbr: bool = False,
    metallic: float = 0.0,
    roughness: float = 0.55,
    ambient: float = 0.34,
    diffuse: float = 0.82,
    specular: float = 0.18,
    specular_power: float = 45.0,
    opacity: float = 1.0,
    lighting: bool = True,
    eye_dome_lighting: bool = False,
    shadows: bool = False,
    anti_aliasing: str = "msaa",
    view: ViewMode = "xz",
    window_size: tuple[int, int] = (1920, 1350),
    screenshot: Path | None = None,
    export_mesh: Path | None = None,
    show: bool = True,
    off_screen: bool = False,
) -> tuple[SWCTree, MeshData]:
    """
    渲染一棵 SWC 血管树。

    这个函数同时支持交互显示、保存截图和导出 mesh。即使不显示窗口，也会
    返回解析后的树和 mesh 数据，方便测试或批处理。
    """
    # 确保 PyVista 可用，后续 plotter 和 PolyData 都依赖它。
    pv = _require_pyvista()
    # 统一输入路径类型，避免字符串和 Path 混用。
    swc_path = Path(swc_path)
    # 可选路径只有在配置了值时才转成 Path。
    screenshot = Path(screenshot) if screenshot else None
    export_mesh = Path(export_mesh) if export_mesh else None

    # 先完成 SWC -> PolyData 的几何转换。
    mesh, tree, mesh_data = swc_to_polydata(
        swc_path,
        radius_style=radius_style,
        sides=sides,
        radius_scale=radius_scale,
        cap_style=cap_style,
        publication_style=publication_style,
    )

    # 如果用户要求导出 mesh，就先保存到磁盘；这一步不要求打开窗口。
    if export_mesh:
        # 确保输出目录存在，避免 save 因目录缺失失败。
        export_mesh.parent.mkdir(parents=True, exist_ok=True)
        # PyVista 根据后缀选择具体 mesh 文件格式。
        mesh.save(export_mesh)

    # 只有需要显示窗口或保存截图时才创建 Plotter。
    needs_render = show or screenshot is not None
    # 如果只是导出 mesh，就直接返回解析数据，节省渲染开销。
    if not needs_render:
        return tree, mesh_data

    # off_screen 可在无显示器环境下截图；有截图但不显示窗口时自动启用。
    # A visible window must never inherit off-screen rendering from a stale or
    # contradictory YAML setting.  When no window is requested, screenshots
    # are rendered off-screen automatically.
    effective_off_screen = (
        False
        if show
        else off_screen or screenshot is not None
    )
    plotter = pv.Plotter(
        off_screen=effective_off_screen,
        window_size=window_size,
    )
    # 设置背景色，默认使用浅色渐变以提高论文截图的空间层次。
    _set_background(plotter, background, background_top)
    _enable_anti_aliasing(plotter, anti_aliasing)
    _configure_lighting(plotter, pv, mesh, lighting)

    mesh_kwargs = {
        "smooth_shading": True,
        "show_edges": show_edges,
        "pbr": bool(pbr),
        "metallic": float(metallic),
        "roughness": float(roughness),
        "ambient": float(ambient),
        "diffuse": float(diffuse),
        "specular": float(specular),
        "specular_power": float(specular_power),
        "opacity": float(opacity),
        "lighting": bool(lighting),
    }

    # color_by_radius=True 时用半径标量映射颜色，便于检查粗细分布。
    if color_by_radius:
        if publication_style:
            scalar_bar_args = {
                "title": "Radius, r (\u00b5m)",
                "title_font_size": font_size,
                "label_font_size": max(int(round(0.82 * font_size)), 12),
                "shadow": False,
                "n_labels": 5,
                "fmt": "%.1f",
                "position_x": 0.9,
                "position_y": 0.25,
                "width": 0.025,
                "height": 0.5,
                "vertical": True,
                "font_family": "arial",
                "color": "#20262d",
            }
        else:
            scalar_bar_args = {
                "title": "radius (um)",
                "title_font_size": 14,
                "label_font_size": 12,
                "shadow": False,
                "n_labels": 4,
                "fmt": "%.2g",
                "position_x": 0.23,
                "position_y": 0.045,
                "width": 0.56,
                "height": 0.045,
            }
        _add_mesh_with_fallback(
            plotter,
            mesh,
            scalars="radius_um",
            cmap=color_map,
            show_scalar_bar=scalar_bar,
            scalar_bar_args=scalar_bar_args,
            interpolate_before_map=True,
            **mesh_kwargs,
        )
    # 否则使用统一颜色，突出网络几何形态。
    else:
        _add_mesh_with_fallback(
            plotter,
            mesh,
            color=color,
            **mesh_kwargs,
        )

    _enable_depth_and_shadow(plotter, eye_dome_lighting=eye_dome_lighting, shadows=shadows)
    if publication_style:
        add_publication_effects(
            plotter,
            pv,
            mesh,
            tree,
            view,
            font_size,
            not show_axes,
        )

    # 加坐标轴帮助判断当前视角和空间方向。
    if show_axes:
        if publication_style:
            add_publication_axes(
                plotter,
                pv,
                mesh,
                view,
                font_size,
            )
        else:
            plotter.add_axes(line_width=2)
    # 左上角文字展示文件名、段数和半径风格，便于截图留档。
    if show_title:
        plotter.add_text(
            f"{swc_path.name} | {tree.edge_count} segments | {radius_style}",
            position="upper_left",
            font_size=font_size,
            color="black",
        )
    # 根据配置设置相机位置。
    camera_scale = 0.68 if publication_style and show_axes else 0.62
    _set_camera(
        plotter,
        mesh,
        view,
        parallel_scale_factor=camera_scale,
    )

    # 有截图路径时，show 会渲染一帧并保存图片。
    if screenshot:
        # 确保截图目录存在。
        screenshot.parent.mkdir(parents=True, exist_ok=True)
        # auto_close=True 防止批处理时窗口资源残留。
        plotter.show(screenshot=str(screenshot), auto_close=True)
    # 没有截图时进入普通交互窗口。
    else:
        plotter.show(auto_close=True)

    # 返回树和 mesh 元数据，命令行入口会用它打印统计信息。
    return tree, mesh_data
