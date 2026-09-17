"""
可视化命令行入口。

读取 visualization_config.yaml，找到要显示的 SWC 文件，然后调用
`render_swc` 完成 mesh 构建、截图保存和交互窗口显示。
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
import re

from .planar_publication import (
    is_xz_planar,
    render_planar_publication_figure,
)
from .pyvista_viewer import render_swc
from .swc_loader import list_swc_files, read_swc
from .volumetric_publication import render_volumetric_publication_figure
from .yaml_config import (
    DEFAULT_CONFIG_PATH,
    DEFAULT_SWC_DIR,
    load_visualization_config,
)


_RESULT_TIMESTAMP = re.compile(
    r"^(?P<date>\d{8})_(?P<time>\d{6})(?:_(?P<microseconds>\d{6}))?(?:_|$)"
)


def _folder_timestamp(directory: Path) -> datetime | None:
    """Parse the generator timestamp prefix from a result-folder name."""
    match = _RESULT_TIMESTAMP.match(directory.name)
    if match is None:
        return None
    microseconds = match.group("microseconds") or "000000"
    try:
        return datetime.strptime(
            f"{match.group('date')}_{match.group('time')}_{microseconds}",
            "%Y%m%d_%H%M%S_%f",
        )
    except ValueError:
        return None


def newest_result_folder(result_root: Path | None = None) -> Path:
    """Return the latest timestamped generator folder containing an SWC file."""
    result_root = (
        DEFAULT_SWC_DIR
        if result_root is None
        else Path(result_root)
    )
    if not result_root.is_dir():
        raise FileNotFoundError(
            f"Vessel result root does not exist: {result_root}."
        )
    candidates = [
        directory
        for directory in result_root.iterdir()
        if directory.is_dir() and any(directory.glob("*.swc"))
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No generated result folder containing an SWC file was found under {result_root}."
        )

    timestamped_candidates = [
        (timestamp, directory)
        for directory in candidates
        if (timestamp := _folder_timestamp(directory)) is not None
    ]
    if timestamped_candidates:
        return max(
            timestamped_candidates,
            key=lambda item: (item[0], item[1].name),
        )[1]

    # Backward compatibility for historical result folders without a
    # timestamp prefix.
    return max(candidates, key=lambda directory: directory.stat().st_mtime)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """解析可视化命令行参数。"""
    # argparse 负责把命令行字符串转换成 Python 对象，后续代码只读取 args。
    parser = argparse.ArgumentParser(description=__doc__)
    # --config 指向可视化 YAML；不给时使用仓库内的默认配置文件。
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help=f"Visualization YAML configuration. Default: {DEFAULT_CONFIG_PATH}",
    )
    # --list 是一个开关参数：出现即为 True，用于只查看目录下有哪些 SWC。
    parser.add_argument(
        "--list",
        action="store_true",
        help="List SWC files from the YAML input.swc_dir and exit.",
    )
    parser.add_argument(
        "--result-folder",
        type=Path,
        default=None,
        help=(
            "Generator result folder to visualize. When omitted, the newest "
            "timestamped result folder is selected automatically."
        ),
    )
    window_group = parser.add_mutually_exclusive_group()
    window_group.add_argument(
        "--show-window",
        dest="show_window",
        action="store_true",
        default=None,
        help=(
            "Open the zoomable PyVista window, overriding the YAML setting "
            "(this is enabled by the default configuration)."
        ),
    )
    window_group.add_argument(
        "--no-window",
        dest="show_window",
        action="store_false",
        help=(
            "Only save image/mesh outputs; do not open an interactive window."
        ),
    )
    # 返回 Namespace，调用者可以通过 args.config、args.list 访问字段。
    return parser.parse_args(argv)


def main(result_folder=None) -> None:
    """可视化主流程。"""
    # 第一步：读取命令行参数，确定配置文件路径和是否只列文件。
    args = parse_args()
    # 第二步：把 YAML 配置解析成结构化对象，并在解析阶段完成路径检查。
    selected_result_folder = (
        args.result_folder
        if args.result_folder is not None
        else result_folder
    )
    if selected_result_folder is None:
        selected_result_folder = newest_result_folder()
        print(
            "RESULT_FOLDER_NAME is None; automatically selected latest result: "
            f"{selected_result_folder}"
        )
    cfg = load_visualization_config(
        args.config,
        result_folder=selected_result_folder,
    )

    # 如果用户只想列出 SWC 文件，就不创建 PyVista 窗口，也不构建 mesh。
    if args.list:
        # 从配置中的输入目录收集所有 .swc 文件。
        files = list_swc_files(cfg.swc_dir)
        # 空目录给出明确提示，避免用户误以为程序卡住或没有输出。
        if not files:
            print(f"No .swc files found under {cfg.swc_dir}")
            return
        # 只打印文件名，便于用户把名字复制到 visualization_config.yaml。
        for path in files:
            print(path.name)
        return

    # 正常模式下读取单个 SWC，并根据配置渲染、截图或导出 mesh。
    source_tree = read_swc(cfg.swc_path)
    use_publication_figure = (
        cfg.publication_style
        and cfg.screenshot_path is not None
    )
    actual_planar = is_xz_planar(source_tree)
    expected_planar = cfg.geometry_mode == "planar_2d"
    if actual_planar != expected_planar:
        raise ValueError(
            "Result geometry metadata/name disagrees with the SWC coordinates: "
            f"geometry_mode={cfg.geometry_mode!r}, planar_coordinates={actual_planar}."
        )
    render_view = "xz" if actual_planar else "iso"
    if use_publication_figure:
        pyvista_screenshot = cfg.result_dir / (
            f"{cfg.swc_path.stem}_{cfg.radius_style}_{render_view}_pyvista_view.png"
        )
    else:
        pyvista_screenshot = cfg.screenshot_path
    show_window = (
        cfg.show_window
        if args.show_window is None
        else bool(args.show_window)
    )
    off_screen = (
        cfg.off_screen
        if args.show_window is None
        else not show_window
    )

    if show_window:
        print(
            "Opening interactive PyVista window: mouse wheel/right drag = zoom, "
            "left drag = rotate, middle drag = pan, R = reset view."
        )

    tree, mesh_data = render_swc(
        cfg.swc_path,
        radius_style=cfg.radius_style,
        sides=cfg.sides,
        radius_scale=cfg.radius_scale,
        cap_style=cfg.cap_style,
        color=cfg.color,
        color_map=cfg.color_map,
        color_by_radius=cfg.color_by_radius,
        background=cfg.background,
        background_top=cfg.background_top,
        show_edges=cfg.show_edges,
        show_axes=cfg.show_axes,
        show_title=cfg.show_title,
        scalar_bar=cfg.scalar_bar,
        publication_style=cfg.publication_style,
        font_size=cfg.font_size,
        pbr=cfg.pbr,
        metallic=cfg.metallic,
        roughness=cfg.roughness,
        ambient=cfg.ambient,
        diffuse=cfg.diffuse,
        specular=cfg.specular,
        specular_power=cfg.specular_power,
        opacity=cfg.opacity,
        lighting=cfg.lighting,
        eye_dome_lighting=cfg.eye_dome_lighting,
        shadows=cfg.shadows,
        anti_aliasing=cfg.anti_aliasing,
        view=render_view,
        window_size=cfg.window_size,
        screenshot=pyvista_screenshot,
        export_mesh=cfg.mesh_path,
        show=show_window,
        off_screen=off_screen,
    )

    publication_files = []
    if use_publication_figure:
        renderer = (
            render_planar_publication_figure
            if actual_planar
            else render_volumetric_publication_figure
        )
        publication_files = renderer(
            tree,
            cfg.screenshot_path,
            color_map=cfg.color_map,
            font_size=cfg.font_size,
        )

    # 下面的日志用于核对实际读取的是哪个配置文件和哪个 SWC 文件。
    print(f"Read config: {Path(args.config).resolve()}")
    print(f"Loaded SWC: {cfg.swc_path}")
    print(f"Geometry mode: {cfg.geometry_mode}")
    # 节点数和边数可以快速判断 SWC 是否完整加载。
    print(f"Nodes: {tree.node_count}, edges: {tree.edge_count}")
    # mesh 点数和面片数用于判断圆柱化后几何体规模是否合理。
    print(
        f"Mesh: {len(mesh_data.points)} points, "
        f"{len(mesh_data.cell_radii_um)} faces, "
        f"skipped_edges={mesh_data.skipped_edges}"
    )
    # 半径风格会显著影响显示效果，因此单独打印。
    print(f"Radius style: {cfg.radius_style}")
    # 只有配置了截图路径时才打印截图结果。
    if cfg.screenshot_path:
        print(f"Screenshot: {cfg.screenshot_path}")
    if use_publication_figure:
        print(f"PyVista popup screenshot: {pyvista_screenshot}")
    for path in publication_files[1:]:
        print(f"Publication file: {path}")
    # 只有配置了 mesh 导出路径时才打印导出结果。
    if cfg.mesh_path:
        print(f"Mesh file: {cfg.mesh_path}")


# 允许用户直接执行 python -m 或 python 文件启动可视化。
if __name__ == "__main__":
    main()
