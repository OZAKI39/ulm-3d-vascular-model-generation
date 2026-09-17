"""最终收敛场的二维面内壁面剪切应力代理 PyVista/VTK 可视化。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .vtk_flow_grid import (
    LUMEN_ARRAY,
    WALL_MASK_ARRAY,
    WALL_SHEAR_ARRAY,
    WALL_SHEAR_DISPLAY_MASK_ARRAY,
)


@dataclass(frozen=True)
class WallShearArtifacts:
    """最终 WSS 代理的离线交互场景和静态预览。"""

    html_path: Path
    preview_path: Path


def render_wall_shear_visualization(
    result_dir: Path,
    *,
    domain=None,
    raster=None,
    flow=None,
) -> WallShearArtifacts:
    """渲染最终 WSS 代理；可直接使用内存场，也可读取已有最终 VTI。"""

    from .pyvista_flow import _customize_exported_html, validate_cfd_flow_dependencies

    validate_cfd_flow_dependencies()
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    artifacts = _artifact_paths(result_dir)
    supplied = (domain is not None, raster is not None, flow is not None)
    if any(supplied) and not all(supplied):
        raise ValueError("domain, raster, and flow must be supplied together for in-memory WSS rendering.")
    if all(supplied):
        from .vtk_flow_grid import build_vtk_stage_grid

        stage_grid = build_vtk_stage_grid(domain, raster, flow, stage="final", include_lic=False)
        plotter, _ = build_wall_shear_plotter_from_image(stage_grid.image_grid)
    else:
        plotter, _ = build_wall_shear_plotter(result_dir)
    try:
        plotter.show(
            auto_close=False,
            interactive=False,
            screenshot=str(artifacts.preview_path),
        )
        plotter.export_html(artifacts.html_path)
        _customize_exported_html(
            artifacts.html_path,
            "Final 2D in-plane wall-shear-stress proxy",
        )
    finally:
        plotter.close()
    return artifacts


def build_wall_shear_plotter(result_dir: Path):
    """从最终 VTI 构造可供离线导出或 trame 使用的 WSS 单面板。"""

    from .pyvista_flow import validate_cfd_flow_dependencies

    validate_cfd_flow_dependencies()
    import pyvista as pv

    result_dir = Path(result_dir)
    field_path = result_dir / "final_flow_field.vti"
    if not field_path.is_file():
        raise ValueError(f"Missing {field_path.name} under {result_dir}.")
    image = pv.read(field_path)
    if not isinstance(image, pv.ImageData):
        raise ValueError("The wall-shear viewer expects final_flow_field.vti to contain vtkImageData.")
    return build_wall_shear_plotter_from_image(image)


def build_wall_shear_plotter_from_image(image):
    """从内存中的最终 vtkImageData 构造 WSS 单面板。"""

    from .pyvista_flow import validate_cfd_flow_dependencies

    validate_cfd_flow_dependencies()
    import pyvista as pv

    if not isinstance(image, pv.ImageData):
        raise ValueError("The wall-shear viewer expects final flow data as vtkImageData.")
    _validate_wall_shear_image(image)
    _ensure_wall_shear_display_mask(image)

    display_mask = np.asarray(image.cell_data[WALL_SHEAR_DISPLAY_MASK_ARRAY], dtype=bool)
    wall_shear = np.asarray(image.cell_data[WALL_SHEAR_ARRAY], dtype=float)
    display_values = wall_shear[display_mask]
    if display_values.size == 0:
        raise ValueError("No closed-wall cells remain after excluding inlet and outlet open faces.")
    if not np.all(np.isfinite(display_values)):
        raise ValueError("The final wall-shear display band contains NaN or Inf values.")
    if np.any(display_values < -1.0e-12):
        raise ValueError("Wall-shear-stress magnitude must be non-negative.")
    upper = float(np.percentile(display_values, 99.5))
    if not np.isfinite(upper) or upper <= np.finfo(float).eps:
        upper = 1.0

    fluid = image.threshold((0.5, 1.5), scalars=LUMEN_ARRAY, preference="cell")
    wall_band = image.threshold(
        (0.5, 1.5),
        scalars=WALL_SHEAR_DISPLAY_MASK_ARRAY,
        preference="cell",
    )
    spacing = float(image.spacing[0])
    wall_band = wall_band.copy(deep=True)
    # VTK.js 的本地 WebGL 深度精度低于桌面 VTK。使用一个网格量级的纯渲染
    # 抬升，避免 wall band 与背景共面后发生 z-fighting；X/Z 物理坐标不变。
    wall_band.points[:, 2] = 1.00 * spacing
    outline = fluid.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    outline = outline.copy(deep=True)
    outline.points[:, 2] = 0.50 * spacing

    plotter = pv.Plotter(
        off_screen=True,
        border=True,
        border_color="#475569",
        window_size=(1500, 900),
    )
    try:
        _populate_wall_shear_plotter(
            plotter,
            image=image,
            wall_band=wall_band,
            outline=outline,
            upper=upper,
        )
    except Exception:
        plotter.close()
        raise
    return plotter, image


def _populate_wall_shear_plotter(plotter, *, image, wall_band, outline, upper: float) -> None:
    from .pyvista_flow import _configure_planar_renderer, _scalar_bar_args

    plotter.set_background("white")
    plotter.add_mesh(
        wall_band,
        scalars=WALL_SHEAR_ARRAY,
        preference="cell",
        cmap="inferno",
        clim=(0.0, upper),
        lighting=False,
        nan_opacity=0.0,
        pickable=True,
        show_edges=False,
        scalar_bar_args=_scalar_bar_args("2D in-plane WSS proxy [Pa]"),
        name="closed_wall_wss_proxy",
    )
    plotter.add_mesh(
        outline,
        color="#334155",
        line_width=1.0,
        lighting=False,
        pickable=False,
        show_scalar_bar=False,
        name="lumen_outline",
    )
    plotter.add_title(
        "Final accepted velocity field\n"
        f"2D in-plane wall-shear-stress proxy — closed-wall band (display clipped at P99.5 = {upper:.4g} Pa)",
        font_size=13,
        color="#111827",
    )

    spacing = float(image.spacing[0])
    nx, nz = int(image.dimensions[0] - 1), int(image.dimensions[1] - 1)
    domain = SimpleNamespace(
        shape=(nx, nz),
        spacing_um=spacing,
        x_coordinates_um=float(image.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * spacing,
        z_coordinates_um=float(image.origin[1]) + (np.arange(nz, dtype=float) + 0.5) * spacing,
    )
    _configure_planar_renderer(plotter, domain)
    plotter.reset_camera()
    plotter.camera.zoom(1.03)


def _ensure_wall_shear_display_mask(image) -> None:
    """为旧 VTI 补建显示掩膜；新 VTI 已直接保存该数组。"""

    if WALL_SHEAR_DISPLAY_MASK_ARRAY in image.cell_data:
        return
    lumen = np.asarray(image.cell_data[LUMEN_ARRAY], dtype=bool)
    wall = np.asarray(image.cell_data[WALL_MASK_ARRAY], dtype=bool)
    inlet = np.asarray(image.cell_data.get("inlet_label", np.zeros(image.n_cells)), dtype=int)
    outlet = np.asarray(image.cell_data.get("outlet_label", np.zeros(image.n_cells)), dtype=int)
    image.cell_data[WALL_SHEAR_DISPLAY_MASK_ARRAY] = (
        lumen & wall & (inlet == 0) & (outlet == 0)
    ).astype(np.uint8)


def _validate_wall_shear_image(image) -> None:
    required = (LUMEN_ARRAY, WALL_MASK_ARRAY, WALL_SHEAR_ARRAY)
    missing = [name for name in required if name not in image.cell_data]
    if missing:
        raise ValueError(
            "final_flow_field.vti is missing wall-shear arrays: "
            + ", ".join(missing)
            + ". Re-run generate_microbubble_trajectories.py with the current renderer."
        )
    stage = str(np.asarray(image.field_data.get("flow_stage", ["unknown"]))[0])
    if stage != "final":
        raise ValueError("Wall-shear visualization is only valid for the final accepted velocity field.")


def _artifact_paths(result_dir: Path) -> WallShearArtifacts:
    return WallShearArtifacts(
        html_path=Path(result_dir) / "final_wall_shear_stress_cfd.html",
        preview_path=Path(result_dir) / "final_wall_shear_stress_cfd.png",
    )
