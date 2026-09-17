"""为 PyVista/VTK CFD 原生产物提供带数值 hover 的 trame 网页查看器。"""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ..vtk.pyvista_flow import _new_cfd_plotter, _populate_cfd_plotter, validate_cfd_flow_dependencies
from ..vtk.vtk_flow_grid import (
    LUMEN_ARRAY,
    PRESSURE_ARRAY,
    SPEED_ARRAY,
    VELOCITY_ARRAY,
    WALL_SHEAR_ARRAY,
    WALL_SHEAR_DISPLAY_MASK_ARRAY,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve a PyVista/VTK CFD scene with a live numerical hover probe.",
    )
    parser.add_argument("--result-dir", type=Path, required=True, help="Directory containing the generated VTI/VTP files.")
    parser.add_argument("--stage", choices=("initial", "final"), default="final", help="Flow stage to inspect.")
    parser.add_argument(
        "--view",
        choices=("flow", "wall-shear"),
        default="flow",
        help="Render the flow dashboard or the final wall-shear-stress proxy.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="trame server host.")
    parser.add_argument("--port", type=int, default=8080, help="trame server port.")
    parser.add_argument("--open-browser", action="store_true", help="Open the probe page in the default browser.")
    return parser.parse_args()


def build_probe_plotter(result_dir: Path, *, stage: str, view: str = "flow"):
    """从原生 VTK 产物重建流场双面板或最终 WSS 单面板。"""

    if view == "wall-shear":
        if stage != "final":
            raise ValueError("The wall-shear-stress proxy is only available for --stage final.")
        from ..vtk.pyvista_wall_shear import build_wall_shear_plotter

        return build_wall_shear_plotter(result_dir)
    if view != "flow":
        raise ValueError("view must be either 'flow' or 'wall-shear'.")

    validate_cfd_flow_dependencies()
    import pyvista as pv

    result_dir = Path(result_dir)
    field_path = result_dir / f"{stage}_flow_field.vti"
    streamline_path = result_dir / f"{stage}_streamlines.vtp"
    if not field_path.is_file() or not streamline_path.is_file():
        raise ValueError(f"Missing {field_path.name} or {streamline_path.name} under {result_dir}.")

    image = pv.read(field_path)
    formal = pv.read(streamline_path)
    if not isinstance(image, pv.ImageData) or not isinstance(formal, pv.PolyData):
        raise ValueError("The CFD probe expects an ImageData VTI and a PolyData VTP.")
    if formal.n_lines == 0:
        raise ValueError("The formal streamline VTP contains no root-to-outlet line.")

    fluid = image.threshold((0.5, 1.5), scalars=LUMEN_ARRAY, preference="cell")
    fluid = fluid.cell_data_to_point_data(pass_cell_data=True)
    fluid.set_active_vectors(VELOCITY_ARRAY, preference="point")
    render_lines = _representative_streamlines(formal, paths_per_outlet=2)

    nx, nz = int(image.dimensions[0] - 1), int(image.dimensions[1] - 1)
    spacing = float(image.spacing[0])
    domain = SimpleNamespace(
        shape=(nx, nz),
        spacing_um=spacing,
        x_coordinates_um=float(image.origin[0]) + (np.arange(nx, dtype=float) + 0.5) * spacing,
        z_coordinates_um=float(image.origin[1]) + (np.arange(nz, dtype=float) + 0.5) * spacing,
    )
    stage_grid = SimpleNamespace(
        stage=stage,
        fluid_grid=fluid,
        speed_um_s=np.asarray(image.cell_data[SPEED_ARRAY]),
        pressure_projection=np.asarray(image.cell_data[PRESSURE_ARRAY]),
    )
    trace = SimpleNamespace(render_lines=render_lines)
    plotter = _new_cfd_plotter()
    try:
        _populate_cfd_plotter(plotter, domain, stage_grid, trace)
    except Exception:
        plotter.close()
        raise
    return plotter, image


def sample_native_field(image_grid, world_position: np.ndarray) -> dict[str, float | int | str | bool] | None:
    """把 hover 世界坐标映射回 cell-centered 物理场并读取数值。"""

    position = np.asarray(world_position, dtype=float)
    if position.shape[0] < 2 or not np.all(np.isfinite(position[:2])):
        return None
    nx, nz = int(image_grid.dimensions[0] - 1), int(image_grid.dimensions[1] - 1)
    spacing_x = float(image_grid.spacing[0])
    spacing_z = float(image_grid.spacing[1])
    ix = int(np.floor((position[0] - float(image_grid.origin[0])) / spacing_x))
    iz = int(np.floor((position[1] - float(image_grid.origin[1])) / spacing_z))
    if ix < 0 or ix >= nx or iz < 0 or iz >= nz:
        return None
    cell_id = ix + nx * iz
    if int(np.asarray(image_grid.cell_data[LUMEN_ARRAY])[cell_id]) <= 0:
        return None

    velocity = np.asarray(image_grid.cell_data[VELOCITY_ARRAY], dtype=float)[cell_id]
    pressure = float(np.asarray(image_grid.cell_data[PRESSURE_ARRAY], dtype=float)[cell_id])
    wall_shear = (
        float(np.asarray(image_grid.cell_data[WALL_SHEAR_ARRAY], dtype=float)[cell_id])
        if WALL_SHEAR_ARRAY in image_grid.cell_data
        else float("nan")
    )
    wall_shear_display_cell = (
        bool(np.asarray(image_grid.cell_data[WALL_SHEAR_DISPLAY_MASK_ARRAY], dtype=bool)[cell_id])
        if WALL_SHEAR_DISPLAY_MASK_ARRAY in image_grid.cell_data
        else False
    )
    return {
        "x_um": float(image_grid.origin[0] + (ix + 0.5) * spacing_x),
        "z_um": float(image_grid.origin[1] + (iz + 0.5) * spacing_z),
        "vx_um_s": float(velocity[0]),
        "vz_um_s": float(velocity[1]),
        "speed_um_s": float(np.asarray(image_grid.cell_data[SPEED_ARRAY], dtype=float)[cell_id]),
        "pressure": pressure,
        "wall_shear_stress_pa": wall_shear,
        "wall_shear_display_cell": wall_shear_display_cell,
        "vessel_id": int(np.asarray(image_grid.cell_data["vessel_id"])[cell_id]),
        "pressure_semantics": str(np.asarray(image_grid.field_data.get("pressure_semantics", ["unknown"]))[0]),
        "wall_shear_semantics": str(
            np.asarray(
                image_grid.field_data.get(
                    "wall_shear_semantics",
                    ["in_plane_wall_shear_stress_magnitude_proxy"],
                )
            )[0]
        ),
    }


def create_probe_server(result_dir: Path, *, stage: str, view: str = "flow"):
    """创建 trame 应用；hover 事件只传坐标，数值始终从原始 VTI 查询。"""

    from pyvista.trame.views import PyVistaLocalView
    from trame.app import get_server
    from trame.ui.vuetify3 import VAppLayout
    from trame.widgets import html, vuetify3

    plotter, image = build_probe_plotter(result_dir, stage=stage, view=view)
    server = get_server(name=f"microvascular_cfd_probe_{view}_{stage}", client_type="vue3")
    state = server.state
    # VtkLocalView 的场景通过 WebSocket 异步同步。若页面刚挂载就启用 hover，
    # 浏览器可能在 renderer 尚未建立时进入 HardwareSelector，触发
    # ``undefined.clear()`` 并让画布保持空白。只在 afterSceneLoaded 后启用拾取。
    state.picking_mode = None
    state.pick_data = None
    state.probe_visible = True
    state.probe_text = "Loading VTK scene..."

    def on_scene_loading(*_args, **_kwargs) -> None:
        state.picking_mode = None
        state.pick_data = None
        state.probe_text = "Loading VTK scene..."

    def on_scene_loaded(*_args, **_kwargs) -> None:
        state.picking_mode = "hover"
        state.probe_text = "Move the mouse over the fluid field to inspect values."

    @state.change("pick_data")
    def on_hover(pick_data, **_kwargs) -> None:
        values = None if not pick_data else sample_native_field(image, pick_data.get("worldPosition", []))
        if values is None:
            state.probe_text = "Solid / outside lumen"
            return
        if view == "wall-shear" and not bool(values["wall_shear_display_cell"]):
            state.probe_text = "Lumen interior / inlet or outlet open face\nNo closed-wall WSS value is displayed here."
            return
        probe_text = (
            f"X = {values['x_um']:.3f} um\n"
            f"Z = {values['z_um']:.3f} um\n"
            f"vx = {values['vx_um_s']:.6g} um/s\n"
            f"vz = {values['vz_um_s']:.6g} um/s\n"
            f"|v| = {values['speed_um_s']:.6g} um/s\n"
            f"pressure = {values['pressure']:.6g} a.u.\n"
            f"vessel_id = {values['vessel_id']}\n"
            f"pressure semantics = {values['pressure_semantics']}"
        )
        if bool(values["wall_shear_display_cell"]) and np.isfinite(float(values["wall_shear_stress_pa"])):
            probe_text += (
                f"\nwall shear stress proxy = {values['wall_shear_stress_pa']:.6g} Pa"
                f"\nwall shear semantics = {values['wall_shear_semantics']}"
            )
        state.probe_text = probe_text

    # 使用最小 VAppLayout，避免扫描环境中无关第三方包的 metadata；某些工程
    # 环境含非 UTF-8 egg-info，标准 SinglePageLayout 的版本提示会因此失败。
    with VAppLayout(server):
        with vuetify3.VMain(style="padding:0;"):
            with html.Div(style="position:relative;width:100%;height:100vh;overflow:hidden;"):
                PyVistaLocalView(
                    plotter,
                    enable_picking=True,
                    picking_modes=("[picking_mode]",),
                    before_scene_loaded=on_scene_loading,
                    after_scene_loaded=on_scene_loaded,
                    hover="pick_data = $event",
                    style="width:100%;height:100%;",
                )
                with vuetify3.VCard(
                    elevation=5,
                    style=(
                        "position:absolute;left:16px;top:16px;z-index:10;"
                        "background:rgba(255,255,255,.90);"
                    ),
                ):
                    html.Div(
                        (
                            "Final 2D wall-shear-stress proxy — live numerical probe"
                            if view == "wall-shear"
                            else f"{stage.capitalize()} microvascular CFD field — live numerical probe"
                        ),
                        style="padding:8px 12px;font:600 15px Arial,sans-serif;",
                    )
                with vuetify3.VCard(
                    elevation=8,
                    style=(
                        "position:absolute;right:16px;top:16px;z-index:10;"
                        "min-width:270px;background:rgba(255,255,255,.92);"
                    ),
                ):
                    vuetify3.VCardTitle("Field probe")
                    with vuetify3.VCardText():
                        html.Pre("{{ probe_text }}", style="margin:0;font:13px/1.45 Consolas,monospace;")

    server.controller.on_server_exited.add(lambda **_: plotter.close())
    # trame 持有回调，但显式保存强引用可避免调用方只保留 server 时 plotter 被回收。
    server._cfd_probe_plotter = plotter
    server._cfd_probe_image = image
    return server


def _representative_streamlines(formal, *, paths_per_outlet: int):
    labels = np.asarray(formal.cell_data["destination_outlet_id"], dtype=int)
    selected: list[int] = []
    for label in sorted(int(value) for value in np.unique(labels)):
        candidates = np.flatnonzero(labels == label)
        count = min(max(1, int(paths_per_outlet)), candidates.size)
        local = np.unique(np.rint(np.linspace(0, candidates.size - 1, count)).astype(int))
        selected.extend(int(candidates[index]) for index in local)
    return formal.extract_cells(selected).extract_surface(algorithm=None)


def main() -> None:
    args = parse_args()
    server = create_probe_server(args.result_dir, stage=args.stage, view=args.view)
    server.start(
        host=args.host,
        port=int(args.port),
        open_browser=bool(args.open_browser),
    )


if __name__ == "__main__":
    try:
        main()
    except (ImportError, ValueError) as exc:
        raise SystemExit(f"CFD probe viewer stopped: {exc}") from None
