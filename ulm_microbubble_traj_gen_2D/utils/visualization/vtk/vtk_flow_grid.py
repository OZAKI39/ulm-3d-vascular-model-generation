"""把二维 X-Z NumPy 流场转换成可供 VTK/PyVista 使用的物理网格。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from ...core.types import FlowField, GridDomain, RasterizedVessels

if TYPE_CHECKING:
    import pyvista as pv


FlowStage = Literal["initial", "final"]

VELOCITY_ARRAY = "velocity_xz_um_s"
SPEED_ARRAY = "speed_um_s"
PRESSURE_ARRAY = "pressure_projection"
WALL_SHEAR_ARRAY = "wall_shear_stress_pa"
LOCAL_SHEAR_ARRAY = "local_shear_stress_pa"
LUMEN_ARRAY = "lumen_mask"
LUMEN_FRACTION_ARRAY = "lumen_fraction"
WALL_MASK_ARRAY = "wall_mask"
WALL_SHEAR_DISPLAY_MASK_ARRAY = "wall_shear_display_mask"
LIC_ARRAY = "lic_intensity"


@dataclass(frozen=True)
class VtkStageGrid:
    """同一求解阶段的完整规则网格和已移除固体的流体网格。"""

    stage: FlowStage
    image_grid: "pv.ImageData"
    fluid_grid: "pv.UnstructuredGrid"
    velocity_xz_um_s: np.ndarray
    speed_um_s: np.ndarray
    pressure_projection: np.ndarray


def build_vtk_stage_grid(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    *,
    stage: FlowStage,
    include_lic: bool = True,
) -> VtkStageGrid:
    """
    构造一个 cell-centered ``vtkImageData``，再提取真实管腔单元。

    项目数组采用 ``(nx, nz)`` 顺序，而 VTK 采用 X 最快变化的线性顺序，
    因此所有 cell data 都必须按 Fortran 顺序展平。当前物理平面映射为
    ``VTK X = physical X``、``VTK Y = physical Z``、``VTK Z = 0``；第三个
    速度分量恒为零，仅用于满足 VTK 三分量向量接口，不虚构出平面流动。
    """

    pv = _require_pyvista()
    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    nx, nz = _validate_domain_and_mask(domain, lumen)
    velocity, speed, pressure = stage_numpy_fields(flow, lumen, stage=stage)
    spacing = float(domain.spacing_um)

    # 求解数组位于 cell center，所以 VTK 点网格需要在每个方向多一个点，
    # 且原点相对第一个 cell center 向外移动半格。
    image = pv.ImageData(
        dimensions=(nx + 1, nz + 1, 1),
        origin=(
            float(domain.x_coordinates_um[0]) - 0.5 * spacing,
            float(domain.z_coordinates_um[0]) - 0.5 * spacing,
            0.0,
        ),
        spacing=(spacing, spacing, spacing),
    )

    image.cell_data[LUMEN_ARRAY] = lumen.astype(np.uint8).ravel(order="F")
    lumen_fraction = np.asarray(raster.lumen_fraction, dtype=np.float32)
    if lumen_fraction.shape != lumen.shape:
        raise ValueError("lumen_fraction must match lumen_mask.shape.")
    image.cell_data[LUMEN_FRACTION_ARRAY] = lumen_fraction.ravel(order="F")
    wall_mask = np.asarray(raster.wall_mask, dtype=bool)
    if wall_mask.shape != lumen.shape:
        raise ValueError("wall_mask must match lumen_mask.shape.")
    image.cell_data[WALL_MASK_ARRAY] = wall_mask.astype(np.uint8).ravel(order="F")
    image.cell_data["vessel_id"] = np.asarray(raster.vessel_id, dtype=np.int32).ravel(order="F")
    image.cell_data["distance_to_wall_um"] = _masked_scalar(
        raster.distance_to_wall_um,
        lumen,
    ).ravel(order="F")
    if stage == "final":
        wall_shear = _masked_scalar(flow.wall_shear_stress_pa, lumen)
    else:
        # WSS 仅在最终接受速度场上由 wall_shear.py 计算；初始场保留同名数组
        # 但全部设为 NaN，避免下游把最终 WSS 错认为 initial WSS。
        wall_shear = np.full(lumen.shape, np.nan, dtype=float)
    image.cell_data[WALL_SHEAR_ARRAY] = wall_shear.ravel(order="F")
    saved_local_shear = getattr(flow, "local_shear_stress_pa", None)
    if stage == "final" and saved_local_shear is not None:
        local_shear = _masked_scalar(saved_local_shear, lumen)
    else:
        local_shear = np.full(lumen.shape, np.nan, dtype=float)
    image.cell_data[LOCAL_SHEAR_ARRAY] = local_shear.ravel(order="F")
    image.cell_data[SPEED_ARRAY] = speed.ravel(order="F")
    image.cell_data[PRESSURE_ARRAY] = pressure.ravel(order="F")
    image.cell_data[VELOCITY_ARRAY] = _vtk_velocity(velocity, lumen)

    inlet = _optional_integer_field(flow.inlet_label, lumen.shape)
    outlet = _optional_integer_field(flow.outlet_label, lumen.shape)
    image.cell_data["inlet_label"] = inlet.ravel(order="F")
    image.cell_data["outlet_label"] = outlet.ravel(order="F")
    # 形态学 wall_mask 也会覆盖开放入口/出口端帽；它们不是固体壁面，WSS 主图
    # 必须排除。原始 WSS 仍完整保留，display mask 只控制可视化有效带。
    wall_shear_display_mask = lumen & wall_mask & (inlet == 0) & (outlet == 0)
    image.cell_data[WALL_SHEAR_DISPLAY_MASK_ARRAY] = wall_shear_display_mask.astype(np.uint8).ravel(order="F")
    if getattr(flow, "boundary_normal_xz", None) is not None:
        boundary_normal = np.asarray(flow.boundary_normal_xz, dtype=float)
        if boundary_normal.shape != (*lumen.shape, 2):
            raise ValueError("boundary_normal_xz must have shape (nx, nz, 2).")
        image.cell_data["boundary_normal_xz"] = _vtk_velocity(boundary_normal, lumen)
    if getattr(flow, "boundary_velocity_xz_um_s", None) is not None:
        boundary_velocity = np.asarray(flow.boundary_velocity_xz_um_s, dtype=float)
        if boundary_velocity.shape != (*lumen.shape, 2):
            raise ValueError("boundary_velocity_xz_um_s must have shape (nx, nz, 2).")
        image.cell_data["boundary_velocity_xz_um_s"] = _vtk_velocity(boundary_velocity, lumen)

    junction = getattr(raster, "junction_core_mask", None)
    if junction is not None:
        image.cell_data["junction_core_mask"] = np.asarray(junction, dtype=np.uint8).ravel(order="F")

    # LIC 是预先烘焙出的普通灰度标量，最终 HTML 不需要特殊 GPU mapper。
    # 计算失败应显式暴露，而不是悄悄退回低质量方向显示。
    if include_lic:
        image.point_data[LIC_ARRAY] = compute_lic_intensity(image)

    # 先 threshold 再把 cell data 插值到 point data。若顺序相反，固体中的
    # NaN/零速度会污染边界点，StreamTracer 便可能错误穿墙或提前停下。
    fluid_cells = image.threshold((0.5, 1.5), scalars=LUMEN_ARRAY, preference="cell")
    fluid = fluid_cells.cell_data_to_point_data(pass_cell_data=True)
    fluid.set_active_vectors(VELOCITY_ARRAY, preference="point")

    image.field_data["flow_stage"] = np.asarray([stage])
    image.field_data["physical_plane"] = np.asarray(["X-Z mapped to VTK X-Y"])
    image.field_data["pressure_unit"] = np.asarray(
        ["not solved" if stage == "initial" else "mmHg"]
    )
    image.field_data["pressure_semantics"] = np.asarray(
        [
            "zero_reference_unsolved"
            if stage == "initial"
            else "DOLFINx gauge pressure relative to one pinned pressure degree of freedom"
        ]
    )
    image.field_data["wall_shear_unit"] = np.asarray(["Pa"])
    image.field_data["wall_shear_semantics"] = np.asarray(
        [
            "not_computed_for_initial_velocity"
            if stage == "initial"
            else "accepted_velocity_in_plane_wall_shear_stress_magnitude_proxy"
        ]
    )
    image.field_data["wall_shear_source_stage"] = np.asarray(
        ["not_available" if stage == "initial" else "final_accepted_velocity"]
    )
    image.field_data["local_shear_unit"] = np.asarray(["Pa"])
    image.field_data["local_shear_semantics"] = np.asarray(
        [
            "not_computed_for_initial_velocity"
            if stage == "initial"
            else "mu*sqrt(2*D:D)_whole_lumen_viscous_shear_magnitude"
        ]
    )
    return VtkStageGrid(
        stage=stage,
        image_grid=image,
        fluid_grid=fluid,
        velocity_xz_um_s=velocity,
        speed_um_s=speed,
        pressure_projection=pressure,
    )


def stage_numpy_fields(
    flow: FlowField,
    lumen_mask: np.ndarray,
    *,
    stage: FlowStage,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """选择初始或最终场，并把所有固体位置严格设置为 NaN。"""

    lumen = np.asarray(lumen_mask, dtype=bool)
    expected_velocity_shape = (*lumen.shape, 2)
    if stage == "initial":
        if flow.initial_velocity_xz_um_s is None:
            raise ValueError("Initial velocity is required for the initial PyVista CFD scene.")
        velocity = np.asarray(flow.initial_velocity_xz_um_s, dtype=float).copy()
        # 初始场尚未做压力投影，因此只能显示明确标注的零参考场，不能误用最终压力。
        pressure = np.zeros(lumen.shape, dtype=float)
    elif stage == "final":
        velocity = np.asarray(flow.velocity_xz_um_s, dtype=float).copy()
        if flow.pressure is None:
            raise ValueError("Final projection pressure is required for the converged PyVista CFD scene.")
        pressure = np.asarray(flow.pressure, dtype=float).copy()
    else:
        raise ValueError("stage must be either 'initial' or 'final'.")

    if velocity.shape != expected_velocity_shape:
        raise ValueError("Velocity must have shape (nx, nz, 2) with components ordered as X, Z.")
    if pressure.shape != lumen.shape:
        raise ValueError("Pressure must have shape (nx, nz).")

    speed = np.linalg.norm(velocity, axis=-1)
    velocity[~lumen] = np.nan
    speed[~lumen] = np.nan
    pressure[~lumen] = np.nan
    return velocity, speed, pressure


def compute_lic_intensity(image_grid: "pv.ImageData", *, steps: int = 16, step_size: float = 0.5) -> np.ndarray:
    """
    用 ``vtkImageDataLIC2D`` 烘焙全分辨率线积分卷积纹理。

    LIC 只表达局部流线轴向，不能区分正反方向，所以渲染时仍会叠加正式
    root-to-outlet Stream Tube 和稀疏箭头。固体速度在这个临时副本中设为零；
    真正显示时仍由 lumen threshold 完全移除固体单元。
    """

    pv = _require_pyvista()
    try:
        from vtkmodules.vtkRenderingCore import vtkRenderWindow
        from vtkmodules.vtkRenderingLICOpenGL2 import vtkImageDataLIC2D
    except ImportError as exc:  # pragma: no cover - 由依赖预检覆盖
        raise ImportError("VTK was built without the OpenGL2 LIC module required by the CFD renderer.") from exc

    source = image_grid.copy(deep=True)
    velocity = np.asarray(source.cell_data[VELOCITY_ARRAY], dtype=float).copy()
    velocity[~np.isfinite(velocity)] = 0.0
    source.cell_data[VELOCITY_ARRAY] = velocity
    point_source = source.cell_data_to_point_data(pass_cell_data=True)
    point_source.set_active_vectors(VELOCITY_ARRAY, preference="point")

    context = vtkRenderWindow()
    context.SetOffScreenRendering(1)
    lic_filter = vtkImageDataLIC2D()
    try:
        lic_filter.SetInputData(point_source)
        lic_filter.SetContext(context)
        lic_filter.SetSteps(max(1, int(steps)))
        lic_filter.SetStepSize(float(step_size))
        lic_filter.Update()
        # 必须在释放 OpenGL context 前深拷贝，否则输出仍引用 filter 内存。
        lic_output = pv.wrap(lic_filter.GetOutput()).copy(deep=True)
        lic_rgb = np.asarray(lic_output.point_data["LIC"], dtype=float)
        if lic_rgb.ndim != 2 or lic_rgb.shape[0] != image_grid.n_points:
            raise RuntimeError("VTK LIC returned an unexpected point-array shape.")
        intensity = np.asarray(lic_rgb[:, 0], dtype=np.float32).copy()
    finally:
        lic_filter.SetContext(None)
        context.Finalize()

    return np.clip(intensity, 0.0, 1.0)


def _validate_domain_and_mask(domain: GridDomain, lumen: np.ndarray) -> tuple[int, int]:
    if lumen.ndim != 2:
        raise ValueError("lumen_mask must be a two-dimensional (nx, nz) array.")
    nx, nz = (int(value) for value in lumen.shape)
    if tuple(domain.shape) != (nx, nz):
        raise ValueError("GridDomain.shape does not match lumen_mask.shape.")
    if float(domain.spacing_um) <= 0.0:
        raise ValueError("Grid spacing must be positive.")
    if np.asarray(domain.x_coordinates_um).shape != (nx,) or np.asarray(domain.z_coordinates_um).shape != (nz,):
        raise ValueError("Physical X/Z coordinate arrays do not match the domain shape.")
    if not np.any(lumen):
        raise ValueError("Cannot visualize an empty lumen mask.")
    return nx, nz


def _vtk_velocity(velocity_xz: np.ndarray, lumen: np.ndarray) -> np.ndarray:
    nx, nz = lumen.shape
    vectors = np.full((nx * nz, 3), np.nan, dtype=np.float32)
    vectors[:, 0] = np.asarray(velocity_xz[..., 0], dtype=np.float32).ravel(order="F")
    vectors[:, 1] = np.asarray(velocity_xz[..., 1], dtype=np.float32).ravel(order="F")
    vectors[:, 2] = 0.0
    return vectors


def _masked_scalar(values: np.ndarray, lumen: np.ndarray) -> np.ndarray:
    scalar = np.asarray(values, dtype=np.float32).copy()
    if scalar.shape != lumen.shape:
        raise ValueError("A scalar field does not match lumen_mask.shape.")
    scalar[~lumen] = np.nan
    return scalar


def _optional_integer_field(values: np.ndarray | None, shape: tuple[int, int]) -> np.ndarray:
    if values is None:
        return np.zeros(shape, dtype=np.int32)
    array = np.asarray(values, dtype=np.int32)
    if array.shape != shape:
        raise ValueError("A boundary-label field does not match lumen_mask.shape.")
    return array


def _require_pyvista():
    try:
        import pyvista as pv
    except ImportError as exc:  # pragma: no cover - 由 facade 给出完整安装提示
        raise ImportError("PyVista is required for VTK CFD visualization.") from exc
    return pv
