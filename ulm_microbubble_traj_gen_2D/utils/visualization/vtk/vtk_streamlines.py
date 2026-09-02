"""使用 VTK 生成并验证连续的 root-to-outlet 二维流线。"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Sequence

import numpy as np

from ...core.types import FlowField, GridDomain, RasterizedVessels
from .vtk_flow_grid import PRESSURE_ARRAY, SPEED_ARRAY, VELOCITY_ARRAY, VtkStageGrid

if TYPE_CHECKING:
    import pyvista as pv

    from ulm_vascular_model_generator.utils.core.models import Vessel

    from ...flow.flow_boundaries import BoundaryFluxFields


@dataclass(frozen=True)
class StreamlineSettings:
    """以物理网格间距定义、且不依赖整幅矩形尺寸的积分参数。"""

    integration_step_cells: float = 0.5
    max_length_factor: float = 2.2
    max_steps: int = 50_000
    boundary_tolerance_cells: float = 2.25
    minimum_direction_alignment: float = 0.98
    render_paths_per_outlet: int = 2


@dataclass(frozen=True)
class StreamlineTraceResult:
    """正式流线、失败诊断流线以及覆盖率审计结果。"""

    formal_lines: "pv.PolyData"
    diagnostic_lines: "pv.PolyData"
    render_lines: "pv.PolyData"
    report_rows: tuple[dict[str, object], ...]
    expected_outlet_labels: tuple[int, ...]
    reached_outlet_labels: tuple[int, ...]
    missing_outlet_labels: tuple[int, ...]
    expected_vessel_ids: tuple[int, ...]
    covered_vessel_ids: tuple[int, ...]
    missing_vessel_ids: tuple[int, ...]


class StreamlineContinuityError(ValueError):
    """无法构造满足出口与血管拓扑覆盖要求的正式流线。"""

    def __init__(self, message: str, result: StreamlineTraceResult) -> None:
        super().__init__(message)
        self.result = result


@dataclass
class _PathRecord:
    points: np.ndarray
    point_data: dict[str, np.ndarray]
    sampled_vessel_id: np.ndarray
    seed_id: int
    root_inlet_id: int
    intended_outlet_id: int
    destination_outlet_id: int
    intended_outlet_vessel_id: int
    destination_outlet_vessel_id: int
    source_vessel_id: int
    vtk_reason: int
    semantic_termination: str
    is_complete: bool
    arc_length_um: float
    start_distance_um: float
    end_distance_um: float
    direction_alignment_fraction: float
    sampled_vessel_ids: tuple[int, ...]
    coverage_vessel_ids: tuple[int, ...]
    failure_code: int
    path_id: int = field(default=-1)


def trace_root_to_outlets(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    stage_grid: VtkStageGrid,
    *,
    settings: StreamlineSettings | None = None,
    open_boundaries: "BoundaryFluxFields | None" = None,
    vessels: "Sequence[Vessel] | None" = None,
    require_complete_coverage: bool = True,
) -> StreamlineTraceResult:
    """
    为每个出口寻找 root 射击点，再从 root 重新向前积分正式流线。

    一般的均匀 root 播种可能完全错过小流量出口对应的极窄吸引域。这里先从
    每个 outlet cell 做 backward 诊断，得到它在 root 入口附近的候选脚点；
    随后使用同一个固定步长 RK4 从这些脚点执行全新的 forward 积分。反向线
    绝不进入正式渲染，只有 forward 结果命中原目标出口才被接纳。
    """

    cfg = settings or StreamlineSettings()
    inlet = _required_labels(flow.inlet_label, domain.shape, "inlet_label")
    outlet = _required_labels(flow.outlet_label, domain.shape, "outlet_label")
    expected_outlets = tuple(int(value) for value in np.unique(outlet[outlet > 0]))
    if not expected_outlets:
        raise ValueError("At least one labeled outlet is required for root-to-outlet streamlines.")
    if not np.any(inlet > 0):
        raise ValueError("At least one labeled root inlet is required for root-to-outlet streamlines.")

    spacing = float(domain.spacing_um)
    max_length = _maximum_integration_length(domain, cfg)
    # 每一个出口标记 cell 都是候选。保留整个候选池能解析分叉 separatrix，
    # 比每个出口只取一个中心点稳定得多。
    outlet_source_points, outlet_source_labels = _boundary_cell_centers(domain, outlet)
    backward = _run_stream_tracer(
        stage_grid.fluid_grid,
        outlet_source_points,
        direction="backward",
        speed=stage_grid.speed_um_s,
        spacing_um=spacing,
        max_length_um=max_length,
        settings=cfg,
    )
    shooting_points, shooting_targets, shooting_roots = _root_shooting_candidates(
        backward,
        outlet_source_labels,
        domain,
        raster,
        inlet,
        cfg,
    )
    if shooting_points.size == 0:
        raise ValueError("Backward VTK diagnostics found no valid root shooting point.")

    forward = _run_stream_tracer(
        stage_grid.fluid_grid,
        shooting_points,
        direction="forward",
        speed=stage_grid.speed_um_s,
        spacing_um=spacing,
        max_length_um=max_length,
        settings=cfg,
    )
    outlet_vessels = _outlet_vessel_map(outlet, raster.vessel_id, open_boundaries)
    ancestor_chains = _ancestor_chains_by_outlet(outlet_vessels, vessels)
    observable_vessels = _globally_observable_vessel_ids(raster, domain)
    records = _classify_forward_paths(
        forward,
        shooting_targets,
        shooting_roots,
        domain,
        raster,
        inlet,
        outlet,
        outlet_vessels,
        getattr(flow, "boundary_normal_xz", None),
        open_boundaries,
        ancestor_chains,
        observable_vessels,
        stage_grid.fluid_grid,
        cfg,
    )

    formal_records = [record for record in records if record.is_complete]
    diagnostic_records = [record for record in records if not record.is_complete]
    for path_id, record in enumerate(formal_records):
        record.path_id = path_id
    for offset, record in enumerate(diagnostic_records, start=len(formal_records)):
        record.path_id = offset

    reached_outlets = tuple(sorted({record.destination_outlet_id for record in formal_records}))
    missing_outlets = tuple(sorted(set(expected_outlets) - set(reached_outlets)))
    if vessels is None:
        expected_vessels = tuple(
            int(value)
            for value in np.unique(np.asarray(raster.vessel_id)[np.asarray(raster.lumen_mask, dtype=bool)])
            if int(value) >= 0
        )
    else:
        expected_vessels = tuple(sorted(int(vessel.vid) for vessel in vessels))
    directly_covered_vessels = {
        value for record in formal_records for value in record.coverage_vessel_ids
    }
    if vessels is not None and ancestor_chains:
        # 只有完全被 junction core 吞没、无法由非 junction 栅格直接观测的极短祖先段，
        # 才允许借助权威父子图补记覆盖；其余血管必须真实、按序出现在正式流线上。
        graph_inferred_vessels = {
            vessel_id
            for outlet_id in reached_outlets
            for vessel_id in ancestor_chains.get(outlet_id, ())
            if vessel_id not in observable_vessels
        }
        directly_covered_vessels.update(graph_inferred_vessels)
    covered_vessels = tuple(sorted(directly_covered_vessels))
    missing_vessels = tuple(sorted(set(expected_vessels) - set(covered_vessels)))

    render_records = _select_render_records(
        formal_records,
        inlet,
        domain,
        expected_vessels,
        cfg.render_paths_per_outlet,
    )
    formal_poly = _records_to_polydata(formal_records)
    diagnostic_poly = _records_to_polydata(diagnostic_records)
    render_poly = _records_to_polydata(render_records)
    report_rows = tuple(_record_to_report_row(record, stage_grid.stage) for record in [*formal_records, *diagnostic_records])

    for mesh in (formal_poly, diagnostic_poly, render_poly):
        mesh.field_data["expected_outlet_labels"] = np.asarray(expected_outlets, dtype=np.int32)
        mesh.field_data["missing_outlet_labels"] = np.asarray(missing_outlets, dtype=np.int32)
        mesh.field_data["expected_vessel_ids"] = np.asarray(expected_vessels, dtype=np.int32)
        mesh.field_data["missing_vessel_ids"] = np.asarray(missing_vessels, dtype=np.int32)

    result = StreamlineTraceResult(
        formal_lines=formal_poly,
        diagnostic_lines=diagnostic_poly,
        render_lines=render_poly,
        report_rows=report_rows,
        expected_outlet_labels=expected_outlets,
        reached_outlet_labels=reached_outlets,
        missing_outlet_labels=missing_outlets,
        expected_vessel_ids=expected_vessels,
        covered_vessel_ids=covered_vessels,
        missing_vessel_ids=missing_vessels,
    )
    if require_complete_coverage and (missing_outlets or missing_vessels):
        raise StreamlineContinuityError(
            "VTK could not build complete root-to-outlet streamlines for "
            f"outlets={list(missing_outlets)} or vessel_ids={list(missing_vessels)}.",
            result,
        )
    return result


def write_streamline_continuity_csv(path: Path, result: StreamlineTraceResult) -> Path:
    """保存每条 forward 路径的入口、出口、连续性和终止原因。"""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "path_id",
        "stage",
        "seed_id",
        "root_inlet_id",
        "intended_outlet_id",
        "destination_outlet_id",
        "intended_outlet_vessel_id",
        "destination_outlet_vessel_id",
        "source_vessel_id",
        "is_complete_root_to_outlet",
        "semantic_termination",
        "vtk_reason_for_termination",
        "n_points",
        "arc_length_um",
        "start_distance_to_root_um",
        "end_distance_to_outlet_um",
        "direction_alignment_fraction",
        "sampled_vessel_ids",
        "coverage_vessel_ids",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(result.report_rows)
    return path


def _run_stream_tracer(
    fluid_grid: "pv.UnstructuredGrid",
    source_points: np.ndarray,
    *,
    direction: str,
    speed: np.ndarray,
    spacing_um: float,
    max_length_um: float,
    settings: StreamlineSettings,
) -> "pv.PolyData":
    pv = _require_pyvista()
    points = np.asarray(source_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Streamline source points must have shape (n, 3).")
    if points.shape[0] == 0:
        return pv.PolyData()
    source = pv.PolyData(points)
    finite_speed = np.asarray(speed, dtype=float)
    finite_speed = finite_speed[np.isfinite(finite_speed)]
    terminal_speed = max(float(np.max(finite_speed)) * 1.0e-10, np.finfo(float).eps)
    step_um = float(settings.integration_step_cells) * float(spacing_um)
    if step_um <= 0.0:
        raise ValueError("Streamline integration step must be positive.")

    # 固定步长 RK4 在分叉 separatrix 上具有可重复的正反向行为；自适应 RK45
    # 会因两次误差控制不同而把 root 射击点送进相邻支路。
    return fluid_grid.streamlines_from_source(
        source,
        vectors=VELOCITY_ARRAY,
        integrator_type=4,
        integration_direction=direction,
        surface_streamlines=True,
        initial_step_length=step_um,
        step_unit="l",
        min_step_length=step_um,
        max_step_length=step_um,
        max_steps=int(settings.max_steps),
        terminal_speed=terminal_speed,
        max_time=None,
        compute_vorticity=False,
        interpolator_type="point",
        max_length=float(max_length_um),
    )


def _root_shooting_candidates(
    backward: "pv.PolyData",
    source_outlet_labels: np.ndarray,
    domain: GridDomain,
    raster: RasterizedVessels,
    inlet_label: np.ndarray,
    settings: StreamlineSettings,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    inlet_points, inlet_values = _boundary_cell_centers(domain, inlet_label, keep_z=False)
    tolerance = float(settings.boundary_tolerance_cells) * float(domain.spacing_um)
    candidates: list[np.ndarray] = []
    targets: list[int] = []
    roots: list[int] = []
    seen: set[tuple[int, int, int]] = set()

    for _, seed_id, _, point_ids in _polyline_cells(backward):
        if seed_id < 0 or seed_id >= source_outlet_labels.size or point_ids.size < 3:
            continue
        points = np.asarray(backward.points[point_ids], dtype=float)
        if not np.all(np.isfinite(points)):
            continue
        root_label, root_distance = _nearest_boundary(points[-1, :2], inlet_points, inlet_values)
        if root_label <= 0 or root_distance > tolerance:
            continue

        # 从仍属于 root 邻域的反向轨迹尾部选择最靠内的原始 VTK 点。相比固定取
        # 倒数第二点，这会避开开放面上的插值混合，同时不横向移动流线或伪造种子。
        candidate = _stable_root_shooting_point(
            points,
            int(root_label),
            inlet_points,
            inlet_values,
            domain,
            raster.lumen_mask,
            tolerance,
        )
        if candidate is None:
            continue
        target = int(source_outlet_labels[seed_id])
        quantization = max(0.02 * float(domain.spacing_um), np.finfo(float).eps)
        key = (target, int(round(candidate[0] / quantization)), int(round(candidate[1] / quantization)))
        if key in seen:
            continue
        seen.add(key)
        candidates.append(candidate)
        targets.append(target)
        roots.append(int(root_label))

    if not candidates:
        return (
            np.empty((0, 3), dtype=float),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
        )
    return (
        np.asarray(candidates, dtype=float),
        np.asarray(targets, dtype=np.int32),
        np.asarray(roots, dtype=np.int32),
    )


def _stable_root_shooting_point(
    backward_points: np.ndarray,
    root_label: int,
    inlet_points: np.ndarray,
    inlet_values: np.ndarray,
    domain: GridDomain,
    lumen_mask: np.ndarray,
    tolerance_um: float,
) -> np.ndarray | None:
    """Select the deepest original backward point that remains in the root neighborhood."""

    points = np.asarray(backward_points, dtype=float)
    maximum_depth = min(float(tolerance_um), float(domain.spacing_um))
    selected: np.ndarray | None = None
    fallback: np.ndarray | None = None
    for index in range(points.shape[0] - 2, -1, -1):
        point = points[index]
        label, distance = _nearest_boundary(point[:2], inlet_points, inlet_values)
        valid_root_point = (
            label == int(root_label)
            and distance <= float(tolerance_um)
            and _nearest_cell_is_lumen(point[:2], domain, lumen_mask)
        )
        if valid_root_point and fallback is None:
            fallback = point.copy()
        if valid_root_point and distance <= maximum_depth:
            selected = point.copy()
            continue
        if selected is not None or fallback is not None:
            break
    return selected if selected is not None else fallback


def _classify_forward_paths(
    forward: "pv.PolyData",
    intended_outlets: np.ndarray,
    intended_roots: np.ndarray,
    domain: GridDomain,
    raster: RasterizedVessels,
    inlet_label: np.ndarray,
    outlet_label: np.ndarray,
    outlet_vessels: dict[int, int],
    boundary_normal_xz: np.ndarray | None,
    open_boundaries: "BoundaryFluxFields | None",
    ancestor_chains: dict[int, tuple[int, ...]],
    observable_vessels: frozenset[int],
    fluid_grid: "pv.UnstructuredGrid",
    settings: StreamlineSettings,
) -> list[_PathRecord]:
    inlet_points, inlet_values = _boundary_cell_centers(domain, inlet_label, keep_z=False)
    outlet_points, outlet_values = _boundary_cell_centers(domain, outlet_label, keep_z=False)
    tolerance = float(settings.boundary_tolerance_cells) * float(domain.spacing_um)
    maximum_segment = 1.75 * float(settings.integration_step_cells) * float(domain.spacing_um)
    records: list[_PathRecord] = []

    for _, seed_id, vtk_reason, point_ids in _polyline_cells(forward):
        if seed_id < 0 or seed_id >= intended_outlets.size or point_ids.size < 2:
            continue
        points = np.asarray(forward.points[point_ids], dtype=float).copy()
        intended_outlet = int(intended_outlets[seed_id])
        intended_root = int(intended_roots[seed_id])
        point_data = _copy_path_point_data(forward, point_ids)
        velocity = np.asarray(point_data.get(VELOCITY_ARRAY, np.empty((0, 3))), dtype=float)
        terminal_point = _target_open_face_intersection(
            points[-2:, :2],
            velocity[-1, :2] if velocity.shape[0] else np.empty(0),
            intended_outlet,
            outlet_label,
            raster.lumen_mask,
            boundary_normal_xz,
            domain,
            open_boundaries,
            maximum_extension_um=1.10 * float(settings.integration_step_cells) * float(domain.spacing_um),
        )
        open_face_hit = terminal_point is not None
        if terminal_point is not None and np.linalg.norm(terminal_point - points[-1, :2]) > 1.0e-8:
            old_endpoint = points[-1].copy()
            points = np.vstack((points, (float(terminal_point[0]), float(terminal_point[1]), 0.0)))
            point_data = _append_terminal_point_data(point_data, old_endpoint, points[-1])
            velocity = np.asarray(point_data.get(VELOCITY_ARRAY, np.empty((0, 3))), dtype=float)

        finite = bool(np.all(np.isfinite(points)))
        deltas = np.diff(points[:, :2], axis=0) if finite else np.empty((0, 2), dtype=float)
        segment_lengths = np.linalg.norm(deltas, axis=1) if deltas.size else np.empty(0, dtype=float)
        arc_length = float(np.sum(segment_lengths))

        root_label, start_distance = _nearest_boundary(points[0, :2], inlet_points, inlet_values)
        nearest_outlet, end_distance = _nearest_boundary(points[-1, :2], outlet_points, outlet_values)
        actual_outlet = int(nearest_outlet) if end_distance <= tolerance else 0
        alignment = _direction_alignment_fraction(deltas, velocity)
        sampled_vessel_id = _sample_integer_grid(points[:, :2], domain, raster.vessel_id, invalid=-1)
        sampled_vessels = tuple(int(value) for value in _compress_ids(sampled_vessel_id) if int(value) >= 0)
        coverage_vessels = _ordered_covered_vessels_outside_junctions(
            points[:, :2],
            domain,
            raster,
        )
        source_vessel = next((value for value in sampled_vessels if value >= 0), -1)
        expected_chain = ancestor_chains.get(intended_outlet)
        if expected_chain is None:
            topology_valid = True
        else:
            expected_observable_chain = tuple(
                vessel_id for vessel_id in expected_chain if vessel_id in observable_vessels
            )
            # Overlapping rasterized cylinders can briefly label a continuous path as a
            # neighboring sibling vessel outside the junction mask. Project the evidence
            # onto the intended ancestor chain, then require every observable ancestor in
            # its exact order. Wrong outlets and discontinuous paths remain separate hard
            # failures, so incidental sibling labels cannot manufacture endpoint coverage.
            matched_expected_chain = _compress_ids(
                np.asarray(
                    [
                        vessel_id
                        for vessel_id in coverage_vessels
                        if vessel_id in expected_chain
                    ],
                    dtype=np.int32,
                )
            )
            matched_observable_chain = tuple(
                vessel_id
                for vessel_id in matched_expected_chain
                if vessel_id in observable_vessels
            )
            topology_valid = (
                _is_ordered_subsequence(matched_expected_chain, expected_chain)
                and matched_observable_chain == expected_observable_chain
            )

        # StreamTracer 在 threshold 后的流体单元上运行，理论上中间点都应有
        # containing cell。首尾点可恰好落在开放边界面，因此只检查内部点。
        if points.shape[0] > 2:
            containing = np.asarray(fluid_grid.find_containing_cell(points[1:-1]), dtype=int)
            inside_fluid = bool(np.all(containing >= 0))
        else:
            inside_fluid = False

        checks = {
            "finite": finite,
            "root": root_label == intended_root and start_distance <= tolerance,
            "outlet": actual_outlet == intended_outlet and end_distance <= tolerance,
            "vtk_exit": int(vtk_reason) == 1,
            "segment": bool(segment_lengths.size and np.max(segment_lengths) <= maximum_segment),
            "direction": alignment >= float(settings.minimum_direction_alignment),
            "open_face": open_face_hit,
            "topology": topology_valid,
            "inside": inside_fluid,
        }
        complete = all(checks.values())
        semantic, failure_code = _semantic_termination(checks, actual_outlet)
        records.append(
            _PathRecord(
                points=points,
                point_data=point_data,
                sampled_vessel_id=sampled_vessel_id,
                seed_id=int(seed_id),
                root_inlet_id=int(root_label),
                intended_outlet_id=intended_outlet,
                destination_outlet_id=int(actual_outlet),
                intended_outlet_vessel_id=int(outlet_vessels.get(intended_outlet, -1)),
                destination_outlet_vessel_id=int(outlet_vessels.get(int(actual_outlet), -1)),
                source_vessel_id=int(source_vessel),
                vtk_reason=int(vtk_reason),
                semantic_termination=semantic,
                is_complete=complete,
                arc_length_um=arc_length,
                start_distance_um=float(start_distance),
                end_distance_um=float(end_distance),
                direction_alignment_fraction=float(alignment),
                sampled_vessel_ids=sampled_vessels,
                coverage_vessel_ids=coverage_vessels,
                failure_code=int(failure_code),
            )
        )
    return records


def _select_render_records(
    records: list[_PathRecord],
    inlet_label: np.ndarray,
    domain: GridDomain,
    expected_vessel_ids: tuple[int, ...],
    paths_per_outlet: int,
) -> list[_PathRecord]:
    if not records:
        return []
    inlet_points, _ = _boundary_cell_centers(domain, inlet_label, keep_z=False)
    centered = inlet_points - np.mean(inlet_points, axis=0, keepdims=True)
    if inlet_points.shape[0] >= 2:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        tangent = vh[0]
    else:
        tangent = np.asarray([1.0, 0.0])

    selected: list[_PathRecord] = []
    for outlet in sorted({record.destination_outlet_id for record in records}):
        group = [record for record in records if record.destination_outlet_id == outlet]
        group.sort(key=lambda record: float(np.dot(record.points[0, :2], tangent)))
        count = min(max(1, int(paths_per_outlet)), len(group))
        indices = np.unique(np.rint(np.linspace(0, len(group) - 1, count)).astype(int))
        selected.extend(group[int(index)] for index in indices)

    # 若每出口的代表线仍漏掉某个栅格 vessel_id，则贪心加入覆盖该段的最短
    # 完整路径。这样可视化 LOD 不会为了限流而删除某段的最后一条拓扑路径。
    selected_ids = {id(record) for record in selected}
    covered = {value for record in selected for value in record.coverage_vessel_ids}
    for vessel_id in sorted(set(expected_vessel_ids) - covered):
        candidates = [record for record in records if vessel_id in record.coverage_vessel_ids]
        if not candidates:
            continue
        candidate = min(candidates, key=lambda record: record.arc_length_um)
        if id(candidate) not in selected_ids:
            selected.append(candidate)
            selected_ids.add(id(candidate))
            covered.update(candidate.coverage_vessel_ids)
    return selected


def _records_to_polydata(records: Iterable[_PathRecord]) -> "pv.PolyData":
    pv = _require_pyvista()
    items = list(records)
    if not items:
        empty = pv.PolyData()
        empty.field_data["note"] = np.asarray(["No paths in this category."])
        return empty

    points: list[np.ndarray] = []
    lines: list[np.ndarray] = []
    offset = 0
    for record in items:
        count = int(record.points.shape[0])
        points.append(record.points)
        lines.append(np.concatenate(([count], np.arange(offset, offset + count, dtype=np.int64))))
        offset += count
    mesh = pv.PolyData(np.vstack(points), lines=np.concatenate(lines))

    common_point_arrays = [VELOCITY_ARRAY, SPEED_ARRAY, PRESSURE_ARRAY, "IntegrationTime", "distance_to_wall_um"]
    for name in common_point_arrays:
        if all(name in record.point_data for record in items):
            mesh.point_data[name] = np.concatenate([record.point_data[name] for record in items], axis=0)
    mesh.point_data["sampled_vessel_id"] = np.concatenate([record.sampled_vessel_id for record in items]).astype(
        np.int32,
        copy=False,
    )
    mesh.point_data["path_id"] = np.concatenate(
        [np.full(record.points.shape[0], record.path_id, dtype=np.int32) for record in items]
    )
    mesh.point_data["path_progress"] = np.concatenate([_normalized_arc_progress(record.points) for record in items])

    mesh.cell_data["path_id"] = np.asarray([record.path_id for record in items], dtype=np.int32)
    mesh.cell_data["SeedId"] = np.asarray([record.seed_id for record in items], dtype=np.int32)
    mesh.cell_data["root_inlet_id"] = np.asarray([record.root_inlet_id for record in items], dtype=np.int32)
    mesh.cell_data["destination_outlet_id"] = np.asarray(
        [record.destination_outlet_id for record in items], dtype=np.int32
    )
    mesh.cell_data["intended_outlet_vessel_id"] = np.asarray(
        [record.intended_outlet_vessel_id for record in items], dtype=np.int32
    )
    mesh.cell_data["destination_outlet_vessel_id"] = np.asarray(
        [record.destination_outlet_vessel_id for record in items], dtype=np.int32
    )
    mesh.cell_data["source_vessel_id"] = np.asarray([record.source_vessel_id for record in items], dtype=np.int32)
    mesh.cell_data["ReasonForTermination"] = np.asarray([record.vtk_reason for record in items], dtype=np.int32)
    mesh.cell_data["semantic_termination_code"] = np.asarray(
        [1 if record.is_complete else record.failure_code for record in items], dtype=np.int32
    )
    mesh.cell_data["is_complete_root_to_outlet"] = np.asarray(
        [record.is_complete for record in items], dtype=np.uint8
    )
    mesh.field_data["semantic_termination_code_1"] = np.asarray(["OUTLET_REACHED"])
    mesh.field_data["trace_provenance"] = np.asarray(["outlet-backward seed discovery; root-forward RK4 verification"])
    return mesh


def _record_to_report_row(record: _PathRecord, stage: str) -> dict[str, object]:
    return {
        "path_id": int(record.path_id),
        "stage": stage,
        "seed_id": int(record.seed_id),
        "root_inlet_id": int(record.root_inlet_id),
        "intended_outlet_id": int(record.intended_outlet_id),
        "destination_outlet_id": int(record.destination_outlet_id),
        "intended_outlet_vessel_id": int(record.intended_outlet_vessel_id),
        "destination_outlet_vessel_id": int(record.destination_outlet_vessel_id),
        "source_vessel_id": int(record.source_vessel_id),
        "is_complete_root_to_outlet": int(record.is_complete),
        "semantic_termination": record.semantic_termination,
        "vtk_reason_for_termination": int(record.vtk_reason),
        "n_points": int(record.points.shape[0]),
        "arc_length_um": f"{record.arc_length_um:.9g}",
        "start_distance_to_root_um": f"{record.start_distance_um:.9g}",
        "end_distance_to_outlet_um": f"{record.end_distance_um:.9g}",
        "direction_alignment_fraction": f"{record.direction_alignment_fraction:.9g}",
        "sampled_vessel_ids": ";".join(str(value) for value in record.sampled_vessel_ids),
        "coverage_vessel_ids": ";".join(str(value) for value in record.coverage_vessel_ids),
    }


def _polyline_cells(mesh: "pv.PolyData"):
    if mesh.n_lines == 0:
        return
    if mesh.n_cells != mesh.n_lines:
        raise ValueError("VTK StreamTracer returned non-line cells unexpectedly.")
    if "SeedIds" not in mesh.cell_data:
        raise ValueError("VTK StreamTracer output is missing the SeedIds cell array.")
    seed_ids = np.asarray(mesh.cell_data["SeedIds"], dtype=int)
    reasons = np.asarray(mesh.cell_data.get("ReasonForTermination", np.full(mesh.n_lines, -1)), dtype=int)
    connectivity = np.asarray(mesh.lines, dtype=np.int64)
    cursor = 0
    for cell_id in range(mesh.n_lines):
        count = int(connectivity[cursor])
        point_ids = connectivity[cursor + 1 : cursor + 1 + count]
        cursor += count + 1
        yield cell_id, int(seed_ids[cell_id]), int(reasons[cell_id]), point_ids


def _copy_path_point_data(mesh: "pv.PolyData", point_ids: np.ndarray) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for name in (VELOCITY_ARRAY, SPEED_ARRAY, PRESSURE_ARRAY, "IntegrationTime", "distance_to_wall_um"):
        if name in mesh.point_data:
            result[name] = np.asarray(mesh.point_data[name][point_ids]).copy()
    return result


def _append_terminal_point_data(
    point_data: dict[str, np.ndarray],
    old_endpoint: np.ndarray,
    new_endpoint: np.ndarray,
) -> dict[str, np.ndarray]:
    """为不足一个积分步的开放面交点复制局部场值，并延长积分时间。"""

    result: dict[str, np.ndarray] = {}
    distance = float(np.linalg.norm(np.asarray(new_endpoint) - np.asarray(old_endpoint)))
    speed_values = np.asarray(point_data.get(SPEED_ARRAY, np.empty(0)), dtype=float)
    terminal_speed = float(speed_values[-1]) if speed_values.size else 0.0
    for name, raw in point_data.items():
        array = np.asarray(raw)
        terminal = array[-1:].copy()
        if name == "IntegrationTime" and terminal_speed > np.finfo(float).eps:
            terminal = terminal + distance / terminal_speed
        result[name] = np.concatenate((array, terminal), axis=0)
    return result


def _boundary_cell_centers(
    domain: GridDomain,
    labels: np.ndarray,
    *,
    keep_z: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    indices = np.argwhere(np.asarray(labels, dtype=int) > 0)
    if indices.size == 0:
        width = 3 if keep_z else 2
        return np.empty((0, width), dtype=float), np.empty(0, dtype=np.int32)
    x = np.asarray(domain.x_coordinates_um, dtype=float)[indices[:, 0]]
    z = np.asarray(domain.z_coordinates_um, dtype=float)[indices[:, 1]]
    points_2d = np.column_stack((x, z))
    values = np.asarray(labels, dtype=np.int32)[indices[:, 0], indices[:, 1]]
    if keep_z:
        return np.column_stack((points_2d, np.zeros(points_2d.shape[0], dtype=float))), values
    return points_2d, values


def _nearest_boundary(
    point_xz: np.ndarray,
    boundary_points: np.ndarray,
    boundary_values: np.ndarray,
) -> tuple[int, float]:
    if boundary_points.size == 0:
        return 0, float("inf")
    distance = np.linalg.norm(boundary_points - np.asarray(point_xz, dtype=float)[None, :], axis=1)
    index = int(np.argmin(distance))
    return int(boundary_values[index]), float(distance[index])


def _sample_integer_grid(
    points_xz: np.ndarray,
    domain: GridDomain,
    values: np.ndarray,
    *,
    invalid: int,
) -> np.ndarray:
    points = np.asarray(points_xz, dtype=float)
    spacing = float(domain.spacing_um)
    ix = np.rint((points[:, 0] - float(domain.x_coordinates_um[0])) / spacing).astype(int)
    iz = np.rint((points[:, 1] - float(domain.z_coordinates_um[0])) / spacing).astype(int)
    result = np.full(points.shape[0], int(invalid), dtype=np.int32)
    inside = (ix >= 0) & (ix < domain.shape[0]) & (iz >= 0) & (iz < domain.shape[1])
    array = np.asarray(values)
    result[inside] = np.asarray(array[ix[inside], iz[inside]], dtype=np.int32)
    return result


def _nearest_cell_is_lumen(point_xz: np.ndarray, domain: GridDomain, lumen_mask: np.ndarray) -> bool:
    return bool(_sample_integer_grid(np.asarray(point_xz, dtype=float)[None, :], domain, lumen_mask, invalid=0)[0])


def _direction_alignment_fraction(deltas: np.ndarray, velocity: np.ndarray) -> float:
    if deltas.shape[0] == 0 or velocity.ndim != 2 or velocity.shape[0] != deltas.shape[0] + 1:
        return 0.0
    midpoint_velocity = 0.5 * (velocity[:-1, :2] + velocity[1:, :2])
    delta_norm = np.linalg.norm(deltas, axis=1)
    velocity_norm = np.linalg.norm(midpoint_velocity, axis=1)
    valid = np.isfinite(delta_norm) & np.isfinite(velocity_norm) & (delta_norm > 0.0) & (velocity_norm > 0.0)
    if not np.any(valid):
        return 0.0
    dot = np.einsum("ij,ij->i", deltas[valid], midpoint_velocity[valid])
    tolerance = -1.0e-8 * delta_norm[valid] * velocity_norm[valid]
    return float(np.mean(dot >= tolerance))


def _semantic_termination(checks: dict[str, bool], actual_outlet: int) -> tuple[str, int]:
    if all(checks.values()):
        return "OUTLET_REACHED", 1
    order = [
        ("finite", "NONFINITE_PATH", -1),
        ("root", "ROOT_NOT_REACHED", -2),
        ("inside", "LEFT_LUMEN_BEFORE_OUTLET", -3),
        ("segment", "DISCONTINUOUS_STEP", -4),
        ("direction", "OPPOSES_FORWARD_VELOCITY", -5),
        ("open_face", "TARGET_OPEN_FACE_NOT_CROSSED", -8),
        ("topology", "VESSEL_TOPOLOGY_MISMATCH", -9),
        ("vtk_exit", "VTK_ABNORMAL_TERMINATION", -6),
        ("outlet", "WRONG_OUTLET" if actual_outlet > 0 else "NO_OUTLET_REACHED", -7),
    ]
    for key, label, code in order:
        if not checks[key]:
            return label, code
    return "UNKNOWN_FAILURE", -99


def _compress_ids(values: np.ndarray) -> tuple[int, ...]:
    sequence: list[int] = []
    for raw in np.asarray(values, dtype=int):
        value = int(raw)
        if not sequence or sequence[-1] != value:
            sequence.append(value)
    return tuple(sequence)


def _target_open_face_intersection(
    final_segment_xz: np.ndarray,
    endpoint_velocity_xz: np.ndarray,
    intended_outlet: int,
    outlet_label: np.ndarray,
    lumen_mask: np.ndarray,
    boundary_normal_xz: np.ndarray | None,
    domain: GridDomain,
    open_boundaries: "BoundaryFluxFields | None",
    *,
    maximum_extension_um: float,
) -> np.ndarray | None:
    """沿最后的 forward 切向补到权威开放面，并返回精确交点。"""

    segment = np.asarray(final_segment_xz, dtype=float)
    velocity = np.asarray(endpoint_velocity_xz, dtype=float)
    if segment.shape != (2, 2) or velocity.shape != (2,):
        return None
    if not np.all(np.isfinite(segment)) or not np.all(np.isfinite(velocity)):
        return None
    faces = _target_open_faces(
        intended_outlet,
        outlet_label,
        lumen_mask,
        boundary_normal_xz,
        domain,
        open_boundaries,
    )
    if not faces:
        return None

    previous, end = segment
    delta = end - previous
    direction_norm = float(np.linalg.norm(delta))
    if direction_norm <= np.finfo(float).eps:
        return None
    direction = delta / direction_norm
    spacing = float(domain.spacing_um)
    candidates: list[tuple[float, np.ndarray]] = []
    for center, normal, length in faces:
        normal = np.asarray(normal, dtype=float)
        normal /= max(float(np.linalg.norm(normal)), np.finfo(float).eps)
        segment_outward = float(np.dot(direction, normal))
        velocity_outward = float(np.dot(velocity, normal))
        if segment_outward <= 1.0e-8:
            continue
        if velocity_outward <= 1.0e-10 * max(float(np.linalg.norm(velocity)), 1.0):
            continue
        distance = float(np.dot(center - end, normal) / segment_outward)
        if distance < -0.05 * spacing or distance > float(maximum_extension_um):
            continue
        distance = max(0.0, distance)
        intersection = end + distance * direction
        tangent = np.asarray([-normal[1], normal[0]], dtype=float)
        if abs(float(np.dot(intersection - center, tangent))) > 0.5 * float(length) + 0.05 * spacing:
            continue
        candidates.append((distance, intersection))
    if not candidates:
        return None
    return min(candidates, key=lambda item: item[0])[1]


def _target_open_faces(
    intended_outlet: int,
    outlet_label: np.ndarray,
    lumen_mask: np.ndarray,
    boundary_normal_xz: np.ndarray | None,
    domain: GridDomain,
    open_boundaries: "BoundaryFluxFields | None",
) -> list[tuple[np.ndarray, np.ndarray, float]]:
    """读取求解器权威 open faces；独立调用时才从标签和 helper normal 重建。"""

    spacing = float(domain.spacing_um)
    if open_boundaries is not None:
        labels = np.asarray(open_boundaries.open_face_label, dtype=int)
        kinds = np.asarray(open_boundaries.open_face_kind, dtype=int)
        selected = np.flatnonzero((labels == int(intended_outlet)) & (kinds > 0))
        faces: list[tuple[np.ndarray, np.ndarray, float]] = []
        for row in selected:
            ix, iz = (int(value) for value in open_boundaries.open_face_cell_ij[row])
            normal = np.asarray(open_boundaries.open_face_normal_xz[row], dtype=float)
            cell_center = np.asarray(
                [domain.x_coordinates_um[ix], domain.z_coordinates_um[iz]],
                dtype=float,
            )
            face_center = cell_center + 0.5 * spacing * normal
            faces.append((face_center, normal, float(open_boundaries.open_face_length_um[row])))
        return faces

    labels = np.asarray(outlet_label, dtype=int)
    lumen = np.asarray(lumen_mask, dtype=bool)
    helper = None if boundary_normal_xz is None else np.asarray(boundary_normal_xz, dtype=float)
    faces = []
    directions = (
        (-1, 0, np.asarray([-1.0, 0.0])),
        (1, 0, np.asarray([1.0, 0.0])),
        (0, -1, np.asarray([0.0, -1.0])),
        (0, 1, np.asarray([0.0, 1.0])),
    )
    for ix, iz in np.argwhere(labels == int(intended_outlet)):
        ix, iz = int(ix), int(iz)
        cell_normal = np.zeros(2, dtype=float) if helper is None else helper[ix, iz]
        for di, dj, normal in directions:
            ni, nj = ix + di, iz + dj
            exposed = ni < 0 or ni >= lumen.shape[0] or nj < 0 or nj >= lumen.shape[1] or not lumen[ni, nj]
            if not exposed:
                continue
            if np.linalg.norm(cell_normal) > np.finfo(float).eps and float(np.dot(normal, cell_normal)) < 0.35:
                continue
            center = np.asarray([domain.x_coordinates_um[ix], domain.z_coordinates_um[iz]], dtype=float)
            faces.append((center + 0.5 * spacing * normal, normal, spacing))
    return faces


def _ordered_covered_vessels_outside_junctions(
    points_xz: np.ndarray,
    domain: GridDomain,
    raster: RasterizedVessels,
) -> tuple[int, ...]:
    """返回沿 root→outlet 顺序出现、且具有有效非 junction 弧长的血管 ID。"""

    points = np.asarray(points_xz, dtype=float)
    if points.shape[0] < 2:
        return ()
    midpoint = 0.5 * (points[:-1] + points[1:])
    segment_length = np.linalg.norm(np.diff(points, axis=0), axis=1)
    midpoint_vessel = _sample_integer_grid(midpoint, domain, raster.vessel_id, invalid=-1)
    junction_mask = getattr(raster, "junction_core_mask", None)
    if junction_mask is None:
        midpoint_junction = np.zeros(midpoint.shape[0], dtype=bool)
    else:
        midpoint_junction = _sample_integer_grid(midpoint, domain, junction_mask, invalid=1).astype(bool)
    length_by_vessel: dict[int, float] = {}
    for vessel_id, length, at_junction in zip(midpoint_vessel, segment_length, midpoint_junction):
        value = int(vessel_id)
        if value < 0 or bool(at_junction):
            continue
        length_by_vessel[value] = length_by_vessel.get(value, 0.0) + float(length)
    minimum_length = 0.50 * float(domain.spacing_um)
    eligible = {
        vessel_id for vessel_id, length in length_by_vessel.items() if length >= minimum_length
    }
    ordered: list[int] = []
    for vessel_id, at_junction in zip(midpoint_vessel, midpoint_junction):
        value = int(vessel_id)
        if value not in eligible or bool(at_junction):
            continue
        if not ordered or ordered[-1] != value:
            ordered.append(value)
    return tuple(ordered)


def _globally_observable_vessel_ids(
    raster: RasterizedVessels,
    domain: GridDomain,
) -> frozenset[int]:
    """列出在非 junction 栅格中具有可观测中心流路的血管。"""

    observable_mask = np.asarray(raster.lumen_mask, dtype=bool).copy()
    junction_mask = getattr(raster, "junction_core_mask", None)
    if junction_mask is not None:
        observable_mask &= ~np.asarray(junction_mask, dtype=bool)
    centerline_distance = getattr(raster, "distance_to_centerline_um", None)
    if centerline_distance is not None:
        # 只剩在 junction 外缘的少量侧壁像素，不代表主流线能真实进入该段；
        # 至少保留一层距理论中心线不超过 0.75 个网格的 cell 才视为可观测。
        observable_mask &= (
            np.asarray(centerline_distance, dtype=float) <= 0.75 * float(domain.spacing_um)
        )
    vessel_ids = np.asarray(raster.vessel_id, dtype=int)[observable_mask]
    return frozenset(int(value) for value in np.unique(vessel_ids) if int(value) >= 0)


def _is_ordered_subsequence(values: tuple[int, ...], expected: tuple[int, ...]) -> bool:
    """检查所有已观测 ID 均沿唯一祖先链前进，拒绝乱序、回跳和兄弟分支。"""

    expected_position = {value: index for index, value in enumerate(expected)}
    positions = [expected_position.get(value, -1) for value in values]
    return all(position >= 0 for position in positions) and all(
        current > previous for previous, current in zip(positions[:-1], positions[1:])
    )


def _ancestor_chains_by_outlet(
    outlet_vessels: dict[int, int],
    vessels: "Sequence[Vessel] | None",
) -> dict[int, tuple[int, ...]]:
    if vessels is None:
        return {}
    by_id = {int(vessel.vid): vessel for vessel in vessels}
    chains: dict[int, tuple[int, ...]] = {}
    for outlet_label, terminal_vessel_id in outlet_vessels.items():
        chain: list[int] = []
        current = int(terminal_vessel_id)
        seen: set[int] = set()
        while current >= 0 and current in by_id and current not in seen:
            seen.add(current)
            chain.append(current)
            current = int(by_id[current].parent_id)
        chains[int(outlet_label)] = tuple(reversed(chain))
    return chains


def _normalized_arc_progress(points: np.ndarray) -> np.ndarray:
    if points.shape[0] <= 1:
        return np.zeros(points.shape[0], dtype=np.float32)
    distance = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(distance)))
    if cumulative[-1] <= np.finfo(float).eps:
        return np.zeros(points.shape[0], dtype=np.float32)
    return np.asarray(cumulative / cumulative[-1], dtype=np.float32)


def _outlet_vessel_map(
    outlet_label: np.ndarray,
    vessel_id: np.ndarray,
    open_boundaries: "BoundaryFluxFields | None" = None,
) -> dict[int, int]:
    result: dict[int, int] = {}
    labels = np.asarray(outlet_label, dtype=int)
    vessels = np.asarray(vessel_id, dtype=int)
    authoritative_ids = tuple(int(value) for value in getattr(open_boundaries, "outlet_ids", ()))
    for label in np.unique(labels[labels > 0]):
        label_value = int(label)
        if 1 <= label_value <= len(authoritative_ids):
            result[label_value] = authoritative_ids[label_value - 1]
            continue
        values = vessels[labels == label]
        values = values[values >= 0]
        if values.size:
            unique, counts = np.unique(values, return_counts=True)
            result[label_value] = int(unique[int(np.argmax(counts))])
    return result


def _required_labels(values: np.ndarray | None, shape: tuple[int, int], name: str) -> np.ndarray:
    if values is None:
        raise ValueError(f"{name} is required for root-to-outlet streamline validation.")
    labels = np.asarray(values, dtype=np.int32)
    if labels.shape != tuple(shape):
        raise ValueError(f"{name} does not match GridDomain.shape.")
    return labels


def _maximum_integration_length(domain: GridDomain, settings: StreamlineSettings) -> float:
    x_extent = float(domain.x_coordinates_um[-1] - domain.x_coordinates_um[0])
    z_extent = float(domain.z_coordinates_um[-1] - domain.z_coordinates_um[0])
    diagonal = float(np.hypot(x_extent, z_extent))
    return max(float(settings.max_length_factor) * diagonal, 100.0 * float(domain.spacing_um))


def _require_pyvista():
    try:
        import pyvista as pv
    except ImportError as exc:  # pragma: no cover - 由 facade 给出完整安装提示
        raise ImportError("PyVista is required for VTK streamline tracing.") from exc
    return pv
