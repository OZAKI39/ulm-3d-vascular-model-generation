"""Identity validation and diagnostics on unchanged arrays.

Area and centroid formulas are independently implemented standard triangle
geometry. Legacy metadata conventions are documented in README_vessel_geometry.
No topology repair, vertex compaction, or face reordering is performed.
"""

import re

import numpy as np

from .io import unit_factor
from .model import BoundaryPatch, GeometryError, VesselGeometry


def integer(value, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not re.fullmatch(r"[+-]?\d+", str(value)):
        raise GeometryError("INVALID_INTEGER", f"{name} 必须为整数，不允许截断：{value!r}")
    return int(value)


def validate_arrays(points, faces, labels) -> None:
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0 or points.dtype.kind not in "fiu":
        raise GeometryError("INVALID_POINTS", f"坐标应为非空 N×3 数值数组：{points.shape}")
    if not np.all(np.isfinite(points)):
        raise GeometryError("NONFINITE_COORDINATES", "点坐标含 NaN / Inf")
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0 or faces.dtype.kind not in "iu":
        raise GeometryError("INVALID_FACES", "连接关系必须是非空 M×3 整数数组")
    if faces.min() < 0 or faces.max() >= len(points):
        raise GeometryError("FACE_INDEX_OUT_OF_RANGE", f"面索引超出 [0, {len(points)-1}]；实际 [{faces.min()}, {faces.max()}]")
    if labels.ndim != 1 or len(labels) != len(faces):
        raise GeometryError("LABEL_LENGTH", f"标签形状 {labels.shape} 与面数 {len(faces)} 不匹配")
    if labels.dtype.kind not in "iu":
        raise GeometryError("NONINTEGER_LABELS", f"CellEntityIds 必须使用整数类型，不能静默截断 {labels.dtype}")
    if labels.dtype.kind == "u" and labels.max() > np.iinfo(np.int64).max:
        raise GeometryError("INVALID_INTEGER", "标签超出有符号 64 位整数范围")


def triangle_geometry(points, faces):
    triangles = np.asarray(points, dtype=np.float64)[faces]
    crosses = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    areas = np.linalg.norm(crosses, axis=1) * 0.5
    return areas, triangles.mean(axis=1), crosses


def topology_diagnostics(points, faces, area_tolerance: float) -> dict:
    """Count edges and edge-connected face components on the supplied copy."""
    directed = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    keys, inverse, counts = np.unique(np.sort(directed, axis=1), axis=0, return_inverse=True, return_counts=True)
    directions = np.where(directed[:, 0] < directed[:, 1], 1, -1)
    sums = np.bincount(inverse, weights=directions, minlength=len(keys))
    owners = np.tile(np.arange(len(faces)), 3)
    order = np.argsort(inverse, kind="stable")
    ordered_group, ordered_owner = inverse[order], owners[order]
    neighbors = np.flatnonzero(ordered_group[1:] == ordered_group[:-1])
    parent = np.arange(len(faces))
    size = np.ones(len(faces), dtype=np.int64)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for index in neighbors:
        a, b = find(ordered_owner[index]), find(ordered_owner[index + 1])
        if a != b:
            if size[a] < size[b]:
                a, b = b, a
            parent[b] = a
            size[a] += size[b]
    areas, _, _ = triangle_geometry(points, faces)
    repeated_vertex = (faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 2] == faces[:, 0])
    duplicate_count = len(faces) - len(np.unique(np.sort(faces, axis=1), axis=0))
    degenerate = int(np.count_nonzero((areas <= area_tolerance) | repeated_vertex))
    boundary_count, nonmanifold = int(np.count_nonzero(counts == 1)), int(np.count_nonzero(counts > 2))
    inconsistent = int(np.count_nonzero((counts == 2) & (sums != 0)))
    referenced = np.unique(faces)
    origin = np.asarray(points, dtype=float)[referenced].mean(0)
    t = np.asarray(points, dtype=float)[faces] - origin
    signed_volume = float(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6)
    components = int(np.count_nonzero(parent == np.arange(len(faces))))
    return {
        "triangle_count": len(faces), "referenced_point_count": len(referenced),
        "boundary_edge_count": boundary_count, "nonmanifold_edge_count": nonmanifold,
        "inconsistently_oriented_shared_edge_count": inconsistent,
        "winding_consistent": inconsistent == 0 and nonmanifold == 0,
        "closed": boundary_count == 0 and nonmanifold == 0,
        "edge_connected_face_components": components,
        "duplicate_triangle_count": duplicate_count, "degenerate_triangle_count": degenerate,
        "degenerate_area_tolerance_m2": area_tolerance,
        "signed_volume_m3": signed_volume,
        "basic_topology_status": "PASS" if not (nonmanifold or inconsistent or duplicate_count or degenerate) else "WARNING",
    }


def build_geometry(points, faces, labels, rows, *, source_length_unit, unit_source,
                   wall_entity_ids, identity_reference=None, cell_data=None,
                   tolerances=None, reference_center_method="area_weighted") -> VesselGeometry:
    points, faces, labels = np.asarray(points), np.asarray(faces), np.asarray(labels)
    validate_arrays(points, faces, labels)
    factor = unit_factor(source_length_unit)
    limits = {"area_absolute_m2": 1e-18, "area_relative": 1e-10,
              "center_absolute_m": 5e-11, "degenerate_area_m2": 1e-28,
              "volume_m3": 1e-30, **(tolerances or {})}
    cell_data = cell_data or {}
    for name, values in cell_data.items():
        if len(values) != len(faces):
            raise GeometryError("LABEL_LENGTH", f"单元数组 {name} 长度与面数不一致")
    walls = [integer(v, "wall_entity_id") for v in wall_entity_ids]
    if not walls or len(walls) != len(set(walls)):
        raise GeometryError("WALL_IDENTITY_MISSING", "壁面身份必须由来源明确指定，且不能重复")
    mapping = {}
    port_ids, boundary_indices = set(), set()
    role_kinds = {"ASSUMED_INLET": "inlet", "ASSUMED_OUTLET": "outlet",
                  "INLET": "inlet", "OUTLET": "outlet", "WALL": "wall",
                  "inlet": "inlet", "outlet": "outlet", "wall": "wall"}
    for row in rows:
        raw_entity = row.get("vmtk_cap_entity_id", row.get("entity_id"))
        entity = integer(raw_entity, "entity_id")
        role = row.get("role")
        if role not in role_kinds:
            raise GeometryError("ROLE_MISSING", f"标签 {entity} 的 role 缺失或无法确认：{role!r}")
        kind = role_kinds[role]
        port = row.get("port_id") or None
        boundary_index = integer(row["boundary_index"], "boundary_index") if row.get("boundary_index") not in (None, "") else None
        if entity in mapping or (port is not None and port in port_ids) or (boundary_index is not None and boundary_index in boundary_indices):
            raise GeometryError("DUPLICATE_MAPPING", f"entity_id / port_id / boundary_index 重复：{row}")
        if (kind == "wall") != (entity in walls):
            raise GeometryError("CONFLICTING_MAPPING", f"清单与已确认壁面身份冲突：entity_id={entity}, role={role}")
        if kind != "wall" and (not port or not row.get("boundary_origin")):
            raise GeometryError("PORT_IDENTITY_MISSING", f"端口缺少 port_id 或 boundary_origin：{row}")
        mapping[entity] = {"row": dict(row), "kind": kind, "role": role, "port": port, "index": boundary_index}
        if port is not None:
            port_ids.add(port)
        if boundary_index is not None:
            boundary_indices.add(boundary_index)
    for wall in walls:
        mapping.setdefault(wall, {"row": {}, "kind": "wall", "role": "WALL", "port": None, "index": None})
    surface_values = set(int(x) for x in np.unique(labels))
    if surface_values != set(mapping):
        raise GeometryError("UNMAPPED_ENTITY", f"标签与已确认分区不一致；未知表面标签={sorted(surface_values-set(mapping))}；没有对应表面的身份={sorted(set(mapping)-surface_values)}")
    reference_by_port = {}
    for row in identity_reference or []:
        if row.get("port_id") in reference_by_port:
            raise GeometryError("DUPLICATE_MAPPING", "来源 QC 存在重复 port_id")
        reference_by_port[row["port_id"]] = row
    if identity_reference is not None and set(reference_by_port) != port_ids:
        raise GeometryError("SOURCE_MISMATCH", "来源 QC 的端口集合与 CSV 不一致")
    meters = points.astype(np.float64) * factor
    if not np.all(np.isfinite(meters)):
        raise GeometryError("NONFINITE_COORDINATES", "单位转换后坐标溢出")
    topology = topology_diagnostics(meters, faces, limits["degenerate_area_m2"])
    topology["complete_surface_status"] = (
        "PASS" if topology["closed"] and topology["edge_connected_face_components"] == 1
        and topology["basic_topology_status"] == "PASS" else "WARNING"
    )
    # Exact duplicate coordinates are merged only in a diagnostic copy (0 m tolerance).
    unique_points, duplicate_map = np.unique(meters, axis=0, return_inverse=True)
    welded_topology = topology_diagnostics(unique_points, duplicate_map[faces], limits["degenerate_area_m2"])
    outward_basis = (topology["closed"] and topology["winding_consistent"]
                     and topology["edge_connected_face_components"] == 1
                     and topology["degenerate_triangle_count"] == 0
                     and topology["duplicate_triangle_count"] == 0
                     and abs(topology["signed_volume_m3"]) > limits["volume_m3"])
    checks = {
        "migration_status": "PASS", "array_validation": "PASS", "identity_mapping": "PASS",
        "tolerances": limits, "boundary_comparisons": [], "original_topology": topology,
        "diagnostic_exact_coordinate_merge": {"tolerance_m": 0.0, "rule": "exact coordinate equality only",
            "unique_point_count_including_unused": len(unique_points), "formal_arrays_changed": False,
            "topology": welded_topology},
        "self_intersections": {"status": "NOT_CHECKED", "reason": "本次未执行复杂三角形自交检查；历史 QC 的 PASS 不作为本次验证"},
        "normal_basis": {"closed_oriented_single_surface": outward_basis,
            "method": "面叉积合量按完整闭合表面的有符号体积确定外向约定；仅调整派生向量，不翻转原三角面",
            "limitation": "自交检查未执行，法向是闭合表面外向约定下的几何法向，不是流速或已测血流方向"},
        "surface_processing": {"translated": False, "rotated": False, "smoothed": False,
            "remeshed": False, "holes_filled": False, "vertices_removed": False,
            "face_winding_changed": False, "only_transform": f"points_m = float64(points_source) * {factor}"},
        "patch_topology": {},
    }
    patches = []
    assignments = np.zeros(len(faces), dtype=np.uint8)
    ordered = sorted(mapping.items(), key=lambda item: (item[1]["kind"] != "wall",
                     item[1]["index"] if item[1]["index"] is not None else 2**63,
                     item[1]["port"] or "", item[0]))
    for entity, identity in ordered:
        ids = np.flatnonzero(labels == entity)
        assignments[ids] += 1
        row, kind, port = identity["row"], identity["kind"], identity["port"]
        areas, centers, cross = triangle_geometry(meters, faces[ids])
        total = float(areas.sum())
        if not np.isfinite(total) or total <= 0:
            raise GeometryError("ZERO_PATCH_AREA", f"分区 {entity} 的面积为零或非有限值")
        center = np.average(centers, axis=0, weights=areas)
        vector = cross.sum(0)
        norm = float(np.linalg.norm(vector))
        normal = None
        normal_status = "NOT_APPLICABLE_WALL" if kind == "wall" else "PENDING"
        if kind != "wall" and outward_basis and norm > total * 1e-12:
            normal = (vector / norm * np.sign(topology["signed_volume_m3"])).tolist()
            normal_status = "OUTWARD_BY_CLOSED_SURFACE_CONVENTION"
        expected_cells = {
            "port_id": port or "", "boundary_origin": row.get("boundary_origin", "WALL"),
            "boundary_index": identity["index"] if identity["index"] is not None else -1,
            "boundary_type_code": {"wall": 0, "inlet": 1, "outlet": 2}[kind],
        }
        for name, expected in expected_cells.items():
            if name in cell_data and (name != "boundary_index" or kind == "wall" or identity["index"] is not None):
                if not np.all(cell_data[name][ids] == expected):
                    raise GeometryError("CELL_IDENTITY_CONFLICT", f"VTP 的 {name} 与清单冲突：entity={entity}, expected={expected!r}")
        comparison = {"entity_id": entity, "port_id": port, "status": "PASS", "checks": []}
        sources = [("boundary_manifest", row)]
        if port in reference_by_port:
            ref = reference_by_port[port]
            for key, value in (("vmtk_cap_entity_id", entity), ("role", identity["role"]), ("boundary_origin", row.get("boundary_origin"))):
                actual = integer(ref[key], key) if key == "vmtk_cap_entity_id" else ref.get(key)
                if actual != value:
                    raise GeometryError("SOURCE_MISMATCH", f"CSV 与 QC 身份不一致：port_id={port}, {key}: {actual!r} != {value!r}")
            if identity["index"] is not None and integer(ref.get("boundary_index"), "boundary_index") != identity["index"]:
                raise GeometryError("SOURCE_MISMATCH", f"CSV 与 QC boundary_index 不一致：{port}")
            sources.append(("identity_qc", ref))
        for source_name, source in sources:
            if source.get("triangle_count") not in (None, ""):
                expected = integer(source["triangle_count"], "triangle_count")
                comparison["checks"].append({"source": source_name, "quantity": "triangle_count", "actual": len(ids), "expected": expected})
                if len(ids) != expected:
                    raise GeometryError("COUNT_MISMATCH", f"端口面数与清单不一致：{port}, {len(ids)} != {expected}", comparison=comparison)
            for area_key, scale in (("area_um2", 1e-12), ("area_m2", 1.0)):
                if source.get(area_key) not in (None, ""):
                    expected = float(source[area_key]) * scale
                    tolerance = limits["area_absolute_m2"] + limits["area_relative"] * abs(expected)
                    delta = abs(total - expected)
                    comparison["checks"].append({"source": source_name, "quantity": "area_m2", "actual": total, "expected": expected, "absolute_error_m2": delta, "allowed_error_m2": tolerance})
                    if not np.isfinite(expected) or expected <= 0 or delta > tolerance:
                        raise GeometryError("AREA_MISMATCH", f"端口面积对照失败：{port}；误差 {delta} m²，容差 {tolerance} m²", comparison=comparison)
            for center_key, scale in (("centroid_um", 1e-6), ("center_m", 1.0), ("centroid_m", 1.0)):
                if center_key in source:
                    expected = np.asarray(source[center_key], dtype=float) * scale
                    if expected.shape != (3,) or not np.all(np.isfinite(expected)):
                        raise GeometryError("INVALID_CENTER", f"来源中心不合法：{port}, {source[center_key]}")
                    # Reproduce only the documented legacy statistic, not its geometry pipeline.
                    source_convention_center = (points[faces[ids]].mean(axis=1).mean(axis=0).astype(float) * factor
                        if reference_center_method == "triangle_centroid_mean_source_dtype" else center)
                    delta = float(np.linalg.norm(source_convention_center - expected))
                    comparison["checks"].append({"source": source_name, "quantity": "center_m",
                        "reference_method": reference_center_method, "actual_same_definition_m": source_convention_center.tolist(),
                        "reference_m": expected.tolist(), "absolute_error_m": delta,
                        "allowed_error_m": limits["center_absolute_m"],
                        "area_weighted_center_m": center.tolist(),
                        "different_definition_distance_m": float(np.linalg.norm(center - expected))})
                    if delta > limits["center_absolute_m"]:
                        raise GeometryError("CENTER_MISMATCH", f"同定义端口中心对照失败：{port}；误差 {delta} m", comparison=comparison)
        checks["boundary_comparisons"].append(comparison)
        checks["patch_topology"][str(entity)] = topology_diagnostics(meters, faces[ids], limits["degenerate_area_m2"])
        checks["patch_topology"][str(entity)]["open_edges_interpretation"] = "wall 分区在端口处有开边属于正常几何分区，不补洞" if kind == "wall" else "独立端口盖片边缘允许开放"
        name = row.get("display_name") or (f"{'假定' if identity['role'].startswith('ASSUMED_') else ''}{'入口' if kind == 'inlet' else '出口'} · {port}" if kind != "wall" else f"管壁 · entity {entity}")
        patches.append(BoundaryPatch(name, entity, port, identity["index"], kind,
            identity["role"], row.get("boundary_origin", "WALL"), ids.tolist(), total,
            center.tolist(), normal, normal_status, {"manifest_row": row,
                "identity_source": "boundary_manifest + CellEntityIds + identity_qc" if port else "identity_qc.wall_entity_id"}))
    if not np.all(assignments == 1):
        raise GeometryError("PARTITION_NOT_EXACT", "每个原始三角面必须恰好属于一个已确认分区")
    checks["partition"] = {"status": "PASS", "total_faces": len(faces),
        "patch_face_count_sum": sum(p.triangle_count for p in patches), "unassigned_faces": 0, "multiply_assigned_faces": 0}
    return VesselGeometry(points.copy(), meters, faces.copy(), np.arange(len(points), dtype=np.int64),
        np.arange(len(faces), dtype=np.int64), labels.copy(), patches, source_length_unit, factor,
        unit_source, {k: v.copy() for k, v in cell_data.items()}, checks)


def assert_same_geometry(a: VesselGeometry, b: VesselGeometry) -> dict:
    from dataclasses import asdict
    names = ("points_source", "points_m", "triangles", "original_point_ids", "original_face_ids", "entity_ids")
    results = {name: bool(np.array_equal(getattr(a, name), getattr(b, name))
                           and getattr(a, name).dtype == getattr(b, name).dtype) for name in names}
    results["patches"] = [asdict(p) for p in a.patches] == [asdict(p) for p in b.patches]
    results["units"] = a.source_length_unit == b.source_length_unit and a.to_meter == b.to_meter and a.unit_source == b.unit_source
    results["source_cell_data"] = a.source_cell_data.keys() == b.source_cell_data.keys() and all(
        np.array_equal(a.source_cell_data[key], b.source_cell_data[key]) and a.source_cell_data[key].dtype == b.source_cell_data[key].dtype for key in a.source_cell_data)
    if not all(results.values()):
        raise GeometryError("ROUNDTRIP_MISMATCH", f"标准化数据包回读不一致：{results}")
    return {"status": "PASS", "exact_array_and_identity_checks": results,
            "maximum_coordinate_difference_m": 0.0, "array_tolerance": 0}
