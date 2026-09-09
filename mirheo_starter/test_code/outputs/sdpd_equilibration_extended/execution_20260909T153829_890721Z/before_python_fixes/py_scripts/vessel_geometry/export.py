"""Exclusive, portable, non-pickle data packages and the import orchestrator."""

from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import shutil
import struct
import sys
import traceback
import uuid

import numpy as np
import yaml

from .io import (git_identity, load_config, read_boundary_manifest, read_json,
                 read_tagged_vtp, require_file, resolve_recorded_path, safe_output,
                 sha256_file, unit_factor, write_json)
from .model import BoundaryPatch, GeometryError, ImportResult, VesselGeometry
from .validation import assert_same_geometry, build_geometry, integer, validate_arrays


def new_run_id() -> str:
    return datetime.now(timezone.utc).strftime("import_%Y%m%dT%H%M%S_%fZ_") + uuid.uuid4().hex[:8]


def copy_exclusive(source: Path, destination: Path) -> None:
    with source.open("rb") as src, destination.open("xb") as dst:
        shutil.copyfileobj(src, dst)


def export_package(geometry: VesselGeometry, directory: Path, source_surface: Path, *, source_project: Path) -> None:
    """Write into an already reserved, empty run directory; never overwrite."""
    directory = safe_output(directory, source_project, [source_surface])
    copy_exclusive(source_surface, directory / "source_surface.vtp")
    arrays = {name: getattr(geometry, name) for name in
              ("points_source", "points_m", "triangles", "original_point_ids", "original_face_ids", "entity_ids")}
    arrays.update({f"cell_{index}": values for index, values in enumerate(geometry.source_cell_data.values())})
    with (directory / "vessel_geometry.npz").open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    write_json(directory / "boundary_patches.json", [asdict(p) for p in geometry.patches])
    write_json(directory / "geometry_metadata.json", {
        "schema_version": 1, "source_length_unit": geometry.source_length_unit,
        "to_meter": geometry.to_meter, "unit_source": geometry.unit_source,
        "source_cell_array_names": list(geometry.source_cell_data),
        "summary": geometry.summary(),
    })


def load_package(directory: Path, *, verify_manifest: bool = True) -> VesselGeometry:
    directory = Path(directory).resolve()
    if verify_manifest:
        if (directory / "migration_failure.json").exists():
            raise GeometryError("PACKAGE_NOT_ACCEPTED", "本次迁移失败，不得沿用历史 PASS")
        if (directory / "migration_intent.json").exists() and not (directory / "migration_manifest.json").exists():
            raise GeometryError("PACKAGE_NOT_ACCEPTED", "本次迁移尚未完成")
        manifest = read_json(directory / "import_manifest.json")
        if manifest.get("schema_version") != 1 or manifest.get("status") != "PASS":
            raise GeometryError("PACKAGE_NOT_ACCEPTED", f"数据包未完成成功导入：{directory}")
        required = {"vessel_geometry.npz", "geometry_metadata.json", "boundary_patches.json",
                    "source_surface.vtp", "import_checks.json", "import_config.yaml"}
        if not required <= manifest.get("output_sha256", {}).keys():
            raise GeometryError("PACKAGE_INVALID", "数据包清单未覆盖必要文件")
        for name, expected in manifest["output_sha256"].items():
            path = (directory / name).resolve()
            if not path.is_relative_to(directory):
                raise GeometryError("PACKAGE_INVALID", f"数据包路径越界：{name}")
            if sha256_file(require_file(path)) != expected:
                raise GeometryError("PACKAGE_HASH_MISMATCH", f"数据包文件哈希不符：{path}")
        if (directory / "migration_checks.json").exists() and not (directory / "migration_manifest.json").is_file():
            raise GeometryError("PACKAGE_NOT_ACCEPTED", "迁移尚未完成，不得将历史 PASS 用作本次 PASS")
        if (directory / "migration_manifest.json").exists():
            migrated = read_json(directory / "migration_manifest.json")
            if migrated.get("schema_version") != 1 or migrated.get("status") != "PASS":
                raise GeometryError("PACKAGE_NOT_ACCEPTED", "迁移清单未通过")
            if migrated.get("source_manifest_sha256") != sha256_file(directory / "import_manifest.json"):
                raise GeometryError("PACKAGE_HASH_MISMATCH", "原始清单与已验收固定哈希不一致")
            required |= {"import_manifest.json", "migration_checks.json", "migration_config.yaml"}
            if not required <= migrated.get("output_sha256", {}).keys():
                raise GeometryError("PACKAGE_INVALID", "迁移清单未覆盖必要文件")
            for name, expected in migrated["output_sha256"].items():
                path = (directory / name).resolve()
                if not path.is_relative_to(directory):
                    raise GeometryError("PACKAGE_INVALID", f"迁移数据包路径越界：{name}")
                if sha256_file(require_file(path)) != expected:
                    raise GeometryError("PACKAGE_HASH_MISMATCH", f"迁移文件哈希不符：{path}")
    metadata = read_json(directory / "geometry_metadata.json")
    if metadata.get("schema_version") != 1 or metadata["to_meter"] != unit_factor(metadata["source_length_unit"]):
        raise GeometryError("PACKAGE_INVALID", "数据包版本或单位换算不合法")
    with np.load(require_file(directory / "vessel_geometry.npz"), allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in
                  ("points_source", "points_m", "triangles", "original_point_ids", "original_face_ids", "entity_ids")}
        cell_data = {name: saved[f"cell_{index}"].copy() for index, name in enumerate(metadata["source_cell_array_names"])}
    validate_arrays(arrays["points_source"], arrays["triangles"], arrays["entity_ids"])
    if not np.array_equal(arrays["points_m"], arrays["points_source"].astype(float) * metadata["to_meter"]):
        raise GeometryError("PACKAGE_INVALID", "回读坐标与显式单位转换不一致")
    for name, count in (("original_point_ids", len(arrays["points_m"])), ("original_face_ids", len(arrays["triangles"]))):
        if not np.array_equal(arrays[name], np.arange(count)):
            raise GeometryError("PACKAGE_INVALID", f"原始编号不一致：{name}")
    patches = [BoundaryPatch(**row) for row in read_json(directory / "boundary_patches.json")]
    assignment = np.zeros(len(arrays["triangles"]), dtype=np.int32)
    entities, ports, boundary_indices = set(), set(), set()
    role_kinds = {"ASSUMED_INLET": "inlet", "ASSUMED_OUTLET": "outlet",
                  "INLET": "inlet", "OUTLET": "outlet", "WALL": "wall",
                  "inlet": "inlet", "outlet": "outlet", "wall": "wall"}
    for patch in patches:
        entity = integer(patch.entity_id, "entity_id")
        if entity in entities or (patch.port_id is not None and patch.port_id in ports):
            raise GeometryError("DUPLICATE_MAPPING", "回读分区 entity_id / port_id 重复")
        entities.add(entity)
        if patch.port_id is not None:
            ports.add(patch.port_id)
        if patch.boundary_index is not None:
            index = integer(patch.boundary_index, "boundary_index")
            if index in boundary_indices:
                raise GeometryError("DUPLICATE_MAPPING", "回读 boundary_index 重复")
            boundary_indices.add(index)
        if role_kinds.get(patch.original_role) != patch.kind or not patch.boundary_origin:
            raise GeometryError("PORT_IDENTITY_MISSING", "回读 role / kind / origin 不一致或缺失")
        if patch.kind != "wall" and not patch.port_id:
            raise GeometryError("PORT_IDENTITY_MISSING", "回读端口缺少 port_id")
        ids = np.asarray(patch.face_ids)
        if ids.dtype.kind not in "iu" or not len(ids) or ids.min() < 0 or ids.max() >= len(assignment):
            raise GeometryError("PACKAGE_INVALID", "回读分区包含非法原始面编号")
        if not np.all(arrays["entity_ids"][ids] == patch.entity_id):
            raise GeometryError("PACKAGE_INVALID", "回读分区与标签不一致")
        np.add.at(assignment, ids, 1)
    if not np.all(assignment == 1):
        raise GeometryError("PACKAGE_INVALID", "回读分区未精确覆盖所有面")
    checks = read_json(directory / "import_checks.json") if (directory / "import_checks.json").exists() else {}
    if verify_manifest and (directory / "migration_manifest.json").exists():
        checks = read_json(directory / "migration_checks.json")
    return VesselGeometry(**arrays, patches=patches, source_length_unit=metadata["source_length_unit"],
        to_meter=metadata["to_meter"], unit_source=metadata["unit_source"], source_cell_data=cell_data, checks=checks)


def _companion_stl_check(path: Path, geometry: VesselGeometry) -> dict:
    raw = require_file(path).read_bytes()
    if len(raw) < 84:
        raise GeometryError("STL_INVALID", f"配套 STL 过短：{path}")
    count = struct.unpack_from("<I", raw, 80)[0]
    if len(raw) != 84 + 50 * count:
        return {"status": "NOT_CHECKED", "reason": "配套 STL 不是严格二进制格式，未比较；身份仍由 VTP 和清单确定"}
    dtype = np.dtype([("normal", "<f4", (3,)), ("points", "<f4", (3, 3)), ("attribute", "<u2")])
    actual = np.frombuffer(raw, dtype=dtype, count=count, offset=84)["points"]
    expected = geometry.points_m[geometry.triangles].astype(np.float32)
    if not np.array_equal(actual, expected):
        raise GeometryError("COMPANION_STL_MISMATCH", f"配套米单位 STL 与 VTP 转换后、按 float32 序列化的逐面坐标不一致：{path}")
    return {"status": "PASS", "triangle_count": count, "units": "m",
            "method": "逐三角形、逐顶点与 float32(points_m[triangles]) 精确比较，不使用 STL 推断标签", "tolerance_after_float32_cast_m": 0}


def _safe_error_details(value):
    if isinstance(value, dict):
        return {str(k): _safe_error_details(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_error_details(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def run_import(config_path: Path, *, progress=print) -> ImportResult:
    """Import once; failures retain a diagnostic record, never a success manifest."""
    config = load_config(config_path)
    if config.get("mode") == "accepted_package":
        from .migration import migrate_package
        return migrate_package(config, progress=progress)
    run_id = new_run_id()
    directory = safe_output(config["output_root"] / run_id, config["source_project"],
                            [config["tagged_surface"], config["boundary_manifest"]])
    directory.mkdir(parents=True, exist_ok=False)
    messages, before_hashes = [], {}
    def log(message):
        line = f"[{datetime.now().astimezone().isoformat(timespec='seconds')}] {message}"
        messages.append(line)
        if progress is not None:
            progress(line)
    try:
        log("[1/5] 检查明确来源、标签清单、单位证据和输入文件")
        keys = ("tagged_surface", "boundary_manifest", "source_unit_config", "identity_qc", "scale_qc", "companion_stl", "source_selection_config")
        paths = [require_file(config[key]) for key in keys if config.get(key)]
        paths.append(require_file(config["config_path"]))
        paths.extend(require_file(p) for p in config["reference_code"])
        before_hashes = {str(path): sha256_file(path) for path in paths}
        rows = read_boundary_manifest(config["boundary_manifest"])
        resolutions = []
        for row in rows:
            if row.get("stl_path"):
                local = resolve_recorded_path(row["stl_path"], config["source_run"])
                paths.append(local)
                resolutions.append({"port_id": row["port_id"], "original_path": row["stl_path"],
                                    "resolved_path": str(local), "rule": "exact selected-run path component; no basename search"})
                if not row.get("display_name"):
                    prefix = {"ASSUMED_INLET": "假定入口", "ASSUMED_OUTLET": "假定出口", "INLET": "入口", "OUTLET": "出口"}.get(row.get("role"), row.get("role", "未确认"))
                    row["display_name"] = f"{prefix} · {local.stem}"
        paths = list(dict.fromkeys(paths))
        for path in paths:
            if str(path) not in before_hashes:
                before_hashes[str(path)] = sha256_file(path)
        if config["source_selection_config"]:
            selection = yaml.safe_load(config["source_selection_config"].read_text(encoding="utf-8-sig"))
            recorded = selection.get("paths", {}).get("source_surface_run")
            if not recorded or (config["source_project"] / recorded).resolve() != config["source_run"]:
                raise GeometryError("SOURCE_AMBIGUOUS", "当前源配置 source_surface_run 与本次选择不一致，阻止真实导入", status="BLOCKED", configured_source_run=str(config["source_run"]), recorded=recorded)
        unit_config = yaml.safe_load(config["source_unit_config"].read_text(encoding="utf-8-sig"))
        recorded_unit = unit_config.get("geometry", {}).get("input_unit")
        if recorded_unit != config["source_length_unit"]:
            raise GeometryError("UNIT_EVIDENCE_CONFLICT", f"明确单位与同次运行配置冲突：{config['source_length_unit']} != {recorded_unit}")
        unit_source = {"configuration": str(config["source_unit_config"]), "key": "geometry.input_unit",
                       "value": recorded_unit, "sha256": before_hashes[str(config["source_unit_config"])],
                       "inferred_from_filename_or_size": False}
        if config["scale_qc"]:
            scale = read_json(config["scale_qc"])
            if scale.get("scale_factor") != unit_factor(recorded_unit):
                raise GeometryError("UNIT_EVIDENCE_CONFLICT", "历史单位 QC 的 scale_factor 与明确单位不一致")
            unit_source["scale_evidence"] = {"path": str(config["scale_qc"]), "scale_factor": scale["scale_factor"], "historical_status_not_reused": True}
        identity_qc = read_json(config["identity_qc"])
        mapping = identity_qc.get("boundary_mapping", identity_qc)
        if "wall_entity_ids" in mapping:
            walls = mapping["wall_entity_ids"]
        elif "wall_entity_id" in mapping:
            walls = [mapping["wall_entity_id"]]
        else:
            raise GeometryError("WALL_IDENTITY_MISSING", f"来源未明确 wall_entity_id，不能将未知标签当作壁面：{config['identity_qc']}")
        log(f"[2/5] 读取原始 VTP：{config['tagged_surface']}（不清理、不三角化、不修复）")
        points, faces, labels, cell_data = read_tagged_vtp(config["tagged_surface"])
        log(f"[3/5] 检查原始数组、分区身份与几何诊断：{len(points)} 点，{len(faces)} 三角面")
        geometry = build_geometry(points, faces, labels, rows, source_length_unit=recorded_unit,
            unit_source=unit_source, wall_entity_ids=walls, identity_reference=mapping.get("boundaries"),
            cell_data=cell_data, tolerances=config["tolerances"], reference_center_method=config["reference_center_method"])
        geometry.checks["unit_evidence"] = {"status": "PASS", **unit_source}
        geometry.checks["companion_stl"] = (_companion_stl_check(config["companion_stl"], geometry)
            if config["companion_stl"] else {"status": "NOT_PROVIDED"})
        log("[4/5] 保存原始 VTP 字节副本、标准化数组及分区；从导出包重新读取")
        export_package(geometry, directory, config["tagged_surface"], source_project=config["source_project"])
        if sha256_file(directory / "source_surface.vtp") != before_hashes[str(config["tagged_surface"])]:
            raise GeometryError("COPY_HASH_MISMATCH", "原始 VTP 副本哈希不一致")
        loaded = load_package(directory, verify_manifest=False)
        geometry.checks["export_roundtrip"] = assert_same_geometry(geometry, loaded)
        source_after = read_tagged_vtp(config["tagged_surface"])
        source_same = all(np.array_equal(a, b) and a.dtype == b.dtype for a, b in zip((points, faces, labels), source_after[:3], strict=True))
        if not source_same:
            raise GeometryError("SOURCE_CHANGED", "导入前后源几何数组发生变化")
        geometry.checks["source_geometry_unchanged"] = {"status": "PASS", "point_order_exact": True,
            "triangle_order_and_connectivity_exact": True, "labels_exact": True, "original_vtp_copy_sha256_exact": True}
        # Original small evidence files are copied byte-for-byte; historical Windows paths stay intact.
        evidence_dir = directory / "source_records"
        evidence_dir.mkdir()
        for key in ("boundary_manifest", "source_unit_config", "identity_qc", "scale_qc", "source_selection_config"):
            if config.get(key):
                copy_exclusive(config[key], evidence_dir / (key + config[key].suffix))
        copy_exclusive(config["config_path"], directory / "import_config.yaml")
        after_hashes = {str(path): sha256_file(path) for path in paths}
        if after_hashes != before_hashes:
            raise GeometryError("SOURCE_CHANGED", "处理前后输入文件 SHA-256 发生变化", before=before_hashes, after=after_hashes)
        geometry.checks["input_hashes_unchanged"] = {"status": "PASS", "before": before_hashes, "after": after_hashes}
        log("[5/5] 数据迁移 PASS；保存独立的几何诊断与来源记录（人工核查仍为 PENDING）")
        write_json(directory / "import_checks.json", geometry.checks)
        with (directory / "import.log").open("x", encoding="utf-8") as stream:
            stream.write("\n".join(messages) + "\n")
        outputs = {str(p.relative_to(directory)): sha256_file(p) for p in sorted(directory.rglob("*")) if p.is_file()}
        manifest = {"schema_version": 1, "run_id": run_id, "status": "PASS",
            "created_utc": datetime.now(timezone.utc).isoformat(), "source_project": str(config["source_project"]),
            "source_run": str(config["source_run"]), "tagged_surface": str(config["tagged_surface"]),
            "boundary_manifest": str(config["boundary_manifest"]), "input_sha256": before_hashes,
            "source_git": git_identity(config["source_project"], paths), "historical_path_resolutions": resolutions,
            "summary": geometry.summary(), "output_sha256": outputs,
            "hash_scope": "output_sha256 covers all package payloads except this self-referential manifest",
            "python_executable": sys.executable, "python_version": sys.version,
            "dependency_versions": {name: version(name) for name in ("numpy", "vtk", "PyYAML")},
            "importer_code_sha256": {str(path): sha256_file(path) for path in sorted(Path(__file__).parent.glob("*.py"))},
            "command_argv": sys.argv, "config_path": str(config["config_path"]),
            "manual_review": "PENDING", "lbm_mesh": "NOT_STARTED", "numerical_boundary_conditions": "NOT_STARTED",
            "lammps_integration": "NOT_STARTED", "source_functions_copied": [],
            "migration_note": "独立实现 VTP 原始数组读取、显式身份映射、三角几何公式和回读。仅使用旧工程输出契约及有记录的统计定义，不复制或运行旧 pipeline。"}
        write_json(directory / "import_manifest.json", manifest)
        return ImportResult(run_id, directory, load_package(directory), manifest)
    except Exception as exc:
        status = exc.status if isinstance(exc, GeometryError) else "FAIL"
        code = exc.code if isinstance(exc, GeometryError) else "IMPORT_INTERNAL_ERROR"
        log(f"[{status}] {exc}")
        failure = {"status": status, "code": code, "message": str(exc), "run_id": run_id,
            "package_path": str(directory), "input_sha256_before": before_hashes,
            "details": _safe_error_details(getattr(exc, "details", {})),
            "null_reason": "错误详情中的非有限数值以 null 表示，不输出 JSON NaN/Infinity",
            "manual_review": "PENDING", "lbm_mesh": "NOT_STARTED", "lammps_integration": "NOT_STARTED"}
        after_failure = {path: sha256_file(Path(path)) if Path(path).is_file() else None for path in before_hashes}
        failure["input_sha256_after"] = after_failure
        failure["input_hashes_unchanged"] = before_hashes == after_failure if before_hashes else None
        write_json(directory / "import_failure.json", failure)
        if not (directory / "import.log").exists():
            with (directory / "import.log").open("x", encoding="utf-8") as stream:
                stream.write("\n".join(messages) + "\n" + traceback.format_exc())
        error = GeometryError(code, str(exc), status=status, run_id=run_id, package_path=str(directory))
        raise error from exc
