"""Migrate an explicitly pinned stage-1 package without rewriting its payloads.

The schema-1 loader, VTP reader, mapping, diagnostics and triangle formulas are
reused from bloodflow_starter. Original manifests remain historical evidence.
See SOURCES.md for provenance; this module never imports a source project.
"""

from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path, PureWindowsPath
import hashlib
import json
import sys

import numpy as np
import yaml

from .io import (PROJECT_ROOT, read_boundary_manifest, read_json, read_tagged_vtp,
                 require_file, safe_output, sha256_file, write_json)
from .model import GeometryError, ImportResult
from .validation import assert_same_geometry, build_geometry


def load_migration_config(path, payload):
    required = {"schema_version", "mode", "source_package", "source_run_id",
                "source_manifest_sha256", "source_acceptance", "output_root", "review_output_root"}
    if required != payload.keys() or payload["schema_version"] != 1:
        raise GeometryError("CONFIG_INVALID", f"迁移配置字段必须为 {sorted(required)}，schema_version=1")
    config = dict(payload, config_path=path)
    for key in required - {"schema_version", "source_acceptance"}:
        if not isinstance(payload[key], str) or not payload[key].strip():
            raise GeometryError("CONFIG_INVALID", f"{key} 必须为非空字符串")
    digest = payload["source_manifest_sha256"]
    if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise GeometryError("CONFIG_INVALID", "source_manifest_sha256 必须为明确的 SHA-256")
    acceptance = payload["source_acceptance"]
    if not isinstance(acceptance, dict) or acceptance.get("status") not in {"ACCEPTED_BY_USER", "SYNTHETIC_TEST_ONLY"} or not acceptance.get("evidence"):
        raise GeometryError("CONFIG_INVALID", "必须记录来源验收状态及依据；合成测试单独标注")
    for key in ("source_package", "output_root", "review_output_root"):
        if PureWindowsPath(payload[key]).drive:
            raise GeometryError("CONFIG_INVALID", f"配置使用本地 WSL 路径：{key}")
        value = Path(payload[key]).expanduser()
        config[key] = (value if value.is_absolute() else PROJECT_ROOT / value).resolve()
    if config["source_package"].name != config["source_run_id"]:
        raise GeometryError("SOURCE_MISMATCH", "source_run_id 与明确选择的数据包目录不一致", status="BLOCKED")
    config["source_project"] = config["source_package"]
    for key in ("output_root", "review_output_root"):
        config[key] = safe_output(config[key], config["source_package"], [config["source_package"], path])
    a, b = config["output_root"], config["review_output_root"]
    if a.is_relative_to(b) or b.is_relative_to(a):
        raise GeometryError("UNSAFE_OUTPUT", "正式数据与核查目录必须相互独立")
    return config


def code_identity():
    paths = sorted(Path(__file__).parent.glob("*.py")) + [
        PROJECT_ROOT / "py_scripts/import_vessel_geometry.py",
        PROJECT_ROOT / "test_code/review_vessel_geometry.py",
        PROJECT_ROOT / "test_code/check_vessel_geometry_browser.cjs"]
    return {str(p.relative_to(PROJECT_ROOT)): sha256_file(p) for p in paths}


def dependency_identity():
    return {name: version(name) for name in ("numpy", "vtk", "PyYAML", "plotly")}


def _package_file(directory, name):
    value = Path(name)
    path = (directory / value).resolve()
    if value.is_absolute() or not path.is_relative_to(directory):
        raise GeometryError("PACKAGE_INVALID", f"清单中的相对路径越界：{name}")
    return require_file(path)


def collect_inputs(config):
    """Match exact manifest hashes; no latest-time selection or basename search."""
    source = config["source_package"]
    manifest_path = require_file(source / "import_manifest.json")
    if sha256_file(manifest_path) != config["source_manifest_sha256"]:
        raise GeometryError("SOURCE_MISMATCH", "已验收来源清单的固定 SHA-256 不符", status="BLOCKED")
    manifest = read_json(manifest_path)
    if manifest.get("run_id") != config["source_run_id"] or manifest.get("status") != "PASS":
        raise GeometryError("SOURCE_MISMATCH", "来源 run_id 或第一阶段自动导入状态不符", status="BLOCKED")
    copies = {"import_manifest.json": manifest_path}
    for name, expected in manifest["output_sha256"].items():
        path = _package_file(source, name)
        if sha256_file(path) != expected:
            raise GeometryError("PACKAGE_HASH_MISMATCH", f"已验收文件哈希不符：{path}", status="BLOCKED")
        copies[name] = path
    # Every originally recorded input is checked in place. These paths are only
    # needed for migration; target loading never dereferences them.
    originals = {}
    for recorded, expected in manifest["input_sha256"].items():
        path = require_file(Path(recorded))
        if sha256_file(path) != expected:
            raise GeometryError("SOURCE_MISMATCH", f"来源文件与已验收记录不符：{recorded}", status="BLOCKED")
        originals[recorded] = path
    # Copy only the referenced complete STL and port STLs, not old flow outputs.
    for row in manifest.get("historical_path_resolutions", []):
        resolved = row["resolved_path"]
        if resolved not in originals:
            raise GeometryError("SOURCE_MISMATCH", "历史端口解析记录未被输入哈希覆盖")
        path = originals[resolved]
        relative = path.relative_to(Path(manifest["source_run"]).resolve())
        copies[str(Path("source_companions") / relative)] = path
    source_config = yaml.safe_load((source / "import_config.yaml").read_text(encoding="utf-8-sig"))
    companion = source_config.get("companion_stl")
    if companion:
        path = (Path(manifest["source_run"]) / companion).resolve()
        if str(path) not in originals or not path.is_relative_to(Path(manifest["source_run"]).resolve()):
            raise GeometryError("SOURCE_MISMATCH", "配套 STL 未由同次运行及原始输入哈希共同确认")
        copies[str(Path("source_companions") / path.relative_to(Path(manifest["source_run"]).resolve()))] = path
    # This exact-run review is evidence, not a new human acceptance assertion.
    historical_review = source.parents[2] / "test_code/outputs/vessel_geometry" / config["source_run_id"]
    report_path = historical_review / "review_report.json"
    if report_path.is_file():
        report = read_json(report_path)
        if report.get("run_id") != config["source_run_id"] or report.get("import_manifest_sha256") != config["source_manifest_sha256"]:
            raise GeometryError("SOURCE_MISMATCH", "历史核查报告不对应所选数据包", status="BLOCKED")
        for name in ("review_report.json", "review_report.md"):
            if (historical_review / name).is_file():
                copies["historical_review/" + name] = require_file(historical_review / name)
    acceptance_path = config["source_acceptance"].get("record_path")
    if acceptance_path:
        path = Path(acceptance_path)
        path = require_file(path if path.is_absolute() else PROJECT_ROOT / path)
        record = yaml.safe_load(path.read_text(encoding="utf-8-sig"))
        if (record.get("geometry_package") != str(source)
                or record.get("geometry_manifest_sha256") != config["source_manifest_sha256"]
                or record.get("stage1_acceptance", {}).get("status") != config["source_acceptance"]["status"]):
            raise GeometryError("SOURCE_MISMATCH", "后续人工验收记录与所选第一阶段包不匹配", status="BLOCKED")
        copies["historical_review/later_acceptance_config.yaml"] = path
    paths = set(copies.values()) | set(originals.values()) | {require_file(config["config_path"])}
    hashes = {str(p): sha256_file(p) for p in sorted(paths)}
    mapping = [{"source_path": str(p), "target_relative_path": name, "sha256": hashes[str(p)]}
               for name, p in sorted(copies.items())]
    # Preserve Windows strings and add relocatable target links in a new record.
    by_source = {str(p): name for name, p in copies.items()}
    for row in manifest.get("historical_path_resolutions", []):
        mapping.append({**row, "target_relative_path": by_source[row["resolved_path"]],
                        "sha256": hashes[row["resolved_path"]]})
    return manifest, copies, hashes, mapping


def validate_accepted_package(directory, geometry):
    """Reuse the existing low-cost checks; never replace any accepted arrays."""
    points, faces, labels, cell_data = read_tagged_vtp(directory / "source_surface.vtp")
    array_checks = {name: bool(np.array_equal(a, b) and a.dtype == b.dtype) for name, a, b in (
        ("points_source", points, geometry.points_source),
        ("triangles", faces, geometry.triangles), ("entity_ids", labels, geometry.entity_ids))}
    array_checks["source_cell_data"] = cell_data.keys() == geometry.source_cell_data.keys() and all(
        np.array_equal(v, geometry.source_cell_data[k]) and v.dtype == geometry.source_cell_data[k].dtype
        for k, v in cell_data.items())
    if not all(array_checks.values()):
        raise GeometryError("SOURCE_MISMATCH", f"VTP 与已验收数组不一致：{array_checks}")
    identity = read_json(directory / "source_records/identity_qc.json")
    mapping = identity.get("boundary_mapping", identity)
    walls = mapping.get("wall_entity_ids", [mapping.get("wall_entity_id")])
    rows = read_boundary_manifest(directory / "source_records/boundary_manifest.csv")
    by_port = {p.port_id: p for p in geometry.patches if p.kind != "wall"}
    for row in rows:
        patch = by_port.get(row.get("port_id"))
        if patch is None or any(patch.provenance["manifest_row"].get(k) != v for k, v in row.items()):
            raise GeometryError("SOURCE_MISMATCH", "原 CSV 行与已验收端口身份不一致")
        if "display_name" in patch.provenance["manifest_row"]:
            row.setdefault("display_name", patch.provenance["manifest_row"]["display_name"])
    source_config = yaml.safe_load((directory / "import_config.yaml").read_text(encoding="utf-8-sig"))
    unit_config = yaml.safe_load((directory / "source_records/source_unit_config.yaml").read_text(encoding="utf-8-sig"))
    if unit_config.get("geometry", {}).get("input_unit") != geometry.source_length_unit:
        raise GeometryError("UNIT_EVIDENCE_CONFLICT", "已验收数组与原单位配置不一致")
    rebuilt = build_geometry(points, faces, labels, rows,
        source_length_unit=geometry.source_length_unit, unit_source=geometry.unit_source,
        wall_entity_ids=walls, identity_reference=mapping.get("boundaries"), cell_data=cell_data,
        tolerances=source_config.get("tolerances"),
        reference_center_method=source_config.get("reference_center_method", "area_weighted"))
    # Geometry and identities are exact. Recomputed floating statistics use the
    # already accepted tolerances, fixed before comparison (NumPy versions differ).
    a = {p.entity_id: p for p in geometry.patches}
    limits = rebuilt.checks["tolerances"]
    comparisons = []
    for patch in rebuilt.patches:
        old = a.get(patch.entity_id)
        if old is None:
            raise GeometryError("SOURCE_MISMATCH", "重核标签集合不一致")
        for key in ("display_name", "entity_id", "port_id", "boundary_index", "kind",
                    "original_role", "boundary_origin", "face_ids", "normal_status", "provenance"):
            if getattr(old, key) != getattr(patch, key):
                raise GeometryError("SOURCE_MISMATCH", f"重核身份不一致：{patch.port_id} / {key}")
        area_error = abs(old.area_m2 - patch.area_m2)
        center_error = float(np.linalg.norm(np.asarray(old.center_m) - patch.center_m))
        area_limit = limits["area_absolute_m2"] + limits["area_relative"] * abs(old.area_m2)
        normal_error = None
        if old.outward_normal is not None and patch.outward_normal is not None:
            normal_error = float(np.linalg.norm(np.asarray(old.outward_normal) - patch.outward_normal))
        elif old.outward_normal != patch.outward_normal:
            raise GeometryError("SOURCE_MISMATCH", "原有法向状态发生变化")
        if area_error > area_limit or center_error > limits["center_absolute_m"] or (normal_error is not None and normal_error > 1e-12):
            raise GeometryError("SOURCE_MISMATCH", f"端口几何统计与已验收数据不符：{patch.port_id}")
        comparisons.append({"entity_id": patch.entity_id, "port_id": patch.port_id,
            "triangle_count": old.triangle_count, "area_error_m2": area_error,
            "area_allowed_error_m2": area_limit, "center_error_m": center_error,
            "center_allowed_error_m": limits["center_absolute_m"],
            "unit_normal_vector_error": normal_error, "normal_allowed_error_dimensionless": 1e-12,
            "null_reason": "管壁没有全局唯一法向" if normal_error is None else None})
    checks = rebuilt.checks
    checks.update(source_vtp_vs_accepted_arrays={"status": "PASS", "exact_checks": array_checks},
                  accepted_patch_statistics={"status": "PASS", "comparisons": comparisons},
                  unit_evidence={"status": "PASS", **geometry.unit_source},
                  source_geometry_unchanged={"status": "PASS", "point_order_exact": True,
                    "triangle_order_and_connectivity_exact": True, "labels_exact": True})
    checks["surface_processing"]["only_transform"] = "NONE: byte-copy existing points_source and points_m; no second unit conversion"
    checks["historical_records"] = {
        "import_checks.json": "HISTORICAL_ONLY; kept byte-for-byte; not this migration's checks",
        "source_records/identity_qc.json": "HISTORICAL_ONLY; identity and reference statistics rechecked",
        "historical_review": "Original PENDING remains unchanged; later user acceptance is separately recorded"}
    return checks


def migrate_package(config, *, progress=print):
    from .export import copy_exclusive, load_package, new_run_id, _safe_error_details
    directory = None
    reserved = False
    before = {}
    def log(message):
        if progress is not None:
            progress(message)
    try:
        source_manifest, copies, before, mappings = collect_inputs(config)
        identity = {"input_sha256": before, "code_sha256": code_identity(),
                    "dependency_versions": dependency_identity()}
        fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
        source_geometry = load_package(config["source_package"])
        for candidate in sorted(config["output_root"].glob("*/migration_manifest.json")):
            safe_output(candidate.parent, config["source_package"])
            try:
                manifest = read_json(candidate)
                if manifest.get("fingerprint") != fingerprint:
                    continue
                loaded = load_package(candidate.parent)
                assert_same_geometry(source_geometry, loaded)
                if {p: sha256_file(require_file(Path(p))) for p in before} != before:
                    raise GeometryError("SOURCE_CHANGED", "复用核查期间输入发生变化")
                log(f"[复用] 输入、配置、相关代码及依赖匹配，目标哈希与源目标数组核查通过：{candidate.parent}")
                return ImportResult(manifest["run_id"], candidate.parent, loaded, manifest)
            except (GeometryError, OSError, ValueError) as exc:
                log(f"[警告] 不复用损坏或不完整的旧结果；保留原目录：{candidate.parent}；{exc}")
        run_id = new_run_id().replace("import_", "migrate_", 1)
        directory = safe_output(config["output_root"] / run_id, config["source_package"], before)
        directory.mkdir(parents=True, exist_ok=False)
        reserved = True
        write_json(directory / "migration_intent.json", {"run_id": run_id,
            "source_run_id": config["source_run_id"], "accept_only_with": "migration_manifest.json status PASS"})
        log("[1/3] 已验收数据包及原始来源哈希一致；复用现有检查器核对数组、身份和端口统计")
        checks = validate_accepted_package(config["source_package"], source_geometry)
        log("[2/3] 按字节复制原数据包、配套 STL 与历史记录；保留原清单及 Windows 路径")
        for relative, source in sorted(copies.items()):
            destination = directory / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            copy_exclusive(source, destination)
            if sha256_file(destination) != before[str(source)]:
                raise GeometryError("COPY_HASH_MISMATCH", f"副本字节不一致：{relative}")
        copy_exclusive(config["config_path"], directory / "migration_config.yaml")
        loaded = load_package(directory, verify_manifest=False)
        checks["export_roundtrip"] = assert_same_geometry(source_geometry, loaded)
        after = {p: sha256_file(require_file(Path(p))) for p in before}
        if before != after:
            raise GeometryError("SOURCE_CHANGED", "处理前后源文件哈希发生变化")
        if identity["code_sha256"] != code_identity():
            raise GeometryError("SOURCE_CHANGED", "处理过程中相关代码发生变化")
        checks["input_hashes_unchanged"] = {"status": "PASS", "before": before, "after": after}
        checks["byte_exact_copies"] = {"status": "PASS", "file_count": len(copies)}
        checks["warnings"] = [
            "ASSUMED_INLET / ASSUMED_OUTLET 是继承的假定身份；未验证真实生理流向。",
            f"全部 {len(loaded.points_m)} 个原始点保留，含 {loaded.summary()['unused_point_count']} 个未引用点。",
            "独立 wall 子集在端口处开放；不补洞、不进行开放壁面的体内外判定。",
            "端口法向沿用闭合表面外向约定；不是流速，复杂自交本次未检查。"]
        checks["not_checked"] = ["complex_self_intersections", "physiological_flow_direction", "user_manual_review"]
        checks["errors"] = []
        write_json(directory / "migration_checks.json", checks)
        states = {"source_model_historical_acceptance": config["source_acceptance"]["status"],
                  "automatic_data_check": "PASS", "human_review": "PENDING",
                  "sdf": "NOT_STARTED", "particles": "NOT_STARTED", "lbm_mesh": "NOT_STARTED",
                  "numerical_boundary_conditions": "NOT_STARTED", "real_blood_flow": "NOT_STARTED",
                  "mirheo_integration": "NOT_STARTED"}
        manifest = {"schema_version": 1, "run_id": run_id, "status": "PASS",
            "created_utc": datetime.now(timezone.utc).isoformat(), "fingerprint": fingerprint,
            **identity, "source_package": str(config["source_package"]),
            "source_run_id": config["source_run_id"], "source_manifest_sha256": config["source_manifest_sha256"],
            "source_acceptance": config["source_acceptance"], "states": states,
            "tagged_surface": source_manifest["tagged_surface"], "summary": loaded.summary(),
            "path_mappings": mappings, "python_executable": sys.executable, "python_version": sys.version,
            "output_sha256": {str(p.relative_to(directory)): sha256_file(p) for p in sorted(directory.rglob("*")) if p.is_file()},
            "hash_scope": "All package payloads except this self-referential migration manifest. Original import_manifest.json unchanged.",
            "loading": "load_package(target) uses only relative files inside target; historical absolute paths are metadata",
            "geometry_copy": "byte-identical NPZ/VTP/JSON; no coordinate or topology transformation"}
        write_json(directory / "migration_manifest.json", manifest)
        final = load_package(directory)
        assert_same_geometry(source_geometry, final)
        log(f"[3/3] 迁移检查 PASS，正式加载器回读精确一致：{directory}")
        return ImportResult(run_id, directory, final, manifest)
    except Exception as exc:
        if directory is not None and not reserved:
            # A collision is never permission to append diagnostics to an old run.
            raise
        if directory is None:
            run_id = new_run_id().replace("import_", "migrate_", 1)
            directory = safe_output(config["output_root"] / run_id, config["source_package"])
            directory.mkdir(parents=True, exist_ok=False)
        failure = {"status": getattr(exc, "status", "FAIL"), "code": getattr(exc, "code", "MIGRATION_FAILED"),
            "message": str(exc), "run_id": directory.name, "package_path": str(directory),
            "input_sha256_before": before, "details": _safe_error_details(getattr(exc, "details", {})),
            "human_review": "PENDING", "sdf": "NOT_STARTED", "real_blood_flow": "NOT_STARTED"}
        write_json(directory / "migration_failure.json", failure)
        raise GeometryError(failure["code"], str(exc), status=failure["status"],
                            run_id=directory.name, package_path=str(directory)) from exc
