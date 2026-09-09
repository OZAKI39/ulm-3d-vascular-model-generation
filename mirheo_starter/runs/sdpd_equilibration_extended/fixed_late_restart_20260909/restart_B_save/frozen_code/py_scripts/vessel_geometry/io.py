"""Explicit source selection and VTP reads, with no geometric processing."""

import csv
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import subprocess

import numpy as np
import yaml

from .model import GeometryError

PROJECT_ROOT = Path(__file__).resolve().parents[2]
READ_ONLY_PROJECTS = (PROJECT_ROOT.parent / "bloodflow_starter", PROJECT_ROOT.parent / "ulm_3D_vascular")
LFS_HEADER = b"version https://git-lfs.github.com/spec/v1"
UNIT_FACTORS = {"m": 1.0, "mm": 1e-3, "um": 1e-6}


def require_file(path: Path) -> Path:
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise GeometryError("INPUT_MISSING", f"输入缺失或为空：{path}", status="BLOCKED", path=str(path))
    with path.open("rb") as stream:
        if stream.read(256).lstrip(b"\xef\xbb\xbf \r\n").startswith(LFS_HEADER):
            raise GeometryError("LFS_POINTER", f"只有 Git LFS 指针，缺少真实文件：{path}", status="BLOCKED", path=str(path))
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    def reject_constant(value):
        raise GeometryError("NONFINITE_JSON", f"JSON 含非法数值 {value}：{path}")
    payload = json.loads(require_file(path).read_text(encoding="utf-8-sig"), parse_constant=reject_constant)
    def finite(value):
        if isinstance(value, float) and not np.isfinite(value):
            reject_constant(value)
        if isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(payload)
    return payload


def write_json(path: Path, payload) -> None:
    # Exclusive writes: a repeated run cannot replace an existing result.
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def safe_output(path: Path, source_project: Path, inputs=()) -> Path:
    resolved = Path(path).resolve()
    source = Path(source_project).resolve()
    if any(resolved == p.resolve() or resolved.is_relative_to(p.resolve())
           for p in (source, *READ_ONLY_PROJECTS)):
        raise GeometryError("UNSAFE_OUTPUT", f"输出指向只读源工程（已解析符号链接）：{resolved}")
    if not resolved.is_relative_to(PROJECT_ROOT) or resolved == PROJECT_ROOT:
        raise GeometryError("UNSAFE_OUTPUT", f"输出必须位于目标工程内部：{resolved}")
    for item in inputs:
        input_path = Path(item).resolve()
        if resolved == input_path or input_path.is_relative_to(resolved):
            raise GeometryError("UNSAFE_OUTPUT", f"输出与输入重叠：{resolved} / {input_path}")
    return resolved


def load_config(path: Path) -> dict:
    path = Path(path).expanduser()
    path = (path if path.is_absolute() else PROJECT_ROOT / path).resolve()
    payload = yaml.safe_load(require_file(path).read_text(encoding="utf-8-sig"))
    if isinstance(payload, dict) and payload.get("mode") == "accepted_package":
        from .migration import load_migration_config
        return load_migration_config(path, payload)
    required = {"source_project", "source_run", "tagged_surface", "boundary_manifest",
                "source_length_unit", "source_unit_config", "identity_qc", "output_root", "review_output_root"}
    optional = {"schema_version", "companion_stl", "scale_qc", "source_selection_config",
                "reference_code", "tolerances", "reference_center_method"}
    if not isinstance(payload, dict) or required - payload.keys():
        raise GeometryError("CONFIG_INVALID", f"配置字段缺失：{path}；需要 {sorted(required)}")
    if payload.keys() - required - optional:
        raise GeometryError("CONFIG_INVALID", f"未知配置字段：{sorted(payload.keys() - required - optional)}")
    if payload.get("schema_version", 1) != 1:
        raise GeometryError("CONFIG_INVALID", "仅支持 schema_version: 1")
    config = dict(payload)
    for key in required:
        if not isinstance(payload[key], str) or not payload[key].strip():
            raise GeometryError("CONFIG_INVALID", f"配置 {key} 必须为非空字符串：{path}")
    if not isinstance(payload.get("tolerances", {}), dict) or not isinstance(payload.get("reference_code", []), list):
        raise GeometryError("CONFIG_INVALID", "tolerances 必须为映射，reference_code 必须为列表")
    config["config_path"] = path
    def resolve(value, base):
        value = Path(str(value)).expanduser()
        if PureWindowsPath(str(value)).drive:
            raise GeometryError("CONFIG_INVALID", f"配置必须使用本地 WSL 路径：{value}")
        return (value if value.is_absolute() else base / value).resolve()
    config["source_project"] = resolve(payload["source_project"], PROJECT_ROOT)
    config["source_run"] = resolve(payload["source_run"], config["source_project"])
    if not config["source_run"].is_relative_to(config["source_project"]):
        raise GeometryError("SOURCE_MISMATCH", "source_run 必须位于 source_project 中")
    for key in ("tagged_surface", "boundary_manifest", "source_unit_config", "identity_qc", "companion_stl", "scale_qc"):
        if payload.get(key):
            config[key] = resolve(payload[key], config["source_run"])
            if not config[key].is_relative_to(config["source_run"]):
                raise GeometryError("SOURCE_MISMATCH", f"{key} 不属于所选同次运行：{config[key]}")
        else:
            config[key] = None
    config["source_selection_config"] = (
        resolve(payload["source_selection_config"], config["source_project"])
        if payload.get("source_selection_config") else None
    )
    config["reference_code"] = [resolve(p, config["source_project"]) for p in payload.get("reference_code", [])]
    inputs = [config[k] for k in ("tagged_surface", "boundary_manifest", "identity_qc", "source_unit_config")]
    for key in ("output_root", "review_output_root"):
        config[key] = safe_output(resolve(payload[key], PROJECT_ROOT), config["source_project"], inputs)
    a, b = config["output_root"], config["review_output_root"]
    if a.is_relative_to(b) or b.is_relative_to(a):
        raise GeometryError("UNSAFE_OUTPUT", "正式数据与核查报告的根目录必须相互独立")
    unit_factor(config["source_length_unit"])
    config["tolerances"] = {"area_absolute_m2": 1e-18, "area_relative": 1e-10,
                            "center_absolute_m": 5e-11, "degenerate_area_m2": 1e-28,
                            "volume_m3": 1e-30, **payload.get("tolerances", {})}
    for name, value in config["tolerances"].items():
        if name not in {"area_absolute_m2", "area_relative", "center_absolute_m", "degenerate_area_m2", "volume_m3"}:
            raise GeometryError("CONFIG_INVALID", f"未知容差：{name}")
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value) or value < 0:
            raise GeometryError("CONFIG_INVALID", f"容差必须是有限非负数：{name}={value}")
    method = config.get("reference_center_method", "area_weighted")
    if method not in {"area_weighted", "triangle_centroid_mean_source_dtype"}:
        raise GeometryError("CONFIG_INVALID", f"未知来源中心定义：{method}")
    config["reference_center_method"] = method
    return config


def unit_factor(unit: str) -> float:
    if not isinstance(unit, str) or unit not in UNIT_FACTORS:
        raise GeometryError("UNKNOWN_UNIT", f"未知长度单位 {unit!r}；必须明确指定 m / mm / um")
    return UNIT_FACTORS[unit]


def resolve_recorded_path(recorded: str, source_run: Path) -> Path:
    """Relocate only a path with an explicit, unique matching run component."""
    win = PureWindowsPath(recorded)
    if win.drive:
        matches = [i for i, part in enumerate(win.parts) if part == source_run.name]
        if len(matches) != 1:
            raise GeometryError("HISTORICAL_PATH_AMBIGUOUS", f"历史路径无法唯一关联本次运行：{recorded}", status="BLOCKED")
        parts = win.parts[matches[0] + 1:]
        if not parts or ".." in parts:
            raise GeometryError("HISTORICAL_PATH_AMBIGUOUS", f"历史路径非法：{recorded}", status="BLOCKED")
        result = source_run.joinpath(*parts).resolve()
    else:
        candidate = Path(recorded)
        result = (candidate if candidate.is_absolute() else source_run / candidate).resolve()
    if not result.is_relative_to(source_run.resolve()):
        raise GeometryError("SOURCE_MISMATCH", f"清单路径不属于同次运行：{recorded} -> {result}")
    return require_file(result)


def read_boundary_manifest(path: Path) -> list[dict]:
    with require_file(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise GeometryError("MANIFEST_INVALID", f"CSV 表头缺失或重复：{path}")
        rows = list(reader)
    if not rows:
        raise GeometryError("MANIFEST_INVALID", f"端口清单为空：{path}")
    return rows


def read_tagged_vtp(path: Path):
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
    reader = vtkXMLPolyDataReader()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *_: errors.append("VTK ErrorEvent"))
    reader.SetFileName(str(require_file(path)))
    reader.Update()
    surface = reader.GetOutput()
    if errors or reader.GetErrorCode() or surface.GetPoints() is None:
        raise GeometryError("VTP_READ_FAILED", f"VTP 读取失败：{path}；{errors}；VTK={reader.GetErrorCode()}")
    if surface.GetNumberOfVerts() or surface.GetNumberOfLines() or surface.GetNumberOfStrips():
        raise GeometryError("NON_TRIANGLE_CELLS", f"仅支持纯三角表面，不支持顶点/线/条带单元：{path}")
    polygons = surface.GetPolys()
    offsets = vtk_to_numpy(polygons.GetOffsetsArray())
    flat = vtk_to_numpy(polygons.GetConnectivityArray())
    if len(offsets) < 2 or offsets[0] != 0 or offsets[-1] != len(flat) or not np.all(np.diff(offsets) == 3):
        raise GeometryError("NON_TRIANGLE_CELLS", f"包含非三角单元或非法连接，不自动三角化：{path}")
    points = vtk_to_numpy(surface.GetPoints().GetData()).copy()
    faces = flat.reshape(-1, 3).copy()
    data = surface.GetCellData()
    label = data.GetArray("CellEntityIds")
    if label is None:
        raise GeometryError("MISSING_LABELS", f"缺少 CellEntityIds 单元标签：{path}")
    labels = vtk_to_numpy(label).copy()
    cell_data = {}
    for i in range(data.GetNumberOfArrays()):
        name = data.GetArrayName(i)
        array = data.GetAbstractArray(i)
        if array.IsA("vtkStringArray"):
            if array.GetNumberOfComponents() != 1:
                raise GeometryError("ARRAY_INVALID", f"不支持多分量字符串数组：{name}")
            values = np.asarray([array.GetValue(j) for j in range(array.GetNumberOfValues())], dtype=str)
        elif data.GetArray(i) is not None:
            values = vtk_to_numpy(data.GetArray(i)).copy()
        else:
            raise GeometryError("ARRAY_INVALID", f"不支持的单元属性：{name}")
        if len(values) != len(faces):
            raise GeometryError("LABEL_LENGTH", f"单元属性 {name} 长度 {len(values)} 与面数 {len(faces)} 不一致")
        cell_data[name] = values
    return points, faces, labels, cell_data


def git_identity(project: Path, paths: list[Path]) -> dict:
    env = {**os.environ, "GIT_OPTIONAL_LOCKS": "0"}
    def git(*args):
        result = subprocess.run(["git", "-C", str(project), *args], env=env,
                                capture_output=True, timeout=20, check=False)
        if result.returncode:
            raise RuntimeError(result.stderr.decode(errors="replace").strip())
        return result.stdout
    try:
        result = {"commit": git("rev-parse", "HEAD").decode().strip(),
                  "branch": git("rev-parse", "--abbrev-ref", "HEAD").decode().strip(), "files": {}}
        for path in paths:
            if not path.is_relative_to(project):
                continue
            relative = str(path.relative_to(project))
            tracked = bool(git("ls-files", "--", relative).strip())
            row = {"tracked": tracked, "actual_sha256": sha256_file(path)}
            if tracked:
                blob = git("show", ":" + relative)
                if blob.startswith(LFS_HEADER):
                    expected = next(line.split(":", 1)[1] for line in blob.decode().splitlines() if line.startswith("oid sha256:"))
                    row.update(index_is_lfs_pointer=True, index_content_sha256=expected,
                               content_differs_from_index=row["actual_sha256"] != expected)
                else:
                    expected = hashlib.sha256(blob).hexdigest()
                    row.update(index_is_lfs_pointer=False, index_content_sha256=expected,
                               content_differs_from_index=row["actual_sha256"] != expected)
            result["files"][relative] = row
        return {"status": "AVAILABLE", **result,
                "comparison": "实际文件字节与 Git 暂存区内容 SHA-256 对照；LFS 使用 oid；未跟踪文件单独标记"}
    except (OSError, RuntimeError, subprocess.TimeoutExpired, StopIteration) as exc:
        return {"status": "UNAVAILABLE", "reason": str(exc)}
