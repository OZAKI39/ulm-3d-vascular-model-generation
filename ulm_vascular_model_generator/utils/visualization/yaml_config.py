"""
可视化 YAML 配置解析。

生成器配置和可视化配置分开，是为了让“生成什么血管”和“怎么显示血管”
互不影响。这个模块把相对路径、默认输出文件名、渲染选项统一解析成
`VisualizationConfig`。
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    # PyYAML 用来读取 visualization_config.yaml。
    import yaml
except ModuleNotFoundError as exc:  # pragma: no cover - import-time guidance
    # 缺少依赖时在导入阶段就给出明确安装提示。
    raise RuntimeError("PyYAML is required. Install it in the pmp environment.") from exc

from .cylinder_mesh import CapStyle, RadiusStyle
from .pyvista_viewer import ViewMode
from .swc_loader import list_swc_files

# PACKAGE_ROOT 指向 ulm_vascular_model_generator 包目录。
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
# CONFIG_DIR 保存默认 YAML 配置。
CONFIG_DIR = PACKAGE_ROOT / "configs"
# 默认可视化配置文件路径。
DEFAULT_CONFIG_PATH = CONFIG_DIR / "visualization_config.yaml"
# 默认从生成器输出目录读取 SWC。
DEFAULT_SWC_DIR = PACKAGE_ROOT / "vessel_swc_models"
# 默认把截图和 mesh 导出到 vis_results。
DEFAULT_RESULT_DIR = DEFAULT_SWC_DIR / "vis_results"
# 只允许这三个顶层 section，防止 YAML 拼写错误被悄悄忽略。
TOP_LEVEL_KEYS = {"input", "output", "render"}


@dataclass(frozen=True)
class VisualizationConfig:
    """已经解析好的可视化设置，路径都转换为绝对路径。"""

    # SWC 搜索目录。
    swc_dir: Path
    # 实际要渲染的 SWC 文件。
    swc_path: Path
    # 结果自身的几何类型；显式结果目录优先读取 generation_validation.json，
    # 旧结果则从 SWC 坐标判断。
    geometry_mode: str
    # 截图和 mesh 的默认输出目录。
    result_dir: Path
    # 截图路径；为 None 表示不保存截图。
    screenshot_path: Path | None
    # mesh 导出路径；为 None 表示不导出 mesh。
    mesh_path: Path | None
    # 半径显示风格：整段同半径或父子端渐变。
    radius_style: RadiusStyle
    # 圆柱截面边数，越大越圆但面片越多。
    sides: int
    # 可视化半径缩放，仅影响显示，不改 SWC 原始值。
    radius_scale: float
    # 圆柱端盖策略。
    cap_style: CapStyle
    # 统一着色时使用的颜色。
    color: str
    # 标量着色时使用的 Matplotlib/PyVista colormap。
    color_map: str
    # 是否按半径标量着色。
    color_by_radius: bool
    # 渲染背景色。
    background: str
    # 渐变背景顶部颜色；为 None 表示纯色背景。
    background_top: str | None
    # 是否显示 mesh 边线。
    show_edges: bool
    # 是否显示坐标轴。
    show_axes: bool
    # 是否在左上角显示标题。
    show_title: bool
    # 是否显示标量条。
    scalar_bar: bool
    publication_style: bool
    font_size: int
    # 是否启用 PyVista/VTK 的 PBR 材质。
    pbr: bool
    # PBR/Phong 材质参数。
    metallic: float
    roughness: float
    ambient: float
    diffuse: float
    specular: float
    specular_power: float
    opacity: float
    # 是否使用自定义三点布光。
    lighting: bool
    # 是否启用 Eye-Dome Lighting 深度增强。
    eye_dome_lighting: bool
    # 是否启用阴影。
    shadows: bool
    # 抗锯齿模式：none/msaa/ssaa/fxaa。
    anti_aliasing: str
    # 相机视角。
    view: ViewMode
    # 窗口尺寸，格式为 (width, height)。
    window_size: tuple[int, int]
    # 是否打开交互窗口。
    show_window: bool
    # 是否使用离屏渲染。
    off_screen: bool


def _read_yaml(path: Path) -> dict[str, Any]:
    """读取可视化 YAML，并检查顶层 key 是否合法。"""
    # safe_load 不执行任意 Python 对象，比普通 load 更适合配置文件。
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    # 顶层必须是 mapping，例如 input/output/render 三个 section。
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping.")
    # 检查未知顶层字段，帮助用户发现拼写错误。
    unknown = set(data) - TOP_LEVEL_KEYS
    if unknown:
        raise ValueError(f"Unknown top-level visualization config key(s): {sorted(unknown)}")
    # 返回原始字典，后续函数会继续解析每个 section。
    return data


def _section(data: dict[str, Any], name: str) -> dict[str, Any]:
    """读取一个 YAML section；缺失时返回空字典。"""
    # 缺少 section 时使用空字典，让每个字段都走默认值。
    raw = data.get(name, {})
    # YAML 里写 section: null 时也按空配置处理。
    if raw is None:
        return {}
    # section 必须是键值对，不能写成列表或字符串。
    if not isinstance(raw, dict):
        raise ValueError(f"'{name}' must be a YAML mapping.")
    # 返回 section 字典供调用方读取字段。
    return raw


def _resolve_path(value: Any, *, default: Path, base: Path = PACKAGE_ROOT) -> Path:
    """把 YAML 中的相对路径转换为包目录下的绝对路径。"""
    # None 或空字符串表示使用调用方提供的默认路径。
    if value in (None, ""):
        return default
    # YAML 中的路径先统一转为 Path。
    path = Path(str(value))
    # 绝对路径原样返回，允许用户把输出放到任意目录。
    if path.is_absolute():
        return path
    # 相对路径默认相对于包目录，避免受当前工作目录影响。
    return base / path


def _newest_swc_file(directory: Path) -> Path:
    """当 swc_file 为 null 时，自动选择目录中最新的 SWC。"""
    # 先列出目录中的所有 SWC。
    files = list_swc_files(directory)
    # 如果目录为空，不能猜测输入文件，直接报错。
    if not files:
        raise FileNotFoundError(f"No .swc files found under {directory}.")
    # 按修改时间选择最新文件，适合刚运行完生成器后的默认可视化。
    return max(files, key=lambda path: path.stat().st_mtime)


def _resolve_swc_file(value: Any, directory: Path) -> Path:
    """解析用户指定的 SWC 文件路径。"""
    # 没有显式指定文件时，使用目录中最新的 SWC。
    if value in (None, ""):
        return _newest_swc_file(directory)

    # 用户可以写绝对路径、当前工作目录相对路径或 swc_dir 内文件名。
    path = Path(str(value))
    # 绝对路径且存在时直接使用。
    if path.is_absolute() and path.exists():
        return path
    # 当前工作目录下能找到时也直接使用。
    if path.exists():
        return path

    # 最常见情况：YAML 里只写文件名，需要拼到 swc_dir 下。
    in_directory = directory / path
    if in_directory.exists():
        return in_directory

    # 三种解析方式都失败时，说明配置指向了不存在的文件。
    raise FileNotFoundError(
        f"Cannot find SWC file {value!r}. Tried {path} and {in_directory}."
    )


def detect_result_geometry_mode(
    result_folder: Path,
    swc_path: Path | None = None,
) -> str:
    """
    Return ``planar_2d`` or ``volumetric_3d`` for one generator result.

    New result folders contain an authoritative ``generation_validation.json``.
    Older folders are still supported by inspecting the SWC coordinate span;
    folder/file naming is used only as a final fallback.
    """
    result_folder = Path(result_folder)
    validation_path = result_folder / "generation_validation.json"
    if validation_path.exists():
        try:
            payload = json.loads(validation_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}
        mode = payload.get("geometry_mode")
        if mode in {"planar_2d", "volumetric_3d"}:
            return str(mode)

    if swc_path is not None and Path(swc_path).exists():
        coordinates = []
        with Path(swc_path).open("r", encoding="ascii") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split()
                if len(fields) >= 5:
                    coordinates.append(
                        [float(fields[2]), float(fields[3]), float(fields[4])]
                    )
        if coordinates:
            xyz = np.asarray(coordinates, dtype=float)
            spans = np.ptp(xyz, axis=0)
            reference_span = max(float(spans[0]), float(spans[2]), 1.0)
            return (
                "planar_2d"
                if float(spans[1]) <= 1.0e-6 * reference_span
                else "volumetric_3d"
            )

    naming_hint = " ".join(
        (
            result_folder.name.lower(),
            "" if swc_path is None else Path(swc_path).name.lower(),
        )
    )
    if "volumetric_3d" in naming_hint or "xyz" in naming_hint:
        return "volumetric_3d"
    if "planar_2d" in naming_hint or "planar" in naming_hint or "xz" in naming_hint:
        return "planar_2d"
    raise ValueError(
        "Cannot determine whether the selected result is planar_2d or "
        f"volumetric_3d: {result_folder}"
    )


def _output_path(
    value: Any,
    *,
    enabled: bool,
    default_name: str,
    result_dir: Path,
) -> Path | None:
    """根据开关和默认文件名解析截图/mesh 输出路径。"""
    # 输出开关关闭时返回 None，渲染层据此跳过保存。
    if not enabled:
        return None
    # 没写文件名时，用 result_dir 加自动生成的默认文件名。
    if value in (None, ""):
        return result_dir / default_name
    # 显式文件名或路径先转成 Path。
    path = Path(str(value))
    # 绝对路径原样使用。
    if path.is_absolute():
        return path
    # 相对路径默认放到 result_dir 下。
    return result_dir / path


def _choice(value: Any, choices: set[str], name: str) -> str:
    """校验字符串枚举值。"""
    # 枚举值必须是字符串，列表/数字通常是配置写错。
    if not isinstance(value, str):
        raise ValueError(f"{name} must be one of {sorted(choices)}.")
    # 字符串必须属于允许集合。
    if value not in choices:
        raise ValueError(f"{name} must be one of {sorted(choices)}, got {value!r}.")
    # 返回原值；调用方可通过 type ignore 标注更窄的 Literal 类型。
    return value


def _window_size(value: Any) -> tuple[int, int]:
    """校验 PyVista 窗口大小。"""
    # window_size 必须写成 [width, height] 两项。
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError("render.window_size must be a two-item list: [width, height].")
    # 转成 int，允许 YAML 中写成 1280.0 这类数值。
    width, height = int(value[0]), int(value[1])
    # 尺寸必须为正，0 或负数会让窗口创建失败。
    if width <= 0 or height <= 0:
        raise ValueError("render.window_size values must be positive.")
    # 返回固定二元组。
    return width, height


def _float_range(value: Any, name: str, lower: float, upper: float) -> float:
    """读取浮点参数，并限制在闭区间内。"""
    result = float(value)
    if not lower <= result <= upper:
        raise ValueError(f"{name} must be in [{lower}, {upper}], got {result}.")
    return result


def _optional_string(value: Any) -> str | None:
    """把 YAML 中的空值规范化为 None，其它值转为字符串。"""
    if value in (None, ""):
        return None
    return str(value)


def load_visualization_config(
    path: Path = DEFAULT_CONFIG_PATH,
    result_folder=None,
) -> VisualizationConfig:
    """读取并解析完整可视化配置。"""
    # 统一配置路径类型。
    path = Path(path)
    # 读取 YAML 原始字典。
    data = _read_yaml(path)
    # 分别取出 input/output/render 三个 section。
    input_cfg = _section(data, "input")
    output_cfg = _section(data, "output")
    render_cfg = _section(data, "render")

    # 解析 SWC 目录；默认指向生成器输出目录。
    swc_dir = _resolve_path(input_cfg.get("swc_dir"), default=DEFAULT_SWC_DIR).resolve()
    # 解析实际 SWC 文件；未指定时自动选择最新文件。
    # 解析输出目录；截图和 mesh 的默认文件名都会放在这里。
    if result_folder is None:
        swc_path = _resolve_swc_file(
            input_cfg.get("swc_file"),
            swc_dir,
        ).resolve()
        result_dir = _resolve_path(
            output_cfg.get("result_dir"),
            default=DEFAULT_RESULT_DIR,
        ).resolve()
    else:
        result_dir = _resolve_path(
            result_folder,
            default=DEFAULT_RESULT_DIR,
        ).resolve()
        swc_dir = result_dir
        swc_path = _newest_swc_file(swc_dir).resolve()

    geometry_mode = detect_result_geometry_mode(result_dir, swc_path)

    # 解析半径显示风格。
    radius_style = _choice(
        render_cfg.get("radius_style", "cylindrical"),
        {"cylindrical", "tapered"},
        "render.radius_style",
    )
    # 解析端盖策略。
    cap_style = _choice(
        render_cfg.get("cap_style", "terminal"),
        {"none", "terminal", "all"},
        "render.cap_style",
    )
    publication_style = bool(render_cfg.get("publication_style", False))
    # 解析相机视角。显式选择一个生成结果并启用论文图时，根据结果自身的
    # 二维/三维类型自动使用正交 X-Z 或三维等轴视图。
    configured_view = _choice(
        render_cfg.get("view", "xz"),
        {"xz", "xy", "yz", "iso"},
        "render.view",
    )
    if result_folder is not None and publication_style:
        view = "xz" if geometry_mode == "planar_2d" else "iso"
    else:
        view = configured_view

    # sides 控制圆柱截面精细程度。
    sides = int(render_cfg.get("sides", 18))
    if sides < 6:
        raise ValueError("render.sides must be at least 6.")
    # radius_scale 只影响可视化显示，不回写 SWC。
    radius_scale = float(render_cfg.get("radius_scale", 1.0))
    if radius_scale <= 0.0:
        raise ValueError("render.radius_scale must be positive.")

    font_size = int(render_cfg.get("font_size", 18))
    if font_size < 8:
        raise ValueError("render.font_size must be at least 8.")
    screenshot_suffix = f"{radius_style}_{view}"
    mesh_suffix = radius_style
    if publication_style:
        if geometry_mode == "planar_2d":
            screenshot_suffix = "planar_xz_publication"
        else:
            screenshot_suffix = "volumetric_3d_publication"
        mesh_suffix += "_publication"

    # 读取输出开关。
    save_screenshot = bool(output_cfg.get("save_screenshot", True))
    export_mesh = bool(output_cfg.get("export_mesh", True))
    # 解析截图路径；默认文件名带上半径风格和视角，避免覆盖不同显示模式。
    screenshot_path = _output_path(
        output_cfg.get("screenshot_name"),
        enabled=save_screenshot,
        default_name=f"{swc_path.stem}_{screenshot_suffix}.png",
        result_dir=result_dir,
    )
    # 解析 mesh 输出路径；默认使用 .vtp，便于 PyVista/ParaView 打开。
    mesh_path = _output_path(
        output_cfg.get("mesh_name"),
        enabled=export_mesh,
        default_name=f"{swc_path.stem}_{mesh_suffix}.vtp",
        result_dir=result_dir,
    )

    # 最终把所有字段放进一个不可变配置对象，后续代码只读这个对象。
    return VisualizationConfig(
        swc_dir=swc_dir,
        swc_path=swc_path,
        geometry_mode=geometry_mode,
        result_dir=result_dir,
        screenshot_path=screenshot_path,
        mesh_path=mesh_path,
        radius_style=radius_style,  # type: ignore[arg-type]
        sides=sides,
        radius_scale=radius_scale,
        cap_style=cap_style,  # type: ignore[arg-type]
        color=str(render_cfg.get("color", "#0b37b7")),
        color_map=str(render_cfg.get("color_map", render_cfg.get("cmap", "turbo"))),
        color_by_radius=bool(render_cfg.get("color_by_radius", False)),
        background=str(render_cfg.get("background", "#d7e4ec")),
        background_top=_optional_string(render_cfg.get("background_top", "#ffffff")),
        show_edges=bool(render_cfg.get("show_edges", False)),
        show_axes=bool(render_cfg.get("show_axes", True)),
        show_title=bool(render_cfg.get("show_title", True)),
        scalar_bar=bool(render_cfg.get("scalar_bar", True)),
        publication_style=publication_style,
        font_size=font_size,
        pbr=bool(render_cfg.get("pbr", False)),
        metallic=_float_range(render_cfg.get("metallic", 0.0), "render.metallic", 0.0, 1.0),
        roughness=_float_range(render_cfg.get("roughness", 0.55), "render.roughness", 0.0, 1.0),
        ambient=_float_range(render_cfg.get("ambient", 0.34), "render.ambient", 0.0, 1.0),
        diffuse=_float_range(render_cfg.get("diffuse", 0.82), "render.diffuse", 0.0, 1.0),
        specular=_float_range(render_cfg.get("specular", 0.18), "render.specular", 0.0, 1.0),
        specular_power=_float_range(render_cfg.get("specular_power", 45.0), "render.specular_power", 1.0, 128.0),
        opacity=_float_range(render_cfg.get("opacity", 1.0), "render.opacity", 0.0, 1.0),
        lighting=bool(render_cfg.get("lighting", True)),
        eye_dome_lighting=bool(render_cfg.get("eye_dome_lighting", False)),
        shadows=bool(render_cfg.get("shadows", False)),
        anti_aliasing=_choice(
            str(render_cfg.get("anti_aliasing", "msaa")).lower(),
            {"none", "msaa", "ssaa", "fxaa"},
            "render.anti_aliasing",
        ),
        view=view,  # type: ignore[arg-type]
        window_size=_window_size(render_cfg.get("window_size", [1920, 1350])),
        show_window=bool(render_cfg.get("show_window", True)),
        off_screen=bool(render_cfg.get("off_screen", False)),
    )
