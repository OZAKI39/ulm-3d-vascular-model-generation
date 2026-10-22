"""YAML-driven SWC preprocessing and representative connected-ROI generation.

With no command-line argument, this entry point reads
``configs/swc_roi_generate.yaml``.  A different YAML file may be supplied as
the sole positional argument.  All scientific and runtime hyperparameters
live in YAML; this file only orchestrates validated processing stages.

Orange arrows encode the SWC parent-to-current relation and do not claim a
measured blood-flow direction.
"""

from __future__ import annotations

import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path

from utils.rodent_vasculature import run_rodent_vasculature_pipeline
from utils.rodent_vasculature.interactive import show_saved_run
from utils.sampling.pipeline import run_sampling_from_rodent_run
from utils.swc_roi_yaml_config import SWCROIRunConfig, load_swc_roi_yaml_config


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_roi_generate.yaml"
USAGE = (
    "Usage:\n"
    "  python s1_swc_roi_generate.py\n"
    "  python s1_swc_roi_generate.py <config.yaml>\n\n"
    "All processing and sampling hyperparameters are defined in YAML.\n"
    f"Default configuration: {DEFAULT_CONFIG}"
)


_PIPELINE_SECTIONS = {
    "inventory": "4.3.1 数据编目与样本选择",
    "preprocess": "4.3.1～4.3.3 数据编目、输入检查、尺度统一与连通性预处理",
    "hierarchical-graph": "4.3.4 从已保存的分析网络构建血管层级图",
    "all": "4.3.1～4.3.4 数据编目、SWC 预处理与血管层级图构建",
}


class _DocumentStageFormatter(logging.Formatter):
    """Label existing console messages without changing the file log records."""

    _sections = (
        (("Starting rodent", "Direction convention:", "Cataloged"),
         "4.3.1 整体处理流程与运行约定"),
        (("SWC-centric preprocessing",), "4.3.3 主血管网络的连通性预处理"),
        (("Preprocessed", "Preprocessing failed"),
         "4.3.2～4.3.3 输入检查、尺度统一与连通性预处理"),
        (("Built directed graph", "Graph construction failed"),
         "4.3.4 血管层级图构建与结构方向表达"),
        (("Wrote",), "4.3.6～4.3.7 可视化图像保存"),
        (("Sampling scope:", "Input models:", "Global node count=", "ROI size um="),
         "4.3.5～4.3.5.2 采样输入 / 候选 ROI、边界语义与特征计算"),
        (("Candidate anchors=",), "4.3.5～4.3.5.2 候选生成与特征筛选结果"),
        (("Feature mode=",), "4.3.5.2 血管尺度与结构特征"),
        (("Clustering method=", "Selection mode="),
         "4.3.5.3 聚类与真实代表区域选择"),
        (("Validation metrics=", "Runtime seconds=", "Output directory=",
          "Run finished:", "Sampling run finished:", "Model loading seconds="),
         "4.3.7 结果保存与功能验收"),
    )

    def format(self, record: logging.LogRecord) -> str:
        message = super().format(record)
        section = next(
            (section for prefixes, section in self._sections if message.startswith(prefixes)),
            None,
        )
        timestamp = self.formatTime(record, "%H:%M:%S")
        label = f"文档 {section}" if section is not None else "运行日志"
        message = message.replace("\n", "\n    ")
        return (
            f"[{timestamp}] [{record.levelname}] [{label}]\n"
            f"    {message}"
        )


class _ConsoleStageLogFilter(logging.Filter):
    """Reapply console formatting when either pipeline recreates its handlers."""

    def __init__(self) -> None:
        super().__init__()
        self.formatter = _DocumentStageFormatter()

    def filter(self, record: logging.LogRecord) -> bool:
        for handler in logging.getLogger("ulm_3d_vascular").handlers:
            if isinstance(handler, logging.StreamHandler) and not isinstance(
                handler, logging.FileHandler
            ):
                handler.setFormatter(self.formatter)
        return True


def _configure_stage_console_logging() -> None:
    logger = logging.getLogger("ulm_3d_vascular")
    if not any(isinstance(item, _ConsoleStageLogFilter) for item in logger.filters):
        logger.addFilter(_ConsoleStageLogFilter())


def _print_stage(section: str, title: str, status: str = "开始") -> None:
    print(
        f"\n{'=' * 72}\n"
        f"[{datetime.now():%H:%M:%S}] [{status}] [文档 {section}] {title}\n"
        f"{'-' * 72}",
        flush=True,
    )


def _print_path(label: str, path: Path, purpose: str) -> None:
    print(f"  {label}：{path}\n    用途：{purpose}", flush=True)


def _print_status(label: str, status: str) -> None:
    meaning = {"PASS": "通过", "WARNING": "有警告，请查看验收详情", "FAIL": "未通过"}
    print(f"  {label}：{status}（{meaning.get(status, status)}）", flush=True)


def _configuration_path(argv: list[str]) -> Path | None:
    if not argv:
        return DEFAULT_CONFIG
    if len(argv) == 1 and argv[0] in {"-h", "--help"}:
        print(USAGE)
        return None
    if len(argv) != 1 or argv[0].startswith("-"):
        raise ValueError("Expected no argument or one YAML configuration path")
    candidate = Path(argv[0]).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (Path.cwd() / candidate).resolve()


def _copy_source_configuration(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def _load_settings(argv: list[str]) -> SWCROIRunConfig | None:
    config_path = _configuration_path(argv)
    if config_path is None:
        return None
    return load_swc_roi_yaml_config(config_path, project_root=PROJECT_ROOT)


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        settings = _load_settings(arguments)
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        print(f"[错误] [运行准备 / YAML 配置读取] {exc}", file=sys.stderr, flush=True)
        print(USAGE, file=sys.stderr)
        return 1
    if settings is None:
        return 0

    config = settings.rodent
    graph_stage = config.stage in {"all", "hierarchical-graph"}
    _configure_stage_console_logging()
    _print_stage("4.3.1", "运行配置与阶段说明", "准备")
    _print_path("输入配置（Configuration）", settings.source_path, "本次实际读取的 YAML 参数文件。")
    print(
        "  文档：references/文本记录 - SWC 预处理与 ROI 生成.md\n"
        "  日志中的章节号对应上述文档；结果保存与验收（4.3.7）在各阶段分别执行。\n"
        f"  全局处理范围：pipeline.stage={config.stage}；{_PIPELINE_SECTIONS[config.stage]}\n"
        f"  ROI 采样开关：{'启用' if settings.sampling_enabled else '关闭'}；"
        f"双视窗开关：{'启用' if config.figure2a_enabled else '关闭'}；"
        f"显示方式：{'交互窗口' if settings.show_gui else '屏幕外预览'}",
        flush=True,
    )
    _print_stage(_PIPELINE_SECTIONS[config.stage], "全局 SWC 处理")
    try:
        run = run_rodent_vasculature_pipeline(config, verbose=settings.verbose)
        copied_config = _copy_source_configuration(
            settings.source_path,
            run.run_root / "source_swc_roi_generate.yaml",
        )
    except Exception as exc:
        print(
            f"[错误] [文档 {_PIPELINE_SECTIONS[config.stage]} / 4.3.7 源配置保存] {exc}",
            file=sys.stderr, flush=True,
        )
        return 1

    _print_stage("4.3.7", "全局 SWC 处理：输出文件与验收结果", "结果")
    _print_path("全局结果目录（Run directory）", run.run_root, "保存本次全局处理的清单、样本数据、图结构和验收文件。")
    _print_path("源配置副本（Saved source configuration）", copied_config, "保存输入 YAML 的原样副本，便于追溯和复现本次参数。")
    _print_path("全局有效配置", run.run_root / "run_config.json", "保存解析后的全局处理参数与实际路径。")
    _print_path("全局验收报告（Acceptance report）", run.html_report, "用浏览器打开，查看本次全局处理的检查项目和通过 / 警告 / 失败原因。")
    _print_path("全局详细日志", run.run_root / "pipeline.log", "查看全局 SWC 处理的详细过程和异常堆栈。")
    _print_status("全局验收状态（Acceptance status）", run.acceptance.overall_status)
    print(f"  验收范围：本次 pipeline.stage={config.stage} 的全局处理；ROI 采样另行验收。", flush=True)

    sampling_run = None
    if run.status != "failed" and graph_stage and settings.sampling_enabled:
        _print_stage("4.3.5", "真实连通局部区域的代表性采样")
        print(
            "  处理顺序：读取 analysis_swc → 生成锚点并裁切连通 ROI（4.3.5）\n"
            "            → 区分 TRUE_TERMINAL / CUT_PORT（4.3.5.1）\n"
            "            → 计算半径与结构特征（4.3.5.2）\n"
            "            → 聚类、选择真实代表 ROI（4.3.5.3）→ 保存与验收（4.3.7）",
            flush=True,
        )
        try:
            sampling_run = run_sampling_from_rodent_run(
                run.run_root,
                settings.sampling,
                verbose=settings.verbose,
            )
            sampling_source_config = _copy_source_configuration(
                settings.source_path,
                sampling_run.run_root / "config" / "source_swc_roi_generate.yaml",
            )
        except Exception as exc:
            print(f"[错误] [文档 4.3.5 / 4.3.7] ROI 采样或源配置保存失败：{exc}", file=sys.stderr, flush=True)
            return 1
        _print_stage("4.3.7", "ROI 采样：输出文件与验收结果", "结果")
        _print_path("ROI 采样结果目录（Sampling run directory）", sampling_run.run_root, "保存基于上述全局结果生成的候选 ROI、特征、聚类和代表样本。")
        _print_path("采样源配置副本（Sampling source configuration）", sampling_source_config, "与全局结果中的源 YAML 副本内容相同，用于追溯本次完整参数。")
        _print_path("采样有效配置", sampling_run.run_root / "config" / "sampling_config.yaml", "保存解析后的 ROI 采样参数。")
        _print_path("采样汇总（Sampling summary）", sampling_run.summary_path, "JSON 文件：候选数、拒绝原因、聚类、代表数量、完整性检查和运行时间。")
        _print_path("最终代表 ROI 清单", sampling_run.run_root / "manifests" / "selected_rois.csv", "列出最终入选的真实代表 ROI；查找后续使用的代表样本时从此处开始。")
        _print_path("ROI 几何库", sampling_run.run_root / "roi_library", "每个 NPZ 保存一个有效候选的几何、半径和局部—全局映射；代表身份以清单为准。")
        _print_path("采样详细日志", sampling_run.run_root / "logs" / "sampling.log", "查看 ROI 候选生成、聚类、代表选择与验收的过程记录。")
        _print_status("ROI 采样验收状态（Sampling status）", sampling_run.status)
        if sampling_run.status == "FAIL":
            print("[结束] ROI 采样验收未通过，退出码 2；请查看上述采样汇总。", flush=True)
            return 2
    else:
        _print_stage("4.3.5", "真实连通局部区域的代表性采样", "跳过")
        print(
            "  原因：" + (
                "全局处理失败。" if run.status == "failed" else
                f"pipeline.stage={config.stage}，当前阶段不执行 ROI 采样。" if not graph_stage else
                "sampling.enabled=false。"
            ),
            flush=True,
        )

    sampling_gui_preview = (
        sampling_run.run_root / "figures" / "interactive_sampling_layer_preview.png"
        if sampling_run is not None
        else None
    )
    if (
        run.status != "failed"
        and graph_stage
        and config.figure2a_enabled
        and not settings.show_gui
        and sampling_run is not None
    ):
        _print_stage("4.3.6", "生成双视窗屏幕外预览（当前不打开交互窗口）")
        try:
            show_saved_run(
                run.run_root,
                sample_id=config.sample_id,
                max_arrows=config.max_direction_arrows,
                volume_opacity=config.figure2a_volume_opacity,
                window_size=config.figure2a_window_size,
                sampling_run_root=sampling_run.run_root,
                screenshot_path=sampling_gui_preview,
                show=False,
            )
            _print_stage("4.3.6～4.3.7", "双视窗预览已保存", "完成")
            _print_path("双视窗预览图（Sampling GUI preview）", sampling_gui_preview, "PNG 图片：左侧展示整体血管与 ROI 位置，右侧展示当前 ROI 的局部结构。")
        except Exception as exc:
            print(f"[错误] [文档 4.3.6] 双视窗预览图保存失败：{exc}", file=sys.stderr, flush=True)
            return 1

    if (
        run.status != "failed"
        and graph_stage
        and config.figure2a_enabled
        and settings.show_gui
    ):
        _print_stage("4.3.6", "创建并打开双视窗交互窗口")
        print("  窗口打开后程序会等待交互操作；关闭窗口后程序才会结束。", flush=True)
        if sampling_run is not None:
            print(
                "  操作：R / S = 已选代表 ROI；A = 全部候选 ROI；C = 下一个聚类。\n"
                "  左键点击左侧 ROI 方框，右侧显示对应局部结构；两个视窗实时水平旋转。",
                flush=True,
            )
        try:
            show_saved_run(
                run.run_root,
                sample_id=config.sample_id,
                max_arrows=config.max_direction_arrows,
                volume_opacity=config.figure2a_volume_opacity,
                window_size=config.figure2a_window_size,
                sampling_run_root=sampling_run.run_root if sampling_run else None,
                screenshot_path=sampling_gui_preview,
            )
        except Exception as exc:
            print(f"[错误] [文档 4.3.6] 双视窗交互显示失败：{exc}", file=sys.stderr, flush=True)
            return 1
        _print_stage("4.3.6", "交互窗口已关闭", "完成")
        if sampling_gui_preview is not None:
            _print_path("双视窗预览图（Sampling GUI preview）", sampling_gui_preview, "本次交互显示保存的 PNG 截图。")
    elif run.status == "failed" or not graph_stage or not config.figure2a_enabled or sampling_run is None:
        _print_stage("4.3.6", "双视窗显示 / 预览", "跳过")
        print(
            "  原因：" + (
                "全局处理失败。" if run.status == "failed" else
                f"pipeline.stage={config.stage}，当前阶段不执行双视窗显示。" if not graph_stage else
                "双视窗可视化未启用（figure2a_enabled=false）。" if not config.figure2a_enabled else
                "当前为屏幕外模式且本次未执行 ROI 采样，无采样交互层预览。"
            ),
            flush=True,
        )
    print(
        "\n[结束] 全局处理验收未通过，退出码 2；请查看上述全局验收报告。"
        if run.status == "failed" else
        "\n[结束] 本次请求的处理已完成，退出码 0；各阶段的验收状态及文件用途见上方结果区。",
        flush=True,
    )
    return 2 if run.status == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
