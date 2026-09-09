"""CLI for geometry import only. --help performs no import or computation."""

import argparse
from pathlib import Path
import sys


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="导入已验收血管数据包并核查；匹配结果直接复用，仅 CPU 几何迁移。")
    parser.add_argument("--config", type=Path, default=Path("py_scripts/vessel_geometry_import.yaml"), help="导入 YAML；相对路径以目标工程根目录为基准")
    args = parser.parse_args(argv)
    from .vessel_geometry.export import run_import
    from .vessel_geometry.model import GeometryError
    try:
        result = run_import(args.config, progress=lambda text: print(text, flush=True))
    except (GeometryError, OSError, ValueError) as exc:
        print(f"[导入未通过] {exc}", file=sys.stderr)
        if isinstance(exc, GeometryError) and exc.details.get("package_path"):
            print(f"失败诊断目录：{exc.details['package_path']}", file=sys.stderr)
        return 2
    print(f"正式数据包：{result.package_path}\n迁移检查：PASS\n人工核查：PENDING\nSDF、粒子、数值边界、真实血流：NOT_STARTED", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
