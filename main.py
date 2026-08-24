"""Command-line entry point for the staged vascular preprocessing workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils.config import (
    HierarchicalGraphConfig,
    MeshCleanupConfig,
    PipelineConfig,
    VisualizationConfig,
    VoxelizationConfig,
)
from utils.hierarchical_graph_pipeline import run_hierarchical_graph_pipeline
from utils.pipeline import run_pipeline


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = PROJECT_ROOT / "unprocess_stl" / "mv03_y.nii.gz_Segment_1.stl"
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Step 1/2 vascular preprocessing and Step 3 hierarchical vascular "
            "representation for a dense vascular STL."
        )
    )
    parser.add_argument(
        "--stage",
        choices=("preprocess", "hierarchical-graph"),
        default="preprocess",
        help=(
            "preprocess runs the existing Step 1/2 pipeline; hierarchical-graph "
            "reads a completed Step 2 run and executes Step 3."
        ),
    )
    parser.add_argument(
        "--source-run",
        type=Path,
        default=None,
        help="Completed Step 2 run directory required by --stage hierarchical-graph.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input STL path.")
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT, help="Root directory for timestamped outputs."
    )
    parser.add_argument(
        "--voxel-size",
        type=float,
        default=2.0,
        metavar="UM",
        help="Isotropic voxel size in micrometers (default: 2.0).",
    )
    parser.add_argument("--padding-voxels", type=int, default=3, help="Empty grid cells around the STL.")
    parser.add_argument(
        "--max-voxel-count",
        type=int,
        default=200_000_000,
        help="Safety limit for total voxel-grid size.",
    )
    parser.add_argument("--min-component-faces", type=int, default=100)
    parser.add_argument("--min-component-area", type=float, default=50.0, metavar="UM2")
    parser.add_argument("--min-component-diagonal", type=float, default=10.0, metavar="UM")
    parser.add_argument(
        "--component-policy",
        choices=("main_network_only", "conservative"),
        default="main_network_only",
        help=(
            "main_network_only keeps one connected STL network (default); "
            "conservative retains large disconnected components for comparison."
        ),
    )
    parser.add_argument(
        "--main-component-id",
        type=int,
        default=None,
        help="Explicit connected-component ID to use as the main network.",
    )
    parser.add_argument(
        "--smoothing-iterations",
        type=int,
        default=0,
        help="Feature-preserving smoothing iterations; 0 preserves the original surface.",
    )
    parser.add_argument("--no-repair", action="store_true", help="Disable repair of retained invalid components.")
    parser.add_argument(
        "--keep-all-voxel-components",
        action="store_true",
        help="Diagnostic mode: do not keep only the largest connected voxel network.",
    )
    parser.add_argument(
        "--voxel-island-warning-fraction",
        type=float,
        default=0.01,
        help="Warn when voxel connectivity filtering removes this fraction (default: 0.01).",
    )
    parser.add_argument(
        "--voxel-island-fail-fraction",
        type=float,
        default=0.10,
        help="Fail acceptance when voxel filtering removes this fraction (default: 0.10).",
    )
    parser.add_argument("--no-visualizations", action="store_true", help="Skip PNG output (not recommended).")
    parser.add_argument(
        "--graph-connectivity",
        type=int,
        choices=(6, 18, 26),
        default=26,
        help="Voxel-neighbor rule used by Step 3 (default: 26).",
    )
    parser.add_argument(
        "--graph-smoothing-window",
        type=int,
        default=5,
        help="Odd smoothing window for derived branch measurements (default: 5).",
    )
    parser.add_argument(
        "--no-graph-smoothing",
        action="store_true",
        help="Keep derived Step 3 measurements unsmoothed; raw geometry is always retained.",
    )
    parser.add_argument(
        "--graph-resample-step",
        type=float,
        default=None,
        metavar="UM",
        help="Physical spacing for the derived smooth sequence; defaults to voxel spacing.",
    )
    parser.add_argument(
        "--junction-direction-distance",
        type=float,
        default=10.0,
        metavar="UM",
        help="Distance from a junction used to estimate outward branch direction.",
    )
    parser.add_argument(
        "--short-branch-warning",
        type=float,
        default=6.0,
        metavar="UM",
        help="Warn, but do not delete, branches shorter than this length.",
    )
    parser.add_argument(
        "--large-junction-warning-voxels",
        type=int,
        default=20,
        help="Warn when one merged junction region contains at least this many voxels.",
    )
    parser.add_argument(
        "--high-degree-warning",
        type=int,
        default=6,
        help="Warn when a graph node has at least this degree.",
    )
    parser.add_argument(
        "--smoothing-deviation-warning",
        type=float,
        default=2.0,
        metavar="UM",
        help="Warn when a derived smooth centerline moves farther than this from resampled raw geometry.",
    )
    parser.add_argument("--verbose", action="store_true", help="Show debug-level console logs.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.stage == "hierarchical-graph":
        if args.source_run is None:
            parser.error("--source-run is required when --stage hierarchical-graph")
        graph_config = HierarchicalGraphConfig(
            neighbor_connectivity=args.graph_connectivity,
            smoothing_enabled=not args.no_graph_smoothing,
            smoothing_window_points=args.graph_smoothing_window,
            resample_step_um=args.graph_resample_step,
            junction_direction_distance_um=args.junction_direction_distance,
            short_branch_warning_um=args.short_branch_warning,
            large_junction_warning_voxels=args.large_junction_warning_voxels,
            high_degree_warning=args.high_degree_warning,
            smoothing_deviation_warning_um=args.smoothing_deviation_warning,
        )
        try:
            graph_result = run_hierarchical_graph_pipeline(
                args.source_run, graph_config, verbose=args.verbose
            )
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1
        print(f"Run directory: {graph_result.run_root}")
        print(f"Acceptance report: {graph_result.html_report}")
        print(f"Acceptance status: {graph_result.acceptance.overall_status}")
        return 2 if graph_result.status == "failed" else 0

    config = PipelineConfig(
        input_stl=args.input,
        output_root=args.output_dir,
        mesh=MeshCleanupConfig(
            component_policy=args.component_policy,
            main_component_id=args.main_component_id,
            min_component_faces=args.min_component_faces,
            min_component_area_um2=args.min_component_area,
            min_component_diagonal_um=args.min_component_diagonal,
            repair_non_manifold=not args.no_repair,
            smoothing_iterations=args.smoothing_iterations,
        ),
        voxel=VoxelizationConfig(
            voxel_size_um=args.voxel_size,
            padding_voxels=args.padding_voxels,
            max_voxel_count=args.max_voxel_count,
            keep_largest_connected_component=not args.keep_all_voxel_components,
            island_warning_fraction=args.voxel_island_warning_fraction,
            island_fail_fraction=args.voxel_island_fail_fraction,
        ),
        visualization=VisualizationConfig(enabled=not args.no_visualizations),
    )
    try:
        result = run_pipeline(config, verbose=args.verbose)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Run directory: {result.run_root}")
    print(f"Acceptance report: {result.html_report}")
    print(f"Acceptance status: {result.acceptance.overall_status}")
    return 2 if result.status == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
