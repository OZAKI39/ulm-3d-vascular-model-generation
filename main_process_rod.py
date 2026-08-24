"""Figure 2(a)-style mouse-brain vascular visualization entry point.

Running this file without arguments processes one representative 192^3 block,
writes acceptance evidence under ``outputs/rodent_vasculature``, and opens a
native PyVista window. The window supports rotate, pan, zoom, reset-view, and
screenshot controls. Orange arrows encode the SWC parent-to-current relation;
they do not claim measured blood-flow direction.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils.rodent_vasculature import RodentVasculatureConfig, run_rodent_vasculature_pipeline
from utils.rodent_vasculature.interactive import show_saved_run
from utils.sampling import SamplingConfig
from utils.sampling.pipeline import run_sampling_from_rodent_run


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "vessel_model"
    / "T - A high-resolution dataset of mouse brain vasculature"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a mouse-brain vascular SWC, select one unchanged analysis component, "
            "render its directed node/branch relationships with optional TIFF context, and "
            "write reproducible acceptance evidence."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--stage",
        choices=("inventory", "preprocess", "hierarchical-graph", "all"),
        default="all",
        help=(
            "inventory catalogs SWC plus optional TIFF context; preprocess validates and "
            "selects the SWC analysis component; "
            "hierarchical-graph reads --source-run; all performs the complete workflow "
            "(default: all)."
        ),
    )
    parser.add_argument(
        "--source-run",
        type=Path,
        default=None,
        help="Completed preprocess/all run required by --stage hierarchical-graph.",
    )
    parser.add_argument(
        "--cohort",
        choices=("raw-analysis", "raw-total", "train", "all"),
        default="raw-analysis",
        help="Dataset subset (default: raw-analysis, the 108 analysis blocks).",
    )
    parser.add_argument(
        "--sample-id",
        default="fMOST_0_5_6_0_0_6_0001_02_01",
        help="Exact catalog sample_id or source stem; eligibility requires only a non-empty SWC.",
    )
    parser.add_argument("--parent-group-id", default=None)
    parser.add_argument("--split", choices=("train", "val", "test"), default=None)
    parser.add_argument(
        "--max-samples",
        type=int,
        default=1,
        help="Number of blocks to process (default: 1 for responsive interactive rendering).",
    )
    parser.add_argument(
        "--spacing-xyz",
        type=float,
        nargs=3,
        default=(1.0, 1.0, 2.0),
        metavar=("SX", "SY", "SZ"),
        help="Voxel spacing in micrometres (default: 1 1 2, as reported by the paper).",
    )
    parser.add_argument(
        "--expected-shape-zyx",
        type=int,
        nargs=3,
        default=(192, 192, 192),
        metavar=("Z", "Y", "X"),
    )
    parser.add_argument(
        "--no-shape-check",
        action="store_true",
        help="Allow TIFF volumes other than the paper's 192 x 192 x 192 block size.",
    )
    parser.add_argument("--smooth-centerlines", action="store_true")
    parser.add_argument("--smoothing-window", type=int, default=5)
    parser.add_argument("--resample-step", type=float, default=1.0, metavar="UM")
    parser.add_argument("--strict-radii", action="store_true")
    parser.add_argument(
        "--analysis-component-id",
        type=int,
        default=None,
        metavar="ID",
        help=(
            "Use this original SWC connected-component ID as analysis_swc. By default, "
            "the component with the greatest physical centerline length is selected."
        ),
    )
    parser.add_argument("--max-direction-arrows", type=int, default=600)
    parser.add_argument(
        "--volume-opacity",
        type=float,
        default=0.32,
        help="Maximum opacity of the raw reference image volume, in (0, 1] (default: 0.32).",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        nargs=2,
        default=(1800, 900),
        metavar=("WIDTH", "HEIGHT"),
    )
    parser.add_argument(
        "--no-sampling",
        action="store_true",
        help="Skip representative connected-ROI sampling after graph construction.",
    )
    parser.add_argument("--sampling-seed", type=int, default=42)
    parser.add_argument(
        "--sampling-anchor-mode",
        choices=("random", "farthest_point", "poisson_disk"),
        default="farthest_point",
    )
    parser.add_argument("--sampling-min-anchor-distance", type=float, default=45.0, metavar="UM")
    parser.add_argument("--sampling-max-candidates", type=int, default=80, metavar="COUNT")
    parser.add_argument(
        "--sampling-roi-size",
        type=float,
        nargs=3,
        default=(80.0, 80.0, 120.0),
        metavar=("SX", "SY", "SZ"),
        help="Physical connected-ROI box size in micrometres.",
    )
    parser.add_argument("--sampling-min-branches", type=int, default=2, metavar="COUNT")
    parser.add_argument(
        "--sampling-feature-mode",
        choices=("radius_only", "radius_plus_structure", "extended_morphology"),
        default="radius_plus_structure",
    )
    parser.add_argument("--sampling-clusters", type=int, default=5, metavar="K")
    parser.add_argument(
        "--sampling-selection-mode",
        choices=("distribution_preserving", "coverage_balanced"),
        default="coverage_balanced",
    )
    parser.add_argument("--sampling-target-count", type=int, default=10, metavar="COUNT")
    parser.add_argument(
        "--sampling-representatives-per-cluster",
        type=int,
        default=2,
        metavar="COUNT",
    )
    parser.add_argument("--sampling-max-overlap", type=float, default=0.25, metavar="FRACTION")
    parser.add_argument("--sampling-min-selected-distance", type=float, default=20.0, metavar="UM")
    parser.add_argument(
        "--no-sampling-comparison",
        action="store_true",
        help="Skip the radius-only versus radius-plus-structure comparison experiment.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open the GUI; still write the same off-screen acceptance preview.",
    )
    parser.add_argument(
        "--no-figure2a",
        action="store_true",
        help="Skip the PyVista Figure 2(a)-style volume preview and GUI.",
    )
    parser.add_argument("--no-visualizations", action="store_true")
    parser.add_argument("--no-graphml", action="store_true")
    parser.add_argument("--no-vtp", action="store_true")
    parser.add_argument("--no-npz", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.stage == "hierarchical-graph" and args.source_run is None:
        parser.error("--source-run is required when --stage hierarchical-graph")

    visualizations_enabled = not args.no_visualizations
    figure2a_enabled = visualizations_enabled and not args.no_figure2a
    graph_stage = args.stage in {"all", "hierarchical-graph"}
    if graph_stage and figure2a_enabled and not args.no_show and args.no_npz:
        parser.error("--no-npz cannot be combined with the interactive window; add --no-show")
    config = RodentVasculatureConfig(
        input_dir=args.input_dir,
        output_root=args.output_dir,
        stage=args.stage,
        source_run=args.source_run,
        cohort=args.cohort,
        sample_id=args.sample_id,
        parent_group_id=args.parent_group_id,
        split=args.split,
        max_samples=args.max_samples,
        spacing_xyz_um=tuple(args.spacing_xyz),
        expected_shape_zyx=None if args.no_shape_check else tuple(args.expected_shape_zyx),
        smoothing_enabled=args.smooth_centerlines,
        smoothing_window_points=args.smoothing_window,
        resample_step_um=args.resample_step,
        strict_nonpositive_radius=args.strict_radii,
        analysis_component_id=args.analysis_component_id,
        visualizations_enabled=visualizations_enabled,
        max_visualization_samples=max(1, args.max_samples),
        max_direction_arrows=args.max_direction_arrows,
        figure2a_enabled=figure2a_enabled,
        figure2a_volume_opacity=args.volume_opacity,
        figure2a_window_size=tuple(args.window_size),
        save_graphml=not args.no_graphml,
        save_vtp=not args.no_vtp,
        save_npz=not args.no_npz,
    )
    try:
        run = run_rodent_vasculature_pipeline(config, verbose=args.verbose)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Run directory: {run.run_root}")
    print(f"Acceptance report: {run.html_report}")
    print(f"Acceptance status: {run.acceptance.overall_status}")

    sampling_run = None
    if run.status != "failed" and graph_stage and not args.no_sampling:
        sampling_config = SamplingConfig(
            output_root=args.output_dir,
            seed=args.sampling_seed,
            anchor_mode=args.sampling_anchor_mode,
            min_anchor_distance_um=args.sampling_min_anchor_distance,
            max_candidate_anchors=args.sampling_max_candidates,
            roi_size_um=tuple(args.sampling_roi_size),
            min_branch_count=args.sampling_min_branches,
            feature_mode=args.sampling_feature_mode,
            n_clusters=args.sampling_clusters,
            selection_mode=args.sampling_selection_mode,
            target_selected_count=args.sampling_target_count,
            representatives_per_cluster=args.sampling_representatives_per_cluster,
            max_selected_overlap=args.sampling_max_overlap,
            min_representative_distance_um=args.sampling_min_selected_distance,
            compare_feature_modes=not args.no_sampling_comparison,
        )
        try:
            sampling_run = run_sampling_from_rodent_run(
                run.run_root,
                sampling_config,
                verbose=args.verbose,
            )
        except Exception as exc:
            print(f"ERROR: connected-ROI sampling failed: {exc}", file=sys.stderr)
            return 1
        print(f"Sampling run directory: {sampling_run.run_root}")
        print(f"Sampling summary: {sampling_run.summary_path}")
        print(f"Sampling status: {sampling_run.status}")
        if sampling_run.status == "FAIL":
            return 2

    sampling_gui_preview = (
        sampling_run.run_root / "figures" / "interactive_sampling_layer_preview.png"
        if sampling_run is not None
        else None
    )
    if (
        run.status != "failed"
        and graph_stage
        and figure2a_enabled
        and args.no_show
        and sampling_run is not None
    ):
        try:
            show_saved_run(
                run.run_root,
                sample_id=args.sample_id,
                max_arrows=args.max_direction_arrows,
                volume_opacity=args.volume_opacity,
                window_size=tuple(args.window_size),
                sampling_run_root=sampling_run.run_root,
                screenshot_path=sampling_gui_preview,
                show=False,
            )
            print(f"Sampling GUI preview: {sampling_gui_preview}")
        except Exception as exc:
            print(f"ERROR: sampling GUI preview could not be written: {exc}", file=sys.stderr)
            return 1

    if run.status != "failed" and graph_stage and figure2a_enabled and not args.no_show:
        print("Opening interactive Figure 2(a)-style window; close it to finish the program.")
        if sampling_run is not None:
            print(
                "GUI controls: R/S=selected representative ROIs, A=all candidate ROIs, "
                "C=next ROI cluster; left-click a box to inspect it."
            )
        try:
            show_saved_run(
                run.run_root,
                sample_id=args.sample_id,
                max_arrows=args.max_direction_arrows,
                volume_opacity=args.volume_opacity,
                window_size=tuple(args.window_size),
                sampling_run_root=sampling_run.run_root if sampling_run else None,
                screenshot_path=sampling_gui_preview,
            )
        except Exception as exc:
            print(f"ERROR: interactive window could not be opened: {exc}", file=sys.stderr)
            return 1
    return 2 if run.status == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
