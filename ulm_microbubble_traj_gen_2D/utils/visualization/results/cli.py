"""Command-line interface for visualizing field-based microbubble flow results."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .plotting import (
    render_animation,
    render_flow_field_map,
    render_hole_map,
    render_local_shear_stress_map,
    render_narrow_cell_map,
    render_snapshot,
    render_wall_shear_stress_map,
)
from .result_loader import DEFAULT_RESULTS_ROOT, find_latest_result, load_result


def parse_args() -> argparse.Namespace:
    """Parse visualization options."""

    parser = argparse.ArgumentParser(
        description="Visualize microbubble movement through the field-based vascular result."
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=None,
        help="A results/<timestamp> folder. If omitted, the newest complete result is used.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List complete result folders and exit.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Animation output path. Defaults to <result-dir>/microbubble_flow.gif.",
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=None,
        help="Optional snapshot output path. Defaults to <result-dir>/microbubble_flow_snapshot.png when set without a value is not supported.",
    )
    parser.add_argument(
        "--snapshot-frame",
        type=int,
        default=-1,
        help="Frame index used for the snapshot. Negative values count from the end.",
    )
    parser.add_argument(
        "--max-bubbles",
        type=int,
        default=0,
        help="Maximum bubbles drawn; 0 (default) draws the full configured per-frame population.",
    )
    parser.add_argument(
        "--tail-length",
        type=int,
        default=20,
        help="Number of previous frames shown as a short bubble trail.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=3,
        help="Use every Nth frame in the animation to keep files compact.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=20,
        help="Animation frames per second.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=130,
        help="Animation DPI.",
    )
    parser.add_argument(
        "--color-mode",
        choices=(
            "speed",
            "wall_shear",
            "local_shear",
            "vessel_id",
            "diameter",
            "clearance",
            "active",
        ),
        default="speed",
        help="Scalar used to color bubble markers.",
    )
    parser.add_argument(
        "--glow-scale",
        type=_glow_scale,
        default=3.0,
        help="Halo scale controlling both physical ratio and overview visibility; default 3.0.",
    )
    parser.add_argument(
        "--no-glow",
        action="store_true",
        help="Disable the visual halo while keeping the physical-size bubble disks.",
    )
    parser.add_argument(
        "--no-animation",
        action="store_true",
        help="Only save a snapshot, not an animation.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm progress bars while rendering animation frames.",
    )
    static_group = parser.add_mutually_exclusive_group()
    static_group.add_argument(
        "--skip-static-maps",
        action="store_true",
        help="Skip the five static diagnostic maps.",
    )
    static_group.add_argument(
        "--force-static-maps",
        action="store_true",
        help="Re-render static maps even when current cached PNG files are reusable.",
    )
    return parser.parse_args()


def _glow_scale(value: str) -> float:
    """Parse a finite halo-to-bubble diameter ratio accepted by the renderer."""

    parsed = float(value)
    if not np.isfinite(parsed) or parsed < 1.0:
        raise argparse.ArgumentTypeError("--glow-scale must be a finite number greater than or equal to 1.0")
    return parsed


def main() -> None:
    """Load a result folder and render requested visualization files."""

    args = parse_args()
    if args.list:
        _list_results()
        return

    data = load_result(args.result_dir)
    display_limit = int(args.max_bubbles)
    display_description = "all active bubbles" if display_limit <= 0 else f"at most {display_limit} bubbles"
    print(
        "Trajectory population: "
        f"frames={data.positions_um.shape[0]}, "
        f"peak active/frame={data.peak_active_bubbles}, "
        f"active/frame={int(data.active_count_per_frame.min())}..{int(data.active_count_per_frame.max())}, "
        f"display={display_description}"
    )
    if args.color_mode == "speed":
        print(
            "Speed color semantics: realized same-permanent-ID center displacement divided by dt; "
            "wall-clamped stationary bubbles therefore display zero translational speed."
        )
    elif args.color_mode == "wall_shear":
        print(
            "Wall-shear color semantics: the lumen background and bubble markers "
            "share one 0..P99.5 WSS scale; values above P99.5 use the top color."
        )
    finite_gap = data.bubble_wall_gap_um[np.isfinite(data.bubble_wall_gap_um)]
    if finite_gap.size:
        penetrations = int(
            np.count_nonzero(finite_gap < float(data.wall_gap_invalid_below_um))
        )
        invalid_label = (
            "invalid accepted negative gaps"
            if data.strict_nonnegative_wall_gap
            else "over-tolerance penetrations (legacy convention)"
        )
        tolerance_label = (
            "contact solve tolerance"
            if data.strict_nonnegative_wall_gap
            else "legacy penetration tolerance"
        )
        print(
            f"Wall clearance: minimum={float(finite_gap.min()):.6g} um, "
            f"contact observations={int(np.count_nonzero(data.bubble_wall_contact))}, "
            f"{invalid_label}={penetrations}, "
            f"{tolerance_label}={float(data.contact_geometry_tolerance_um):.6g} um"
        )
    else:
        print("Wall clearance: unavailable in this legacy trajectory; regenerate particles for contact validation.")
    numerical_wall_locks = int(np.count_nonzero(data.bubble_numerical_wall_lock))
    if data.numerical_wall_lock_source != "unavailable":
        print(
            "Numerical wall-lock diagnostic: "
            f"records={numerical_wall_locks}, "
            f"source={data.numerical_wall_lock_source}"
        )
    if data.molecular_target.visible:
        density_um2 = data.molecular_target.target_density_molecules_per_m2 * 1.0e-12
        print(
            "Molecular target overlay: "
            f"wall sites={data.molecular_target.target_wall_site_count}, "
            f"density={density_um2:.6g} molecules/um^2, "
            f"source={data.molecular_target.source_path}"
        )
    else:
        print("Molecular target overlay: unavailable or empty for this result directory.")
    animation_path = args.output or (data.result_dir / "microbubble_flow.gif")
    holes_path = data.result_dir / "lumen_holes.png"
    narrow_cells_path = data.result_dir / "narrow_lumen_cells.png"
    initial_flow_path = data.result_dir / "initial_flow_field.png"
    final_flow_path = data.result_dir / "final_flow_field.png"
    wall_shear_path = data.result_dir / "wall_shear_stress_distribution.png"
    local_shear_path = data.result_dir / "local_shear_stress_distribution.png"
    snapshot_path = args.snapshot
    if snapshot_path is None and args.no_animation:
        snapshot_path = data.result_dir / "microbubble_flow_snapshot.png"

    if not args.skip_static_maps:
        dependencies = _static_map_dependencies(data)
        static_outputs = (
            (
                "Holes map",
                holes_path,
                max(args.dpi, 240),
                lambda path, dpi: render_hole_map(data, path, dpi=dpi),
            ),
            (
                "Narrow-cells map",
                narrow_cells_path,
                max(args.dpi, 240),
                lambda path, dpi: render_narrow_cell_map(data, path, dpi=dpi),
            ),
            (
                "Initial flow-field map",
                initial_flow_path,
                max(args.dpi, 220),
                lambda path, dpi: render_flow_field_map(data, path, stage="initial", dpi=dpi),
            ),
            (
                "Final flow-field map",
                final_flow_path,
                max(args.dpi, 220),
                lambda path, dpi: render_flow_field_map(data, path, stage="final", dpi=dpi),
            ),
            (
                "Wall-shear-stress distribution map",
                wall_shear_path,
                max(args.dpi, 220),
                lambda path, dpi: render_wall_shear_stress_map(data, path, dpi=dpi),
            ),
            (
                "Whole-lumen local-shear-stress distribution map",
                local_shear_path,
                max(args.dpi, 220),
                lambda path, dpi: render_local_shear_stress_map(data, path, dpi=dpi),
            ),
        )
        for label, path, output_dpi, renderer in static_outputs:
            reusable = not args.force_static_maps and _is_reusable_png(
                path,
                dpi=output_dpi,
                dependencies=dependencies,
            )
            if not reusable:
                renderer(path, output_dpi)
            status = "reused" if reusable else "saved"
            print(f"{label} {status}: {path}")
    else:
        print("Static diagnostic maps skipped.")

    if snapshot_path is not None:
        render_snapshot(
            data,
            snapshot_path,
            frame_index=args.snapshot_frame,
            max_bubbles=args.max_bubbles,
            tail_length=args.tail_length,
            color_mode=args.color_mode,
            dpi=max(args.dpi, 150),
            glow_enabled=not args.no_glow,
            glow_scale=args.glow_scale,
        )
        print(f"Snapshot saved: {snapshot_path}")

    if not args.no_animation:
        render_animation(
            data,
            animation_path,
            max_bubbles=args.max_bubbles,
            tail_length=args.tail_length,
            stride=args.stride,
            fps=args.fps,
            color_mode=args.color_mode,
            dpi=args.dpi,
            show_progress=not args.no_progress,
            glow_enabled=not args.no_glow,
            glow_scale=args.glow_scale,
        )
        print(f"Animation saved: {animation_path}")

    print(f"Visualized result: {data.result_dir}")


def _static_map_dependencies(data) -> tuple[Path, ...]:
    """Return source files whose changes invalidate cached static maps."""

    dependencies = [
        data.result_dir / "velocity_and_wall_shear_field.npz",
        Path(__file__).resolve(),
        Path(__file__).with_name("plotting.py").resolve(),
        Path(__file__).with_name("result_loader.py").resolve(),
    ]
    if data.molecular_target.source_path is not None:
        dependencies.append(Path(data.molecular_target.source_path))
    return tuple(path for path in dependencies if path.is_file())


def _is_reusable_png(output_path: Path, *, dpi: int, dependencies: tuple[Path, ...]) -> bool:
    """Return whether a PNG is current and has the requested resolution metadata."""

    path = Path(output_path)
    if not path.is_file():
        return False
    newest_source_time = max((source.stat().st_mtime_ns for source in dependencies), default=0)
    if path.stat().st_mtime_ns < newest_source_time:
        return False
    try:
        from PIL import Image

        with Image.open(path) as image:
            stored_dpi = image.info.get("dpi")
    except (OSError, ValueError):
        return False
    if not stored_dpi or len(stored_dpi) < 2:
        return False
    return all(abs(float(value) - float(dpi)) <= 1.0 for value in stored_dpi[:2])


def _list_results() -> None:
    """Print result folders that contain both field and trajectory files."""

    root = DEFAULT_RESULTS_ROOT
    if not root.exists():
        print(f"No results directory found: {root}")
        return
    for item in sorted(root.iterdir()):
        if not item.is_dir():
            continue
        if (item / "velocity_and_wall_shear_field.npz").exists() and (item / "microbubble_field_trajectories.npz").exists():
            print(item)
    latest = find_latest_result(root)
    print(f"Latest: {latest}")
