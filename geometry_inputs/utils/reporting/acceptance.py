"""Automatic, explainable acceptance checks for generated artifacts."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import numpy as np

if TYPE_CHECKING:
    from ..config import VoxelizationConfig
    from ..mesh.cleanup import MeshCleanupResult
    from ..voxel.skeleton import SkeletonResult
    from ..voxel.voxelize import VoxelizationResult


@dataclass(frozen=True, slots=True)
class AcceptanceCheck:
    name: str
    status: str
    message: str

    def to_dict(self) -> dict[str, str]:
        return asdict(self)


@dataclass(slots=True)
class AcceptanceResult:
    overall_status: str
    checks: list[AcceptanceCheck]

    def report(self) -> dict[str, Any]:
        return {
            "overall_status": self.overall_status,
            "counts": {
                "PASS": sum(item.status == "PASS" for item in self.checks),
                "WARNING": sum(item.status == "WARNING" for item in self.checks),
                "FAIL": sum(item.status == "FAIL" for item in self.checks),
            },
            "checks": [item.to_dict() for item in self.checks],
        }


def _relative_change(before: float, after: float) -> float:
    return abs(after - before) / max(abs(before), np.finfo(float).eps)


def _overall(checks: Iterable[AcceptanceCheck]) -> str:
    statuses = {item.status for item in checks}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARNING" in statuses:
        return "WARNING"
    return "PASS"


def evaluate_acceptance(
    cleanup: MeshCleanupResult,
    voxel: VoxelizationResult,
    skeleton: SkeletonResult,
    required_files: list[Path],
    voxel_config: VoxelizationConfig,
) -> AcceptanceResult:
    import nibabel as nib

    after = cleanup.cleaned_quality
    checks: list[AcceptanceCheck] = []

    for name, value in (
        ("No degenerate triangles", after.degenerate_triangle_count),
        ("No duplicate triangles", after.duplicate_triangle_count),
        ("Closed surface (no boundary edges)", after.boundary_edge_count),
        ("No non-manifold edges", after.non_manifold_edge_count),
    ):
        checks.append(
            AcceptanceCheck(
                name=name,
                status="PASS" if value == 0 else "FAIL",
                message=f"Detected count: {value}",
            )
        )

    if cleanup.repair_failures:
        checks.append(
            AcceptanceCheck(
                "Requested component repairs",
                "FAIL",
                f"Repair failed for component IDs: {cleanup.repair_failures}",
            )
        )
    else:
        checks.append(AcceptanceCheck("Requested component repairs", "PASS", "No repair failed."))

    summary = cleanup.summary()
    expected_mesh_components = (
        1
        if cleanup.main_network_selection.policy == "main_network_only"
        else summary["kept_component_count"]
    )
    checks.append(
        AcceptanceCheck(
            "Final STL has the expected connected network count",
            "PASS"
            if after.connected_component_count == expected_mesh_components
            else "FAIL",
            f"Expected {expected_mesh_components}, detected {after.connected_component_count}.",
        )
    )

    selection = cleanup.main_network_selection
    checks.append(
        AcceptanceCheck(
            "Main-network ranking agrees",
            "PASS" if selection.ranking_agrees else "WARNING",
            "Surface area, triangle count, and spatial span select the same component."
            if selection.ranking_agrees
            else (
                "The three size rankings disagree; visually confirm selected component "
                f"{selection.selected_component_id}."
            ),
        )
    )

    main_record = next(
        item for item in cleanup.components if item.component_type == "main_network"
    )
    if selection.policy == "main_network_only":
        main_area_delta = _relative_change(
            main_record.surface_area_um2, after.surface_area_um2
        )
        main_surface_status = (
            "PASS"
            if main_area_delta <= 0.02
            else ("WARNING" if main_area_delta <= 0.10 else "FAIL")
        )
        main_surface_message = f"Selected-main surface-area change: {main_area_delta:.3%}."
    else:
        main_surface_status = "PASS"
        main_surface_message = "Conservative comparison mode keeps additional components."
    checks.append(
        AcceptanceCheck(
            "Selected main surface is preserved",
            main_surface_status,
            main_surface_message,
        )
    )
    checks.append(
        AcceptanceCheck(
            "Disconnected STL parts were classified",
            "PASS",
            f"Removed {summary['small_fragment_count']} small fragments and "
            f"{summary['island_network_count']} island networks; total removed surface "
            f"fraction={summary['removed_surface_area_fraction']:.3%}.",
        )
    )

    checks.append(
        AcceptanceCheck(
            "Voxel mask is non-empty",
            "PASS" if voxel.foreground_voxel_count > 0 else "FAIL",
            f"Foreground voxels: {voxel.foreground_voxel_count:,}.",
        )
    )
    fill_status = "PASS" if 0 < voxel.foreground_fraction < 0.5 else "FAIL"
    checks.append(
        AcceptanceCheck(
            "Voxel mask fill is plausible",
            fill_status,
            f"Foreground fraction: {voxel.foreground_fraction:.3%}.",
        )
    )
    voxel_volume_delta = _relative_change(after.enclosed_volume_um3, voxel.mask_volume_um3)
    checks.append(
        AcceptanceCheck(
            "Voxelized volume agrees with mesh",
            "PASS"
            if voxel_volume_delta <= 0.10
            else ("WARNING" if voxel_volume_delta <= 0.25 else "FAIL"),
            f"Relative volume difference: {voxel_volume_delta:.3%}.",
        )
    )

    checks.append(
        AcceptanceCheck(
            "Final voxel mask is one connected network",
            "PASS" if voxel.connected_component_count == 1 else "FAIL",
            f"Before filtering={voxel.initial_connected_component_count}, "
            f"after filtering={voxel.connected_component_count}.",
        )
    )

    removed_fraction = voxel.removed_island_fraction
    if removed_fraction >= voxel_config.island_fail_fraction:
        voxel_island_status = "FAIL"
    elif removed_fraction >= voxel_config.island_warning_fraction:
        voxel_island_status = "WARNING"
    else:
        voxel_island_status = "PASS"
    checks.append(
        AcceptanceCheck(
            "Voxel-island removal is within the safety limit",
            voxel_island_status,
            f"Removed {voxel.removed_island_voxel_count:,} voxels "
            f"({removed_fraction:.3%}); warning threshold="
            f"{voxel_config.island_warning_fraction:.1%}, fail threshold="
            f"{voxel_config.island_fail_fraction:.1%}.",
        )
    )

    checks.append(
        AcceptanceCheck(
            "Coarse skeleton is non-empty",
            "PASS" if skeleton.skeleton_voxel_count > 0 else "FAIL",
            f"Skeleton voxels: {skeleton.skeleton_voxel_count:,}.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Skeleton remains inside mask",
            "PASS" if skeleton.voxels_outside_mask == 0 else "FAIL",
            f"Skeleton voxels outside mask: {skeleton.voxels_outside_mask}.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Coarse skeleton is one connected network",
            "PASS" if skeleton.connected_component_count == 1 else "FAIL",
            f"Skeleton connected components: {skeleton.connected_component_count}.",
        )
    )

    missing = [str(path) for path in required_files if not path.is_file() or path.stat().st_size == 0]
    checks.append(
        AcceptanceCheck(
            "Required output files exist",
            "PASS" if not missing else "FAIL",
            "All required files were written." if not missing else f"Missing or empty: {missing}",
        )
    )
    nifti_files = [path for path in required_files if path.name.endswith(".nii.gz")]
    unreadable: list[str] = []
    for path in nifti_files:
        try:
            image = nib.load(str(path))
            if np.prod(image.shape) == 0:
                unreadable.append(str(path))
        except Exception:
            unreadable.append(str(path))
    checks.append(
        AcceptanceCheck(
            "NIfTI outputs can be reopened",
            "PASS" if not unreadable else "FAIL",
            "All NIfTI files reopened successfully."
            if not unreadable
            else f"Unreadable NIfTI files: {unreadable}",
        )
    )
    return AcceptanceResult(overall_status=_overall(checks), checks=checks)
