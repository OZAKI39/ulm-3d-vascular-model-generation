"""Create a locally extended, tagged CFD surface from saved PASS inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from utils.cfd_surface_prepare.config import load_surface_prepare_config
from utils.cfd_surface_prepare.vmtk_pipeline import (
    PASS_STATUSES,
    print_vmtk_experiment_header,
    print_vmtk_result,
    run_vmtk_surface_prepare,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create one official VMTK TPS boundary-normal candidate, remesh one "
            "active region spanning the local CORE collar and all extensions, "
            "keep FAR_CORE excluded, and cap without global remeshing."
        )
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG,
        type=Path,
        help=f"strict YAML configuration (default: {DEFAULT_CONFIG})",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        config = load_surface_prepare_config(args.config, project_root=PROJECT_ROOT)
        print_vmtk_experiment_header()
        result = run_vmtk_surface_prepare(config, project_root=PROJECT_ROOT)
        print_vmtk_result(result)
        return 0 if result.status in PASS_STATUSES else 2
    except Exception as error:
        allowed_failures = {
            "VMTK_ENVIRONMENT_BLOCKED",
            "VMTK_TPS_EXTENSION_FAILED",
            "VMTK_EXTENSION_GEOMETRY_FAILED",
            "VMTK_REMESH_CORE_FIDELITY_FAILED",
            "VMTK_SURFACE_QC_FAILED",
            "INVALID_VMTK_EXTENSION_MODE",
            "BOUNDARY_NORMAL_INPUT_PLANE_MISMATCH",
            "VMTK_BOUNDARY_NORMAL_RAW_GEOMETRY_FAILED",
            "VMTK_BOUNDARY_NORMAL_FINAL_SURFACE_FAILED",
            "INVALID_VMTK_POSTPROCESS_CONFIGURATION",
            "VMTK_RAW_CORE_NOT_EXACT_COPY",
            "VMTK_RAW_EXTENSION_MESH_QUALITY_FAILED",
            "VMTK_RAW_DIRECT_CAP_FAILED",
            "VMTK_CAPONLY_TOPOLOGY_FAILED",
            "VMTK_CAPONLY_CORE_PRESERVATION_FAILED",
            "VMTK_CAPONLY_RADIUS_FIDELITY_FAILED",
            "VMTK_CAPONLY_BOUNDARY_MAPPING_FAILED",
            "ORIGINAL_ULTRALISER_GEOMETRY_MODIFIED",
            "VMTK_ENTITY_EXCLUSION_NOT_SAFE",
            "VMTK_ENTITY_ASSIGNMENT_FAILED",
            "VMTK_ENTITY_REMESH_CORE_MODIFIED",
            "VMTK_ENTITY_REMESH_NO_EFFECT",
            "VMTK_ENTITY_REMESH_GEOMETRY_FAILED",
            "VMTK_ENTITY_REMESH_TOPOLOGY_FAILED",
            "VMTK_ENTITY_REMESH_BOUNDARY_MAPPING_FAILED",
            "VMTK_ENTITY_REMESH_RADIUS_FAILED",
            "VMTK_METER_SCALE_SERIALIZATION_FAILED",
            "ENTITY_REMESH_INTERSECTION_DETECTOR_REVIEW_REQUIRED",
            "VMTK_GUARDED_ENTITY_EXCLUSION_NOT_SAFE",
            "VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED",
            "GUARD_REGION_CLASSIFICATION_SUSPICIOUS",
            "VMTK_GUARDED_ENTITY_REMESH_CORE_MODIFIED",
            "VMTK_GUARDED_ENTITY_REMESH_GUARD_MODIFIED",
            "VMTK_GUARDED_ENTITY_REMESH_NO_EFFECT",
            "VMTK_GUARDED_ENTITY_REMESH_GEOMETRY_FAILED",
            "VMTK_GUARDED_ENTITY_REMESH_TOPOLOGY_FAILED",
            "VMTK_GUARDED_ENTITY_REMESH_BOUNDARY_MAPPING_FAILED",
            "VMTK_GUARDED_ENTITY_REMESH_RADIUS_FAILED",
            "VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED",
            "VMTK_CROSS_SEAM_FAR_CORE_MODIFIED",
            "VMTK_CROSS_SEAM_COLLAR_GEOMETRY_FAILED",
            "VMTK_CROSS_SEAM_RING_TOPOLOGY_PRESERVED",
            "VMTK_CROSS_SEAM_TOPOLOGY_FAILED",
            "VMTK_CROSS_SEAM_GEOMETRY_FAILED",
        }
        failure = str(error).split(":", maxsplit=1)[0]
        if failure not in allowed_failures:
            failure = "VMTK_TPS_EXTENSION_FAILED"
        print(f"CFD surface preparation failed: {error}")
        print(f"Final status: {failure}")
        print("NEXT: REVIEW CROSS-SEAM REMESH FAILURE")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
