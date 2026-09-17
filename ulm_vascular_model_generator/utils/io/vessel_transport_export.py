from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
from ..core.models import Vessel


def vessel_transport_path_from_swc(swc_path: Path) -> Path:
    """
    Given a SWC file path, return the corresponding machine-readable vessel information path.
    """

    return Path(swc_path).with_suffix(".vessels.npz")


def write_vessel_transport_npz(vessels: list[Vessel], path: Path, metadata: dict[str, Any] | None = None) -> None:
    """
    把 Vessel 列表保存成供微泡轨迹生成器直接读取的 npz 文件。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(vessels)
    parent_id = np.asarray([v.parent_id for v in vessels], dtype=int)
    children_count = np.asarray([len(v.children) for v in vessels], dtype=int)
    max_children = int(children_count.max(initial=0))
    children = np.full((n, max_children), -1, dtype=int)
    for row, vessel in enumerate(vessels):
        if vessel.children:
            children[row, : len(vessel.children)] = np.asarray(vessel.children, dtype=int)

    metadata = {} if metadata is None else dict(metadata)
    geometry_mode = str(metadata.get("geometry_mode", "")).strip().lower()
    if geometry_mode not in {"planar_2d", "volumetric_3d"}:
        raise ValueError(
            "Transport format version 3 requires metadata geometry_mode to be "
            "'planar_2d' or 'volumetric_3d'."
        )
    flow_quantity = str(metadata.get("flow_quantity", "")).strip().lower()
    required_quantity = (
        "planar_flux_per_unit_depth"
        if geometry_mode == "planar_2d"
        else "volume_flow"
    )
    if flow_quantity != required_quantity:
        raise ValueError(
            f"{geometry_mode} transport data must declare "
            f"flow_quantity={required_quantity!r}; found {flow_quantity!r}."
        )
    payload = {
        "format_version": np.asarray([3], dtype=int),
        "vessel_id": np.asarray([v.vid for v in vessels], dtype=int),
        "parent_id": parent_id,
        "children": children,
        "children_count": children_count,
        "x_p": np.vstack([np.asarray(v.x_p, dtype=float) for v in vessels])
        if vessels
        else np.empty((0, 3)),
        "x_d": np.vstack([np.asarray(v.x_d, dtype=float) for v in vessels])
        if vessels
        else np.empty((0, 3)),
        "branching_mode": np.asarray(
            [v.branching_mode for v in vessels], dtype=str
        ),
        "role": np.asarray([v.role for v in vessels], dtype=str),
        "is_main_trunk": np.asarray(
            [v.is_main_trunk for v in vessels], dtype=bool
        ),
        "mean_velocity_um_s": np.asarray(
            [v.mean_velocity for v in vessels], dtype=float
        ),
        "flow_conservation_residual": np.asarray(
            [v.flow_conservation_residual for v in vessels], dtype=float
        ),
        "murray_residual": np.asarray(
            [v.murray_residual for v in vessels], dtype=float
        ),
        "metadata_keys": np.asarray(list(metadata.keys()), dtype=str),
        "metadata_values": np.asarray(
            [str(value) for value in metadata.values()], dtype=str
        ),
    }
    if geometry_mode == "planar_2d":
        payload.update(
            {
                "half_width_um": np.asarray(
                    [v.radius for v in vessels], dtype=float
                ),
                "prescribed_outflux_um2_s": np.asarray(
                    [v.prescribed_outflow for v in vessels], dtype=float
                ),
                "planar_flux_um2_s": np.asarray(
                    [v.flow_rate for v in vessels], dtype=float
                ),
            }
        )
    else:
        payload.update(
            {
                "radius_um": np.asarray(
                    [v.radius for v in vessels], dtype=float
                ),
                "prescribed_outflow_um3_s": np.asarray(
                    [v.prescribed_outflow for v in vessels], dtype=float
                ),
                "volume_flow_um3_s": np.asarray(
                    [v.flow_rate for v in vessels], dtype=float
                ),
            }
        )
    np.savez_compressed(path, **payload)


def read_vessel_transport_npz(path: Path) -> list[Vessel]:
    """从 npz 文件恢复 Vessel 列表，供下游仿真直接使用。"""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"找不到血管信息文件: {path}")

    with np.load(path, allow_pickle=False) as data:
        if "format_version" not in data.files:
            raise ValueError("Vessel transport data must declare format_version.")
        version = int(np.asarray(data["format_version"], dtype=int)[0])
        if version not in {2, 3}:
            raise ValueError(
                "Vessel transport format_version must be 2 or 3; "
                f"found {version}."
            )
        if version == 3:
            keys = np.asarray(data["metadata_keys"], dtype=str)
            values = np.asarray(data["metadata_values"], dtype=str)
            metadata = dict(zip(keys.tolist(), values.tolist(), strict=True))
            geometry_mode = str(metadata.get("geometry_mode", "")).strip().lower()
            flow_quantity = str(
                metadata.get("flow_quantity", "")
            ).strip().lower()
            if geometry_mode == "planar_2d":
                if flow_quantity != "planar_flux_per_unit_depth":
                    raise ValueError(
                        "planar_2d format-v3 data must declare "
                        "flow_quantity='planar_flux_per_unit_depth'."
                    )
                size_values = data["half_width_um"]
                prescribed_values = data["prescribed_outflux_um2_s"]
                flow_values = data["planar_flux_um2_s"]
            elif geometry_mode == "volumetric_3d":
                if flow_quantity != "volume_flow":
                    raise ValueError(
                        "volumetric_3d format-v3 data must declare "
                        "flow_quantity='volume_flow'."
                    )
                size_values = data["radius_um"]
                prescribed_values = data["prescribed_outflow_um3_s"]
                flow_values = data["volume_flow_um3_s"]
            else:
                raise ValueError(
                    "Transport format version 3 requires valid geometry_mode metadata."
                )
        else:
            size_values = data["radius_um"]
            prescribed_values = data["prescribed_outflow"]
            flow_values = data["flow_rate"]
        vessel_ids = np.asarray(data["vessel_id"], dtype=int)
        vessels: list[Vessel] = []
        for row, vid in enumerate(vessel_ids):
            vessel = Vessel(
                vid=int(vid),
                parent_id=int(data["parent_id"][row]),
                children=[],
                x_p=np.asarray(data["x_p"][row], dtype=float),
                x_d=np.asarray(data["x_d"][row], dtype=float),
                radius=float(size_values[row]),
                branching_mode=str(data["branching_mode"][row]),
                role=str(data["role"][row]),
                is_main_trunk=bool(data["is_main_trunk"][row]),
                prescribed_outflow=float(prescribed_values[row]),
                flow_rate=float(flow_values[row]),
                mean_velocity=float(data["mean_velocity_um_s"][row]),
                flow_conservation_residual=float(
                    data["flow_conservation_residual"][row]
                ),
                murray_residual=float(data["murray_residual"][row]),
            )
            vessels.append(vessel)

        child_counts = np.asarray(data["children_count"], dtype=int)
        children = np.asarray(data["children"], dtype=int)
        for row, vessel in enumerate(vessels):
            count = int(child_counts[row])
            vessel.children = [int(child) for child in children[row, :count] if int(child) >= 0]

    return vessels


def read_vessel_transport_metadata(path: Path) -> dict[str, str]:
    """Read the scalar provenance metadata stored beside the vessel arrays."""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"找不到血管信息文件: {path}")
    with np.load(path, allow_pickle=False) as data:
        keys = np.asarray(data.get("metadata_keys", np.empty(0, dtype=str)), dtype=str)
        values = np.asarray(
            data.get("metadata_values", np.empty(0, dtype=str)), dtype=str
        )
    if keys.ndim != 1 or values.ndim != 1 or keys.shape != values.shape:
        raise ValueError("Vessel transport metadata keys and values must be equal-length vectors.")
    if len(set(keys.tolist())) != int(keys.size):
        raise ValueError("Vessel transport metadata contains duplicate keys.")
    return dict(zip(keys.tolist(), values.tolist(), strict=True))
