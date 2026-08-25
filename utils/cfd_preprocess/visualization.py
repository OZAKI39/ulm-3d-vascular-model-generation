"""The two required diagnostic figures for CFD preprocessing."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from utils.sampling.sampling_types import GlobalVascularModel

from .one_d_flow import GlobalFlowResult
from .port_transfer import PortTransfer


def global_pressure_figure(
    path: Path,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
) -> Path:
    segments = np.asarray(
        [
            [edge.upstream_position_um, edge.downstream_position_um]
            for edge in model.edges
        ],
        dtype=float,
    )
    values = np.asarray(
        [
            0.5
            * (
                flow.pressures_pa[model.node_index_by_id[edge.upstream_node_id]]
                + flow.pressures_pa[model.node_index_by_id[edge.downstream_node_id]]
            )
            for edge in model.edges
        ]
    )
    figure = plt.figure(figsize=(9, 7), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    collection = Line3DCollection(segments, cmap="viridis", linewidths=0.65)
    collection.set_array(values)
    axis.add_collection3d(collection)
    positions = model.node_positions_um
    axis.set_xlim(float(positions[:, 0].min()), float(positions[:, 0].max()))
    axis.set_ylim(float(positions[:, 1].min()), float(positions[:, 1].max()))
    axis.set_zlim(float(positions[:, 2].min()), float(positions[:, 2].max()))
    root = positions[model.node_index_by_id[flow.root_node_id]]
    axis.scatter(*root, color="red", s=38, label="ASSUMED_GLOBAL_INLET")
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    axis.set_title("Global 1D pressure (SWC parent→current simulation direction)")
    axis.legend(loc="upper right", fontsize=8)
    figure.colorbar(collection, ax=axis, shrink=0.68, label="Pressure (Pa)")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _circle(center: np.ndarray, normal: np.ndarray, radius: float) -> np.ndarray:
    helper = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(helper, normal))) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    second = np.cross(normal, first)
    angles = np.linspace(0, 2 * np.pi, 49)
    return center + radius * (
        np.cos(angles)[:, None] * first + np.sin(angles)[:, None] * second
    )


def roi_boundary_figure(
    path: Path,
    surface_path: Path,
    transfers: list[PortTransfer],
    *,
    final_status: str,
) -> Path:
    surface = pv.read(surface_path).triangulate()
    triangles = surface.faces.reshape((-1, 4))[:, 1:]
    stride = max(1, len(triangles) // 12000)
    polygons = surface.points[triangles[::stride]]
    figure = plt.figure(figsize=(10, 8), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(
        polygons,
        facecolor=(0.68, 0.76, 0.84, 0.20),
        edgecolor=(0.35, 0.42, 0.48, 0.10),
        linewidth=0.15,
    )
    axis.add_collection3d(collection)
    for item in transfers:
        color = "#d62728" if item.role == "ASSUMED_INLET" else "#1f77b4"
        circle = _circle(
            item.center_um, item.geometry.outward_normal, item.source_radius_um
        )
        axis.plot(circle[:, 0], circle[:, 1], circle[:, 2], color=color, linewidth=2.2)
        length = max(6.0, item.geometry.extension_length_um * 0.35)
        vector = item.geometry.outward_normal * length
        axis.quiver(
            *item.center_um,
            *vector,
            color=color,
            arrow_length_ratio=0.25,
            linewidth=1.8,
        )
        label = (
            f"{'INLET' if item.role == 'ASSUMED_INLET' else 'OUTLET'}\n"
            f"Q={item.flow_pl_s:.3g} pL/s\nP={item.pressure_pa:.3g} Pa"
        )
        axis.text(*(item.center_um + vector), label, color=color, fontsize=7)
    bounds = np.asarray(surface.bounds).reshape((3, 2))
    axis.set_xlim(*bounds[0])
    axis.set_ylim(*bounds[1])
    axis.set_zlim(*bounds[2])
    axis.set_box_aspect(bounds[:, 1] - bounds[:, 0])
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    axis.set_title(f"ROI boundary transfer diagnostic — {final_status}")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path
