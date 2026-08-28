"""Exactly three concise post-solution flow-review figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from .geometry import SurfacePartition


def _scatter_review(grid: pv.DataSet, values: np.ndarray, path: Path, title: str, label: str) -> None:
    centers = np.asarray(grid.cell_centers().points)
    stride = max(1, len(centers) // 80_000)
    points = centers[::stride]
    scalar = np.asarray(values)[::stride]
    figure = plt.figure(figsize=(10, 8), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    artist = axis.scatter(points[:, 0] * 1e6, points[:, 1] * 1e6, points[:, 2] * 1e6, c=scalar, s=1.0, cmap="viridis")
    figure.colorbar(artist, ax=axis, shrink=0.7, label=label)
    axis.set_title(title)
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def create_flow_figures(grid: pv.UnstructuredGrid, partition: SurfacePartition, output: Path) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    velocity = np.asarray(grid.cell_data["velocity_phy"])
    speed = np.linalg.norm(velocity, axis=1)
    pressure = np.asarray(grid.cell_data["pressure_gauge_pa"])
    velocity_path = output / "velocity_magnitude_review.png"
    pressure_path = output / "gauge_pressure_review.png"
    streamline_path = output / "streamlines_from_inlet.png"
    _scatter_review(grid, speed, velocity_path, "Steady velocity magnitude", "|u| (m/s)")
    _scatter_review(grid, pressure, pressure_path, "Steady gauge pressure", "p_gauge (Pa)")

    inlet = partition.patch("inlet")
    center = inlet.center_um * 1.0e-6 - inlet.outward_normal * 2.0e-7
    inward = -inlet.outward_normal
    reference = np.array([1.0, 0.0, 0.0]) if abs(inward[0]) < 0.8 else np.array([0.0, 1.0, 0.0])
    basis_1 = np.cross(inward, reference)
    basis_1 /= np.linalg.norm(basis_1)
    basis_2 = np.cross(inward, basis_1)
    radius = inlet.equivalent_radius_um * 0.75e-6
    angles = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    seeds = np.vstack([center, *(center + radius * (np.cos(angle) * basis_1 + np.sin(angle) * basis_2) for angle in angles)])
    source = pv.PolyData(seeds)
    cell_grid = grid.cell_data_to_point_data(pass_cell_data=True)
    stream = cell_grid.streamlines_from_source(source, vectors="velocity_phy", max_time=2.0e-4, initial_step_length=2.0e-7)
    figure = plt.figure(figsize=(10, 8), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    if stream.n_cells:
        points = np.asarray(stream.points) * 1.0e6
        lines = np.asarray(stream.lines, dtype=np.int64)
        offset = 0
        while offset < len(lines):
            point_count = int(lines[offset])
            line = lines[offset + 1 : offset + 1 + point_count]
            offset += point_count + 1
            segment = points[line]
            axis.plot(segment[:, 0], segment[:, 1], segment[:, 2], linewidth=0.7)
    axis.scatter(seeds[:, 0] * 1e6, seeds[:, 1] * 1e6, seeds[:, 2] * 1e6, c="red", s=10)
    axis.set_title("Streamlines seeded from inlet")
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_zlabel("z (µm)")
    figure.savefig(streamline_path, dpi=180)
    plt.close(figure)
    return [velocity_path, pressure_path, streamline_path]
