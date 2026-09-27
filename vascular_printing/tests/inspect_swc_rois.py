"""Inspect three saved ROI selections using ONLY original SWC geometry.

Run with the project's .venv Python; see README_swc_roi_inspection.md.
No production renderer, processed geometry, clipping nodes or spline is used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pyvista as pv
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TUBE_SIDES = 48
MODES = ("1: SWC tubes (straight segments)", "2: SWC centerlines + nodes", "3: SWC radius spheres")


@dataclass(frozen=True)
class Node:
    id: int
    kind: int
    xyz: tuple[float, float, float]
    radius: float
    parent: int
    line: int


@dataclass
class Region:
    anchor: int
    roi_id: str
    nodes: list[Node]
    edges: list[tuple[int, int]]  # original parent ID, original child ID


def read_swc(path: Path) -> dict[int, Node]:
    """Keep source numbers unchanged, including IDs and physical line numbers."""
    nodes = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        fields = raw.partition("#")[0].split()
        if not fields:
            continue
        if len(fields) != 7:
            raise ValueError(f"{path}:{line_number}: expected seven SWC columns")
        values = np.asarray([float(value) for value in fields])
        if not np.isfinite(values).all() or any(values[i] != int(values[i]) for i in (0, 1, 6)):
            raise ValueError(f"{path}:{line_number}: nonfinite values or noninteger IDs")
        node = Node(int(values[0]), int(values[1]), tuple(values[2:5]), float(values[5]), int(values[6]), line_number)
        if node.id in nodes or node.id == node.parent:
            raise ValueError(f"{path}:{line_number}: duplicate ID or self parent")
        if node.radius <= 0:
            raise ValueError(f"{path}:{line_number}: radius <= 0; inspector will not replace source radii")
        nodes[node.id] = node
    if not nodes:
        raise ValueError(f"Empty SWC: {path}")
    for node in nodes.values():
        if node.parent != -1 and node.parent not in nodes:
            raise ValueError(f"Node {node.id}: missing parent {node.parent}")
    return nodes


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def latest_run(root: Path = PROJECT_ROOT) -> Path:
    for run in sorted((root / "outputs/human_brava/sampling").glob("*"), reverse=True):
        summary = run / "report/sampling_summary.json"
        manifest = run / "manifests/selected_rois.csv"
        if summary.is_file() and manifest.is_file():
            info = json.loads(summary.read_text())
            if info.get("status") == "PASS" and len(read_rows(manifest)) >= 3:
                return run
    raise ValueError("No successful human sampling run with >= 3 ROIs; pass --run-dir")


def choose_rows(run: Path, anchors: list[int] | None) -> list[dict[str, str]]:
    if anchors is None:
        rows = sorted(read_rows(run / "manifests/selected_rois.csv"), key=lambda row: int(row["selection_rank"]))[:3]
    else:
        if len(anchors) != 3 or len(set(anchors)) != 3:
            raise ValueError("Choose exactly three distinct anchor IDs")
        available = read_rows(run / "manifests/candidate_rois.csv")
        rows = []
        for anchor in anchors:
            matches = [row for row in available if int(row["anchor_id"]) == anchor]
            if len(matches) != 1:
                raise ValueError(f"Anchor {anchor}: expected one candidate ROI, found {len(matches)}")
            rows.append(matches[0])
    if len(rows) != 3 or len({row["source_model_id"] for row in rows}) != 1:
        raise ValueError("Three ROIs from the same SWC sample are required")
    return rows


def source_for(run: Path, row: dict[str, str]) -> tuple[Path, str]:
    config = yaml.safe_load((run / "config/source_swc_roi_generate.yaml").read_text())
    directory = Path(config["paths"]["input_dir"])
    if not directory.is_absolute():
        directory = PROJECT_ROOT / directory
    unit = config["coordinates"]["swc_units"]
    if unit not in ("mm", "um"):
        raise ValueError("Inspector requires an SWC in physical mm or um units")
    cohort = config["pipeline"]["cohort"]
    prefix = f"{cohort}__"
    model = row["source_model_id"]
    if not model.startswith(prefix):
        raise ValueError("ROI source model does not match the saved source configuration")
    return directory / f"{model.removeprefix(prefix)}.swc", unit


def extract_region(row: dict[str, str], source: dict[int, Node], unit: str) -> Region:
    ids = [int(value) for value in row["global_node_ids"].split(";") if value]
    if not ids or len(set(ids)) != len(ids) or any(value not in source for value in ids):
        raise ValueError(f"Invalid original SWC membership in {row['roi_id']}")
    nodes = [source[value] for value in ids]
    xyz = np.asarray([node.xyz for node in nodes]) * (1000 if unit == "mm" else 1)
    lower = np.asarray([float(row[f"bbox_min_{axis}_um"]) for axis in "xyz"])
    upper = np.asarray([float(row[f"bbox_max_{axis}_um"]) for axis in "xyz"])
    if not ((xyz >= lower - 1e-7) & (xyz <= upper + 1e-7)).all():
        raise ValueError("Original node coordinates disagree with saved ROI bounds")
    members = set(ids)
    edges = [(node.parent, node.id) for node in nodes if node.parent in members]
    if not edges:
        raise ValueError(f"No complete source SWC edges in {row['roi_id']}")
    return Region(int(row["anchor_id"]), row["roi_id"], nodes, edges)


def branch_paths(region: Region) -> list[list[int]]:
    """Join degree-two source edges without adding or moving a single node."""
    adjacent = {node.id: [] for node in region.nodes}
    for parent, child in region.edges:
        adjacent[parent].append(child)
        adjacent[child].append(parent)
    visited = set()
    paths = []

    def walk(start, following):
        path = [start, following]
        visited.add(tuple(sorted((start, following))))
        while len(adjacent[path[-1]]) == 2:
            nxt = next(node for node in adjacent[path[-1]] if node != path[-2])
            edge = tuple(sorted((path[-1], nxt)))
            if edge in visited:
                break
            visited.add(edge)
            path.append(nxt)
        paths.append(path)

    for start, neighbors in adjacent.items():
        if len(neighbors) != 2:
            for neighbor in neighbors:
                if tuple(sorted((start, neighbor))) not in visited:
                    walk(start, neighbor)
    for parent, child in region.edges:  # also handle a closed component
        if tuple(sorted((parent, child))) not in visited:
            walk(parent, child)
    return paths


def geometry(region: Region, panel: int) -> tuple[pv.PolyData, pv.PolyData, pv.PolyData, pv.PolyData]:
    nodes = {node.id: node for node in region.nodes}
    source_ids, cells = [], []
    for path in branch_paths(region):
        offset = len(source_ids)
        source_ids.extend(path)
        cells.extend([len(path), *range(offset, offset + len(path))])
    line = pv.PolyData(np.asarray([nodes[node].xyz for node in source_ids], dtype=float))
    line.lines = cells
    line["source_node_id"] = np.asarray(source_ids, dtype=np.int64)
    line["radius"] = np.asarray([nodes[node].radius for node in source_ids])
    line["roi_panel"] = np.full(len(source_ids), panel, dtype=np.int32)
    # VTK joins straight segments; rings at all original samples have exact r.
    # No smoothing, interpolation knots, internal caps or implicit-surface union.
    tube = line.tube(scalars="radius", absolute=True, n_sides=TUBE_SIDES, capping=False)
    points = pv.PolyData(np.asarray([node.xyz for node in region.nodes], dtype=float))
    points["source_node_id"] = np.asarray([node.id for node in region.nodes], dtype=np.int64)
    points["radius"] = np.asarray([node.radius for node in region.nodes])
    points["roi_panel"] = np.full(len(region.nodes), panel, dtype=np.int32)
    # Double precision glyphs limit float32 coordinate error in displayed spheres.
    sphere = pv.Sphere(radius=1, theta_resolution=32, phi_resolution=24)
    sphere.points = sphere.points.astype(np.float64)
    spheres = points.glyph(geom=sphere, scale="radius", orient=False, factor=1)
    return line, tube, points, spheres


def write_audit(output: Path, run: Path, source_path: Path, unit: str, regions: list[Region]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    report = {"source_swc": str(source_path.resolve()), "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
              "sampling_run": str(run.resolve()), "unit": unit, "geometry": "original SWC nodes and parent edges; no smoothing",
              "roi_membership": "saved retained original node IDs; synthetic cut nodes and border-crossing edges omitted", "regions": []}
    for region in regions:
        with (output / f"roi_{region.anchor:06d}_source_nodes.csv").open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            writer.writerow(["swc_node_id", "type", "x", "y", "z", "radius", "parent_id", "source_line", "unit"])
            writer.writerows((n.id, n.kind, *n.xyz, n.radius, n.parent, n.line, unit) for n in region.nodes)
        nodes = {node.id: node for node in region.nodes}
        changes = sorted([(parent, child, nodes[parent].radius, nodes[child].radius,
                           nodes[child].radius / nodes[parent].radius,
                           float(np.linalg.norm(np.asarray(nodes[child].xyz) - nodes[parent].xyz)))
                          for parent, child in region.edges], key=lambda row: abs(np.log(row[4])), reverse=True)
        with (output / f"roi_{region.anchor:06d}_radius_changes.csv").open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            writer.writerow(["parent_id", "child_id", "parent_radius", "child_radius", "child_parent_ratio", "distance"])
            writer.writerows(changes)
        report["regions"].append({"anchor": region.anchor, "roi_id": region.roi_id, "original_nodes": len(region.nodes),
                                  "original_edges": len(region.edges), "largest_adjacent_change": changes[0]})
    (output / "source_audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")


class Inspector:
    def __init__(self, source: dict[int, Node], regions: list[Region], source_path: Path, unit: str,
                 output: Path, off_screen: bool = False):
        self.source, self.regions, self.unit, self.output = source, regions, unit, output
        self.mode, self.colored = 0, False
        self.plotter = pv.Plotter(shape=(2, 3), row_weights=(0.7, 0.3), window_size=(1800, 1000), border=False, off_screen=off_screen,
                                 title="Original SWC - three ROI inspector")
        self.actors, self.cameras = [], []
        self.children = {node_id: [] for node_id in source}
        for node in source.values():
            if node.parent in self.children:
                self.children[node.parent].append(node.id)
        radii = [node.radius for region in regions for node in region.nodes]
        self.clim = (min(radii), max(radii))
        if self.clim[0] == self.clim[1]:
            self.clim = (0, self.clim[1])
        for panel, region in enumerate(regions):
            self.plotter.subplot(0, panel)
            self.plotter.set_background("#111923")
            line, tube, points, spheres = geometry(region, panel)
            shared = dict(scalars="radius", clim=self.clim, cmap="viridis", show_scalar_bar=False,
                          opacity=1, pickable=True, ambient=0.35, diffuse=0.65)
            tube_actor = self.plotter.add_mesh(tube, **shared)
            line_actor = self.plotter.add_mesh(line, color="#a5b5c4", line_width=2, pickable=False)
            node_actor = self.plotter.add_mesh(points, point_size=7, render_points_as_spheres=True, **shared)
            sphere_actor = self.plotter.add_mesh(spheres, **shared)
            self.actors.append((tube_actor, line_actor, node_actor, sphere_actor))
            self.plotter.add_text(f"ROI {panel + 1} / anchor {region.anchor}\n{source_path.name} / {unit}\n"
                                  f"{len(region.nodes)} source nodes / {len(region.edges)} edges",
                                  position="upper_left", font_size=11, color="white")
            self.plotter.view_isometric()
            self.plotter.reset_camera()
            self.cameras.append(self.plotter.camera_position)
            self.plotter.subplot(1, panel)
            self.plotter.set_background("#111923")
            self.plotter.add_text("1 Tubes   2 Nodes   3 Spheres\nC Color   R Reset   S Save\n"
                                  "Left drag: rotate / Wheel: zoom\nRight click / P: inspect source node",
                                  position=(12, 145), font_size=10, color="#b5c5d5")
            self.show_node(panel, region.anchor)
        self.update_view()
        self.plotter.enable_point_picking(callback=self.picked, picker="point", use_picker=True,
                                          show_message=False, show_point=False, tolerance=0.018)
        for index in range(3):
            self.plotter.add_key_event(str(index + 1), lambda index=index: self.set_mode(index))
        self.plotter.add_key_event("c", self.toggle_color)
        self.plotter.add_key_event("r", self.reset)
        self.plotter.add_key_event("s", self.screenshot)

    def update_view(self):
        for panel, actors in enumerate(self.actors):
            tube, line, points, spheres = actors
            for actor, visible in zip(actors, (self.mode == 0, self.mode == 1, self.mode == 1, self.mode == 2)):
                actor.SetVisibility(visible)
            for actor in (tube, points, spheres):
                actor.mapper.SetScalarVisibility(self.colored)
                actor.prop.color = "#61c9bc"
            self.plotter.subplot(1, panel)
            color_text = f"radius color: {self.clim[0]:g} - {self.clim[1]:g} {self.unit}" if self.colored else "uniform color / original radii"
            self.plotter.add_text(MODES[self.mode] + "\n" + color_text, position=(12, 245),
                                  font_size=10, color="#74ddcc", name="mode", render=False)
        self.plotter.render()

    def set_mode(self, mode: int):
        self.mode = mode
        self.update_view()

    def toggle_color(self):
        self.colored = not self.colored
        self.update_view()

    def reset(self):
        for panel, camera in enumerate(self.cameras):
            self.plotter.subplot(0, panel)
            self.plotter.camera_position = camera
        self.plotter.render()

    def show_node(self, panel: int, node_id: int):
        node = self.source[node_id]
        parent = self.source.get(node.parent)
        parent_text = "root" if parent is None else f"{parent.id}, r={parent.radius:g}, ratio={node.radius / parent.radius:.3g}"
        children = ", ".join(f"{child}:{self.source[child].radius:g}" for child in self.children[node_id]) or "none"
        message = (f"SWC node {node.id} / file line {node.line}\n"
                   f"r = {node.radius:g} {self.unit} / d = {2 * node.radius:g} {self.unit}\n"
                   f"Parent: {parent_text}\nChildren (id:r): {children}\n"
                   f"xyz = ({node.xyz[0]:g}, {node.xyz[1]:g}, {node.xyz[2]:g})")
        self.plotter.subplot(1, panel)
        self.plotter.add_text(message, position=(12, 15), font_size=10, color="white", name="source_readout", render=False)

    def picked(self, point, picker):
        dataset, point_id = picker.GetDataSet(), picker.GetPointId()
        if dataset is None or point_id < 0:
            return
        data = pv.wrap(dataset)
        if "source_node_id" not in data.point_data or "roi_panel" not in data.point_data:
            return
        node_id, panel = int(data["source_node_id"][point_id]), int(data["roi_panel"][point_id])
        self.show_node(panel, node_id)
        self.plotter.subplot(0, panel)
        self.plotter.add_points(np.asarray([self.source[node_id].xyz]), color="#ffb34d", point_size=13,
                                render_points_as_spheres=True, pickable=False, name="picked_source_node")
        self.plotter.render()

    def screenshot(self):
        path = self.output / f"inspection_mode{self.mode + 1}_{datetime.now():%Y%m%d_%H%M%S_%f}.png"
        self.plotter.screenshot(str(path))
        print(f"Screenshot: {path}", flush=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, help="Saved human sampling run; default: latest successful run")
    parser.add_argument("--anchors", nargs=3, type=int, help="Three original anchor IDs; default: selection ranks 1-3")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--off-screen", action="store_true", help="Render all three modes and exit")
    args = parser.parse_args(argv)
    run = args.run_dir.resolve() if args.run_dir else latest_run()
    rows = choose_rows(run, args.anchors)
    source_path, unit = source_for(run, rows[0])
    source = read_swc(source_path)
    regions = [extract_region(row, source, unit) for row in rows]
    output = args.output_dir or PROJECT_ROOT / "tests/swc_roi_inspection_output" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    write_audit(output, run, source_path, unit, regions)
    viewer = Inspector(source, regions, source_path, unit, output, args.off_screen)
    print(f"原始 SWC: {source_path}\nROI 锚点: {[region.anchor for region in regions]}\n"
          f"单位: {unit}；节点数据与半径审计: {output}\n"
          "1 管面，2 中心线与节点，3 原半径球；右键选择节点；C 颜色，R 复位，S 截图。", flush=True)
    if args.off_screen:
        for mode in range(3):
            viewer.set_mode(mode)
            viewer.plotter.screenshot(str(output / f"mode_{mode + 1}.png"))
        viewer.plotter.close()
    else:
        viewer.plotter.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
