"""Create a stable, non-overwriting layout for pipeline outputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


@dataclass(frozen=True, slots=True)
class OutputLayout:
    sample_root: Path
    run_root: Path
    meshes: Path
    volumes: Path
    reports: Path
    visualizations: Path
    log_file: Path
    status_file: Path
    config_file: Path
    html_report: Path
    latest_file: Path


@dataclass(frozen=True, slots=True)
class HierarchicalGraphOutputLayout:
    sample_root: Path
    run_root: Path
    graphs: Path
    tables: Path
    volumes: Path
    reports: Path
    visualizations: Path
    log_file: Path
    status_file: Path
    config_file: Path
    html_report: Path
    latest_file: Path


def safe_sample_name(path: Path) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem).strip("._")
    return name or "vascular_sample"


def create_output_layout(output_root: Path, input_stl: Path, now: datetime | None = None) -> OutputLayout:
    timestamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    sample_root = output_root / safe_sample_name(input_stl)
    run_root = sample_root / f"run_{timestamp}"
    counter = 1
    while run_root.exists():
        run_root = sample_root / f"run_{timestamp}_{counter:02d}"
        counter += 1

    meshes = run_root / "meshes"
    volumes = run_root / "volumes"
    reports = run_root / "reports"
    visualizations = run_root / "visualizations"
    for folder in (meshes, volumes, reports, visualizations):
        folder.mkdir(parents=True, exist_ok=False)

    return OutputLayout(
        sample_root=sample_root,
        run_root=run_root,
        meshes=meshes,
        volumes=volumes,
        reports=reports,
        visualizations=visualizations,
        log_file=run_root / "pipeline.log",
        status_file=run_root / "run_status.json",
        config_file=run_root / "run_config.json",
        html_report=run_root / "acceptance_report.html",
        latest_file=sample_root / "latest_run.json",
    )


def create_hierarchical_graph_output_layout(
    source_run: Path, now: datetime | None = None
) -> HierarchicalGraphOutputLayout:
    source_run = Path(source_run).resolve()
    sample_root = source_run.parent
    timestamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    run_root = sample_root / f"hierarchical_graph_run_{timestamp}"
    counter = 1
    while run_root.exists():
        run_root = sample_root / f"hierarchical_graph_run_{timestamp}_{counter:02d}"
        counter += 1

    graphs = run_root / "graphs"
    tables = run_root / "tables"
    volumes = run_root / "volumes"
    reports = run_root / "reports"
    visualizations = run_root / "visualizations"
    for folder in (graphs, tables, volumes, reports, visualizations):
        folder.mkdir(parents=True, exist_ok=False)

    return HierarchicalGraphOutputLayout(
        sample_root=sample_root,
        run_root=run_root,
        graphs=graphs,
        tables=tables,
        volumes=volumes,
        reports=reports,
        visualizations=visualizations,
        log_file=run_root / "hierarchical_graph.log",
        status_file=run_root / "run_status.json",
        config_file=run_root / "run_config.json",
        html_report=run_root / "acceptance_report.html",
        latest_file=sample_root / "latest_hierarchical_graph_run.json",
    )
