"""Non-overwriting output layout for Schmid directed-graph runs."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from ..reporting.output_layout import safe_sample_name


@dataclass(frozen=True, slots=True)
class SchmidOutputLayout:
    sample_root: Path
    run_root: Path
    normalized: Path
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


def create_schmid_output_layout(
    output_root: Path, input_dir: Path, now: datetime | None = None
) -> SchmidOutputLayout:
    sample_name = safe_sample_name(input_dir)
    timestamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    sample_root = output_root / sample_name
    run_root = sample_root / f"directed_graph_run_{timestamp}"
    counter = 1
    while run_root.exists():
        run_root = sample_root / f"directed_graph_run_{timestamp}_{counter:02d}"
        counter += 1
    folders = {
        name: run_root / name
        for name in ("normalized", "graphs", "tables", "volumes", "reports", "visualizations")
    }
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    return SchmidOutputLayout(
        sample_root=sample_root,
        run_root=run_root,
        normalized=folders["normalized"],
        graphs=folders["graphs"],
        tables=folders["tables"],
        volumes=folders["volumes"],
        reports=folders["reports"],
        visualizations=folders["visualizations"],
        log_file=run_root / "directed_pkl_pipeline.log",
        status_file=run_root / "run_status.json",
        config_file=run_root / "run_config.json",
        html_report=run_root / "acceptance_report.html",
        latest_file=sample_root / "latest_directed_graph_run.json",
    )
