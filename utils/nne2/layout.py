"""Non-overwriting output layout for NNE2 batch runs."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from ..reporting.output_layout import safe_sample_name


@dataclass(frozen=True, slots=True)
class NNE2OutputLayout:
    dataset_root: Path
    run_root: Path
    inventory: Path
    stacks: Path
    trees: Path
    reports: Path
    visualizations: Path
    cache: Path
    log_file: Path
    config_file: Path
    status_file: Path
    html_report: Path
    latest_file: Path


def create_nne2_output_layout(
    output_root: Path,
    input_dir: Path,
    *,
    stage: str = "all",
    now: datetime | None = None,
) -> NNE2OutputLayout:
    timestamp = (now or datetime.now()).strftime("%Y%m%d_%H%M%S")
    dataset_root = output_root / "NNE2_HDbase_v1.0"
    prefix = {
        "inventory": "inventory_run",
        "preprocess": "preprocess_run",
        "hierarchical-graph": "hierarchical_graph_run",
        "all": "nne2_hierarchy_run",
    }[stage]
    run_root = dataset_root / f"{prefix}_{timestamp}"
    counter = 1
    while run_root.exists():
        run_root = dataset_root / f"{prefix}_{timestamp}_{counter:02d}"
        counter += 1
    folders = {
        name: run_root / name
        for name in ("inventory", "stacks", "trees", "reports", "visualizations")
    }
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    cache = dataset_root / "_cache"
    cache.mkdir(parents=True, exist_ok=True)
    return NNE2OutputLayout(
        dataset_root=dataset_root,
        run_root=run_root,
        inventory=folders["inventory"],
        stacks=folders["stacks"],
        trees=folders["trees"],
        reports=folders["reports"],
        visualizations=folders["visualizations"],
        cache=cache,
        log_file=run_root / "nne2_pipeline.log",
        config_file=run_root / "run_config.json",
        status_file=run_root / "run_status.json",
        html_report=run_root / "acceptance_report.html",
        latest_file=dataset_root / f"latest_{prefix}.json",
    )


def stack_output_dir(layout: NNE2OutputLayout, stack_name: str) -> Path:
    path = layout.stacks / safe_sample_name(Path(stack_name))
    for child in ("volumes", "graphs", "tables", "reports", "visualizations"):
        (path / child).mkdir(parents=True, exist_ok=True)
    return path


def tree_output_dir(
    layout: NNE2OutputLayout, tree_key: str, stack_name: str
) -> Path:
    name = f"{tree_key}__{safe_sample_name(Path(stack_name))}"
    path = layout.trees / name
    for child in ("volumes", "graphs", "tables", "reports", "visualizations"):
        (path / child).mkdir(parents=True, exist_ok=True)
    return path
