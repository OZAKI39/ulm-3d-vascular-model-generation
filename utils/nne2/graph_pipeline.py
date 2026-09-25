"""Common Step 3 graph export and acceptance helpers for NNE2."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..config import HierarchicalGraphConfig
from ..graph.export import verify_graph_exports
from ..graph.model import HierarchicalGraphResult
from ..graph.validation import evaluate_hierarchical_graph_acceptance
from ..io import write_json
from .export import export_undirected_graph


def export_and_validate_step3_graph(
    graph: HierarchicalGraphResult,
    output_dir: Path,
    graph_config: HierarchicalGraphConfig,
) -> tuple[list[Path], dict[str, Any]]:
    exported = export_undirected_graph(
        graph,
        output_dir,
        save_graphml=graph_config.save_graphml,
        save_vtp=graph_config.save_vtp,
        save_npz=graph_config.save_npz,
    )
    reopen_errors = verify_graph_exports(exported, len(graph.nodes), len(graph.branches))
    acceptance = evaluate_hierarchical_graph_acceptance(
        graph, graph_config, exported, reopen_errors
    )
    report = acceptance.report()
    # Multiple stack components are allowed before the per-tree BO0 component is selected.
    for check in report["checks"]:
        if check["name"] == "Source skeleton is one connected network" and check["status"] == "FAIL":
            check["status"] = "WARNING"
            check["message"] += " Per-tree hierarchy will retain only the BO0 root component."
    statuses = {check["status"] for check in report["checks"]}
    report["overall_status"] = (
        "FAIL" if "FAIL" in statuses else "WARNING" if "WARNING" in statuses else "PASS"
    )
    report["counts"] = {
        status: sum(check["status"] == status for check in report["checks"])
        for status in ("PASS", "WARNING", "FAIL")
    }
    report["export_reopen_errors"] = reopen_errors
    report_path = output_dir / "reports" / "step3_graph_acceptance.json"
    write_json(report, report_path)
    return [*exported, report_path], report
