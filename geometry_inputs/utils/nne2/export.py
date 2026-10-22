"""Portable NNE2 volume, graph, hierarchy and HTML exports."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import networkx as nx
import nibabel as nib
import numpy as np

from ..graph.export import export_hierarchical_graph
from ..graph.model import HierarchicalGraphResult
from ..io import write_csv, write_json
from .model import NNE2HierarchyResult


def write_stack_nifti(
    data_zyx: np.ndarray,
    path: Path,
    spacing_xyz_um: tuple[float, float, float],
    description: str,
    *,
    dtype: np.dtype[Any] | type | None = None,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    data_xyz = np.transpose(data_zyx, (2, 1, 0))
    affine = np.diag([*spacing_xyz_um, 1.0]).astype(np.float64)
    if dtype is None:
        dtype = np.uint8 if np.asarray(data_xyz).dtype == bool else np.asarray(data_xyz).dtype
    image = nib.Nifti1Image(np.asarray(data_xyz, dtype=dtype), affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("micron")
    image.header["descrip"] = description[:79]
    nib.save(image, str(path))
    return path


# UTF-8-safe acceptance page used by current runs.
def write_acceptance_html(
    path: Path,
    summary: dict[str, Any],
    images: list[Path],
    run_root: Path,
) -> Path:
    rows = "\n".join(
        f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in summary.items()
    )
    image_html = "\n".join(
        f'<figure><img src="{html.escape(image.relative_to(run_root).as_posix())}" '
        f'alt="{html.escape(image.stem)}"><figcaption>{html.escape(image.stem)}</figcaption></figure>'
        for image in images
        if image.is_file()
    )
    path.write_text(
        f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>NNE2 hierarchy acceptance</title>
<style>body{{font-family:Arial,sans-serif;max-width:1200px;margin:2rem auto;color:#172033}}
table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd3dd;padding:.5rem;text-align:left}}
th{{background:#eef3f8;width:34%}}figure{{margin:1.5rem 0}}img{{max-width:100%;border:1px solid #ccd3dd}}
.note{{padding:1rem;background:#fff4d6;border-left:5px solid #e8a400}}</style></head>
<body><h1>NNE2 Step 1-3 与有向 hierarchy 验收报告</h1>
<p class="note">方向来自 diving trunk、Branching Order 和三维连接关系的推断，
不是实测血流方向。缺失关键数据的记录不会进入图像处理。</p>
<table>{rows}</table><h2>可视化结果</h2>{image_html}</body></html>\n""",
        encoding="utf-8",
    )
    return path


def _graphml_ready(graph: nx.Graph) -> nx.Graph:
    output = graph.__class__()
    output.graph.update(
        {
            key: value if isinstance(value, (str, int, float, bool)) else json.dumps(value)
            for key, value in graph.graph.items()
        }
    )
    for node, data in graph.nodes(data=True):
        output.add_node(
            node,
            **{
                key: value if isinstance(value, (str, int, float, bool)) else json.dumps(value)
                for key, value in data.items()
                if value is not None
            },
        )
    if graph.is_multigraph():
        for left, right, key, data in graph.edges(keys=True, data=True):
            output.add_edge(
                left,
                right,
                key=key,
                **{
                    name: value
                    if isinstance(value, (str, int, float, bool))
                    else json.dumps(value)
                    for name, value in data.items()
                    if value is not None
                },
            )
    else:
        for left, right, data in graph.edges(data=True):
            output.add_edge(
                left,
                right,
                **{
                    name: value
                    if isinstance(value, (str, int, float, bool))
                    else json.dumps(value)
                    for name, value in data.items()
                    if value is not None
                },
            )
    return output


def export_undirected_graph(
    graph: HierarchicalGraphResult,
    output_dir: Path,
    *,
    save_graphml: bool = True,
    save_vtp: bool = True,
    save_npz: bool = True,
) -> list[Path]:
    """Export the complete common Step 3 representation plus an NNE2 report alias."""
    graphs = output_dir / "graphs"
    tables = output_dir / "tables"
    reports = output_dir / "reports"
    exported = export_hierarchical_graph(
        graph,
        graphs,
        tables,
        save_graphml=save_graphml,
        save_vtp=save_vtp,
        save_npz=save_npz,
    )
    report_json = reports / "undirected_branch_graph_report.json"
    write_json(graph.report(), report_json)
    return [*exported, report_json]


def export_hierarchy(result: NNE2HierarchyResult, output_dir: Path) -> list[Path]:
    graphs = output_dir / "graphs"
    tables = output_dir / "tables"
    reports = output_dir / "reports"
    directed_graphml = graphs / "directed_hierarchy.graphml"
    primary_graphml = graphs / "primary_parent_tree.graphml"
    nx.write_graphml(_graphml_ready(result.directed_graph), directed_graphml)
    nx.write_graphml(_graphml_ready(result.primary_tree), primary_graphml)
    hierarchy_json = graphs / "directed_hierarchy.json"
    write_json(
        {
            "schema_version": "1.0",
            "representation": {
                "name": "NNE2 directed anatomical hierarchy",
                "tree_key": result.tree_key,
                "stack_name": result.stack_name,
                "directed": True,
                "measured_flow_direction": False,
                "is_strict_tree": False,
                "ambiguous_links_preserved": True,
            },
            "summary": result.report(),
            "anchors": [item.to_dict() for item in result.anchors],
            "branches": [item.to_dict() for item in result.branches],
            "parent_child_relations": [item.to_dict() for item in result.relations],
            "excluded_disconnected_branch_ids": result.excluded_branch_ids,
        },
        hierarchy_json,
    )
    branches_csv = tables / "directed_branches.csv"
    relations_csv = tables / "parent_child_relations.csv"
    anchors_csv = tables / "measurement_anchor_matches.csv"
    unresolved_csv = tables / "unresolved_and_cross_links.csv"
    write_csv([item.to_dict() for item in result.branches], branches_csv)
    write_csv([item.to_dict() for item in result.relations], relations_csv)
    write_csv([item.to_dict() for item in result.anchors], anchors_csv)
    write_csv(
        [
            {
                "branch_id": item.branch_id,
                "direction_status": item.direction_status,
                "confidence": item.confidence,
                "is_cross_link": item.is_cross_link,
                "order_source": item.order_source,
            }
            for item in result.branches
            if item.branch_id in result.unresolved_branch_ids
            or item.branch_id in result.cross_link_branch_ids
            or item.branch_id in result.order_conflict_branch_ids
        ],
        unresolved_csv,
    )
    report_json = reports / "hierarchy_report.json"
    write_json(result.report(), report_json)
    return [
        directed_graphml,
        primary_graphml,
        hierarchy_json,
        branches_csv,
        relations_csv,
        anchors_csv,
        unresolved_csv,
        report_json,
    ]


def _write_acceptance_html_mojibake(
    path: Path,
    summary: dict[str, Any],
    images: list[Path],
    run_root: Path,
) -> Path:
    rows = "\n".join(
        f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in summary.items()
    )
    image_html = "\n".join(
        f'<figure><img src="{html.escape(image.relative_to(run_root).as_posix())}" '
        f'alt="{html.escape(image.stem)}"><figcaption>{html.escape(image.stem)}</figcaption></figure>'
        for image in images
        if image.is_file()
    )
    path.write_text(
        f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>NNE2 hierarchy acceptance</title>
<style>body{{font-family:Arial,sans-serif;max-width:1200px;margin:2rem auto;color:#172033}}
table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd3dd;padding:.5rem;text-align:left}}
th{{background:#eef3f8;width:34%}}figure{{margin:1.5rem 0}}img{{max-width:100%;border:1px solid #ccd3dd}}
.note{{padding:1rem;background:#fff4d6;border-left:5px solid #e8a400}}</style></head>
<body><h1>NNE2 Step 1-3 与有向 hierarchy 验收报告</h1>
<p class="note">方向来自 diving trunk、Branching Order 和三维连接关系的推断，不是实测血流方向。缺失关键数据的记录未参与处理。</p>
<table>{rows}</table><h2>可视化结果</h2>{image_html}</body></html>\n""",
        encoding="utf-8",
    )
    return path
