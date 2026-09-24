"""Standalone HTML acceptance report for Step 3."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from jinja2 import Template

from .acceptance import AcceptanceResult


_TEMPLATE = Template(
    """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Hierarchical vascular representation acceptance report</title>
  <style>
    :root { color-scheme: light; font-family: Arial, "Microsoft YaHei", sans-serif; }
    body { margin: 0; background: #f4f6f8; color: #17212b; }
    main { max-width: 1280px; margin: 0 auto; padding: 28px; }
    h1, h2 { margin-top: 0; }
    .card { background: white; border-radius: 10px; padding: 20px; margin: 0 0 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,.06); }
    .notice { border-left: 5px solid #d89b26; background: #fff8e7; }
    .status { display: inline-block; padding: 7px 13px; border-radius: 999px; font-weight: 700; }
    .PASS { color: #12633d; background: #dff4e8; }
    .WARNING { color: #7a4b00; background: #fff0c7; }
    .FAIL { color: #8c2020; background: #ffe0e0; }
    table { border-collapse: collapse; width: 100%; }
    th, td { text-align: left; vertical-align: top; padding: 10px; border-bottom: 1px solid #e6e9ed; }
    th { background: #f7f8fa; }
    .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 16px; }
    figure { margin: 0; background: white; border: 1px solid #e4e8ec; border-radius: 8px; overflow: hidden; }
    figure img { display: block; width: 100%; height: auto; }
    figcaption { padding: 9px 12px; color: #4a5662; }
    code { background: #eef1f4; padding: 2px 5px; border-radius: 4px; word-break: break-all; }
    .small { color: #52606d; font-size: .92rem; }
  </style>
</head>
<body><main>
  <section class="card">
    <h1>Step 3 - Hierarchical Vascular Representation</h1>
    <p><span class="status {{ acceptance.overall_status }}">{{ acceptance.overall_status }}</span></p>
    <p><b>Source run:</b> <code>{{ source_run }}</code></p>
    <p><b>Output run:</b> <code>{{ run_root }}</code></p>
    <p class="small">PASS passed automatically. WARNING needs visual review. FAIL should not be used downstream.</p>
  </section>

  <section class="card notice">
    <h2>Important scope</h2>
    <p>This is a whole-network coarse navigation graph. Raw centerline voxels are retained, but radius and curvature are coarse voxel-derived estimates. They are not final CFD geometry or final Geometry Generator training truth. Flow direction is unknown, so parent/daughter fields are not invented.</p>
  </section>

  <section class="card">
    <h2>Automatic checks</h2>
    <table><thead><tr><th>Status</th><th>Check</th><th>Result</th></tr></thead><tbody>
    {% for item in acceptance.checks %}
      <tr><td><span class="status {{ item.status }}">{{ item.status }}</span></td><td>{{ item.name }}</td><td>{{ item.message }}</td></tr>
    {% endfor %}
    </tbody></table>
  </section>

  <section class="card">
    <h2>Three representation scales</h2>
    <table><thead><tr><th>Scale</th><th>Question</th><th>Stored data</th></tr></thead><tbody>
      <tr><td>1 - Topology</td><td>What connects to what?</td><td>terminals, junctions, branches, cycles, branch-as-node relations</td></tr>
      <tr><td>2 - Branch morphology</td><td>What does a whole branch look like?</td><td>length, tortuosity, direction, coarse radius/taper, coarse curvature</td></tr>
      <tr><td>3 - Dense geometry</td><td>What happens along every position?</td><td>raw P(s), smoothed P(s), coarse r(s), local direction and curvature sequences</td></tr>
    </tbody></table>
  </section>

  <section class="card">
    <h2>Key measurements</h2>
    <table><tbody>
      <tr><th>Skeleton voxels</th><td>{{ report.skeleton_voxel_count }}</td><th>Represented voxels</th><td>{{ report.represented_voxel_count }}</td></tr>
      <tr><th>Missing / extra</th><td>{{ report.missing_voxel_count }} / {{ report.extra_voxel_count }}</td><th>Duplicated interiors</th><td>{{ report.duplicate_interior_voxel_count }}</td></tr>
      <tr><th>Nodes</th><td>{{ report.node_count }}</td><th>Branches</th><td>{{ report.branch_count }}</td></tr>
      <tr><th>Terminal nodes</th><td>{{ report.terminal_node_count }}</td><th>Junction nodes</th><td>{{ report.junction_node_count + report.complex_junction_node_count }}</td></tr>
      <tr><th>Cycle rank</th><td>{{ report.cycle_rank }}</td><th>Stored cycle basis</th><td>{{ report.cycle_basis_count }}</td></tr>
      <tr><th>Branch-as-node nodes</th><td>{{ report.branch_as_node_count }}</td><th>Relations</th><td>{{ report.branch_relation_count }}</td></tr>
      <tr><th>Median branch length</th><td>{{ "%.3f"|format(report.branch_length_um.median) }} um</td><th>Total raw branch length</th><td>{{ "%.3f"|format(report.branch_length_um.total) }} um</td></tr>
      <tr><th>Coordinate system</th><td>{{ report.coordinate_system }}</td><th>Voxel spacing</th><td>{{ report.spacing_um }} um</td></tr>
    </tbody></table>
  </section>

  <section class="card">
    <h2>Visual review</h2>
    <div class="grid">
    {% for image in images %}
      <figure><img src="{{ image.path }}" alt="{{ image.label }}"><figcaption>{{ image.label }}</figcaption></figure>
    {% endfor %}
    </div>
  </section>

  <section class="card">
    <h2>Files</h2>
    <ul>{% for file in files %}<li><a href="{{ file.path }}">{{ file.label }}</a></li>{% endfor %}</ul>
  </section>
</main></body></html>"""
)


def write_hierarchical_graph_html_report(
    path: Path,
    *,
    source_run: Path,
    run_root: Path,
    acceptance: AcceptanceResult,
    report: dict[str, Any],
    images: list[Path],
    files: list[Path],
) -> None:
    image_items = [
        {
            "path": item.relative_to(run_root).as_posix(),
            "label": item.stem.replace("_", " ").title(),
        }
        for item in images
        if item.is_file()
    ]
    file_items = [
        {
            "path": item.relative_to(run_root).as_posix(),
            "label": item.relative_to(run_root).as_posix(),
        }
        for item in files
        if item.is_file()
    ]
    path.write_text(
        _TEMPLATE.render(
            source_run=source_run,
            run_root=run_root,
            acceptance=acceptance,
            report=report,
            images=image_items,
            files=file_items,
        ),
        encoding="utf-8",
    )
