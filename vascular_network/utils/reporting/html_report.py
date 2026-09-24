"""Generate a standalone, browser-readable acceptance report."""

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
  <title>ULM vascular preprocessing acceptance report</title>
  <style>
    :root { color-scheme: light; font-family: Arial, "Microsoft YaHei", sans-serif; }
    body { margin: 0; background: #f4f6f8; color: #17212b; }
    main { max-width: 1240px; margin: 0 auto; padding: 28px; }
    h1, h2 { margin-top: 0; }
    .card { background: white; border-radius: 10px; padding: 20px; margin: 0 0 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,.06); }
    .status { display: inline-block; padding: 7px 13px; border-radius: 999px; font-weight: 700; }
    .PASS { color: #12633d; background: #dff4e8; }
    .WARNING { color: #7a4b00; background: #fff0c7; }
    .FAIL { color: #8c2020; background: #ffe0e0; }
    table { border-collapse: collapse; width: 100%; }
    th, td { text-align: left; vertical-align: top; padding: 10px; border-bottom: 1px solid #e6e9ed; }
    th { background: #f7f8fa; }
    .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 16px; }
    figure { margin: 0; background: white; border: 1px solid #e4e8ec; border-radius: 8px; overflow: hidden; }
    figure img { display: block; width: 100%; height: auto; }
    figcaption { padding: 9px 12px; color: #4a5662; }
    code { background: #eef1f4; padding: 2px 5px; border-radius: 4px; }
    .small { color: #52606d; font-size: .92rem; }
  </style>
</head>
<body><main>
  <section class="card">
    <h1>ULM 3-D vascular preprocessing</h1>
    <p><span class="status {{ acceptance.overall_status }}">{{ acceptance.overall_status }}</span></p>
    <p><b>Input:</b> <code>{{ input_stl }}</code></p>
    <p><b>Run directory:</b> <code>{{ run_root }}</code></p>
    <p class="small">PASS means the automatic check passed. WARNING requires visual review. FAIL means the result should not be used downstream.</p>
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
    <h2>Key measurements</h2>
    <table><tbody>
      <tr><th>Input triangles</th><td>{{ cleanup.input_quality.triangle_count }}</td><th>Cleaned triangles</th><td>{{ cleanup.cleaned_quality.triangle_count }}</td></tr>
      <tr><th>Input components</th><td>{{ cleanup.input_quality.connected_component_count }}</td><th>Final STL components</th><td>{{ cleanup.cleaned_quality.connected_component_count }}</td></tr>
      <tr><th>Main component ID</th><td>{{ cleanup.main_network_component_id }}</td><th>Main surface fraction</th><td>{{ "%.3f%%"|format(cleanup.main_network_surface_area_fraction * 100) }}</td></tr>
      <tr><th>Small fragments removed</th><td>{{ cleanup.small_fragment_count }}</td><th>Island networks removed</th><td>{{ cleanup.island_network_count }}</td></tr>
      <tr><th>Total removed components</th><td>{{ cleanup.removed_component_count }}</td><th>Removed surface fraction</th><td>{{ "%.3f%%"|format(cleanup.removed_surface_area_fraction * 100) }}</td></tr>
      <tr><th>Voxel dimensions</th><td>{{ voxel.dimensions_xyz }}</td><th>Voxel spacing</th><td>{{ voxel.spacing_um }} um</td></tr>
      <tr><th>Voxel components before</th><td>{{ voxel.initial_connected_component_count }}</td><th>Voxel components after</th><td>{{ voxel.connected_component_count }}</td></tr>
      <tr><th>Voxel islands removed</th><td>{{ voxel.removed_island_voxel_count }}</td><th>Voxel removal fraction</th><td>{{ "%.3f%%"|format(voxel.removed_island_fraction * 100) }}</td></tr>
      <tr><th>Final mask voxels</th><td>{{ voxel.foreground_voxel_count }}</td><th>Skeleton components</th><td>{{ skeleton.connected_component_count }}</td></tr>
      <tr><th>Skeleton voxels</th><td>{{ skeleton.skeleton_voxel_count }}</td><th>Skeleton outside mask</th><td>{{ skeleton.voxels_outside_mask }}</td></tr>
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


def write_html_report(
    path: Path,
    *,
    input_stl: Path,
    run_root: Path,
    acceptance: AcceptanceResult,
    cleanup_summary: dict[str, Any],
    voxel_report: dict[str, Any],
    skeleton_report: dict[str, Any],
    images: list[Path],
    files: list[Path],
) -> None:
    image_items = [
        {"path": item.relative_to(run_root).as_posix(), "label": item.stem.replace("_", " ").title()}
        for item in images
        if item.is_file()
    ]
    file_items = [
        {"path": item.relative_to(run_root).as_posix(), "label": item.relative_to(run_root).as_posix()}
        for item in files
        if item.is_file()
    ]
    path.write_text(
        _TEMPLATE.render(
            input_stl=input_stl,
            run_root=run_root,
            acceptance=acceptance,
            cleanup=cleanup_summary,
            voxel=voxel_report,
            skeleton=skeleton_report,
            images=image_items,
            files=file_items,
        ),
        encoding="utf-8",
    )
