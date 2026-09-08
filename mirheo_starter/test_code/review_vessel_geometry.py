"""Offline review of the exported-and-reloaded geometry, using formal APIs."""

import argparse
from html import escape
from html.parser import HTMLParser
import json
from pathlib import Path
import sys

# Direct file execution puts test_code, rather than the project root, on sys.path.
if __name__ == "__main__" and __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

TITLE = "本次仅核查血管几何和标签，尚未生成 SDF 或运行流体。"


def _json(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")


def source_package_comparison(geometry, directory: Path) -> dict:
    """Independently reread the byte-exact source VTP copy for the review table."""
    import numpy as np
    from py_scripts.vessel_geometry.io import read_tagged_vtp
    from py_scripts.vessel_geometry.model import GeometryError
    points, faces, labels, _ = read_tagged_vtp(directory / "source_surface.vtp")
    equal = {"points": np.array_equal(points, geometry.points_source),
             "triangles": np.array_equal(faces, geometry.triangles),
             "labels": np.array_equal(labels, geometry.entity_ids)}
    if not all(equal.values()):
        raise GeometryError("REVIEW_SOURCE_MISMATCH", f"原始 VTP 副本与重新载入的数据包不一致：{equal}")
    meters = points.astype(float) * geometry.to_meter
    return {"status": "PASS", "method": "独立读取 source_surface.vtp，与 NPZ 回读数组逐项比较",
            "point_count": len(points), "triangle_count": len(faces),
            "entity_ids": [int(x) for x in np.unique(labels)],
            "extent_m": np.ptp(meters, axis=0).tolist(),
            "referenced_extent_m": np.ptp(meters[np.unique(faces)], axis=0).tolist(),
            "exact_comparisons": equal}


def make_figure(geometry):
    import numpy as np
    import plotly.graph_objects as go
    from plotly.colors import qualitative
    figure = go.Figure()
    controls = {"walls": [], "ports": {}}
    points = geometry.points_m * 1e6
    arrow_length = float(np.max(np.ptp(points[np.unique(geometry.triangles)], axis=0))) * 0.035
    for order, patch in enumerate(geometry.patches):
        faces = geometry.triangles[patch.face_ids]
        used, local = np.unique(faces.ravel(), return_inverse=True)
        local = local.reshape(-1, 3)
        xyz = points[used]
        color = "#91a4b7" if patch.kind == "wall" else qualitative.Dark24[(order - 1) % len(qualitative.Dark24)]
        index = len(figure.data)
        group = str(patch.entity_id)
        if patch.kind == "wall":
            controls["walls"].append(index)
        else:
            controls["ports"][str(patch.entity_id)] = [index]
        figure.add_trace(go.Mesh3d(x=xyz[:, 0].tolist(), y=xyz[:, 1].tolist(), z=xyz[:, 2].tolist(),
            i=local[:, 0].tolist(), j=local[:, 1].tolist(), k=local[:, 2].tolist(),
            customdata=used.tolist(), name=patch.display_name, legendgroup=group, showlegend=True,
            color=color, opacity=0.35 if patch.kind == "wall" else 1.0,
            flatshading=patch.kind != "wall", meta={"entity_id": patch.entity_id, "source_face_ids": patch.face_ids},
            hovertemplate="原始点编号 %{customdata}<br>x=%{x:.4f} µm<br>y=%{y:.4f} µm<br>z=%{z:.4f} µm<extra>%{fullData.name}</extra>"))
        if patch.kind == "wall":
            continue
        center = np.asarray(patch.center_m) * 1e6
        controls["ports"][str(patch.entity_id)].append(len(figure.data))
        figure.add_trace(go.Scatter3d(x=[float(center[0])], y=[float(center[1])], z=[float(center[2])],
            mode="markers+text", text=[f"{patch.display_name} / ID {patch.entity_id}"], textposition="top center",
            marker={"size": 5, "color": color}, name=f"{patch.display_name} 面积加权中心", legendgroup=group, showlegend=False))
        if patch.outward_normal is not None:
            vector = np.asarray(patch.outward_normal)
            tip = center + arrow_length * vector
            controls["ports"][str(patch.entity_id)].append(len(figure.data))
            figure.add_trace(go.Scatter3d(x=[float(center[0]), float(tip[0])], y=[float(center[1]), float(tip[1])],
                z=[float(center[2]), float(tip[2])], mode="lines", line={"width": 5, "color": color},
                legendgroup=group, showlegend=False, name="几何法向（非流速）"))
            controls["ports"][str(patch.entity_id)].append(len(figure.data))
            figure.add_trace(go.Cone(x=[float(tip[0])], y=[float(tip[1])], z=[float(tip[2])],
                u=[float(vector[0])], v=[float(vector[1])], w=[float(vector[2])],
                sizemode="absolute", sizeref=arrow_length * 0.3, anchor="tip", showscale=False,
                colorscale=[[0, color], [1, color]], legendgroup=group, showlegend=False, name="几何法向（非流速）"))
    figure.update_layout(template="plotly_white", margin={"l": 0, "r": 0, "t": 10, "b": 0},
        scene={"aspectmode": "data", "xaxis_title": "x / µm", "yaxis_title": "y / µm", "zaxis_title": "z / µm",
               "camera": {"eye": {"x": 1.4, "y": 1.5, "z": 1.1}}},
        legend={"orientation": "h", "y": -0.03, "groupclick": "togglegroup"}, height=740)
    # Display-only local indices above do not change points, faces, or numbering in the package.
    return json.loads(figure.to_json()), controls


def write_review_html(geometry, manifest, destination: Path, report: dict) -> None:
    from plotly.offline import get_plotlyjs
    figure, controls = make_figure(geometry)
    rows = []
    options = ['<option value="">全部端口</option>']
    for patch in geometry.patches:
        if patch.kind == "wall":
            continue
        options.append(f'<option value="{patch.entity_id}">{escape(patch.display_name)} · entity {patch.entity_id}</option>')
        center = ", ".join(f"{x * 1e6:.6f}" for x in patch.center_m)
        normal = ", ".join(f"{x:.6f}" for x in patch.outward_normal) if patch.outward_normal is not None else "待核查"
        columns = [patch.display_name, patch.port_id, str(patch.entity_id), str(patch.boundary_index),
                   patch.original_role, patch.boundary_origin, str(patch.triangle_count),
                   f"{patch.area_m2 * 1e12:.9f}", center, normal + "；" + patch.normal_status]
        rows.append(f'<tr data-entity="{patch.entity_id}">' + "".join(f"<td>{escape(str(value))}</td>" for value in columns) + "</tr>")
    summary = geometry.summary()
    source = report.get("source_vs_package") or source_package_comparison(geometry, Path(report["package_path"]))
    source_rows = [
        ("点数（含未引用点）", source["point_count"], summary["point_count"]),
        ("三角面数", source["triangle_count"], summary["triangle_count"]),
        ("原始标签集合", source["entity_ids"], summary["entity_ids"]),
        ("全部点尺寸 / µm", [float(x * 1e6) for x in source["extent_m"]], [float(x * 1e6) for x in summary["extent_m"]]),
        ("实际三角表面尺寸 / µm", [float(x * 1e6) for x in source["referenced_extent_m"]], [float(x * 1e6) for x in summary["referenced_extent_m"]]),
        ("坐标单位", f"原始 VTP: {geometry.source_length_unit}；已验收数组: m", "m；显示时换算为 µm，未二次换算正式数组"),
        ("坐标关系、连接、编号、标签和分区", "原始输入", geometry.checks["export_roundtrip"]["status"] + "，逐项精确一致"),
    ]
    comparison = "".join("<tr>" + "".join(f"<td>{escape(str(x))}</td>" for x in row) + "</tr>" for row in source_rows)
    hash_rows = "".join(f"<tr><td>{escape(path)}</td><td><code>{sha}</code></td></tr>" for path, sha in manifest["input_sha256"].items())
    user_checks = ["血管形状是否与旧工程一致？", "各个开口是否位于预期位置？", "入口出口身份是否正确？",
                   "管壁有没有被错误标成开口？", "单位与尺寸是否合理？", "模型是否被压扁、拉伸或意外平移？"]
    checklist = "".join(f"<label class=check><input type=checkbox> {text}</label>" for text in user_checks)
    topology = geometry.checks["original_topology"]
    diagnostic_rows = [
        ("自动数据迁移", geometry.checks["migration_status"]),
        ("完整表面基础诊断（本次复算）", topology["complete_surface_status"]),
        ("来源模型历史人工验收", report["states"].get("source_model_historical_acceptance", "SYNTHETIC_TEST_ONLY")),
        ("完整表面边界边 / 非流形边", f"{topology['boundary_edge_count']} / {topology['nonmanifold_edge_count']}"),
        ("重复 / 退化三角形", f"{topology['duplicate_triangle_count']} / {topology['degenerate_triangle_count']}"),
        ("通过共享边连接的面分量", topology["edge_connected_face_components"]),
        ("共享边绕向一致", str(topology["winding_consistent"])),
        ("输入哈希前后对照", geometry.checks["input_hashes_unchanged"]["status"]),
        ("原始 VTP 副本与 NPZ 独立对照", source["status"]),
        ("自交复杂检查", "NOT_CHECKED；历史结果未冒充本次检查"),
    ]
    diagnostic_table = "<table><tbody>" + "".join(f"<tr><td>{escape(str(a))}</td><td>{escape(str(b))}</td></tr>" for a, b in diagnostic_rows) + "</tbody></table>"
    payload = {"summary": summary, "figure": figure, "controls": controls,
               "package_path": report["package_path"], "source_tagged_surface": manifest["tagged_surface"],
               "states": report["states"], "all_source_face_ids_rendered": True,
               "visualization_source": "exported package loaded with verified hashes; complete triangle mesh"}
    document = """<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>__TITLE__</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#f3f6fa;color:#183044;font:15px/1.7 system-ui,"Microsoft YaHei",sans-serif}
header{padding:30px max(3vw,20px);background:#102c40;color:white}h1{font-size:26px;margin:0 0 8px}h2{font-size:21px;margin:0 0 14px}
main{max-width:1480px;margin:22px auto;padding:0 18px}section{background:white;border:1px solid #dce5ed;border-radius:12px;padding:22px;margin-bottom:20px}
.states{display:flex;flex-wrap:wrap;gap:10px;margin-top:16px}.state{background:#edf5fc;color:#17354b;padding:8px 14px;border-radius:7px}.pass{background:#d9f5e4;color:#075632}
.controls{display:flex;flex-wrap:wrap;gap:20px;align-items:center;background:#f6f8fa;padding:12px}.note{color:#526779}.attention{border-left:4px solid #d59124;padding:10px 16px;background:#fff7e5}
select,button{font:inherit;padding:6px;max-width:500px}#vessel-view{width:100%;min-height:740px}#render-status{font-weight:600}
.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #dce5ed;padding:9px;vertical-align:top;min-width:90px}th{background:#eef3f8;text-align:left}
td{overflow-wrap:anywhere}tr.selected{background:#fff0c2}code{font-family:monospace;overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f6f8fa;padding:14px;font-size:12px}
.check{display:block;padding:7px}summary{cursor:pointer;font-weight:600}a{color:#14699c}.stats{display:flex;gap:32px;flex-wrap:wrap;font-size:17px}
</style></head><body><header><h1>__TITLE__</h1>
<div>已有形状 · 原始三角形与身份保留 · 正式数据以米保存 · 页面以微米显示</div>
<div class="states"><span class="state pass">自动迁移检查：PASS</span><span class="state">可视化文件：GENERATED</span>
<span class="state">浏览器交互测试：__BROWSER_STATUS__</span><span class="state">本次人工核查：PENDING</span>
<span class="state">来源模型历史人工验收：__SOURCE_ACCEPTANCE__</span>
<span class="state">SDF、粒子、数值边界、真实血流：NOT_STARTED</span></div></header>
<main><section><h2>1. 查看实际导入的三维血管</h2><div class=stats>__STATS__</div>
<p class=note>显示完整三角表面；数据来自导出后重新载入的数据包。拖动旋转，滚轮缩放；图例点击可隐藏分区。</p>
<div class=controls><label><input id=wall-visible type=checkbox checked> 显示管壁</label>
<label>壁面透明度 <input id=wall-opacity type=range min=0 max=1 step=0.05 value=0.35></label>
<label>突出端口 <select id=port-select>__OPTIONS__</select></label><button id=reset-view>恢复视角</button></div>
<p id=render-status>当前浏览器渲染：等待加载（不代表人工核查通过）</p><div id=vessel-view></div>
<p class=attention>ASSUMED_INLET / ASSUMED_OUTLET 保留“假定”身份，不是已证实的生理入口出口。箭头仅是基于闭合、单分量、绕向一致表面的几何外向约定，不是流速。自交检查本次为 NOT_CHECKED。</p>
</section><section><h2>2. 端口身份和实际几何</h2><div class=scroll><table><thead><tr>
<th>显示名</th><th>port_id</th><th>原始 entity_id</th><th>boundary_index</th><th>原始 role</th><th>boundary_origin</th><th>三角面数</th><th>面积 / µm²</th><th>面积加权中心 (x,y,z) / µm</th><th>几何单位法向与状态</th>
</tr></thead><tbody>__PORT_ROWS__</tbody></table></div></section>
<section><h2>3. 源数据与迁移数据对照</h2><p><b>正式数据包：</b><code>__PACKAGE__</code></p>
<table><thead><tr><th>项目</th><th>原始数据</th><th>迁移数据</th></tr></thead><tbody>__COMPARISON__</tbody></table>
<p class=attention>保留全部 __UNUSED__ 个未引用点；“全部点尺寸”与“实际三角表面尺寸”分别列出。原始点未删除。旧 QC 的中心按其原统计定义复核，与面积加权中心的差值单独记录，不混用定义。</p>
<details><summary>输入文件 SHA-256（处理前后相同）</summary><div class=scroll><table><tbody>__HASHES__</tbody></table></div></details>
</section><section><h2>4. 已执行检查与未执行项目</h2>
<p>数据迁移检查与几何诊断分别报告。单独 wall 分区的端口开边是正常现象；本次没有补洞、清理、合并或翻转正式网格。</p>
__DIAGNOSTICS__<details><summary>展开完整记录：检查容差、实际误差、单位证据与前后哈希</summary><pre>__CHECKS__</pre></details>
<p><b>本次用户人工核查仍为 PENDING。</b> __BROWSER_NOTE__ 下方勾选仅辅助观察，不改变报告中的 PENDING。</p>
</section><section><h2>5. 请逐项人工核查</h2>__CHECKLIST__
<p>本次仅迁移和核查已验收血管几何与标签。尚未生成 SDF、粒子或数值边界，尚未运行真实血流。</p>
</section></main><script id="geometry-audit-data" type="application/json">__PAYLOAD__</script>
<script>__PLOTLY__</script><script>
const audit=JSON.parse(document.getElementById('geometry-audit-data').textContent);
const view=document.getElementById('vessel-view');const status=document.getElementById('render-status');
Plotly.newPlot(view,audit.figure.data,audit.figure.layout,{responsive:true,displaylogo:false}).then(()=>{
 status.textContent='当前浏览器：三维场景已加载。人工核查仍为 PENDING。';
 const walls=audit.controls.walls;
 document.getElementById('wall-visible').addEventListener('change',e=>Plotly.restyle(view,{visible:e.target.checked},walls));
 document.getElementById('wall-opacity').addEventListener('input',e=>Plotly.restyle(view,{opacity:Number(e.target.value)},walls));
 document.getElementById('port-select').addEventListener('change',e=>{
   for(const [id,indices] of Object.entries(audit.controls.ports)){Plotly.restyle(view,{opacity:(!e.target.value||id===e.target.value)?1:0.12},indices);}
   for(const row of document.querySelectorAll('tr[data-entity]')){row.classList.toggle('selected',row.dataset.entity===e.target.value);}
 });
 document.getElementById('reset-view').addEventListener('click',()=>Plotly.relayout(view,{'scene.camera':{eye:{x:1.4,y:1.5,z:1.1}}}));
}).catch(error=>{status.textContent='三维渲染失败：'+String(error)+'。请使用支持 WebGL 的 Edge / Chrome，并保留此错误信息。';status.style.color='#b91c1c';});
</script></body></html>"""
    substitutions = {"TITLE": escape(TITLE), "STATS": f"<span>{summary['point_count']:,} 个原始点</span><span>{summary['triangle_count']:,} 个三角面</span><span>{len(options)-1} 个端口</span>",
        "OPTIONS": "".join(options), "PORT_ROWS": "".join(rows), "PACKAGE": escape(report["package_path"]),
        "COMPARISON": comparison, "UNUSED": str(summary["unused_point_count"]), "HASHES": hash_rows,
        "CHECKS": escape(json.dumps(geometry.checks, indent=2, ensure_ascii=False, allow_nan=False)),
        "SOURCE_ACCEPTANCE": escape(report["states"].get("source_model_historical_acceptance", "SYNTHETIC_TEST_ONLY")),
        "BROWSER_STATUS": escape(report["states"].get("browser_verification", "NOT_TESTED")),
        "BROWSER_NOTE": escape(report.get("browser_note", "浏览器交互未执行，NOT_TESTED。")),
        "CHECKLIST": checklist, "DIAGNOSTICS": diagnostic_table, "PAYLOAD": _json(payload), "PLOTLY": get_plotlyjs()}
    for key, value in substitutions.items():
        document = document.replace("__" + key + "__", value)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(document)


def check_html(path: Path, geometry) -> dict:
    import numpy as np
    class AuditParser(HTMLParser):
        def __init__(self):
            super().__init__(); self.collect = False; self.payload = []; self.remote_scripts = []; self.ids = set()
        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if "id" in attrs:
                self.ids.add(attrs["id"])
            if tag == "script":
                if "src" in attrs:
                    self.remote_scripts.append(attrs["src"])
                self.collect = attrs.get("id") == "geometry-audit-data"
        def handle_data(self, data):
            if self.collect:
                self.payload.append(data)
        def handle_endtag(self, tag):
            if tag == "script":
                self.collect = False
    parser = AuditParser()
    document = path.read_text(encoding="utf-8")
    parser.feed(document)
    payload = json.loads("".join(parser.payload))
    faces = []
    coordinates_match, connectivity_match, labels_match = True, True, True
    for trace in payload["figure"]["data"]:
        if trace["type"] == "mesh3d":
            assert len(trace["i"]) == len(trace["j"]) == len(trace["k"]) == len(trace["meta"]["source_face_ids"])
            face_ids = trace["meta"]["source_face_ids"]
            point_ids = np.asarray(trace["customdata"])
            coordinates_match &= np.array_equal(np.column_stack([trace[key] for key in ("x", "y", "z")]), geometry.points_m[point_ids] * 1e6)
            connectivity_match &= np.array_equal(point_ids[np.column_stack([trace[key] for key in ("i", "j", "k")])], geometry.triangles[face_ids])
            labels_match &= bool(np.all(geometry.entity_ids[face_ids] == trace["meta"]["entity_id"]))
            faces.extend(face_ids)
    checks = {"title_present": TITLE in document, "no_external_script_resources": not parser.remote_scripts,
        "all_original_triangles_once": sorted(faces) == list(range(len(geometry.triangles))),
        "display_coordinates_equal_package_um": bool(coordinates_match),
        "display_connectivity_equal_package": bool(connectivity_match),
        "display_labels_equal_package": bool(labels_match),
        "controls_present": {"wall-visible", "wall-opacity", "port-select", "reset-view"} <= parser.ids,
        "true_aspect_ratio": payload["figure"]["layout"]["scene"]["aspectmode"] == "data",
        "manual_review_pending": payload["states"]["human_review"] == "PENDING",
        "exported_data_source": payload["visualization_source"].startswith("exported package"),
        "inline_plotly_present": "plotly.js" in document and "Plotly.newPlot" in document}
    if not all(checks.values()):
        raise ValueError(f"HTML 结构或完整性检查失败：{checks}")
    return {"status": "PASS", "checks": checks, "rendered_triangle_count": len(faces),
            "browser_execution": "NOT_TESTED", "html_size_bytes": path.stat().st_size}


def verify_browser(html_path, output_dir):
    import subprocess
    from py_scripts.vessel_geometry.io import read_json
    script = Path(__file__).with_name("check_vessel_geometry_browser.cjs")
    arguments = [subprocess.check_output(["wslpath", "-w", str(p)], text=True).strip()
                 for p in (script, html_path, output_dir)]
    result = subprocess.run(["/mnt/d/Program Files/nodejs/node.exe", *arguments],
                            capture_output=True, text=True, timeout=180)
    report_path = output_dir / "browser_checks.json"
    if not report_path.is_file():
        raise RuntimeError(f"浏览器检查未产生报告：{result.stdout}\n{result.stderr}")
    report = read_json(report_path)
    if result.returncode or report["status"] != "PASS":
        raise RuntimeError(f"浏览器检查失败，详见 {report_path}：{report.get('error')}")
    return report


def ensure_review(config_path, *, progress=print, browser_test=False):
    from importlib.metadata import version
    from py_scripts.vessel_geometry.export import load_package, run_import
    from py_scripts.vessel_geometry.io import load_config, read_json, safe_output, sha256_file, write_json
    config = load_config(config_path)
    result = run_import(config_path, progress=progress)
    geometry = load_package(result.package_path)
    report_dir = safe_output(config["review_output_root"] / result.run_id, config["source_project"])
    report_path = report_dir / "review_report.json"
    html_path = report_dir / "geometry_review.html"
    if report_path.exists():
        report = read_json(report_path)
        if (report.get("review_code_sha256") == sha256_file(Path(__file__))
                and report.get("html_sha256") == sha256_file(html_path)
                and report.get("package_path") == str(result.package_path)):
            check_html(html_path, geometry)
            if progress is not None:
                progress(f"[复用] 已有完整网格核查页面，文件保持不变：{html_path}")
            return result, report
        raise ValueError(f"既有核查页不匹配，未覆盖：{report_dir}")
    report_dir.mkdir(parents=True, exist_ok=False)
    report = {"run_id": result.run_id, "package_path": str(result.package_path),
        "html_path": str(html_path), "command_argv": sys.argv,
        "states": {**result.manifest.get("states", {}),
                   "automatic_data_check": geometry.checks["migration_status"],
                   "basic_geometry_diagnostics": geometry.checks["original_topology"]["complete_surface_status"],
                   "visualization_file": "GENERATED", "browser_verification": "NOT_TESTED", "human_review": "PENDING",
                   "sdf": "NOT_STARTED", "particles": "NOT_STARTED", "lbm_mesh": "NOT_STARTED",
                   "numerical_boundary_conditions": "NOT_STARTED", "real_blood_flow": "NOT_STARTED",
                   "mirheo_integration": "NOT_STARTED"},
        "browser_note": "仅 HTML 生成不能证明浏览器交互成功；没有执行记录时 NOT_TESTED。人工核查始终 PENDING。",
        "summary": geometry.summary(), "import_manifest_sha256": sha256_file(result.package_path / "import_manifest.json"),
        "plotly_version": version("plotly"), "geometry_checks": geometry.checks}
    if (result.package_path / "migration_manifest.json").exists():
        report["migration_manifest_sha256"] = sha256_file(result.package_path / "migration_manifest.json")
    report["source_vs_package"] = source_package_comparison(geometry, result.package_path)
    report["review_code_sha256"] = sha256_file(Path(__file__))
    if browser_test:
        # The candidate uses this same target package, template, JavaScript and
        # full mesh. Only the evidence/status text changes in the final page.
        preview = report_dir / "geometry_review_preview.html"
        write_review_html(geometry, result.manifest, preview, report)
        check_html(preview, geometry)
        report["browser_preview"] = verify_browser(preview, report_dir / "browser_preview")
        report["states"]["browser_verification"] = "PASS"
        report["browser_note"] = "Windows Chrome 软件渲染已实际检查完整网格、管壁隐藏与透明度、端口突出、鼠标旋转、滚轮缩放和恢复视角；无远程请求。最终 HTML 的执行记录及精确 SHA-256 见同目录 browser_final/browser_checks.json。此结果不代表用户人工验收。"
    write_review_html(geometry, result.manifest, html_path, report)
    report["html_structure_checks"] = check_html(html_path, geometry)
    report["html_sha256"] = sha256_file(html_path)
    if browser_test:
        report["browser_final"] = verify_browser(html_path, report_dir / "browser_final")
    write_json(report_path, report)
    with (report_dir / "review_report.md").open("x", encoding="utf-8") as stream:
        stream.write(f"# {TITLE}\n\n正式数据包：`{result.package_path}`\n\nHTML：`{html_path}`\n\n"
            f"{len(geometry.points_m)} 个原始点，{len(geometry.triangles)} 个三角面。完整显示全部三角面，无降采样。\n\n"
            "本次迁移自动检查 PASS；正式数组、身份与字节副本精确一致。源文件前后哈希不变。\n\n"
            "低成本拓扑与端口统计已复算；复杂自交 NOT_CHECKED，wall 子集开边保留，ASSUMED 身份未改。\n\n"
            f"HTML GENERATED；结构、网格与离线资源检查 PASS；浏览器交互 {report['states']['browser_verification']}；本次人工核查 PENDING。\n\n"
            "历史记录保持原字节；后续来源人工验收见迁移清单。SDF、粒子、数值边界、真实血流 NOT_STARTED。\n")
    return result, report


def open_windows_file(path):
    import subprocess
    # The caller explicitly requested --open. No wsl.exe, shell or server.
    converted = subprocess.run(["wslpath", "-w", str(path)], check=True, capture_output=True, text=True, timeout=10)
    subprocess.run(["explorer.exe", converted.stdout.strip()], check=False, timeout=15)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="一键核查已验收血管迁移；匹配包与离线 HTML 直接复用，仅 CPU 几何。")
    parser.add_argument("--config", type=Path, default=Path("py_scripts/vessel_geometry_import.yaml"), help="YAML 相对路径及其中相对路径均以目标项目根目录为基准")
    parser.add_argument("--open", action="store_true", help="完成后用 Windows 默认浏览器打开 HTML；不启动服务器")
    parser.add_argument("--browser-test", action="store_true", help="首次生成时额外用本地 Windows Chrome 和 Node 执行有限离线交互测试并截图")
    args = parser.parse_args(argv)
    try:
        result, report = ensure_review(args.config, progress=lambda text: print(text, flush=True), browser_test=args.browser_test)
        print(f"\n正式数据包：{result.package_path}\n离线三维 HTML：{report['html_path']}\n"
              f"迁移自动检查：PASS；HTML：GENERATED；浏览器交互：{report['states']['browser_verification']}；人工核查：PENDING", flush=True)
        if args.open:
            open_windows_file(Path(report["html_path"]))
        return 0
    except Exception as exc:
        print(f"[{getattr(exc, 'status', 'FAIL')}] {exc}", file=sys.stderr)
        if getattr(exc, "details", {}).get("package_path"):
            print(f"诊断目录：{exc.details['package_path']}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
