#!/usr/bin/env python3
"""只读查看 hello_rbc_flow.py 的结果：真实 PLY 动画 + CSV 位置/速度图。

安装查看依赖（不需要重新编译 Mirheo）：
    python -m pip install "plotly==6.5.2" "trimesh==4.11.1"

使用：
    python py_scripts/view_rbc_flow.py --open
    python py_scripts/view_rbc_flow.py --run /absolute/path/to/run

本程序不 import mirheo，不启动 MPI，不推进仿真，不修改原始输出。
生成一个包含 JavaScript 的离线 HTML。不读取/渲染 HDF5 体速度场。
PLY 文件编号只表示保存帧，不擅自当作物理时间；CSV 时间单独显示。
适用于单个 RBC 的上述 Hello World；多对象或其他文件格式需另适配。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric_key(path: Path) -> tuple[int, str]:
    match = re.search(r"_(\d+)\.ply$", path.name, re.I)
    if match is None:
        raise ValueError(f"不认识 PLY 帧命名：{path.name}；预期 rbc_数字.ply")
    return int(match.group(1)), path.name


def select_run(requested: Path | None, project: Path) -> Path:
    if requested is not None:
        chosen = requested.expanduser().resolve()
        if not chosen.is_dir():
            raise ValueError(f"结果目录不存在：{chosen}")
        return chosen
    root = project / "runs"
    candidates = sorted(
        (p for p in root.iterdir()
         if p.is_dir() and (p / "rbc_track.csv").is_file()
         and any((p / "ply").glob("rbc_*.ply"))),
        key=lambda p: p.name,
    ) if root.is_dir() else []
    if not candidates:
        raise ValueError(
            f"在 {root} 没找到同时有 rbc_track.csv 和 ply/rbc_*.ply 的结果。\n"
            "请先确认仿真输出，或用 --run 指定实际运行目录。"
        )
    print("找到以下结果；请选择要查看的那一次，程序不会替你选择最新目录：")
    for index, candidate in enumerate(candidates, 1):
        print(f"  [{index}] {candidate.name}")
    if not sys.stdin.isatty():
        raise ValueError("非交互运行请用 --run 明确指定结果目录。")
    value = input("请输入方括号中的编号，然后按 Enter：").strip()
    if not value.isdigit() or not 1 <= int(value) <= len(candidates):
        raise ValueError("编号无效，请重新运行并输入列表中的数字。")
    return candidates[int(value) - 1].resolve()


def read_mesh(path: Path, np: Any) -> tuple[Any, Any]:
    # 使用成熟的 PLY 读取器；直接使用数组，不做 clean/平滑/重网格化。
    from trimesh.exchange.ply import load_ply
    with path.open("rb") as stream:
        data = load_ply(stream, fix_texture=False, skip_materials=True)
    points = np.asarray(data.get("vertices"), dtype=np.float64)
    faces = np.asarray(data.get("faces"))
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError(f"顶点数据无效：{path}")
    if not np.isfinite(points).all():
        raise ValueError(f"顶点含 NaN/Inf：{path}。不隐藏异常，请检查仿真。")
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0:
        raise ValueError(f"需要非空三角面 PLY：{path}")
    if not np.issubdtype(faces.dtype, np.integer):
        raise ValueError(f"面索引不是整数：{path}")
    if int(faces.min()) < 0 or int(faces.max()) >= len(points):
        raise ValueError(f"面索引越界：{path}")
    return points, faces.astype(np.int64, copy=False)


def read_track(path: Path, np: Any) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"轨迹 CSV 没有表头：{path}")
        reader.fieldnames = [key.strip().lstrip("#").strip() for key in reader.fieldnames]
        required = ["objId", "time", "comx", "comy", "comz", "vx", "vy", "vz"]
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"轨迹 CSV 缺少列：{sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("rbc_track.csv 只有表头，没有数据。")
    identifiers = {row["objId"].strip() for row in rows}
    if len(identifiers) != 1:
        raise ValueError("检测到多个对象；本查看器仅用于单 RBC Hello World。")
    columns = {key: np.asarray([float(row[key]) for row in rows], dtype=float)
               for key in required if key != "objId"}
    if any(not np.isfinite(arr).all() for arr in columns.values()):
        raise ValueError("轨迹或速度存在 NaN/Inf；不删除异常行伪装成正常结果。")
    if len(rows) > 1 and not np.all(np.diff(columns["time"]) > 0):
        raise ValueError("CSV 时间没有严格递增，可能是不同运行混在同一个文件中。")
    return {"id": next(iter(identifiers)), "nrows": len(rows), **columns}


def split_at_jumps(x: Any, ys: Any, jumps: Any) -> tuple[list[Any], list[list[Any]]]:
    """只断开绘图线段，不移动/解包裹原数据；避免画穿过周期盒的假直线。"""
    xx: list[Any] = []
    out: list[list[Any]] = [[] for _ in ys]
    for k in range(len(x)):
        if k and bool(jumps[k - 1]):
            xx.append(None)
            for values in out:
                values.append(None)
        xx.append(float(x[k]))
        for values, source in zip(out, ys):
            values.append(float(source[k]))
    return xx, out


def box_edges(domain: Any) -> tuple[list[Any], list[Any], list[Any]]:
    corners = [(i * domain[0], j * domain[1], k * domain[2])
               for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    components: list[list[Any]] = [[], [], []]
    for i, p in enumerate(corners):
        for j in range(i + 1, len(corners)):
            q = corners[j]
            if sum(p[d] != q[d] for d in range(3)) == 1:
                for d in range(3):
                    components[d].extend([float(p[d]), float(q[d]), None])
    return components[0], components[1], components[2]


def build_report(run: Path, max_frames: int, frame_ms: int) -> Path:
    try:
        import numpy as np
        import plotly
        import plotly.graph_objects as go
        import plotly.io as pio
        import trimesh
    except ImportError as exc:
        raise ValueError(
            '缺少查看依赖。请在当前 .venv 执行：\n'
            'python -m pip install "plotly==6.5.2" "trimesh==4.11.1"'
        ) from exc

    files = sorted((run / "ply").glob("rbc_*.ply"), key=numeric_key)
    if not files:
        raise ValueError(f"未找到 {run}/ply/rbc_*.ply。请检查运行是否完成并正确输出。")
    file_ids = [numeric_key(path)[0] for path in files]
    if len(set(file_ids)) != len(file_ids):
        raise ValueError("PLY 存在重复的保存编号，不能建立唯一播放顺序。")
    indices = np.unique(np.linspace(0, len(files) - 1, min(len(files), max_frames), dtype=int))
    selected = [files[int(i)] for i in indices]
    messages = [
        "仅播放保存的真实 PLY 帧，不生成中间形变；播放快慢不是仿真性能。",
        "所有坐标、速度和 CSV 时间均按原示例模拟单位显示，未经 SI 标定。",
        "PLY 的文件编号不是时间戳。本页不假设 PLY 与 CSV 逐行精确同步。",
        "周期盒不是血管壁；本 Hello World 没有真实血管、微泡或入口出口。",
        "本页未读取或渲染 HDF5 体速度场。轨迹曲线不是液体流线。",
        "这是可视化读取检查，不代表膜力学、耦合或物理精度验证通过。",
    ]
    if len(selected) < len(files):
        messages.append(f"为控制页面大小，仅按确定性间隔显示 {len(selected)}/{len(files)} 个文件；"
                        "每个被显示帧的全部三角面保留。CSV 曲线不抽样。")
    if len(files) == 1:
        messages.append("只有一帧 PLY，只能查看形状，不能显示随时间的运动。")

    info_path = run / "case_info.json"
    info: dict[str, Any] = {}
    domain = None
    if info_path.is_file():
        info = json.loads(info_path.read_text(encoding="utf-8"))
        domain = np.asarray(info.get("domain"), dtype=float)
        if domain.shape != (3,) or not np.isfinite(domain).all() or (domain <= 0).any():
            raise ValueError("case_info.json 的 domain 不合法。")
    else:
        messages.append("未找到 case_info.json；坐标范围取显示数据，未画计算盒，未判断周期跳变。")

    track_path = run / "rbc_track.csv"
    track = read_track(track_path, np) if track_path.is_file() else None
    if track is None:
        messages.append("缺少 rbc_track.csv：只生成膜动画，不伪造中心轨迹或速度。")

    meshes = []
    manifest_meshes = []
    total_vertices = 0
    for index, path in enumerate(selected, 1):
        points, faces = read_mesh(path, np)
        total_vertices += len(points)
        if total_vertices > 2_000_000:
            raise ValueError("显示顶点总量超过本查看器的保守上限。使用 --max-frames 50 减少显示帧。")
        meshes.append((points, faces))
        manifest_meshes.append({"file": path.name, "sha256": sha256(path),
                                "vertices": len(points), "faces": len(faces)})
        if index == 1 or index % 25 == 0 or index == len(selected):
            print(f"读取膜网格：{index}/{len(selected)}", flush=True)

    points0, faces0 = meshes[0]
    def surface(points: Any, faces: Any, **kwargs: Any) -> Any:
        return go.Mesh3d(x=points[:, 0].tolist(), y=points[:, 1].tolist(), z=points[:, 2].tolist(),
                         i=faces[:, 0].tolist(), j=faces[:, 1].tolist(), k=faces[:, 2].tolist(),
                         **kwargs)
    scene = go.Figure()
    scene.add_trace(surface(points0, faces0, name="当前 RBC 膜", showlegend=True,
                            opacity=1.0, flatshading=False,
                            hovertemplate="x=%{x:.5g}<br>y=%{y:.5g}<br>z=%{z:.5g}<extra>真实膜节点</extra>"))
    scene.add_trace(surface(points0, faces0, name="第一保存帧（半透明参照，可隐藏）",
                            showlegend=True, opacity=0.18, hoverinfo="skip"))
    csv_jumps = None
    if track is not None:
        com = np.column_stack([track[k] for k in ("comx", "comy", "comz")])
        csv_jumps = (np.abs(np.diff(com, axis=0)) > 0.5 * domain).any(axis=1) \
            if domain is not None else np.zeros(max(0, len(com)-1), dtype=bool)
        ts, pos = split_at_jumps(track["time"], com.T, csv_jumps)
        scene.add_trace(go.Scatter3d(x=pos[0], y=pos[1], z=pos[2], mode="lines+markers",
                                    name="CSV 全程中心轨迹（静态参照）", line=dict(width=4),
                                    marker=dict(size=2), customdata=ts,
                                    hovertemplate="CSV t=%{customdata}<br>x=%{x:.5g}<br>y=%{y:.5g}<br>z=%{z:.5g}<extra></extra>"))
        if bool(csv_jumps.any()):
            messages.append(f"CSV 中有 {int(csv_jumps.sum())} 处相邻坐标跳变超过半个盒长："
                            "绘图在此断线，原坐标没有修改。可能为周期绕回，也可能采样过稀；"
                            "本程序不能据此恢复真实跨盒次数。")
    if domain is not None:
        bx, by, bz = box_edges(domain)
        scene.add_trace(go.Scatter3d(x=bx, y=by, z=bz, mode="lines", line=dict(width=2, dash="dot"),
                                    name="周期计算盒轮廓（不是管壁）", hoverinfo="skip"))
    lower = np.min([p.min(axis=0) for p, _ in meshes], axis=0)
    upper = np.max([p.max(axis=0) for p, _ in meshes], axis=0)
    if domain is not None:
        lower = np.minimum(lower, 0.0)
        upper = np.maximum(upper, domain)
    if track is not None:
        lower = np.minimum(lower, com.min(axis=0))
        upper = np.maximum(upper, com.max(axis=0))
    pad = np.maximum((upper-lower)*0.04, 1e-9)
    axes = {f"{axis}axis": dict(title=f"{axis.upper()}（模拟长度）", autorange=False,
                              range=[float(lower[d]-pad[d]), float(upper[d]+pad[d])])
            for d, axis in enumerate("xyz")}
    animation_frames = []
    slider_steps = []
    for path, (points, faces) in zip(selected, meshes):
        label = numeric_key(path)[0]
        animation_frames.append(go.Frame(name=path.name, traces=[0],
            data=[surface(points, faces)],
            layout=go.Layout(title=dict(text=f"RBC 保存帧 {label} · {html.escape(path.name)}"))))
        slider_steps.append(dict(method="animate", label=str(label),
                                 args=[[path.name], {"mode":"immediate", "frame":{"duration":0,"redraw":True},
                                                     "transition":{"duration":0}}]))
    scene.frames = animation_frames
    scene.update_layout(
        title=dict(text=f"RBC 保存帧 {file_ids[0]} · {html.escape(selected[0].name)}"),
        height=690, margin=dict(l=10,r=10,t=60,b=100),
        scene=dict(**axes, aspectmode="data", uirevision="keep_camera",
                   camera=dict(eye=dict(x=1.4,y=1.2,z=0.9))),
        uirevision="keep_camera", legend=dict(orientation="h", y=1.10),
        sliders=[dict(active=0, currentvalue=dict(prefix="保存帧编号："), steps=slider_steps,
                      x=0.04, len=0.92, y=-0.01, pad=dict(t=35))],
        updatemenus=[dict(type="buttons", direction="left", x=0.02, y=-0.18,
                         buttons=[dict(label="▶ 播放", method="animate",
                                       args=[None,{"fromcurrent":True,"mode":"immediate",
                                                   "frame":{"duration":frame_ms,"redraw":True},
                                                   "transition":{"duration":0}}]),
                                  dict(label="暂停", method="animate",
                                       args=[[None],{"mode":"immediate","frame":{"duration":0,"redraw":False},
                                                     "transition":{"duration":0}}]),
                                  dict(label="回到第一帧", method="animate",
                                       args=[[selected[0].name],{"mode":"immediate","frame":{"duration":0,"redraw":True},
                                                                 "transition":{"duration":0}}])])])
    config = dict(displaylogo=False, responsive=True, scrollZoom=True)
    fragments = [pio.to_html(scene, include_plotlyjs=True, full_html=False,
                             auto_play=False, config=config, div_id="rbc_animation")]
    if track is not None:
        fig_pos = go.Figure()
        for key in ("comx", "comy", "comz"):
            xs, vals = split_at_jumps(track["time"], [track[key]], csv_jumps)
            fig_pos.add_trace(go.Scatter(x=xs,y=vals[0],mode="lines",name=key))
        fig_pos.update_layout(title="中心位置随模拟时间变化（CSV 原始坐标）",height=370,
                              xaxis_title="CSV time（模拟时间；不是秒）",yaxis_title="中心坐标（模拟长度）")
        fragments.append(pio.to_html(fig_pos, include_plotlyjs=False, full_html=False,
                                     config=config, div_id="rbc_position"))
        fig_vel = go.Figure()
        for key in ("vx", "vy", "vz"):
            fig_vel.add_trace(go.Scatter(x=track["time"].tolist(),y=track[key].tolist(),mode="lines",name=key))
        fig_vel.update_layout(title="中心速度分量随模拟时间变化",height=370,
                              xaxis_title="CSV time（模拟时间；不是秒）",yaxis_title="中心速度（模拟长度/模拟时间）")
        fragments.append(pio.to_html(fig_vel, include_plotlyjs=False, full_html=False,
                                     config=config, div_id="rbc_velocity"))

    field_files = sorted([p.name for p in (run/"h5").glob("*") if p.suffix.lower() in {".h5",".xmf",".xdmf"}])
    summary = {
        "run_directory": str(run), "ply_files_found": len(files), "ply_files_displayed": len(selected),
        "selected_frames": manifest_meshes, "csv_rows": None if track is None else track["nrows"],
        "csv_sha256": sha256(track_path) if track is not None else None,
        "case_info_sha256": sha256(info_path) if info_path.is_file() else None,
        "playback_frame_milliseconds": frame_ms, "units": info.get("units", "未标定模拟单位"),
        "h5_xmf_files_present": field_files, "field_rendered": False,
        "warnings_and_limits": messages, "visual_validation": "USER_REVIEW_PENDING",
        "dependencies": {"numpy":np.__version__,"trimesh":trimesh.__version__,"plotly":plotly.__version__}}
    rows_html = "".join(f"<tr><th>{html.escape(k)}</th><td>{html.escape(str(v))}</td></tr>"
                        for k,v in [("运行目录",run),("PLY 文件 / 显示帧数",f"{len(files)} / {len(selected)}"),
                                    ("CSV 行数",'缺失' if track is None else track['nrows']),
                                    ("模拟参数 dt",info.get('dt','未知')),
                                    ("播放间隔",f"{frame_ms} ms/保存帧，仅影响观看"),
                                    ("状态","可视化读取完成；物理验证未执行")])
    notes = "".join(f"<li>{html.escape(m)}</li>" for m in messages)
    page = """<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Mirheo RBC 结果查看</title>
<style>body{margin:0;font-family:system-ui,'Microsoft YaHei',sans-serif;background:#f5f5f5;color:#202020}
main{max-width:1250px;margin:24px auto;padding:0 18px}section{background:white;border:1px solid #ddd;
border-radius:10px;margin:18px 0;padding:20px}h1{font-size:27px}p,li,td,th{line-height:1.7}
table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left;word-break:break-all}
th{width:200px}code{background:#f2f2f2;padding:2px 5px}li{margin:6px 0}.plot{padding:10px}</style></head><body><main>
<h1>Mirheo 单红细胞 · 已有计算结果查看</h1>
<p>这不是重新仿真：播放的是磁盘中已经保存的膜网格，曲线来自 rbc_track.csv。</p>
<section><h2>1. 本次读取了什么</h2><table>""" + rows_html + """</table></section>
<section><h2>2. 怎样观看</h2><p>先点击“▶ 播放”，也可以拖动保存帧滑块。
按住鼠标左键拖动以旋转，滚轮缩放；点击图例可隐藏第一帧参照、轨迹或计算盒。
视角与坐标范围固定，不会每帧追着细胞自动居中。</p>
<p>半透明的第一帧只是位置参照，不是第二个红细胞。CSV 全程轨迹是静态参照，
不是与每帧同步的当前位置标记。右上角相机图标可保存当前图像。</p></section>
""" + "".join(f'<section class="plot">{fragment}</section>' for fragment in fragments) + \
        f'<section><h2>3. 必须保留的说明</h2><ul>{notes}</ul></section>' + \
        f'<section><h2>4. 体流场文件</h2><p>检测到 {len(field_files)} 个 H5/XMF/XDMF 文件。'
    page += "本页没有绘制它们；查看液体箭头或切片，请在 ParaView 中打开配套 XMF/XDMF，保留其 H5 文件。</p></section>"
    page += "</main></body></html>"
    outdir = run / "visualization"
    outdir.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output = outdir / f"rbc_view_{stamp}.html"
    output.write_text(page,encoding="utf-8")
    output.with_suffix(".json").write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    print(f"\nPLY：发现 {len(files)} 帧，显示 {len(selected)} 帧。")
    print(f"CSV：{'未读取' if track is None else str(track['nrows'])+' 行'}。")
    print("可视化页面：",output)
    print("只写入 visualization 子目录，没有修改原始仿真文件。")
    print("Windows 打开结果文件夹的命令：")
    print(f'explorer.exe "$(wslpath -w {shlex.quote(str(outdir))})"')
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run",type=Path,help="包含 ply/ 和 rbc_track.csv 的实际运行目录；省略则列出选项")
    parser.add_argument("--max-frames",type=int,default=200,help="页面最多显示多少个保存帧，默认200，包含首尾")
    parser.add_argument("--frame-ms",type=int,default=100,help="每个显示帧的播放间隔毫秒数，默认100；不改变模拟时间")
    parser.add_argument("--open",action="store_true",help="完成后在 Windows 文件资源管理器中打开页面所在文件夹")
    args = parser.parse_args()
    if args.max_frames < 2 or args.frame_ms < 1:
        parser.error("--max-frames 至少为2；--frame-ms 必须为正整数。")
    try:
        run = select_run(args.run,Path(__file__).resolve().parents[1])
        result = build_report(run,args.max_frames,args.frame_ms)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"查看失败：{exc}",file=sys.stderr)
        return 1
    if args.open:
        if shutil.which("wslpath") and shutil.which("explorer.exe"):
            try:
                winpath = subprocess.check_output(["wslpath","-w",str(result.parent)],text=True).strip()
                subprocess.run(["explorer.exe",winpath],check=False)
            except (OSError, subprocess.SubprocessError) as exc:
                print(f"页面已生成，但打开资源管理器失败：{exc}")
        else:
            print("当前环境未发现 WSL 文件资源管理器接口，请手动打开输出 HTML。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
