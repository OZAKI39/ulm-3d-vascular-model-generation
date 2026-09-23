#!/usr/bin/env python3
"""Reproducible CPU-only Particle-0 plots/data. Run --stage 0..6 or --stage all.

All plotted measurements are retained in data/. Scientific surfaces are read
from hash-verified frozen inputs. No image synthesis or field smoothing.
"""
from pathlib import Path
import argparse
import csv
from dataclasses import fields
import json
import os
import sys
import time
import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField, FlowSample
from particle_3d.validation_cases import (SEED, AFFINE_ATOL, AFFINE_RTOL, affine_case, real_node_case,
    shared_faces, shared_face_case, classification_case, real_flow_case, interior_positions)

REPORT = PACKAGE / "reports/particle0"
DATA, FIGURES = REPORT / "data", REPORT / "figures"
COLORS = dict(WALL="#9ab2bd", INLET="#d83d38", OUTLET_01="#1c81c5", OUTLET_02="#25964a", OUTLET_03="#9d50b7")
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 9,
                     "figure.dpi": 110, "savefig.dpi": 180})


def json_write(path, value):
    def convert(x):
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, np.ndarray): return x.tolist()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False, default=convert) + "\n")


def csv_write(name, rows):
    with (DATA / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def save(fig, name, footer):
    fig.text(.5, .025, footer, ha="center", va="bottom", fontsize=10,
             bbox=dict(facecolor="#f1f4f6", edgecolor="none", pad=8))
    fig.savefig(FIGURES / name, bbox_inches="tight")
    plt.close(fig)


def note(stage, title, text, filenames):
    (REPORT / f"{stage:02d}_{title}.md").write_text(f"# P0.{stage} 检查说明\n\n{text}\n\n" +
        "\n".join(f"![检查图](figures/{f})" for f in filenames) + "\n")


def vessel(ax, boundaries, alpha=.14, labels=True):
    for role, surface in boundaries.items():
        triangles = surface.faces.reshape(-1, 4)[:, 1:]
        xyz = np.asarray(surface.points)[triangles] * 1e6  # visualization only
        artist = Poly3DCollection(xyz, facecolor=COLORS[role], edgecolor="none",
                                   alpha=alpha if role == "WALL" else .9, rasterized=True)
        ax.add_collection3d(artist)
        if labels and role != "WALL":
            p = np.asarray(surface.points).mean(0) * 1e6
            ax.text(*p, "  " + role, color=COLORS[role], fontsize=10, fontweight="bold")
    xyz = np.asarray(boundaries["WALL"].points) * 1e6
    lo, hi = xyz.min(0), xyz.max(0)
    margin = .06 * (hi - lo)
    ax.set(xlim=(lo[0] - margin[0], hi[0] + margin[0]), ylim=(lo[1] - margin[1], hi[1] + margin[1]),
           zlim=(lo[2] - margin[2], hi[2] + margin[2]), xlabel="x (µm)", ylabel="y (µm)", zlabel="z (µm)")
    ax.set_box_aspect(hi - lo)
    ax.view_init(elev=23, azim=-57)


def stage0(summary, mesh, flow, boundaries, field):
    summary = dict(summary)
    summary.update(seed=SEED, condition_inf_max=float(field.geometry.condition_inf.max()),
                   condition_inf_median=float(np.median(field.geometry.condition_inf)),
                   weight_tolerance_max=float(field.geometry.weight_tolerance.max()),
                   candidate_padding_m=field.geometry.candidate_padding_m,
                   figure_geometry="All original boundary triangles; no remeshing, smoothing or scalar modification.")
    json_write(DATA / "00_frozen_input_summary.json", summary)
    fig = plt.figure(figsize=(14, 9)); ax = fig.add_subplot(121, projection="3d")
    vessel(ax, boundaries, .55)
    ax.set_title("Check: frozen vessel and five boundary roles")
    text_ax = fig.add_subplot(122); text_ax.axis("off")
    text_ax.text(.03, .9, "FROZEN INPUT AUDIT", fontsize=19, weight="bold")
    details = [f"Nodes: {summary['node_count']:,}", f"Tetrahedra: {summary['tetra_count']:,}",
               "POINT Velocity: float64[3], m/s", "POINT Pressure: float64, Pa", "Coordinates / core calculations: SI (m)",
               f"Mesh SHA: {summary['mesh_sha256'][:20]}…", f"Flow SHA: {summary['flow_sha256'][:20]}…",
               "Point coordinates AND point order: identical", "Canonical tetra rows: corresponding",
               "Flow local vertex permutation: [1, 0, 2, 3]", "Interpolation uses canonical volume connectivity"]
    text_ax.text(.03, .80, "\n\n".join(details), fontsize=11, va="top")
    ax.legend(handles=[Line2D([0], [0], color=c, lw=5, label=f"{role} (ID {summary['boundary_face_ids'][role]})") for role, c in COLORS.items()],
              loc="lower left", bbox_to_anchor=(-.1, -.04))
    fig.subplots_adjust(bottom=.14, wspace=.08)
    save(fig, "00_frozen_input_overview.png", f"PASS: 31 scientific hashes checked; {summary['node_count']:,} points / {summary['tetra_count']:,} tetra rows agree.\nAll {summary['permuted_tetra_count']:,} rows have the documented local permutation. Five boundary IDs agree.")
    note(0, "frozen_input_audit", f"读取的是独立 worktree 内的冻结 handoff。网格 {summary['node_count']} 个节点、{summary['tetra_count']} 个四面体；节点坐标和顺序逐项完全一致。\n\n全部单元的局部顶点顺序为 [1,0,2,3]，正式编号和插值顶点一律取体网格。文件哈希、单位和五个边界编号均已核对。旧 WSL FEM 开发目录缺少 handoff 文件且有历史改动，所以保持原样，未拿它替代冻结输入。\n\n图形数据入口为 data/00_frozen_input_summary.json，其中保存各表面路径和 SHA；画图读取全部原始三角面。\n\nSonoVue 仓库只登记存在，没有接入；2D 代码仅阅读结构，没有复制。", ["00_frozen_input_overview.png"])
    return summary


def stage1():
    case = affine_case(); b = case["batch"]
    uerr = np.max(np.abs(b.velocity_m_s - case["velocity"]), axis=1)
    perr = np.abs(b.pressure_pa - case["pressure"])
    gerr = np.max(np.abs(b.velocity_gradient_s_inv - case["gradient"]), axis=(1, 2))
    rows = []
    for k, p in enumerate(case["positions"]):
        row = dict(sample_id=k, kind=case["kinds"][k], x_m=p[0], y_m=p[1], z_m=p[2],
                   analytic_pressure_pa=case["pressure"][k], sampled_pressure_pa=b.pressure_pa[k],
                   velocity_error_m_s=uerr[k], pressure_error_pa=perr[k], gradient_error_s_inv=gerr[k])
        actual_weights = case["field"].geometry.weights(p, 0)
        for i in range(4): row.update({f"analytic_weight_{i}": case["weights"][k, i], f"sampled_weight_{i}": actual_weights[i]})
        for i, axis in enumerate("xyz"):
            row.update({f"analytic_velocity_{axis}_m_s": case["velocity"][k, i], f"sampled_velocity_{axis}_m_s": b.velocity_m_s[k, i],
                        f"analytic_vorticity_{axis}_s_inv": case["vorticity"][i], f"sampled_vorticity_{axis}_s_inv": b.vorticity_s_inv[k, i]})
        for i in range(3):
            for j in range(3):
                row.update({f"analytic_gradient_{i}{j}_s_inv": case["gradient"][i, j], f"sampled_gradient_{i}{j}_s_inv": b.velocity_gradient_s_inv[k, i, j],
                            f"analytic_strain_{i}{j}_s_inv": case["strain"][i, j], f"sampled_strain_{i}{j}_s_inv": b.strain_rate_s_inv[k, i, j]})
        rows.append(row)
    csv_write("01_affine_field_validation.csv", rows)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for col, (name, unit, truth, actual, error) in enumerate([
        ("Velocity components", "m/s", case["velocity"], b.velocity_m_s, uerr),
        ("Pressure", "Pa", case["pressure"], b.pressure_pa, perr),
        ("Velocity gradient components", "1/s", np.broadcast_to(case["gradient"], b.velocity_gradient_s_inv.shape), b.velocity_gradient_s_inv, gerr)]):
        ax = axes[0, col]
        ax.scatter(truth.ravel(), actual.ravel(), s=14, alpha=.55, label="sampled")
        bounds = [float(truth.min()), float(truth.max())]
        ax.plot(bounds, bounds, "k--", lw=1, label="exact identity")
        ax.set(xlabel=f"Analytic ({unit})", ylabel=f"Sampled ({unit})", title=f"Check: {name}")
        ax.legend(); ax.grid(alpha=.2)
        ax = axes[1, col]
        ax.plot(np.arange(len(error)), error, ".", label="absolute error")
        ax.set(xlabel="Query index (dimensionless)", ylabel=f"Absolute error ({unit})", title=f"Max error = {error.max():.3e} {unit}")
        ax.legend(); ax.grid(alpha=.2)
    fig.suptitle("P0.1 — known affine answers: interior, faces, edges and vertices", fontsize=16)
    fig.subplots_adjust(bottom=.17, top=.9, hspace=.4, wspace=.35)
    passed = all(np.allclose(actual, truth, atol=AFFINE_ATOL, rtol=AFFINE_RTOL) for actual, truth in [
        (b.velocity_m_s, case["velocity"]), (b.pressure_pa, case["pressure"]), (b.velocity_gradient_s_inv, case["gradient"]),
        (b.vorticity_s_inv, case["vorticity"]), (b.strain_rate_s_inv, case["strain"])])
    save(fig, "01_affine_field_validation.png", f"{'PASS' if passed else 'FAIL'}: N={len(rows)}; max velocity={uerr.max():.3e} m/s; pressure={perr.max():.3e} Pa; gradient={gerr.max():.3e} 1/s.\nFixed atol=rtol=256×float64 eps={AFFINE_ATOL:.3e}; all five point categories included.")
    metric = dict(passed=passed, count=len(rows), max_affine_velocity_error=float(uerr.max()), max_affine_pressure_error=float(perr.max()),
                  max_affine_gradient_error=float(gerr.max()), atol=AFFINE_ATOL, rtol=AFFINE_RTOL)
    note(1, "affine_validation", f"人工设置一个每个方向都线性变化的速度和压力场，答案事先知道。barycentric 就是位置对四个顶点的四个贡献权重；gradient 是速度沿各方向变化的快慢。\n\n{len(rows)} 个位置覆盖中心、随机内部、面、边和顶点。最大速度误差 {uerr.max():.3e} m/s，压力误差 {perr.max():.3e} Pa，梯度误差 {gerr.max():.3e} 1/s；保持预先固定的容差。", ["01_affine_field_validation.png"])
    return metric


def stage2(field, boundaries):
    rows = real_node_case(field, boundaries); csv_write("02_node_back_sampling_error.csv", rows)
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    metric = dict(count=len(rows), passed=all(r["inside_lumen"] and r["velocity_error_m_s"] <= r["velocity_bound_m_s"] and r["pressure_error_pa"] <= r["pressure_bound_pa"] for r in rows))
    for ax, key, name, unit in zip(axes, ["velocity_error_m_s", "pressure_error_pa"], ["Velocity", "Pressure"], ["m/s", "Pa"]):
        errors = np.array([r[key] for r in rows]); ordered = np.argsort(errors)
        bounds = np.array([r[key.replace("error", "bound")] for r in rows])
        ax.scatter(np.arange(len(rows)), errors[ordered], s=9, label="measured error")
        ax.plot(bounds[ordered], alpha=.45, lw=.8, label="local roundoff bound")
        ax.set_yscale("symlog", linthresh=max(float(errors.max()) * 1e-4, 1e-30))
        ax.set(title=f"Check: recover nodal {name.lower()}", xlabel="Sorted query index (dimensionless)", ylabel=f"Absolute error ({unit})")
        ax.text(.04, .95, f"N={len(rows)}\nmax={errors.max():.3e}\nmedian={np.median(errors):.3e}\np95={np.percentile(errors, 95):.3e}",
                transform=ax.transAxes, va="top", bbox=dict(facecolor="white", alpha=.9))
        ax.legend(loc="lower right"); ax.grid(alpha=.2)
        metric[f"max_node_{name.lower()}_error"] = float(errors.max())
        metric[f"median_{name.lower()}_error"] = float(np.median(errors))
        metric[f"p95_{name.lower()}_error"] = float(np.percentile(errors, 95))
    fig.subplots_adjust(bottom=.18, wspace=.3)
    save(fig, "02_node_back_sampling_error.png", f"{'PASS' if metric['passed'] else 'FAIL'}: all {len(rows)} nodes classified inside; each error <= its geometry-derived bound (CSV).\nIncludes fixed/random nodes, all ports, wall and interior. Symlog preserves exact zeros.")
    note(2, "node_back_sampling", f"直接把原网格节点坐标送回查询接口，与文件内节点值逐项比较，共 {len(rows)} 次。包括四个进出口、壁面和血管内部。\n\n最大速度误差 {metric['max_node_velocity_error']:.3e} m/s，压力误差 {metric['max_node_pressure_error']:.3e} Pa；图中同时画出每个点由局部数值精度推得的上界。零误差没有被替换成非零数据。", ["02_node_back_sampling_error.png"])
    return metric


def stage3(field):
    faces = shared_faces(field); rows, checks = shared_face_case(field, faces)
    csv_write("03_shared_face_continuity.csv", rows); json_write(DATA / "03_shared_face_metrics.json", checks)
    fig, axes = plt.subplots(3, 4, figsize=(18, 11))
    for col in range(4):
        subset = [r for r in rows if r["face_id"] == col]
        for i, (key, label) in enumerate([("speed_m_s", "Speed (m/s)"), ("pressure_pa", "Pressure (Pa)"), ("gradient_norm_s_inv", "Gradient norm (1/s)")]):
            ax = axes[i, col]
            ax.plot([r["signed_distance_m"] * 1e6 for r in subset], [r[key] for r in subset], "o-", ms=3, label="queried field")
            ax.axvline(0, color="gray", ls="--", lw=.8, label="shared face")
            ax.set(xlabel="Signed distance (µm)", ylabel=label,
                   xticks=np.array([-1., 0., 1.]) * faces[col]["distance"] * 1e6)
            ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3g"))
            ax.ticklabel_format(axis="y", useOffset=False); ax.grid(alpha=.2)
            if i == 0: ax.set_title(f"Face {col}: cells {checks[col]['left']} / {checks[col]['right']}")
            if col == 0: ax.legend()
    fig.suptitle("Check: velocity/pressure continuous across faces; gradient may jump (4 of 8 tested faces shown)", fontsize=15)
    fig.subplots_adjust(bottom=.13, top=.92, hspace=.48, wspace=.45)
    umax = max(c["velocity_trace_error"] for c in checks); pmax = max(c["pressure_trace_error"] for c in checks)
    jumps = [c["gradient_jump_s_inv"] for c in checks]
    passed = all(c["velocity_trace_error"] <= c["velocity_bound"] and c["pressure_trace_error"] <= c["pressure_bound"] for c in checks)
    save(fig, "03_shared_face_continuity.png", f"{'PASS' if passed else 'FAIL'}: 8 faces; on-face trace max error: velocity {umax:.3e} m/s, pressure {pmax:.3e} Pa (bounds in JSON).\nMeasured gradient jump range {min(jumps):.3e}–{max(jumps):.3e} 1/s; no jump threshold and no smoothing.")
    note(3, "shared_face_continuity", f"每个共享面都沿法线取两侧位置和面上位置。两侧在面上的速度误差最大 {umax:.3e} m/s，压力误差最大 {pmax:.3e} Pa。\n\n速度和压力连续不意味着它们必须处处相等：离开面后会随位置变化。梯度在四面体内保持常数，跨面可以跳变；这里测得跳变量范围 {min(jumps):.3e}–{max(jumps):.3e} 1/s，没有平滑，也没有把跳变判为失败。左右全部九个梯度分量均已保存。\n\n顶点排序回归还覆盖人工四面体全部 24 种 flow 局部排列，以及真实网格另一种局部排列。", ["03_shared_face_continuity.png"])
    return dict(passed=passed, face_count=len(faces), max_velocity_trace_error=umax, max_pressure_trace_error=pmax,
                gradient_jump_min=min(jumps), gradient_jump_max=max(jumps))


def stage4(field, boundaries):
    rows = classification_case(field, boundaries); csv_write("04_inside_outside_classification.csv", rows)
    errors = sum(not r["correct"] for r in rows)
    fig = plt.figure(figsize=(17, 9)); ax = fig.add_subplot(121, projection="3d", computed_zorder=False)
    vessel(ax, boundaries, .10, labels=False)
    styles = [("inside", lambda r: r["expected_inside"] and r["kind"] != "boundary", "#167bb5", "o"),
              ("outside (near surface)", lambda r: not r["expected_inside"] and r["kind"] != "bbox_outside", "#d13937", "x"),
              ("boundary", lambda r: r["kind"] == "boundary", "#e3a200", "^"),
              ("WRONG classification", lambda r: not r["correct"], "#ff00dc", "*")]
    for name, select, color, marker in styles:
        subset = [r for r in rows if select(r)]
        if subset:
            xyz = np.array([[r[k] for k in ["x_m", "y_m", "z_m"]] for r in subset]) * 1e6
            ax.scatter(*xyz.T, color=color, marker=marker, s=10 if marker == "o" else 28 if marker != "*" else 90,
                       label=f"{name}: {len(subset)}", depthshade=False, zorder=10 + len(name))
        else: ax.scatter([], [], [], color=color, marker=marker, label=f"{name}: 0")
    ax.set_title("Check: inside / boundary / near-outside locations")
    ax.legend(loc="lower center", bbox_to_anchor=(.5, -.10), ncol=2)
    ax = fig.add_subplot(122)
    near = [r for r in rows if r["role"]]
    for role_index, (role, color) in enumerate(COLORS.items()):
        subset = [r for r in near if r["role"] == role]
        for inside, marker in [(True, "o"), (False, "x")]:
            selected = [r for r in subset if r["inside_lumen"] == inside]
            xs = np.array([r["offset_m"] for r in selected]) * 1e9
            ax.scatter(xs, np.full(len(xs), role_index), color=color, marker=marker, s=24, alpha=.6)
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set(xlabel="Signed outward offset from surface (nm)", ylabel="Boundary role (category)",
           yticks=range(5), yticklabels=list(COLORS), ylim=(-.5, 4.5), title="Check: even small outward offsets are rejected")
    ax.set_xticks([-.1, 0., 1e-4, .1])
    ax.axvline(0, color="gray", ls="--"); ax.grid(alpha=.2)
    ax.legend(handles=[Line2D([], [], marker="o", color="black", ls="", label="inside / boundary: True"),
                       Line2D([], [], marker="x", color="black", ls="", label="outside: False")], loc="upper left")
    fig.subplots_adjust(bottom=.19, wspace=.65)
    distances = [r["offset_m"] for r in rows if r["kind"] == "near_outside"]
    save(fig, "04_inside_outside_classification.png", f"{'PASS' if errors == 0 else 'FAIL'}: {len(rows)} queries; {errors} misclassified; six additional bounding-box outside points tested (CSV).\nSmallest outward test offset={min(distances):.3e} m; locator candidate padding={field.geometry.candidate_padding_m:.3e} m. Outside fields are all NaN.")
    note(4, "inside_outside_contract", f"真实几何共检查 {len(rows)} 个位置，错误分类 {errors} 个。精确表面点视为场存在；稍向域外移动后必须返回域外，最小测试偏移为 {min(distances):.3e} m。\n\n左图中很近的内外点会重叠，右图用实际有符号距离放大显示分类。所有域外量都是 NaN，不是零速度。这个点的定义不表示有半径的微泡或 RBC 能合法贴在墙上；本阶段没有 wall gap。", ["04_inside_outside_classification.png"])
    return dict(passed=errors == 0, count=len(rows), errors=errors, minimum_outside_offset_m=min(distances),
                boundary_count=sum(r["kind"] == "boundary" for r in rows), outside_count=sum(not r["expected_inside"] for r in rows))


def stage5(field, boundaries):
    rows, positions, ids, weights, batch = real_flow_case(field); csv_write("05_real_flow_samples.csv", rows)
    physical = np.concatenate([getattr(batch, name).ravel() for name in ["velocity_m_s", "pressure_pa", "velocity_gradient_s_inv", "vorticity_s_inv", "strain_rate_s_inv"]])
    nan_count, inf_count = int(np.isnan(physical).sum()), int(np.isinf(physical).sum())
    speed = np.linalg.norm(batch.velocity_m_s, axis=1)
    fig = plt.figure(figsize=(12, 10)); ax = fig.add_subplot(111, projection="3d")
    vessel(ax, boundaries, .18)
    xyz = positions * 1e6
    colors = plt.cm.viridis(plt.Normalize(speed.min(), speed.max())(speed))
    # Normalize only display lengths so slow branches retain visible directions.
    # Actual speeds remain in the color scale and CSV. No core values are changed.
    vectors = batch.velocity_m_s / speed[:, None] * 1.5
    ax.quiver(*xyz.T, *vectors.T, colors=colors, length=1, normalize=False, linewidth=.8, arrow_length_ratio=.3)
    ax.scatter(*xyz.T, c=speed, cmap="viridis", s=2)
    scalar = plt.cm.ScalarMappable(norm=plt.Normalize(speed.min(), speed.max()), cmap="viridis")
    fig.colorbar(scalar, ax=ax, shrink=.58, pad=.06, label="Actual speed (m/s)")
    ax.set_title("Check: 400 frozen-flow velocity directions inside the real vessel", pad=20)
    ax.legend(handles=[Line2D([0], [0], color=c, lw=4, label=role) for role, c in COLORS.items()], loc="lower left")
    fig.subplots_adjust(bottom=.16)
    passed = nan_count == inf_count == 0 and bool(batch.inside_lumen.all())
    save(fig, "05_real_flow_velocity_vectors.png", f"{'PASS' if passed else 'FAIL'}: 400/400 valid; NaN={nan_count}, Inf={inf_count}; speed range {speed.min():.3e}–{speed.max():.3e} m/s.\nSeed={SEED + 5}; volume-uniform samples. Display arrows normalized to 1.5 µm; color = actual speed; no field changes.")
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for ax, key, label, unit in zip(axes.ravel(), ["speed_m_s", "pressure_pa", "vorticity_norm_s_inv", "strain_norm_s_inv"],
                                  ["Speed", "Pressure", "Vorticity magnitude", "Strain-rate Frobenius norm"], ["m/s", "Pa", "1/s", "1/s"]):
        values = np.array([r[key] for r in rows])
        artist = ax.scatter(xyz[:, 0], xyz[:, 1], c=values, s=17, cmap="viridis")
        ax.set(xlabel="x (µm)", ylabel="y (µm)", title=f"Check: {label}\nrange {values.min():.3e}–{values.max():.3e} {unit}")
        ax.set_aspect("equal"); ax.grid(alpha=.15)
        fig.colorbar(artist, ax=ax, label=f"{label} ({unit})", shrink=.85)
    fig.suptitle("P0.5 — scalar diagnostics, x-y projection (all original z coordinates retained in CSV)")
    fig.subplots_adjust(bottom=.14, top=.88, hspace=.42, wspace=.3)
    save(fig, "06_real_flow_scalar_diagnostics.png", f"{'PASS' if passed else 'FAIL'}: 400 valid point samples; NaN={nan_count}, Inf={inf_count} across velocity, pressure and all derivative components.\nThese are field diagnostics; no new FEM convergence or physics acceptance is implied.")
    note(5, "real_flow_visual_check", f"固定 seed={SEED + 5}，按四面体体积加权选取单元，再在单元内均匀取点，共 400 个位置。速度范围 {speed.min():.3e}–{speed.max():.3e} m/s，有效样本 NaN={nan_count}、Inf={inf_count}。\n\nvorticity 是速度场在局部的旋转趋势；strain rate 是流动在局部拉伸和剪切的快慢。图中应检查箭头与血管方向是否协调、颜色是否出现孤立异常；箭头长度只作显示缩放，原值在 CSV。\n\n这些图仍等待用户人工审核，自动检查没有代替视觉判断。", ["05_real_flow_velocity_vectors.png", "06_real_flow_scalar_diagnostics.png"])
    return dict(passed=passed, count=len(rows), nan_count_valid_samples=nan_count, inf_count_valid_samples=inf_count,
                speed_min=float(speed.min()), speed_max=float(speed.max()))


def stage6(field):
    faces = shared_faces(field)
    positions, _, _ = interior_positions(field, 400)
    special = np.vstack((positions, [f["center"] for f in faces], field.points_m[[0, 100, 1000]], [[0., 0., 0.], [1., 1., 1.]]))
    batch = field.sample_many(special)
    differences = 0
    for i, p in enumerate(special):
        scalar = field.sample(p)
        for f in fields(FlowSample):
            if not np.array_equal(getattr(batch, f.name)[i], getattr(scalar, f.name), equal_nan=True): differences += 1
    json_write(DATA / "07_scalar_batch_equivalence.json", dict(count=len(special), differing_fields=differences,
        categories=["random_inside", "shared_face", "vertex", "outside"], comparison="exact, equal_nan=True"))
    rows = []
    field.sample_many(positions[:10])  # explicit warm-up, construction excluded
    for count in [100, 1000, 10000]:
        points, _, _ = interior_positions(field, count, SEED + 6)
        start = time.perf_counter(); scalar = [field.sample(p) for p in points]; scalar_time = time.perf_counter() - start
        start = time.perf_counter(); batch = field.sample_many(points); batch_time = time.perf_counter() - start
        rows.append(dict(point_count=count, scalar_wall_time_s=scalar_time, batch_wall_time_s=batch_time,
                         valid_count=int(batch.inside_lumen.sum()), seed=SEED + 6))
        print(f"Timing {count}: scalar={scalar_time:.6f}s batch={batch_time:.6f}s", flush=True)
    csv_write("07_sampling_runtime_smoke.csv", rows)
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot([r["point_count"] for r in rows], [r["scalar_wall_time_s"] for r in rows], "o-", label="loop over sample()")
    ax.plot([r["point_count"] for r in rows], [r["batch_wall_time_s"] for r in rows], "s-", label="sample_many()")
    for r in rows: ax.annotate(f"{r['batch_wall_time_s']:.3f} s", (r["point_count"], r["batch_wall_time_s"]), xytext=(5, 8), textcoords="offset points")
    ax.set(xscale="log", yscale="log", xlabel="Number of query points (count)", ylabel="Wall time (s)",
           title="Check: scalar/batch consistency and CPU runtime smoke\nINFORMATIONAL PERFORMANCE SMOKE\nNOT A PRODUCTION BENCHMARK")
    ax.legend(); ax.grid(alpha=.3)
    fig.subplots_adjust(bottom=.2)
    save(fig, "07_sampling_runtime_smoke.png", f"{'PASS' if differences == 0 else 'FAIL'}: exact equality for all seven fields at {len(special)} mixed queries; differing fields={differences}.\nTiming: one repetition after 10-point warm-up; includes result allocation, excludes mesh load. No performance gate.")
    note(6, "scalar_batch_and_runtime", f"逐点接口与批量接口在 {len(special)} 个混合位置的七项返回值逐项完全一致，差异项数 {differences}，域外 NaN 按相同无效值比较。\n\n计时仅做 WSL CPU 冒烟检查：100、1000、10000 点，每组一次，排除载入网格和构建 locator 时间，包含返回值分配。性能不设通过门槛，也不是生产基准。", ["07_sampling_runtime_smoke.png"])
    return dict(passed=differences == 0, count=len(special), differing_fields=differences, timings=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["all", "0", "1", "2", "3", "4", "5", "6"], default="all")
    parser.add_argument("--fem-root", type=Path, default=Path(os.environ.get("PARTICLE0_FEM_ROOT", PACKAGE.parent / "formal_3D_flow_solver/FEM_SimVascular")))
    args = parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    if args.stage == "1":
        json_write(DATA / "metrics_p01.json", stage1()); return
    summary, mesh, flow, boundaries = read_frozen(args.fem_root)
    field = FrozenFEMField.from_grids(mesh, flow)
    functions = {0: lambda: stage0(summary, mesh, flow, boundaries, field), 1: stage1,
                 2: lambda: stage2(field, boundaries), 3: lambda: stage3(field),
                 4: lambda: stage4(field, boundaries), 5: lambda: stage5(field, boundaries), 6: lambda: stage6(field)}
    for stage in range(7) if args.stage == "all" else [int(args.stage)]:
        print(f"Generating P0.{stage}", flush=True)
        metrics = functions[stage]()
        json_write(DATA / f"metrics_p0{stage}.json", metrics)
        print(json.dumps({k: v for k, v in metrics.items() if k in ["passed", "count", "errors", "max_node_velocity_error", "max_affine_gradient_error"]}), flush=True)


if __name__ == "__main__": main()
