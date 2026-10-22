"""Build the review from closed result records; never infer completed simulations."""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.io.artifacts import write_json
from mevo_mb.visualization.plots import save


def main():
    root = Path.cwd()
    result_dir = root/"results"
    demo = json.loads((result_dir/"DEMO_CAMPAIGN.json").read_text())
    assert demo["status"] == "COMPLETE"
    convergence = json.loads((result_dir/"CONVERGENCE_STATE.json").read_text())
    assert convergence["status"] == "COMPLETE"
    lines = ["# PZT → acoustic field → microbubble demo：实际验证结果", "",
             "本报告只读取已关闭run。基础单位SI；当前为线性声场、径向响应替代模型、单向耦合，不是经生理标定的SonoVue预测。", "",
             "## 声场", "", "| Frequency | Basis/direct relative L2 | Masked phase RMS (rad) |", "|---|---:|---:|"]
    figures = []
    resources = []
    for frequency in [2, 1]:
        state = json.loads((result_dir/f"acoustic{frequency}_STATE.json").read_text())
        assert state["status"] == "COMPLETE" and state["validation"]["status"] == "PASS"
        val = state["validation"]
        lines += [f"| {frequency} MHz | {val['relative_L2']:.6g} | {val['phase_RMS_rad_masked']:.6g} |"]
        figures.append(Path(state["validation_run"])/"figures/superposition.png")
        for meta in sorted((root/state["run"]).glob("acoustic/basis_*/metadata.json")):
            resources.append(json.loads(meta.read_text()))
        resources.append(json.loads((result_dir/f"grid{frequency}_STATE.json").read_text())["solver_metadata"])
    lines += ["", "CUDA binary的HDF5依赖失败保留；本次正式声学demo使用CuPy GPU。整个域不保存完整时间历史，末4周期仅记录ROI。", "",
              "## 迁移与聚集", "", "| Case | N | Time (ms) | Mean lateral displacement (µm) | NN initial → final (µm) | Largest cluster initial → final | Near-contact final |", "|---|---:|---:|---:|---:|---:|---:|"]
    for case in demo["scenarios"]:
        run = Path(case["run"])
        summary = case["summary"]
        a, b = summary["initial"], summary["final"]
        lines.append(f"| {case['scenario']}/{case['pattern']} | {summary['particle_count']} | {summary['duration_s']*1e3:g} | {b['mean_lateral_displacement_m']*1e6:.6g} | {a['mean_nearest_neighbor_m']*1e6:.5g} → {b['mean_nearest_neighbor_m']*1e6:.5g} | {a['largest_cluster_size']} → {b['largest_cluster_size']} | {b['near_contact_fraction']:.4g} |")
        figures.extend(run/"figures"/name for name in ["acoustic_and_primary_force.png", "trajectories_and_distribution.png", "metrics.png", "pair_distance_distribution.png"])
    figures.insert(0, Path(demo["scenarios"][0]["run"])/"figures/geometry_and_field.png")
    control = json.loads((result_dir/"SECONDARY_OFF_STATE.json").read_text())
    assert control["status"] == "COMPLETE"
    active = next(x for x in demo["scenarios"] if x["scenario"] == "aggregation")
    on, off = active["summary"], control["summary"]
    with h5py.File(Path(active["run"])/"particles/trajectories.h5") as a, h5py.File(Path(control["run"])/"particles/trajectories.h5") as b:
        assert np.array_equal(a["positions_m"][0], b["positions_m"][0])
        assert np.array_equal(a["radii_m"][:], b["radii_m"][:])
    on_control = json.loads((Path(active["run"])/"acoustic/CONTROL.json").read_text())
    off_control = json.loads((Path(control["run"])/"acoustic/CONTROL.json").read_text())
    for key in ["source_control_real_Pa", "source_control_imag_Pa", "basis_sha256"]:
        assert on_control[key] == off_control[key]
    lines += ["", "### 相同初值的泡间作用对照", "", "| Secondary force | Final NN (µm) | Largest cluster | Near-contact fraction |", "|---|---:|---:|---:|",
              f"| ON | {on['final']['mean_nearest_neighbor_m']*1e6:.6g} | {on['final']['largest_cluster_size']} | {on['final']['near_contact_fraction']:.6g} |",
              f"| OFF | {off['final']['mean_nearest_neighbor_m']*1e6:.6g} | {off['final']['largest_cluster_size']} | {off['final']['near_contact_fraction']:.6g} |", "",
              "两次初始位置、粒径、源控制量一致。此对照只支持当前数值模型内部归因，不代表真实微泡近场聚集定律已验证。", ""]
    migrations = [x for x in demo["scenarios"] if x["scenario"] == "migration"]
    with h5py.File(Path(migrations[0]["run"])/"particles/trajectories.h5") as a, h5py.File(Path(migrations[1]["run"])/"particles/trajectories.h5") as b:
        assert np.array_equal(a["positions_m"][0], b["positions_m"][0])
        x0, x1, x2 = a["positions_m"][0], a["positions_m"][-1], b["positions_m"][-1]
        endpoint_difference = np.linalg.norm(x1-x2, axis=1)
    lines += ["", f"相同初始ID，两个相位图样的末位置向量差：均值 {endpoint_difference.mean()*1e6:.6g} µm，最大 {endpoint_difference.max()*1e6:.6g} µm。", ""]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, x, label in zip(axes, [x1, x2], ["Baseline phase", "Shifted phase"], strict=True):
        delta = (x-x0)*1e6
        ax.scatter(delta[:, 1], delta[:, 2], s=12)
        ax.set(xlabel="Lateral displacement Δy (µm)", ylabel="Lateral displacement Δz (µm)", title=label)
        ax.grid(alpha=0.25)
    (root/"results/review_figures").mkdir(exist_ok=True)
    save(fig, root/"results/review_figures/phase_pattern_displacements")
    figures.append(Path("results/review_figures/phase_pattern_displacements.png"))
    lines += ["## 两档数值敏感性", "", "| Scenario | Pressure L2 dx comparison | Gradient L2 | Primary-force L2 | dt endpoint RMS (µm) | dt NN relative difference |", "|---|---:|---:|---:|---:|---:|"]
    aggregate_dt = convergence["results"]["aggregation"]["timestep"]
    lines[-2:-2] = ["结论：尚未证明聚集轨迹收敛。20 ms聚集对照的同ID末位置RMS差为 "
                    f"{aggregate_dt['endpoint_RMS_difference_m']*1e6:.4g} µm，平均位移相差 "
                    f"{100*aggregate_dt['metrics']['mean_displacement_m']['relative_difference']:.4g}%。"
                    "最近邻距离已接近接触直径，其稳定不能证明单泡轨迹准确。", ""]
    for name, data in convergence["results"].items():
        g, dt = data["grid_field"], data["timestep"]
        lines.append(f"| {name} | {g['complex_pressure_relative_L2']:.5g} | {g['complex_gradient_relative_L2']:.5g} | {g['primary_force_relative_L2']:.5g} | {dt['endpoint_RMS_difference_m']*1e6:.5g} | {dt['metrics']['mean_nearest_neighbor_m']['relative_difference']:.5g} |")
        figures.append(Path(data["grid_report"])/"figures/grid_convergence.png")
    figures.append(Path(convergence["run"])/"figures/timestep_convergence.png")
    lines += ["", "| Scenario | Peak coarse / fine (kPa) | Center coarse / fine (kPa) | Mean displacement dt / dt/2 (µm) | NN dt / dt/2 (µm) | Grid particle endpoint RMS (µm) |", "|---|---:|---:|---:|---:|---:|"]
    for name, data in convergence["results"].items():
        g, dt = data["grid_field"], data["timestep"]["metrics"]
        lines.append(f"| {name} | {g['coarse_peak_pressure_Pa']/1e3:.5g} / {g['fine_at_common_points_peak_pressure_Pa']/1e3:.5g} | {g['coarse_center_pressure_amplitude_Pa']/1e3:.5g} / {g['fine_center_pressure_amplitude_Pa']/1e3:.5g} | {dt['mean_displacement_m']['coarse']*1e6:.5g} / {dt['mean_displacement_m']['fine']*1e6:.5g} | {dt['mean_nearest_neighbor_m']['coarse']*1e6:.5g} / {dt['mean_nearest_neighbor_m']['fine']*1e6:.5g} | {data['grid_particle']['endpoint_RMS_difference_m']*1e6:.5g} |")
    lines += ["", "压力峰值和中心值使用相同每元2 kPa源控制；聚集粒子另按记录的scale放大到初始云中位PNP 75 kPa。", ""]
    lines += ["", "这些是0.1/0.08 mm及25/12.5 µs的实测差异，不把基底叠加PASS当作网格无关。几何位移限制子步仍生效，具体计数/最小步长见JSON。", "",
              "## 模型与失败证据", "", "首次聚集试算在逐对接触速度迭代达到上限后停止；部分状态及源文件保留。最终采用相同摩擦无关非穿透约束的整体NNLS投影，30泡接触链测试通过。", "",
              "100泡开发探针因密集接触计算成本于1.941884 ms主动中断，7935状态保留。最终32泡是独立的完整20 ms演示，不能混称100泡结果。", "",
              "近场secondary Bjerknes经过显式数值正则化；接触不代表融合，cluster阈值不代表黏附。Marmottant响应尚未实现，不伪造壳参数。", ""]
    effective_path = result_dir/"EFFECTIVE_STEP_STATE.json"
    if effective_path.exists():
        effective = json.loads(effective_path.read_text())
        assert effective["status"] == "COMPLETE"
        lines += ["### 实际子步补充检查", "",
                  "迁移的25/12.5 µs名义对照均产生800子步，不能独立证明时间收敛。补跑更小有效步长；聚集补充检查只覆盖0–2 ms形成阶段。", "",
                  "| Scenario | Horizon (ms) | Coarse / finer substeps | Endpoint RMS difference (µm) | NN coarse / finer (µm) |", "|---|---:|---:|---:|---:|"]
        for name, data in effective["results"].items():
            lines.append(f"| {name} | {data['horizon_s']*1e3:g} | {data['coarse_steps_to_horizon']} / {data['refined_steps_to_horizon']} | {data['endpoint_RMS_difference_m']*1e6:.6g} | {data['coarse_mean_NN_m']*1e6:.6g} / {data['refined_mean_NN_m']*1e6:.6g} |")
        lines += ["", "此补充检查不改变原始主算例，也不推翻20 ms密集接触轨迹尚未收敛的结论。", ""]
    for case in demo["scenarios"]:
        lines.append(f"- {case['scenario']}/{case['pattern']}: max predicted |ξ|/R0 = {case['summary']['maximum_linear_radial_amplitude_over_R0']:.5g}；该数值不认证线性替代模型的真实生理有效性。")
    lines += ["", "## 实际计算资源", "", "| f (MHz) | Grid | Nt | ROI sensors | Solver seconds | Sampled device VRAM MiB |", "|---|---|---:|---:|---:|---:|"]
    for meta in resources:
        lines.append(f"| {meta['frequency_Hz']/1e6:g} | {meta['grid_dimensions']} | {meta['Nt']} | {meta['sensor_count']} | {meta['wall_clock_s']:.3f} | {meta['sampled_device_peak_VRAM_MiB']} |")
    lines += ["", "显存是整设备0.5 s采样峰值，不是进程私有分配上界。生产声场全部实际在RTX4090执行。", "",
              "外置PML后FFT网格为288×160×160；细档为360×200×200。这是由metadata物理网格与PML宽度计算的实际求解网格。", "",
              "## 文件入口", "", "- 配置、接口、安装与限制：../README.md", "- 测试证据：TEST_RESULTS.txt", "- 全帧数据检查：RESULT_INTEGRITY.json", "- 源码与运行输出身份：DELIVERY_MANIFEST.json", ""]
    (root/"docs/VALIDATION_REPORT.md").write_text("\n".join(lines))
    # Self-contained HTML uses exact tables and local artifact links; no external web dependencies.
    body = ["<!doctype html><meta charset='utf-8'><title>Microbubble acoustics review</title><style>body{font:17px system-ui;max-width:1250px;margin:36px auto;padding:0 24px;color:#172638;background:#fafbfc}img{width:100%;background:white}table{border-collapse:collapse;background:white;font-size:14px;width:100%;margin:20px 0}th,td{border:1px solid #ddd;padding:9px;text-align:left}th{background:#e5edf4}p{line-height:1.65}a{color:#125c91}section{margin:32px 0}</style>",
            "<h1>PZT → 3D acoustics → microbubble transport</h1><p>实际结果与数值敏感性。Linear surrogate · one-way coupling · no coalescence</p>",
            "<p><a href='README.md'>README</a> · <a href='docs/VALIDATION_REPORT.md'>中文验证报告</a> · <a href='results/DEMO_CAMPAIGN.json'>Run index</a> · <a href='results/CONVERGENCE_STATE.json'>Convergence</a></p>",
            render_text(lines)]
    for path in figures:
        if not (root/path).is_file():
            raise FileNotFoundError(path)
        body.append(f"<section><h2>{html.escape(path.parent.parent.name)} / {html.escape(path.stem)}</h2><a href='{path.with_suffix('.pdf')}'>PDF</a><img src='{path}' loading='lazy'></section>")
    body.append("<h2>ParaView</h2><ul>")
    for case in demo["scenarios"]:
        path = Path(case["run"])/"vtk/particles.pvd"
        body.append(f"<li><a href='{path}'>{case['scenario']} / {case['pattern']} particles.pvd</a></li>")
    body.append("</ul>")
    (root/"OPEN_RESULTS.html").write_text("\n".join(body))
    write_json(root/"docs/REVIEW_FIGURES.json", {str(p): hashlib.file_digest((root/p).open("rb"), "sha256").hexdigest() for p in figures})


def render_text(lines):
    """Render our small generated report format without network or dependencies."""
    result, in_table = [], False
    for line in lines:
        if line.startswith("|"):
            if set(line.replace("|", "").replace(":", "").strip()) <= {"-"}:
                continue
            tag = "td" if in_table else "th"
            if not in_table:
                result.append("<table>")
                in_table = True
            result.append("<tr>"+"".join(f"<{tag}>{html.escape(cell.strip())}</{tag}>" for cell in line.strip("|").split("|"))+"</tr>")
        else:
            if in_table:
                result.append("</table>")
                in_table = False
            if line.startswith("#"):
                level = len(line)-len(line.lstrip("#"))
                result.append(f"<h{level}>"+html.escape(line.lstrip("# "))+f"</h{level}>")
            elif line:
                result.append("<p>"+html.escape(line)+"</p>")
    if in_table:
        result.append("</table>")
    return "\n".join(result)


if __name__ == "__main__":
    main()
