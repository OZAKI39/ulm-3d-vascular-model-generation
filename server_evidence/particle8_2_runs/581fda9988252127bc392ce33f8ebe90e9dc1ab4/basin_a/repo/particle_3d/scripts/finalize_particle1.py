#!/usr/bin/env python3
"""Offline P1 final checks and Chinese review; never grants manual approval."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

PACKAGE = Path(__file__).resolve().parents[1]
REPO = PACKAGE.parent
REPORT = PACKAGE / "reports/particle1"
DATA = REPORT / "data"
LOGS = REPORT / "logs"
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen, sha256, FEM_BRANCH
from particle_3d.particle1_audit import check_dependency, P0_COMMIT, FROZEN_COMMIT, P1_BRANCH
from particle_3d.sonovue_adapter import read_sonovue


def git(*args, cwd=REPO):
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def read_data(name):
    return json.loads((DATA / name).read_text())


def run(name, args, cwd):
    print(f"Checking {name}", flush=True)
    result = subprocess.run(args, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log = LOGS / f"{name}.log"
    log.write_text(result.stdout)
    return dict(name=name, command=args, cwd=str(cwd), returncode=result.returncode,
                log=str(log.relative_to(REPO)))


def junit(path):
    cases = list(ET.parse(path).iter("testcase"))
    failures = sum(c.find("failure") is not None or c.find("error") is not None for c in cases)
    skipped = sum(c.find("skipped") is not None for c in cases)
    return dict(total=len(cases), passed=len(cases) - failures - skipped, failed=failures, skipped=skipped)


def review(v, initial, uniform, rotation, comparison):
    status = v["automated_checks"]
    case_table = "\n".join(
        f"| {c['case_id']} | {c['validation_dt_s']:.16g} | {c['row_count']} | {c['exit_boundary']} | "
        f"{c['exit_time_s']:.12f} | {c['trajectory_length_m']:.12e} | {c['max_trial_step_length_m']:.6e} |"
        for c in v["real_trajectory_cases"])
    diffs = comparison["exit_time_adjacent_differences_s"]
    text = f"""# Particle-1 人工审核报告

## 1. 这一阶段做了什么

这一阶段第一次让一个球形微泡在冻结的三维 FEM 流场里真正运动。
微泡每到一个位置，就向 Particle-0 查询当地的流体速度和旋转信息，
再用显式 Euler 更新位置。这里只有一个微泡；三个时间步结果是同一微泡的三次独立验证重放。
每一步开发都永久保留测试、PNG、CSV/JSON 源数据和简短中文说明。
本轮在 WSL CPU 完成，没有运行 CFD，也没有进入 Particle-2。

开发目录：{REPO}

Particle-1 branch：{v['git_branch']}

被测试的实现 commit：{v['git_commit']}

机器记录中的 git_commit 指向已提交的源码版本；随后可有仅保存最终报告和日志的证据提交。
代码、测试、脚本、图片和数据的逐文件 SHA 保存在 PARTICLE1_VALIDATION.json。

## 2. 沿用了什么

| 依赖 | 本轮处理 |
|---|---|
| Particle-0 implementation | {P0_COMMIT}，源码、数学、容差和原测试保持字节一致 |
| P0 人工审核 | 用户已审核现有报告、JSON 和 8 图；记录 PASS / USER_CHAT_REVIEW |
| P0 evidence-only commit | {v['particle0_acceptance_evidence_commit']}，只更新原报告和 JSON |
| Frozen FEM base | {FROZEN_COMMIT}，只读 |
| Frozen branch | {FEM_BRANCH}，未改动 |
| mesh SHA256 | {v['mesh_sha256']} |
| flow SHA256 | {v['flow_sha256']} |
| SonoVue histogram SHA256 | {v['sonovue_histogram_sha256']} |
| SonoVue version | {v['sonovue_sampler_version']} |

P0 人工审核日期和授权证据完整保存在本次 JSON 的 particle0_manual_review_evidence，
包括 reviewer=USER_CHAT_REVIEW、blocking issue=NONE、不要求增强 P0 测试、允许开始 P1。
本次没有增加全节点回采样等 P0 增强项目，只重复原有 46 项回归测试。

SonoVue 原目录为 {v['sonovue_root']}，全部 {v['sonovue_manifest_file_count']} 个 manifest 文件重新核验一致。
适配层只调用原函数抽取一个样本，并完成 µm→m→半径的换算，没有复制或修改 sampler。
seed={v['single_mb_seed']}，直径 {v['single_mb_diameter_um']:.15g} µm，
半径 {v['single_mb_radius_m']:.15g} m；formal_simulation_population=false。
样本 CSV 旁边有绑定文件 SHA 的 metadata。
本次环境 NumPy {v['environment']['packages']['numpy']}、Python {v['environment']['python_version']}，
核对的是同一环境中原始 sampler 的精确输出，不声称不同 NumPy 版本逐位一致。

旧 2D 参考代码只用于了解 field 与 state 分开的组织方式，没有复制二维坐标、µm 核心单位、
旧网格 cache、RBC drift 或二维边界行为。FEM 没有重新求解。

## 3. 单微泡现在有哪些状态

| 状态 | 含义 | 单位 / 形状 |
|---|---|---|
| particle_id | 当前微泡的编号 | 非负整数 |
| position_m | 中心的位置 | m，float64，3 分量 |
| radius_m | 微泡半径 | m，有限正 float64 |
| velocity_m_s | 平移速度 | m/s，float64，3 分量 |
| angular_velocity_s_inv | 旋转速度 | 1/s，float64，3 分量 |

状态拒绝非有限输入，数组不与外部输入共享可写内存。
返回新状态，不原地修改旧状态。没有质量、线加速度、角加速度、quaternion 或惯性积分。
球形外观不随朝向改变，因此当前只保存旋转速度。

位置更新为 x_new=x_old+dt*u(x_old)，dt 由调用者显式传入，必须有限且大于零，没有默认值。
更新后在新位置重新查询一次，只刷新输出 V 和 Ω，不再次移动位置。
因此每行速度与该行位置一致，位置积分仍是显式 Euler。
普通 motion API 的域外起点/终点会报错；逐段边界检查由独立的验证程序完成。

## 4. 为什么 V = u_inf

当前假设小球很快达到流体作用下的平衡速度。保留的阻力公式是：

F_hydro = −6πμa(V−u_inf)

μ 是动力黏度，a 是半径，V 是粒子速度，u_inf 是背景流速。
因为当前没有其它力，平衡要求 F_hydro=0，所以 V=u_inf。
这一步没有求解 m dV/dt=F，也没有偷偷加入重力、浮力或其它未批准的物理。
半径仍保存在状态和来源记录中，但在本阶段无其它力的平移速度公式中会消去。

## 5. 为什么 Omega = 0.5*vorticity

Particle-0 给出的 vorticity 是 curl(u)。在整体刚体旋转的人工流场中，
curl(u) 恰好是流体实际旋转速度的两倍，所以球形微泡取一半：

Ω_particle = 0.5 curl(u_inf)

测试专门使用三个非零旋转分量；若误写成 Ω=vorticity，永久回归测试会失败。

## 6. 自动测试

Particle-0：**{v['particle0_regression_tests']['passed']} passed / {v['particle0_regression_tests']['failed']} failed**。
Particle-1：**{v['particle1_tests']['passed']} passed / {v['particle1_tests']['failed']} failed**。
Frozen handoff：**{v['frozen_tests']['passed']} passed / {v['frozen_tests']['failed']} failed**。
上述测试 skipped 均为 {v['particle1_tests']['skipped']}。

| 检查 | 结果 | 误差 / 证据 | test file |
|---|---|---|---|
| P0 依赖与范围 | {status} | P0 原 26 个 Python 文件 SHA 一致；原 46 项测试不变 | test_particle0_dependency_locked.py、test_particle1_scope.py |
| 状态、单位和输入保护 | {status} | float64；3 分量；拒绝别名改写和非正半径 | test_microbubble_state.py、test_microbubble_units.py |
| SonoVue 原函数与只读 | {status} | seed 重现；原函数同环境精确匹配；24 个文件 SHA 一致 | test_sonovue_single_mb_adapter.py、test_sonovue_sampler_unchanged.py |
| Stokes 平衡与零旋转 | {status} | V−u 最大误差 0 m/s；Ω=0 | test_stokes_following_uniform_flow.py、test_zero_vorticity_zero_rotation.py |
| 均匀流解析直线 | {status} | 最大位置分量误差 {uniform['max_position_error_m']:.16e} m | test_uniform_flow_trajectory.py |
| 纯旋转及 1/2 | {status} | 最大旋转分量误差 {rotation['max_angular_velocity_error_s_inv']:.3e} 1/s | test_pure_rotation_angular_velocity.py、test_vorticity_half_factor.py |
| 显式一步与无惯性状态 | {status} | 使用旧点速度；输入保持不变；非法 dt 拒绝 | test_single_mb_step_api.py、test_step_does_not_modify_input_state.py、test_invalid_dt.py、test_no_inertial_acceleration_state.py |
| 确定性真实起点 | {status} | 全部 {initial['candidate_count']} 个候选；tetra_id={initial['initial_tetra_id']} | test_real_single_mb_initialization.py |
| 真实 V=u | {status} | 全部 4016 行重新查询 P0；最大误差 {v['real_max_velocity_relation_error_m_s']:.3e} m/s | test_real_single_mb_velocity_matches_fem.py |
| 真实 Ω=curl/2 | {status} | 全部 4016 行各分量；最大误差 {v['real_max_angular_relation_error_s_inv']:.3e} 1/s | test_real_single_mb_rotation_matches_half_vorticity.py |
| 逐段边界与出口 | {status} | 4013 条线段；无 WALL/INLET；均为 OUTLET_02 | test_boundary_segment_classification.py、test_real_single_mb_no_wall_crossing.py、test_real_single_mb_outlet_classification.py |
| 时间步对比、图与源数据 | {status} | 3 个预先确定的减半 dt；8 图及数据齐全；可重新作图 | test_validation_timestep_comparison.py、test_particle1_review_artifacts.py、test_particle1_plot_generation.py |

测试文件均在 [particle1 tests](../../tests/particle1/)。Frozen integrity={v['frozen_integrity']}：
当前 worktree 的 validate_frozen.py 和 verify_manifest.py 全部通过，4602 个 payload 文件和 31 个科学 SHA 保持一致。
原有 18 项 handoff pytest 锁定 frozen branch 名称，因此在原干净 frozen checkout 执行，没有修改历史测试。

均匀流预设舍入界是 16*(step_count+1)*eps*max(abs(x0)+abs(u)*duration)，
因为累加次数更多可能产生更多浮点舍入；这里没有把极小舍入误差误称为 Euler 截断误差。
纯旋转预设 atol=rtol=256*eps={rotation['atol']:.16e}，没有根据结果放宽阈值。

边界检查覆盖全部五个原始表面，使用 VTK 搜索候选后计算线段与三角形首次交点。
误差预算由 float64 eps、实际坐标尺度及三角形条件数决定。
同一边角同时命中时优先 WALL，再 INLET，再出口。专门的人工穿壁用例确认首次 WALL 就记录 FAIL 并停止，
保存原始试探端点和交点，没有反弹、推回、继续走或重调时间步。

真实起点是入口相邻 tetra 的中心，选择中心向内流速最大者，精确并列取最小 canonical tetra_id。
位置为 {initial['initial_position_m']} m，速度为 {initial['initial_velocity_m_s']} m/s。
最短局部边长为 {initial['local_shortest_edge_m']:.12e} m，起点速率为 {initial['initial_speed_m_s']:.12e} m/s。
首次运行前按 dt0=0.25*最短边长/起点速率确定以下三个验证步长，没有根据是否穿壁试调参数。

| case | VALIDATION_ONLY dt (s) | 保存行数 | 首次边界 | 出口时间 (s) | 轨迹长度 (m) | 最大试探步长 (m) |
|---|---|---|---|---|---|---|
{case_table}

三条轨迹的全部 active 与终止状态有限；每条线段重新检查，没有 WALL 或 INLET 穿出。
出口时间按 Euler 线段内的交点比例计算，事件发生后立即停止。
最大一步位移只作信息记录，没有用任意跳跃阈值宣布人工审核通过。
相邻出口时间差为 {diffs[0]:.12e}、{diffs[1]:.12e} s，差值略减小，
但不能据此宣布一阶收敛、验证真实连续流的误差或确定正式时间步。

完整命令和日志见 [final_check_commands.json](logs/final_check_commands.json)，
[P0 测试](logs/particle0_pytest.log)、[P1 测试](logs/particle1_pytest.log) 和 [Frozen 测试](logs/frozen_pytest.log)。
开发过程曾发生作图 JSON 的 NumPy 布尔值序列化错误，使后续图未生成，相关测试报 7 failed / 105 passed。
已修复序列化并重新生成全部图；原失败日志保留在 combined_development_tests.log，
修复后完整回归为 113 passed / 0 failed，最终独立回归结果如上。没有改变科学模型或容差来修复该问题。

## 7. 人工审核图片

![00 范围与来源](figures/00_particle1_scope_and_provenance.png)

- 应该看什么：P0 已接受、FEM 只读、只有一个球，依赖和 SHA 有来源。
- 实际看到了什么：流程从冻结场查询到 V=u、Ω=curl/2，再到位置更新。
- 有没有异常：自动来源核验通过；不代表用户已完成本轮图审。

![01 尺寸来源](figures/01_single_mb_size_provenance.png)

- 应该看什么：原 histogram / CDF 上 THIS PARTICLE 的位置、直径、半径和 seed。
- 实际看到了什么：直径 {v['single_mb_diameter_um']:.9f} µm，半径 {v['single_mb_radius_m'] * 1e6:.9f} µm，仅一个验证样本。
- 有没有异常：单位换算和原函数一致，没有重拟合分布或正式粒子群。

![02 均匀流轨迹](figures/02_uniform_flow_trajectory.png)

- 应该看什么：三套 dt 的点是否在解析直线上，误差轴的数量级。
- 实际看到了什么：轨迹重合，最大位置误差约 {uniform['max_position_error_m']:.3e} m，速度误差为零。
- 有没有异常：均在预设舍入界内；较多步数产生较多舍入是预期现象。

![03 纯旋转](figures/03_pure_rotation_validation.png)

- 应该看什么：左侧流向，以及右侧解析与微泡角速度的三个分量。
- 实际看到了什么：49 个位置的角速度误差为零，Ω=0.5*vorticity。
- 有没有异常：自动检查未发现漏乘 1/2 的错误；箭头长度仅为显示比例。

![04 单步推进](figures/04_single_step_dataflow.png)

- 应该看什么：从旧位置沿旧点速度指向新位置，旁边的中文步骤。
- 实际看到了什么：位移等于 dt*u_old，新状态查场刷新 V/Ω，输入状态不变。
- 有没有异常：自动检查通过；图上的 marker 只作示意，不代表真实半径。

![05 真实血管轨迹](figures/05_real_vessel_single_mb_trajectory.png)

- 应该看什么：真实表面、入口和三个出口，以及微泡起终点、速度箭头、是否出现明显一步跳跃。
- 实际看到了什么：展示最细验证重放，轨迹由 INLET 附近走向 OUTLET_02；方向箭头叠加显示以免被透明表面遮挡。
- 有没有异常：逐段自动检查无 WALL/INLET 穿出；突然跳跃的最终视觉审核仍待用户。marker 已标明视觉放大，不证明有限半径壁面间隙。

![06 真实场诊断](figures/06_real_vessel_single_mb_diagnostics.png)

- 应该看什么：四个独立单位的曲线，旋转速度与 vorticity 的一半关系。
- 实际看到了什么：速度和压力三套结果接近；旋转曲线细碎变化明显，Ω 始终满足逐分量的一半关系。
- 有没有异常：数值全部有限。导数跨 tetra 跳变来自原始离散场，未平滑，也未当作连续流体真解。

![07 验证步长对比](figures/07_validation_timestep_comparison.png)

- 应该看什么：轨迹重合程度、出口分类、出口时间和路径长度随验证 dt 的变化。
- 实际看到了什么：三条轨迹接近、出口相同；出口时间差略减小，最细路径只是比较参考，不是解析真值。
- 有没有异常：没有自动判定异常边界事件；趋势只作开发验证，不确定 production dt，也不是 FEM timestep study。

## 8. 当前限制

- 只有一个 MB，且 MB 是球；没有第二个微泡或 RBC。
- 没有质量、线加速度、角加速度、惯性积分或 NVE。
- 没有 wall force、wall correction、wall reaction、contact、lubrication、adhesion 或 molecular binding。
- 没有 particle-particle hydrodynamics、MB-MB / RBC-MB interaction、collision 或多体阻力系统。
- 没有 Brownian、acoustic force、gravity、buoyancy、lift、added mass 或 Basset history。
- 没有 LAMMPS、正式 continuous injection、粒子删除生命周期、flux schedule、浓度或 hematocrit。
- 不穿 WALL 只说明本次中心轨迹在三个验证步长下未穿过表面，不证明有限半径球的 clearance。
- validation dt 不是 production dt；production_particle_timestep_frozen=false。API 仍要求调用者传入 dt。
- FEM 仍然 frozen，未运行新的 CFD、网格收敛或 FEM 时间步研究。
- Particle-0 gradient 仍然是 tetra 内常数，跨单元可跳变；未改变数学、容差或性能实现。
- 一次单微泡验证不等于正式 simulation population，不代表全部入口位置或全部尺寸已验证。
- 轨迹视觉上是否有突然跳跃仍需用户审核，不由任意硬阈值代替。

## 9. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

PRODUCTION_PARTICLE_TIMESTEP_FROZEN = false

FEM_MODIFIED = NO

PARTICLE2_STARTED = false

Particle-1 在此停止。只做本地提交，没有 push 或 merge main。
全部新增文件见 [FILES_CHANGED.txt](FILES_CHANGED.txt)，接口和复现命令见 [PARTICLE1_README.md](../../PARTICLE1_README.md)。
"""
    (REPORT / "PARTICLE1_REVIEW.md").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff-root", type=Path, required=True)
    args = parser.parse_args()
    handoff = args.handoff_root.resolve()
    fem = REPO / "formal_3D_flow_solver/FEM_SimVascular"
    LOGS.mkdir(parents=True, exist_ok=True)
    source_paths = ["particle_3d/src", "particle_3d/tests", "particle_3d/scripts",
                    "particle_3d/PARTICLE1_README.md", "particle_3d/pyproject.toml"]
    if git("status", "--porcelain", "--", *source_paths):
        raise RuntimeError("Commit implementation first; git_commit must identify the tested sources")
    if git("branch", "--show-current") != P1_BRANCH:
        raise RuntimeError("Finalization requires the dedicated Particle-1 branch")
    dependency = check_dependency(REPO)
    if git("branch", "--show-current", cwd=handoff) != FEM_BRANCH:
        raise RuntimeError("Handoff tests require the unchanged original frozen branch checkout")
    if git("rev-parse", "HEAD", cwd=handoff) != FROZEN_COMMIT or git("status", "--porcelain", cwd=handoff):
        raise RuntimeError("Original frozen checkout changed; do not repair it here")
    summary, _, _, _ = read_frozen(fem)
    scope = read_data("00_particle1_scope_and_provenance.json")
    contract, _, inventory = read_sonovue(scope["sonovue_root"])
    sonovue_unchanged = (inventory == scope["sonovue_inventory"] and
                        sha256(Path(scope["sonovue_root"]) / "SHA256SUMS") == scope["sonovue_sha256sums_sha256"])
    commands = []
    for name, script in [("frozen_validate", "validate_frozen.py"), ("frozen_manifest", "verify_manifest.py")]:
        commands.append(run(name, [sys.executable, "-B", "scripts/fem_freeze_sync/" + script], fem))
        if commands[-1]["returncode"]:
            write_json(LOGS / "final_check_commands.json", commands)
            raise RuntimeError("Frozen integrity blocker; preserve inputs and inspect logs")
    frozen_tests = ["tests/test_fem_sync_manifest.py"] + [
        str(p.relative_to(handoff)) for p in sorted((handoff / "tests").glob("test_frozen_*.py"))
    ] + ["tests/test_particle_handoff_contract.py", "tests/test_no_particle_implementation_yet.py", "tests/test_main_history_preserved.py"]
    suites = [("frozen", handoff, frozen_tests),
              ("particle0", REPO, ["particle_3d/tests/particle0"]),
              ("particle1", REPO, ["particle_3d/tests/particle1"])]
    tests = {}
    for name, cwd, paths in suites:
        xml = LOGS / f"{name}_junit.xml"
        commands.append(run(name + "_pytest", [sys.executable, "-B", "-m", "pytest", "-q", "-p", "no:cacheprovider",
                                               f"--junitxml={xml}", *paths], cwd))
        tests[name] = junit(xml)
    write_json(LOGS / "final_check_commands.json", commands)
    # Post-test checks include scientific contents, not just Git status.
    check_dependency(REPO)
    _, _, final_inventory = read_sonovue(scope["sonovue_root"])
    sonovue_unchanged = sonovue_unchanged and final_inventory == inventory
    frozen_diff = git("diff", FROZEN_COMMIT, "--", "formal_3D_flow_solver/FEM_SimVascular")
    p0_diff = git("diff", P0_COMMIT, "--", *dependency["source_sha256"])
    uniform = read_data("02_uniform_metrics.json")
    rotation = read_data("03_rotation_metrics.json")
    initial = read_data("05_real_initialization.json")
    cases = read_data("05_real_trajectory_sweep.json")
    size = read_data("01_single_mb_size_provenance.json")
    comparison = read_data("07_comparison_metrics.json")
    figures = sorted(str(p.relative_to(REPO)) for p in (REPORT / "figures").glob("*.png"))
    test_files = sorted(str(p.relative_to(REPO)) for p in (PACKAGE / "tests/particle1").glob("test_*.py"))
    passed = (all(c["returncode"] == 0 for c in commands) and
              all(t["failed"] == t["skipped"] == 0 for t in tests.values()) and
              tests["particle0"]["passed"] == 46 and tests["frozen"]["passed"] == 18 and tests["particle1"]["passed"] > 0 and
              not frozen_diff and not p0_diff and sonovue_unchanged and len(figures) == 8 and
              uniform["passed"] and rotation["passed"] and len(cases) == 3 and all(c["passed"] for c in cases))
    sources = sorted(p for top in [PACKAGE / "src", PACKAGE / "tests", PACKAGE / "scripts"] for p in top.rglob("*.py"))
    sources.append(PACKAGE / "PARTICLE1_README.md")
    acceptance_commit = git("log", "-1", "--format=%H", "--", "particle_3d/reports/particle0/PARTICLE0_VALIDATION.json")
    v = dict(
        stage="Particle-1", automated_checks="PASS" if passed else "FAIL",
        git_branch=git("branch", "--show-current"), git_commit=git("rev-parse", "HEAD"),
        git_commit_semantics="Tested committed implementation. A later evidence-only commit may contain this record.",
        particle0_dependency_commit=P0_COMMIT, particle0_regression_tests=tests["particle0"],
        particle0_acceptance_evidence_commit=acceptance_commit,
        particle0_manual_review_evidence=dependency["manual_review_evidence"],
        particle0_source_unchanged=not bool(p0_diff), particle0_enhanced_tests_added=False,
        particle1_tests=tests["particle1"], frozen_tests=tests["frozen"],
        frozen_base_commit=FROZEN_COMMIT, frozen_branch=FEM_BRANCH,
        frozen_integrity="PASS" if not frozen_diff and all(c["returncode"] == 0 for c in commands[:3]) else "FAIL",
        frozen_payload_modified=bool(frozen_diff),
        mesh_sha256=summary["mesh_sha256"], flow_sha256=summary["flow_sha256"],
        mesh_path=summary["mesh_path"], flow_path=summary["flow_path"],
        sonovue_root=scope["sonovue_root"], sonovue_histogram_sha256=contract["histogram_sha256"],
        sonovue_sampler_version=contract["sampler_version"], sonovue_sampler_unchanged=sonovue_unchanged,
        sonovue_manifest_file_count=len(inventory), sonovue_inventory=inventory,
        single_mb_seed=size["seed"], single_mb_diameter_um=size["diameter_um"],
        single_mb_diameter_m=size["diameter_m"], single_mb_radius_m=size["radius_m"],
        formal_simulation_population=False,
        uniform_flow_test="PASS" if uniform["passed"] else "FAIL",
        uniform_flow_max_velocity_error=uniform["max_velocity_error_m_s"],
        uniform_flow_max_position_error=uniform["max_position_error_m"],
        uniform_flow_error_units={"velocity": "m/s", "position": "m"},
        pure_rotation_test="PASS" if rotation["passed"] else "FAIL",
        pure_rotation_max_angular_velocity_error=rotation["max_angular_velocity_error_s_inv"],
        pure_rotation_error_unit="s^-1", pure_rotation_atol=rotation["atol"], pure_rotation_rtol=rotation["rtol"],
        real_trajectory_start=initial["initial_position_m"], real_trajectory_start_units="m",
        real_initial_tetra_id=initial["initial_tetra_id"], real_initialization=initial,
        validation_timesteps_s=initial["validation_timesteps_s"], timestep_role="VALIDATION_ONLY",
        real_trajectory_finite=all(c["active_and_terminal_finite"] for c in cases),
        wall_crossing=any(c["wall_crossing"] for c in cases), inlet_crossing=any(c["inlet_crossing"] for c in cases),
        exit_boundary=[c["exit_boundary"] for c in cases], exit_time_s=[c["exit_time_s"] for c in cases],
        case_order="All list fields follow validation_timesteps_s from coarse to fine; no production case selected.",
        real_trajectory_cases=cases, every_segment_checked=all(c["every_segment_checked"] for c in cases),
        checked_segment_count=sum(c["checked_segment_count"] for c in cases),
        real_max_velocity_relation_error_m_s=max(c["max_velocity_relation_error_m_s"] for c in cases),
        real_max_angular_relation_error_s_inv=max(c["max_angular_relation_error_s_inv"] for c in cases),
        velocity_relation_check="PASS" if all(c["max_velocity_relation_error_m_s"] == 0 for c in cases) else "FAIL",
        angular_velocity_relation_check="PASS" if all(c["max_angular_relation_error_s_inv"] == 0 for c in cases) else "FAIL",
        validation_timestep_comparison=comparison, figures=figures, test_files=test_files,
        source_sha256={str(p.relative_to(REPO)): sha256(p) for p in sources},
        figure_sha256={p: sha256(REPO / p) for p in figures},
        data_sha256={str(p.relative_to(REPO)): sha256(p) for p in sorted(DATA.iterdir()) if p.is_file()},
        log_sha256={str(p.relative_to(REPO)): sha256(p) for p in sorted(LOGS.iterdir()) if p.is_file() and p.suffix in [".log", ".xml", ".json"]},
        production_particle_timestep_frozen=False, manual_visual_review="PENDING_USER_REVIEW",
        visual_step_jump_review="PENDING_USER_REVIEW", particle2_started=False, no_cfd_executed=True,
        no_gpu_server_started=True, pushed=False, merged_main=False,
        recorded_at_utc=datetime.now(timezone.utc).isoformat(),
        environment=dict(python=sys.version, python_version=platform.python_version(), executable=sys.executable,
                         platform=platform.platform(), packages={n: version(n) for n in ["numpy", "vtk", "pyvista", "matplotlib", "pytest", "Pillow"]}))
    write_json(REPORT / "PARTICLE1_VALIDATION.json", v)
    review(v, initial, uniform, rotation, comparison)
    (LOGS / "TEST_SUMMARY.md").write_text(
        f"# Particle-1 最终检查\n\nAUTOMATED_CHECKS = {v['automated_checks']}\n\n"
        f"P0：{tests['particle0']['passed']} passed / {tests['particle0']['failed']} failed。\n\n"
        f"P1：{tests['particle1']['passed']} passed / {tests['particle1']['failed']} failed。\n\n"
        f"Frozen handoff：{tests['frozen']['passed']} passed / {tests['frozen']['failed']} failed。\n\n"
        f"Frozen integrity：{v['frozen_integrity']}。原 P0 / SonoVue / FEM 保持只读。\n\n"
        "开发失败日志保留；以最终独立 suite 日志与 JUnit 为最终验收。\n\n"
        "MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW；没有冻结 production dt，没有 CFD，没有 Particle-2。\n")
    changed = set(git("diff", "--name-only", acceptance_commit).splitlines())
    changed.update(git("ls-files", "--others", "--exclude-standard", "particle_3d").splitlines())
    changed.add("particle_3d/reports/particle1/FILES_CHANGED.txt")
    (REPORT / "FILES_CHANGED.txt").write_text("\n".join(sorted(changed)) + "\n")
    print(json.dumps({k: v[k] for k in ["automated_checks", "git_commit", "particle0_regression_tests",
                                       "particle1_tests", "frozen_tests", "manual_visual_review"]}), flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
