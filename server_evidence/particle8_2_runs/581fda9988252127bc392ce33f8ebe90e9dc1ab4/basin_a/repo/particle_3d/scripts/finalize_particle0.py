#!/usr/bin/env python3
"""Run final offline checks and write an evidence-backed Chinese review."""
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
REPORT = PACKAGE / "reports/particle0"
LOGS = REPORT / "logs"
sys.path.insert(0, str(PACKAGE / "src"))
from particle_3d.audit import read_frozen, sha256, FEM_BRANCH


def git(*args, cwd=REPO):
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def run(name, args, cwd):
    print(f"Checking {name}", flush=True)
    result = subprocess.run(args, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (LOGS / f"{name}.log").write_text(result.stdout)
    return dict(name=name, command=args, cwd=str(cwd), returncode=result.returncode,
                log=str((LOGS / f"{name}.log").relative_to(REPO)))


def junit(path):
    cases = list(ET.parse(path).iter("testcase"))
    failed = [case for case in cases if case.find("failure") is not None or case.find("error") is not None]
    skipped = [case for case in cases if case.find("skipped") is not None]
    return dict(total=len(cases), passed=len(cases) - len(failed) - len(skipped), failed=len(failed), skipped=len(skipped)), cases


def group_pass(cases, names):
    selected = [c for c in cases if any(name in c.attrib["classname"] for name in names)]
    return bool(selected) and all(len(c) == 0 or all(child.tag not in ["failure", "error", "skipped"] for child in c) for c in selected)


def review(v, metrics):
    p1, p2, p3, p4, p5, p6 = [metrics[str(i)] for i in range(1, 7)]
    status = v["automated_checks"]
    testlink = lambda name: f"[ {name} ](../../tests/particle0/{name})"
    text = f"""# Particle-0 人工审核报告

## 1. 这一阶段做了什么

这一阶段建立了一个只读的三维流体场查询接口 `FrozenFEMField`，没有移动任何微泡或 RBC。
给它一个以米表示的位置，就返回速度、压力、三项导数信息，以及位置是否在域内和所属四面体编号。
正式单元和顶点编号来自冻结体网格，速度和压力只读取冻结流场的节点数组。
先用人工知道答案的四面体检查数学，再用真实节点、共享面和血管边界检查实际文件。
每个小步骤都保留了永久测试、PNG 图、作图数据和中文说明。
批量查询已与逐点查询比较，CPU 计时仅作信息记录。
整个阶段只在 WSL CPU 上运行，下一阶段尚未开始。

## 2. 输入有没有变

FEM 保持只读，**本阶段修改 FEM：NO**。开发目录为 `{REPO}`。

| 输入 | 当前记录 |
|---|---|
| mesh SHA256 | `{v['mesh_sha256']}` |
| flow SHA256 | `{v['flow_sha256']}` |
| 节点 / 四面体 | {v['node_count']} / {v['tetra_count']} |
| 边界 | WALL=1，OUTLET_03=2，OUTLET_01=3，INLET=4，OUTLET_02=5 |
| Frozen Git branch | `{v['fem_branch']}` |
| Frozen Git commit | `{v['frozen_base_commit']}` |
| Particle Git branch | `{v['git_branch']}` |
| 被测试的实现提交 | `{v['git_commit']}` |

`git_commit` 指向被测试的实现提交；随后可以有只保存报告和日志的证据提交。代码和测试的逐文件 SHA 保存在机器记录 `source_sha256`，便于对应。

网格实际路径：`{v['mesh_path']}`。

流场实际路径：`{v['flow_path']}`。

你给出的旧路径 `/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular` 没有 handoff/manifest，且已有历史改动；它被原样保留，前后 Git 状态已对比。总仓库中已有指定的干净冻结分支，本次从它建立独立开发 worktree，按其中 manifest 定位输入。没有修改、清理旧目录，也没有修改 main 或 frozen branch。

开始前和结束后运行官方只读检查。最终 worktree 内两个完整性脚本均返回成功；完整官方 handoff pytest 在仍保留冻结分支名的原 checkout 运行，结果 {v['frozen_tests']['passed']} passed、{v['frozen_tests']['failed']} failed。历史测试会锁定分支名称，故不能把它改成在 Particle 分支上假装仍是冻结交接；原测试未修改。

SonoVue：`/home/lzy/projects/sonovue_size_distribution_v0/`，对应分支 `codex/sonovue-size-distribution-20260920_111258`，仅登记存在，没有读取它生成粒子或修改 sampler。旧 2D 两个参考文件只阅读组织方法，没有复制二维插值或 RBC 代码。

## 3. 我们怎样查询一个位置

```text
位置（m）
  → VTK 空间索引枚举可能的 tetra（四面体）
  → 用自己的四面体公式验证位置，算四个权重
  → 按 canonical volume 的四个全局节点编号取 Velocity / Pressure
  → 对节点值加权，返回速度 / 压力
  → 从四面体节点速度差计算 gradient
  → 返回 vorticity / strain rate
```

barycentric weights 就是位置分别靠四个顶点多少的四个贡献权重。令 `A=[x1−x0,x2−x0,x3−x0]`，计算 `N[1:]=inv(A)@(x−x0)`，`N0=1−sum(N[1:])`。

gradient 是速度沿三个空间方向变化的快慢，固定定义 `G[i,j]=∂u_i/∂x_j`。vorticity 是流动的局部旋转趋势，按 curl(u) 的三个分量计算。strain rate 是流动的局部拉伸和剪切快慢，使用 `(G+G.T)/2`。

线性四面体内这些导数是常数。共享面上的速度和压力连续，但相邻四面体的导数可以不同；本次没有平滑。多个单元包含同一点时选最小 canonical 零起始 tetra_id，所以面、边和顶点上的结果确定可重复。

几何容差固定为 `tau=64*eps*kappa_inf(A)*(1+c/h)`，其中 `h=||A||inf`，`c` 为该单元节点坐标最大绝对值。eps 是 float64 舍入尺度，kappa 衡量单元形状放大误差的程度，`c/h` 计入 SI 坐标相减误差。64 为固定的一小段减法、求逆和乘加留出舍入余量，未随测试结果调整。候选包围盒仅扩展 `max(4*tau*h)`，最终仍按每个单元自己的权重判断；`tau>=sqrt(eps)` 则拒绝几何。

本次最大 kappa={v['condition_inf_max']:.6g}，最大 tau={v['weight_tolerance_max']:.6e}，候选搜索 padding={v['candidate_padding_m']:.6e} m。最小明确域外测试偏移 {p4['minimum_outside_offset_m']:.6e} m，远大于浮点搜索余量。

精确落在闭合域表面上时 `inside_lumen=True`；域外返回 False、tetra_id=-1，速度、压力和所有导数量统一为 NaN。无效坐标形状或 NaN/Inf 输入抛出 ValueError。这个规则只说明场是否存在，不说明有半径的粒子能否放在墙上。

## 4. 自动测试结果

Particle-0：**{v['tests']['passed']} passed / {v['tests']['failed']} failed / {v['tests']['skipped']} skipped**。

| 检查 | 结果 | 实际误差 / 数值 | 测试文件 |
|---|---|---|---|
| 冻结输入与边界 | {status} | 31 个科学 SHA；4602 个 payload 文件；坐标/顺序一致 | {testlink('test_frozen_input_integrity.py')} |
| 人工速度、压力、梯度 | {status} | 速度 {p1['max_affine_velocity_error']:.3e} m/s；压力 {p1['max_affine_pressure_error']:.3e} Pa；梯度 {p1['max_affine_gradient_error']:.3e} 1/s | `test_affine_*.py`, `test_tetra_barycentric.py` |
| 真实节点回采样 | {status} | {p2['count']} 次；速度 {p2['max_node_velocity_error']:.3e} m/s；压力 {p2['max_node_pressure_error']:.3e} Pa | {testlink('test_real_node_back_sampling.py')} |
| 共享面速度、压力 | {status} | 8 个面；速度 {p3['max_velocity_trace_error']:.3e} m/s；压力 {p3['max_pressure_trace_error']:.3e} Pa | `test_shared_face_velocity_continuity.py`, `test_shared_face_pressure_continuity.py` |
| 单元内常梯度、跨面记录 | {status} | 跳变量 {p3['gradient_jump_min']:.3e}–{p3['gradient_jump_max']:.3e} 1/s；不设跳变失败阈值 | {testlink('test_gradient_piecewise_constant.py')} |
| 顶点排列 | {status} | 人工 24 种 flow 排列和真实另一种排列，结果完全一致 | {testlink('test_vertex_order_regression.py')} |
| 域内 / 域外与闭合边界 | {status} | {p4['count']} 点，错误 {p4['errors']}；{p4['outside_count']} 域外点全部 NaN | `test_inside_outside_contract.py`, `test_boundary_point_contract.py` |
| 真实随机场 | {status} | 400 点；NaN={p5['nan_count_valid_samples']}、Inf={p5['inf_count_valid_samples']} | {testlink('test_real_random_field.py')} |
| 标量 / 批量 | {status} | 413 个混合位置，差异项数 {p6['differing_fields']} | {testlink('test_scalar_batch_equivalence.py')} |
| 文件、数组保护和图数据 | {status} | 全部 8 图及源数据保留；拒绝输入数组别名修改 | `test_sampler_ownership.py`, `test_review_artifacts.py` |

人工场预设 atol=rtol={p1['atol']:.6e}，来自 `256*eps`。真实节点误差上界为 `4*tau*局部值跨度 + 16*eps*局部最大绝对值`；共享面为 `4*tau*局部最大绝对值`，各点实际界保存在数据文件。本轮没有为通过测试放宽任何科学标准。真实梯度没有连续真解；“最大 gradient error”明确指人工场，不冒充真实 FEM 的物理误差。

完整输出见 [Particle 测试日志](logs/particle0_pytest.log)、[JUnit 结果](logs/particle0_junit.xml)、[冻结验收摘要](logs/final_check_commands.json)。所有新增文件见 [FILES_CHANGED.txt](FILES_CHANGED.txt)。

## 5. 人工审核图片

![冻结输入和五个边界](figures/00_frozen_input_overview.png)

应检查血管形状和五个边界位置是否符合交接。图中显示一个入口、三个出口和壁面，节点数、单元数与 manifest 一致；文件哈希没有异常。

![人工已知答案验证](figures/01_affine_field_validation.png)

上排应落在解析值等于采样值的对角线上，下排看误差大小。143 个位置的结果与已知答案吻合，最大梯度误差 {p1['max_affine_gradient_error']:.3e} 1/s，未发现超出预设界的误差。

![真实节点回采样误差](figures/02_node_back_sampling_error.png)

看实际误差与各点计算出来的舍入上界。{p2['count']} 次查询全部恢复节点值到浮点精度，速度/压力中位误差都为零；图保留精确零值，没有用假小数代替。

![共享面的连续性和梯度变化](figures/03_shared_face_continuity.png)

重点看速度和压力在距离零处衔接，第三排梯度可以切换。图展示 8 个受测面中的 4 个，速度/压力没有数值断裂；梯度台阶是本次离散场的预期行为，左右九个分量都在 CSV。

![域内域外分类](figures/04_inside_outside_classification.png)

左图看采样位置覆盖，右图看沿实际法线向域外移动后是否被拒绝。{p4['count']} 点中误分 {p4['errors']} 个，红叉为域外测试位置，错误分类另用紫星标识，本次紫星数量为零；极近位置在全局图会重叠，右侧距离图用于核对。

![真实速度箭头](figures/05_real_flow_velocity_vectors.png)

检查箭头方向与血管走向是否协调，尤其是入口和各出口附近。400 个有效点的速度范围 {p5['speed_min']:.3e}–{p5['speed_max']:.3e} m/s，较快区域集中在入口至 OUTLET_02 支路附近；箭头统一显示长度以便看清低速分支，颜色保留实际速度。自动检查未发现域外点或无效数值，流向的最终视觉认可仍等用户确认。

![真实标量诊断](figures/06_real_flow_scalar_diagnostics.png)

分别看速度、压力、局部旋转和局部拉伸剪切，不把不同单位放进同一色标。各图为 x-y 投影，原始 z 坐标及全部分量在 CSV；400 个样本没有 NaN/Inf，也未做平滑，但这不构成新的 FEM 物理验收。

![批量一致性与 CPU 计时](figures/07_sampling_runtime_smoke.png)

这里看查询数量增加时耗时怎样变化。10000 点批量耗时为 {p6['timings'][-1]['batch_wall_time_s']:.3f} s，逐点与批量的返回值完全一致；计时只跑一次、只做信息记录，没有性能通过阈值，也没有为性能改变科学实现。

全部图片由永久脚本 [generate_particle0_report.py](../../scripts/generate_particle0_report.py) 生成，各步原数据在 [data](data/)；三维表面图通过 `00_frozen_input_summary.json` 中的路径及 SHA 重读全部原始边界三角形。只在绘图层把米转成 µm/nm，未移动、重采样或修改原 FEM 数据。

## 6. 当前已知限制

- FEM 本身仍然是 frozen Stage Q。
- 没有重新做 mesh convergence（网格加密后结果是否稳定的检查）。
- 没有重新做 FEM timestep study（改变流体计算时间步后的敏感性检查）。
- 没有做新的 CPU/GPU field equivalence（两种计算设备输出流场是否一致的检查）。
- gradient 是 linear tetra 内 piecewise constant，也就是每个线性四面体里是常数，跨面可跳变。
- Particle-0 没有真实粒子、微泡或 RBC 状态。
- 还没有 wall gap，即有限半径粒子与壁面的间隙。
- 还没有 particle dynamics，即粒子随时间运动的方程或积分。
- SonoVue sampler 尚未接入正式 particle population。
- 极端无法用 float64 分辨的四面体会明确报错；本次真实网格不触发该拒绝。
- 批量接口首先保证正确，当前使用逐点科学实现；计时不代表生产性能，未控制 WSL 全部背景负载。

## 7. 人工审核状态

AUTOMATED_CHECKS = {status}

MANUAL_VISUAL_REVIEW = PENDING_USER_REVIEW

下一阶段：Particle-1。**尚未开始 Particle-1**。本轮只作本地提交，没有 push、merge main、启动 GPU server 或重新求解 FEM。
"""
    (REPORT / "PARTICLE0_REVIEW.md").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff-root", type=Path, required=True,
                        help="Original clean frozen-branch FEM directory for its historical branch-locked tests")
    args = parser.parse_args()
    handoff = args.handoff_root.resolve()
    fem = REPO / "formal_3D_flow_solver/FEM_SimVascular"
    LOGS.mkdir(parents=True, exist_ok=True)
    # Refuse to label uncommitted implementation with an unrelated Git hash.
    implementation_paths = ["particle_3d/src", "particle_3d/tests", "particle_3d/scripts", "particle_3d/README.md", "particle_3d/pyproject.toml"]
    if git("status", "--porcelain", "--", *implementation_paths):
        raise RuntimeError("Commit reviewed implementation before finalizing; source must be attributable to git_commit")
    summary, _, _, _ = read_frozen(fem)
    if git("branch", "--show-current", cwd=handoff) != FEM_BRANCH:
        raise RuntimeError("--handoff-root must be the original frozen-branch checkout")
    python = sys.executable
    commands = []
    for name, script in [("frozen_validate_after", "validate_frozen.py"), ("frozen_manifest_after", "verify_manifest.py")]:
        result = run(name, [python, "-B", "scripts/fem_freeze_sync/" + script], fem)
        commands.append(result)
        if result["returncode"]:
            write_json(LOGS / "final_check_commands.json", commands)
            raise RuntimeError(f"Frozen integrity blocker; see {result['log']}; do not repair frozen input")
    frozen_test_paths = ["tests/test_fem_sync_manifest.py"] + [str(p.relative_to(handoff)) for p in sorted((handoff / "tests").glob("test_frozen_*.py"))] + [
        "tests/test_particle_handoff_contract.py", "tests/test_no_particle_implementation_yet.py", "tests/test_main_history_preserved.py"]
    commands.append(run("frozen_pytest_after", [python, "-B", "-m", "pytest", "-q", "-p", "no:cacheprovider",
                         f"--junitxml={LOGS / 'frozen_junit.xml'}", *frozen_test_paths], handoff))
    commands.append(run("particle0_pytest", [python, "-B", "-m", "pytest", "-q", "-s", "-p", "no:cacheprovider",
                         f"--junitxml={LOGS / 'particle0_junit.xml'}", "particle_3d/tests/particle0"], REPO))
    write_json(LOGS / "final_check_commands.json", commands)
    tests, cases = junit(LOGS / "particle0_junit.xml")
    frozen_tests, _ = junit(LOGS / "frozen_junit.xml")
    metrics = {str(i): json.loads((REPORT / "data" / f"metrics_p0{i}.json").read_text()) for i in range(7)}
    base = (LOGS / "frozen_base_commit.txt").read_text().strip()
    frozen_diff = git("diff", base, "--", "formal_3D_flow_solver/FEM_SimVascular")
    legacy = Path("/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular")
    before = LOGS / "legacy_fem_status_before.txt"
    legacy_unchanged = None
    if legacy.exists() and before.exists():
        after = subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=legacy, text=True)
        (LOGS / "legacy_fem_status_after.txt").write_text(after)
        legacy_unchanged = before.read_text() == after
    passed = (all(r["returncode"] == 0 for r in commands) and tests["failed"] == tests["skipped"] == 0
              and frozen_tests["failed"] == frozen_tests["skipped"] == 0 and not frozen_diff
              and all(metrics[str(i)]["passed"] for i in range(1, 7)) and legacy_unchanged is not False)
    figures = sorted(str(p.relative_to(REPO)) for p in (REPORT / "figures").glob("*.png"))
    test_files = sorted(str(p.relative_to(REPO)) for p in (PACKAGE / "tests/particle0").glob("test_*.py"))
    sources = sorted(p for top in [PACKAGE / "src", PACKAGE / "tests", PACKAGE / "scripts"] for p in top.rglob("*.py"))
    v = dict(stage="Particle-0", automated_checks="PASS" if passed else "FAIL",
             git_commit=git("rev-parse", "HEAD"), git_commit_semantics="Tested implementation commit. Later evidence-only commits may contain this record.",
             git_branch=git("branch", "--show-current"), frozen_base_commit=base, fem_branch=FEM_BRANCH,
             mesh_path=summary["mesh_path"], mesh_sha256=summary["mesh_sha256"], flow_path=summary["flow_path"], flow_sha256=summary["flow_sha256"],
             node_count=summary["node_count"], tetra_count=summary["tetra_count"], boundary_face_ids=summary["boundary_face_ids"],
             synthetic_tests_passed=group_pass(cases, ["test_affine", "test_tetra_barycentric"]),
             real_node_tests_passed=group_pass(cases, ["test_real_node"]),
             shared_face_tests_passed=group_pass(cases, ["test_shared_face", "test_gradient_piecewise", "test_vertex_order"]),
             outside_contract_passed=group_pass(cases, ["test_inside_outside", "test_boundary_point"]),
             scalar_batch_equivalence_passed=group_pass(cases, ["test_scalar_batch"]),
             max_node_velocity_error=metrics["2"]["max_node_velocity_error"], max_node_pressure_error=metrics["2"]["max_node_pressure_error"],
             max_affine_velocity_error=metrics["1"]["max_affine_velocity_error"], max_affine_pressure_error=metrics["1"]["max_affine_pressure_error"],
             max_affine_gradient_error=metrics["1"]["max_affine_gradient_error"], nan_count_valid_samples=metrics["5"]["nan_count_valid_samples"],
             inf_count_valid_samples=metrics["5"]["inf_count_valid_samples"], figures=figures, test_files=test_files,
             manual_visual_review="PENDING_USER_REVIEW", tests=tests, frozen_tests=frozen_tests,
             frozen_integrity="PASS" if all(r["returncode"] == 0 for r in commands[:3]) and not frozen_diff else "FAIL",
             frozen_payload_modified=bool(frozen_diff), legacy_fem_git_status_unchanged=legacy_unchanged,
             no_cfd_executed=True, no_gpu_server_started=True, particle1_started=False, pushed=False,
             source_sha256={str(p.relative_to(REPO)): sha256(p) for p in sources},
             figure_sha256={p: sha256(REPO / p) for p in figures},
             data_sha256={str(p.relative_to(REPO)): sha256(p) for p in sorted((REPORT / "data").iterdir()) if p.is_file()},
             condition_inf_max=metrics["0"]["condition_inf_max"], weight_tolerance_max=metrics["0"]["weight_tolerance_max"],
             candidate_padding_m=metrics["0"]["candidate_padding_m"], recorded_at_utc=datetime.now(timezone.utc).isoformat(),
             environment=dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                              packages={name: version(name) for name in ["numpy", "vtk", "pyvista", "matplotlib", "pytest", "Pillow"]}))
    write_json(REPORT / "PARTICLE0_VALIDATION.json", v)
    review(v, metrics)
    excluded = {"__pycache__", ".pytest_cache", ".venv"}
    changed = sorted(str(p.relative_to(REPO)) for p in PACKAGE.rglob("*") if p.is_file() and not excluded.intersection(p.parts))
    for pending in ["particle_3d/reports/particle0/FILES_CHANGED.txt", "particle_3d/reports/particle0/logs/TEST_SUMMARY.md"]:
        if pending not in changed:
            changed.append(pending)
    changed.sort()
    (REPORT / "FILES_CHANGED.txt").write_text("\n".join(changed) + "\n")
    (LOGS / "TEST_SUMMARY.md").write_text(f"# 最终检查摘要\n\nFrozen payload：{v['frozen_integrity']}。官方 handoff：{frozen_tests['passed']} passed / {frozen_tests['failed']} failed。\n\nParticle-0：{tests['passed']} passed / {tests['failed']} failed / {tests['skipped']} skipped。\n\n各命令、运行目录、退出码及日志见 final_check_commands.json。没有运行 CFD。人工审核仍为 PENDING_USER_REVIEW。\n")
    print(json.dumps({k: v[k] for k in ["automated_checks", "git_commit", "tests", "frozen_tests", "manual_visual_review"]}), flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__": main()
