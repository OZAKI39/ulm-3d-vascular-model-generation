"""Data-derived Chinese report and byte-preserving progress-document insertion."""

from pathlib import Path
import hashlib
import difflib
import re
from .common import ROOT, BASE, BASELINE, js, write, sha


def main():
    T = js(BASE / "analysis/LUBRICATION_TRAJECTORIES.json")
    M = js(BASE / "analysis/MORPHOLOGY.json")["comparisons"]
    B = js(BASE / "analysis/B2_GEOMETRY_TRAJECTORIES.json")
    K = js(BASE / "analysis/COST_AND_YIELD.json")
    E = js(BASE / "analysis/EDGE_AUDIT.json")
    D = js(BASE / "geometry_v1/DENSE_STATIC.json")
    S = js(BASE / "geometry_v1/SMALL_SYSTEMS.json")
    bg = B["comparisons"][0]
    v = bg["maxima"]
    med = K["matched_diagnostic_medians"]
    a = med["M0_LUB2_PBD12_DT10"]
    b = med["M0_LUB4_PBD12_DT10"]
    cost_ratio = b["integration_quality_wall_s"] / a["integration_quality_wall_s"]
    full_ratio = b["resident_body_s"] / a["resident_body_s"]
    remote = js(BASE / "REMOTE_VERIFICATION.json")
    report = [
        "# F10H：局部润滑、几何软接触与逐泡形态收敛\n",
        "实际开发、GPU科学计算、CPU独立参考、逐步审计和科研图已完成。四轮润滑通过原两个短窗口，但O/1102在6.370 ms再次拒绝，不能完成10 ms。几何修正版B2在同模型10/5 µs比较中通过位置与形状门槛，却仍未通过连边时间和顺序门槛。优先继续有界B2研究，不升级生产默认；EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false。\n",
        "## 身份、数据与实际调用\n",
        f"完整F10G基准 {BASELINE}；独立分支 dev/fast-pzt-f10h-lubrication-contact-20261012。WSL工程：{ROOT}；服务器：/workspace/mevo_acoustic_microbubbles/f10h。103项F10G归档、13项实际输入以及已归档执行提交均核对，原F10G工作区和历史源码未修改。\n",
        "只新增f10h_lubrication.py、f10h_contact.py和f10h_diagnostic.py。原Physics已支持四轮，GPU迭代核直接复用；F3在默认参数下分配，然后显式选择候选轮数。B2接触为CPU有界非线性原型，GPU保留原声力、CFD和接受检查，没有GPU生产软接触或新全局QP。实际顺序包括预测后的壁面检查、每三轮PBD后的壁面检查、最终CFD与几何检查，详见数值模型文档。\n",
        "全部方案使用原始坐标和完整控制历史确定性重放。没有使用失败提案作为重启、没有新种子、压力或声场表；没有重算CFD。控制事件、Graph/诊断一致性、两轮回归、故意拒绝和时钟不推进均通过。新模式B>1明确拒绝，不能声称已验证新模式批量成员隔离。\n",
        "执行提交与远程快照：\n\n"
        + "\n".join(
            f"- `{c}`：{n}项源码文件核对。" for c, n in remote["source_files_against_git_archive"].items()
        )
        + "\n",
        "## 按任务要求回答的16个问题\n",
        "1. **四轮是否真实通过两个窗口？**是。O/1102从原始初态到5.000 ms，准备1104到2.350 ms，分别500和235个接受步，原17.8 nm与178 nm预算不变。\n",
        "2. **完整10 ms是否通过？**否。只尝试O一次，在第637步、6.360→6.370 ms触发2048。预测侵入82.601 nm，PBD12后18.011186 nm，最后合法时间6.360 ms；该步最大单轮修正27.795 nm，墙面修正零。失败立即停止，未增加轮数或继续非法积分。\n",
        "实际运行清单（各变体不是独立初态）：\n\n|案例|模式|目标/ms|接受步|实际合法终点/ms|状态|\n|---|---|---:|---:|---:|---|",
    ]
    for t in T["trials"] + B["trials"]:
        dt = 5e-6 if t["mode"].endswith("DT5") else 1e-5
        report.append(
            f"|{t['case']}|{t['mode']}|{t['target_s'] * 1e3:.3f}|{t['accepted_steps']}|{t['accepted_steps'] * dt * 1e3:.3f}|{t['status']}|"
        )
    report += [
        f"\n3. **四轮比两轮贵多少？**相同O初态、控制、0.2 ms、50 µs保存和逐步诊断，三次暖重复：积分＋质量中位数{a['integration_quality_wall_s']:.6f}→{b['integration_quality_wall_s']:.6f} s（增加{(cost_ratio - 1) * 100:.2f}%）；含输出主体{a['resident_body_s']:.6f}→{b['resident_body_s']:.6f} s（比值{full_ratio:.4f}）。输出占主导且存在计时波动，不能把后者解释为四轮生产加速。CUDA跨度含CPU/同步间隙。\n",
        "4. **四轮是否充分收敛？**否。复用F10G同状态参考，四轮最大逐泡速度差约为自由速度尺度的23.52%／33.58%。本轮还将每个真实积分步的速度代回同一阻力算子，分别保存原始量纲、块预条件残差、离线参考差与最差法向泡对；这些不是实际流体误差。\n",
        "|真实四轮窗口|最大预条件残差/自由速度|最大逐泡速度差/自由速度|最大接受重叠/nm|最大单轮修正/nm|\n|---|---:|---:|---:|---:|",
    ]
    for r in T["rows"]:
        report.append(
            f"|{'/'.join(r['path'].split('/')[-2:])}|{100 * r['max_normalized_preconditioned']:.3f}%|{100 * r['max_normalized_velocity_error']:.3f}%|{r['max_accepted_overlap_m'] * 1e9:.6f}|{r['max_single_pass_m'] * 1e9:.6f}|"
        )
    report += [
        "\n5. **四轮10/5 µs质心和形态差异？**使用相同模式、初态和控制，按真实共同接受时间比较；没有借用历史两轮C1。\n",
        "|案例|共同终点/ms|质心差/nm|固定ID差/nm|Kabsch/nm|Rg相对差/%|长短比/%|张量/%|连边数10/5 µs|\n|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for c in M:
        z = c["maxima"]
        report.append(
            f"|{c['case']}|{c['common_end_s'] * 1e3:.3f}|{z['centroid_norm_m'] * 1e9:.6f}|{z['max_fixed_ID_m'] * 1e9:.6f}|{z['Kabsch_RMSD_m'] * 1e9:.6f}|{z['Rg_relative'] * 100:.6f}|{z['aspect_relative'] * 100:.6f}|{z['tensor_relative'] * 100:.6f}|{c['event_counts']}|"
        )
    report += [
        "\n6. **四轮的内部排列和事件收敛吗？**否。O的位置、Kabsch、Rg、长短比、张量和边事件均未过门槛；准备阶段的Rg、长短比和张量通过，但固定ID、Kabsch、成员和边事件未通过。准备过程存在持续合并与右删失，未观察到分裂；右删失不当作PASS。O一直为32泡连通团，不能把这一点当作内部排列稳定。没有继续缩小dt。\n",
        "7. **B2还有分离后排斥吗？**旧F10G有17个泡对步，固定ID、时间、真实间隙、线性间隙、新旧法向、距离和反力均逐项追溯。新算法在每次迭代和最终位置上重算真实球心距离、法向和接触集合；17次相同输入重解以及新1500个接受步均未出现实际分离后的几何弹簧反力。这里是几何罚力，不是预测性不可穿透约束。\n",
        f"8. **修复后的小系统和32泡怎样？**固定F10G的2kc={js(ROOT / 'configs/fast_pzt/f10h/protocol.json')['B2_kc_N_m']:.17g} N/m。12个两泡与6个四泡案例通过独立非线性参考，另含ID交换与可用的共线解析式；最大参考位置差{max(r['reference_error_m'] for r in S['rows']) * 1e12:.6f} pm。两真实32泡静态残余为{max(D['rows'][0]['result']['geometric_overlap_m']) * 1e9:.6f}／{max(D['rows'][1]['result']['geometric_overlap_m']) * 1e9:.6f} nm，完整GPU环境检查均通过。\n",
        "新B2 O/1102短窗：\n\n|dt/µs|接受步|最大重叠/nm|最大单次修正/nm|平均/最大外层迭代|最大方程残差/m|错误分离排斥|\n|---|---:|---:|---:|---|---:|---:|",
    ]
    for r in B["rows"]:
        report.append(
            f"|{5 if r['mode'].endswith('DT5') else 10}|{r['accepted_steps']}|{r['max_accepted_overlap_m'] * 1e9:.6f}|{r['max_single_iteration_m'] * 1e9:.6f}|{r['mean_outer_iterations']:.3f}/{r['max_outer_iterations']}|{r['max_equation_residual_m']:.3e}|{len(r['false_reactions'])}|"
        )
    report += [
        f"\n9. **B2自己的步长收敛？**同一几何B2、刚度、阻力和控制的10/5 µs均完成5 ms。质心差{v['centroid_norm_m'] * 1e9:.6f} nm，固定ID差{v['max_fixed_ID_m'] * 1e9:.6f} nm，Kabsch {v['Kabsch_RMSD_m'] * 1e9:.6f} nm；Rg／长短比／张量相对差{v['Rg_relative'] * 100:.6f}%／{v['aspect_relative'] * 100:.6f}%／{v['tensor_relative'] * 100:.6f}%，这些全部通过登记阈值，成员相同。\n",
        f"但29／29条边事件虽有相同的泡对和开关序列，最大时间差仍为{bg['matched_edge_max_time_difference_s'] * 1e6:.3f} µs，且部分跨泡对顺序不同。最差为ID 7–24断边：1.660／1.790 ms；按所有接受位置独立重构七条试验的每一条GPU边事件均完全一致，排除了事件记录时钟错位。该对相同时间最大间隙差约{E['worst_pair_max_gap_difference_m'] * 1e9:.3f} nm，说明阈值附近的小空间差异仍会改变事件时刻，不能据此放宽10 µs门槛。连通阈值表面间隙178 nm与实际接触的零间隙是两回事。完整形态仍MORPHOLOGY_NOT_CONVERGED。\n",
        "10. **PBD与B2形态差异如何解释？**PBD使用有限轮局部润滑和几何投影；B2使用同一离散阻力的充分求解和隐式弹簧。二者的迁移率近似与接触规则同时不同，逐泡形态不应逐位相同。四轮PBD与几何B2的同5 ms比较：最大ID差3.764 µm、Kabsch 1.553 µm；旧冻结法向B2与几何B2仅35.501 nm／15.569 nm。可以归因于离散方案差异，但现有比较不能分离每项的全部贡献，也不能证明哪个更符合真实实验。\n",
        "11. **优先继续哪个候选？**几何B2更值得有限研究：它修复单边性，且同模型固定ID差仅约10.41 nm，四轮PBD对应差为7.81 µm，外层最多3次；不过接触事件仍未过门槛，CPU原型更贵。下一轮优先事件敏感性、求解误差与时间离散误差分离，再考虑小规模GPU可行性，不直接生产化。\n",
        "12. **数值通过但形态尚未全面可信的方案？**四轮10/5 µs短窗均通过接触预算，却未通过位置与连边形态；几何B2两短窗通过接触、位置与张量等门槛，但边时间和顺序失败。B2的O始终一团，无持续分裂合并事件，不能推广成已验证丰富拓扑变化。\n",
        "13. **有正式GPU批量生产证据吗？**没有。原M0默认不变；新模式B>1明确拒绝，未跑B64或新完整OFF/S/O，没有AI训练。\n",
        f"14. **四组图路径？**{BASE / 'figures'}。Fig_F10H_1_lubrication_contact、Fig_F10H_2_failure_window_replay、Fig_F10H_3_geometry_morphology、Fig_F10H_4_accuracy_cost_admission，均有英文PNG/PDF。真实半径1.78 µm，全部图实际打开审查；来源、执行提交、输入SHA和VISUAL_REVIEW在同目录。\n",
        "15. **科研文档已备份与更新？**在microbubble_aggregation.md最小插入本轮证据；备份为microbubble_aggregation.md.bak_F10H_20261012。before/after、diff、SHA和原标题结构核对在progress_document目录，删除插入段即可逐字节恢复原文。\n",
        "16. **F10I优先什么？**先围绕ID 7–24等已识别的真实阈值跨越检查事件定位、接触求解容差和时间离散贡献，保持原科学阈值并报告阈值定义的敏感性；再以CPU参考为依据测试少量残差受控的GPU润滑／接触实现。不要直接延长B2到10 ms、加大刚度、扩大种子或用成员不变冒充事件收敛。\n",
        "## 实际成本与数据边界\n",
        "|模式与窗口|积分＋质量/s|详细诊断输出/s|保存读回/s|HDF5/s|主体/s|\n|---|---:|---:|---:|---:|---:|",
    ]
    for row in T["trials"] + B["trials"]:
        t = row["timing"]
        report.append(
            f"|{row['case']} {row['mode']} {row['target_s'] * 1e3:g} ms|{t['integration_quality_wall_s']:.6f}|{t['observer_output_s']:.6f}|{t['sync_readback_wall_s']:.6f}|{t['hdf5_output_s']:.6f}|{t['resident_body_s']:.6f}|"
        )
    report += [
        f"\n九个远程进程总墙时{K['total_process_wall_s']:.6f} s，包括CPU参考、工程检查、失败和诊断输出；不是全部研发时间或纯GPU时间。科学试验共4341个接受步、1个科学拒绝（另有工程故意拒绝）；完整10 ms合格样本0，新独立初态0，完整同模式OFF/S/O三联0。CPU接触初始化、迭代、阻力装配和提交读回分别保存在B2_GEOMETRY_TRAJECTORIES。\n",
        "模型仍没有校准脂质壳局部接触刚度、近场声学、润滑下限和壁面水动力。几何单边性精确到离散求解容差，并不等于真实壳接触本构。非线性能量未证明全局凸或解唯一。所有变体保留原split_group，通用记录器短窗samples/hour不作为生产吞吐。\n",
        "一次本地分析曾在rsync未完成时找不到SOFT_STEPS而停止；没有修改科学数据，等待回传和哈希完成后仅重跑离线分析。记录在OFFLINE_ATTEMPTS.json。服务器科学和工程进程均正常退出。\n",
        "复现：先确认回传完整，依次执行scripts.fast_pzt_f10h.analyze、events、figures、delivery、verify；重画图不重新跑轨迹。所有科学执行使用独立Supervisor标签和不可变Git快照，源代码按git archive核对，170项远程结果与13项输入哈希通过。\n",
    ]
    (ROOT / "docs/FAST_PZT_F10H_DEVELOPMENT.md").write_text("\n".join(report))
    insertion = f"""
**F10H把四轮润滑真正放进连续轨迹。**F10G中同一失败状态的四轮预测能过门槛，这次保留原声力、CFD、PBD12和10 µs，从原始初态与完整控制重放：O/1102完成5 ms，准备1104完成2.35 ms。但唯一一次O的10 ms尝试在6.370 ms又被18.011 nm残余重叠拒绝，最后合法时间6.360 ms，墙面修正为零。四轮延迟了失败，没有解决全程问题，也没有成为充分收敛的润滑解。

**步长对照把内部排列的限制看得更清楚。**同一四轮模式10/5 µs均通过两个短窗，O的质心最大差约1.107 nm，固定ID差却达7.813 µm、Kabsch 2.547 µm，Rg和长短比相对差约1.863%／11.863%；准备阶段固定ID与Kabsch为394.122／109.361 nm，也超过事先登记的178／89 nm。成员和连边另行比较，准备过程的右删失没有算通过。不能用质心稳定替代逐泡收敛，也没有继续无界缩小步长或加轮数。

**几何软接触先修复一件明确的事。**保持F10G的2kc约0.0464132 N/m和原离散润滑关系，改为每次迭代更新真实球心距离、法向和接触集合，实际分离时弹簧力为零，不叠加PBD。12项两泡与6项四泡独立非线性参考通过，两真实32泡静态残余17.498／17.317 nm；旧17个错误分离排斥泡对步的同输入重解不再出现该问题。随后只对O做5 ms的10/5 µs真实对照，共1500个接受步，最大重叠约15.260／15.325 nm，未放宽原预算。

**位置形态更一致，连边时刻还没有过关。**几何B2的同模型最大固定ID差约10.412 nm、Kabsch 4.638 nm，Rg／长短比／张量相对差约0.00636%／0.03765%／0.05983%，均通过登记门槛，成员相同，也没有分离后排斥。但29／29条边事件的最大时间差仍为130 µs，超过10 µs门槛；ID 7–24在1.660／1.790 ms断边，部分跨泡对顺序也不同。全部七条实际试验的连边已按接受位置独立重构，与GPU记录完全一致，因此不能归因于保存时钟错位，更不能凭最终连通图相同就宣布完整形态收敛。

**更值得继续研究，还没有生产准入。**匹配0.2 ms的三次诊断计时，四轮相对两轮积分与质量约增加{(cost_ratio - 1) * 100:.2f}%，完整主体由0.613到0.606 s，差异被输出成本和波动主导，不是生产加速。几何B2的5 ms积分与质量约4.053／7.404 s，详细输出另约10.540／20.563 s；九个远程进程共{K["total_process_wall_s"]:.3f} s，不是全部研发时间或纯GPU时间。下一轮优先B2事件敏感性和残差受控低成本实现，不直接扩大长轨迹或批量。近场声学、润滑下限和接触刚度仍未物理校准，原M0默认不变，继续 **EXPLORATORY / PHYSICS_UNCALIBRATED / MORPHOLOGY_NOT_CONVERGED，training_ready=false**。

> **F10H证据与下一步：**docs/FAST_PZT_F10H_DEVELOPMENT.md、docs/FAST_PZT_F10H_NUMERICAL_MODEL.md与results/fast_pzt/f10h/保留预登记、原始拒绝、同模型逐泡与连边对照、几何反力、独立参考、真实成本、四组英文PNG/PDF及源码和回传哈希。F10I应先解释已识别泡对的阈值跨越误差，再判断GPU实现是否同时带来数值依据和计算收益；不把几何罚力的单边性当作实测壳接触校准。

"""
    target = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = target.with_name(target.name + ".bak_F10H_20261012")
    raw = target.read_bytes()
    anchor = b"\n---\n\n### " + "尾声：故事里还没讲完的部分".encode()
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    inc = insertion.encode()
    if backup.exists():
        before = backup.read_bytes()
        assert raw == before.replace(anchor, inc + anchor, 1), (
            "Progress file changed; inspect before updating"
        )
    else:
        assert (
            hashlib.sha256(raw).hexdigest()
            == "92b8631a56b32a2ea8290d8632d3ebf6f4577271e2ea06e2c8aa63c39e659d1e"
        ), "Progress file changed since reviewed F10G document"
        assert raw.count(anchor) == 1
        before = raw
        backup.write_bytes(raw)
        target.write_bytes(raw.replace(anchor, inc + anchor, 1))
    after = target.read_bytes()
    assert after.replace(inc, b"", 1) == before

    def titles(b):
        return [s for s in b.decode().splitlines() if re.match(r"^#{1,6} ", s)]

    assert titles(before) == titles(after)
    for name, data in [("before.md", before), ("after.md", after), ("insertion.md", inc)]:
        (out / name).write_bytes(data)
    (out / "changes.diff").write_text(
        "".join(
            difflib.unified_diff(
                before.decode().splitlines(True),
                after.decode().splitlines(True),
                fromfile="before",
                tofile="after",
            )
        )
    )
    write(
        out / "VERIFICATION.json",
        dict(
            status="PASS",
            before_sha256=sha(backup),
            after_sha256=sha(target),
            insertion_sha256=sha(out / "insertion.md"),
            backup=str(backup),
            original_bytes_preserved=True,
            headings_preserved=True,
            headings=titles(before),
        ),
    )
    print("Report and minimal progress insertion complete")


if __name__ == "__main__":
    main()
