"""Create a data-derived report and a minimal, byte-preserving progress insertion."""

from pathlib import Path
import difflib
import hashlib
import re
from .common import BASE, ROOT, BASELINE, CASES, js, sha, write


def main():
    morph = js(BASE / "analysis/MORPHOLOGY.json")
    traj = js(BASE / "analysis/TRAJECTORY_AUDIT.json")["rows"]
    cost = js(BASE / "analysis/COST_AND_YIELD.json")
    lub = js(BASE / "lubrication_v1/RESULTS.json")["cases"]
    perf = js(BASE / "performance_v2/RESULTS.json")["rows"]
    report = [
        "# F10G：过阻尼接触与局部润滑可行性\n",
        "实际开发、服务器验证与图表已完成。离散润滑参考、小系统隐式接触和两条有界短窗通过；逐泡形态与物理校准没有因此通过。原默认没有改变，training_ready=false。\n",
        "## 身份与范围\n",
        f"完整 F10F 基准 {BASELINE}，原工作区干净、90 项交付清单及 13 项实际输入 SHA-256 已核对。新分支 dev/fast-pzt-f10g-overdamped-contact-20261011；WSL 根目录 {ROOT}；服务器 /workspace/mevo_acoustic_microbubbles/f10g。\n",
        "固定执行提交：润滑 6254996fb；两泡／四泡／32泡 fd10ae32a；完整环境静态及两短窗 308578ccd；补充适配器检查 5152d418a；当前 active-mask 核计时 d67cfd2ed。早期不带 active-mask 的性能对照 0bafe95c8 单独保留并计费，不作当前生产核对照。\n",
        "新模块 f10g_contact.py 为 CPU 独立参考和有界隐式原型；f10g_diagnostic.py 为 B=1 的 GPU 原受力／CPU 接触诊断适配器。不改任何 F1–F10F 源文件、原轨迹、默认与预算；没有新增种子、场表、CFD、PZT、压力或 AI 训练。\n",
        "## 按任务要求回答的 16 个问题\n",
        "1. **原两轮与参考差多少？**均值补偿前的最大块预条件残差为自由速度尺度的 18.564%／15.192%；补偿后为 18.532%／15.152%。两轮与离散参考的最大逐泡速度差分别为自由速度尺度的 32.531%／35.266%。这是同一个近似阻力系统的代数误差，不是实验流速误差。参考矩阵经对称性、正定性、稀疏／Cholesky 交叉验证；条件数约 1514／1204，间隙下限激活 42／35 对。\n",
        "2. **增加轮数是否改善接触？**是，限于两例同状态、同外力、同 CFD、同 10 µs 的只读预测。没有运行四轮润滑的真实独立轨迹，也没有把瞬时结果推广到全程。\n",
        "|状态|润滑|预测侵入/nm|PBD12残余/nm|原预算|\n|---|---|---:|---:|---|",
    ]
    for c in lub:
        for r in c["rows"]:
            report.append(
                f"|{c['case']}|{r['mode']}|{r['predicted_overlap_m'] * 1e9:.6f}|{r['final_overlap_m'] * 1e9:.6f}|{'PASS' if r['contact_budget_pass'] else 'FAIL'}|"
            )
    report += [
        "\n四轮已使两例越过原17.8 nm门槛，但八轮的最大速度误差仍为自由速度尺度的约19.74%／30.34%，不能当作收敛解。LUBRICATION_CONVERGENCE_NOT_DECISIVE=false；接触相关性明确，未做无限迭代。\n",
        "3. **两泡是否独立验证？**通过90项固定法向测试，覆盖无接触、接近、轻微重叠、分离、切向速度、三档dt与三档kc，并逐项交换ID；另做硬接触极限。解析后向欧拉式与独立力空间活跃集枚举一致。没有把惯性XPBD的dt²和质量矩阵直接搬入Stokes迁移率。\n",
        "4. **四泡是否稳定？**三种非对称共享节点构型、B1/B2各三档刚度共18项均收敛并匹配独立小系统参考；反力非负、合力平衡、位移与迁移率方程一致。这里的PASS是指定线性化模型的求解验证；某些人为构型本身仍有数十nm重叠，不是全部通过32泡原预算或获得物理校准。\n",
        "5. **刚度减半／加倍改变什么？**kc=0.02320659759189329 N/m在软接触试验前由实际法向压紧力／8.9 nm登记，不使用Marmottant的chi代替，也没有终点反调。静态残余如下。\n",
        "|状态|模式|0.5kc / kc / 2kc 残余/nm|0.5kc至2kc最大ID差/nm|\n|---|---|---|---:|",
    ]
    for s in morph["stiffness"]:
        vals = " / ".join(f"{v['max_overlap_m'] * 1e9:.3f}" for v in s["rows"])
        report.append(f"|{s['case']}|{s['mode']}|{vals}|{s['half_vs_double']['max_fixed_ID_m'] * 1e9:.3f}|")
    report += [
        "\n6. **是否有自洽的润滑＋接触？**B2使用同一A的逆，U=A⁻¹U_free+A⁻¹F_contact/ζ；声力与接触力各经过一次迁移率，不重复润滑，不叠加PBD，不加RPY。这只是原离散近似的一致组合。\n",
        "7. **真实32泡是否通过？**B1三档未能同时通过两例。B2只有2kc同时通过两个完整环境静态检查：残余17.501／17.323 nm，最大单次修正1.538／1.008 nm，原压力、CFD、场域和壁面质量通过，墙面修正为零。相同输入GPU重新求得的力／自由速度与保存值一致，候选位置与CPU参考一致。\n",
        "8. **减少了原失败吗？**仅对这两条候选短窗，是：从共同原始初态和控制重新推进，O达到5 ms，准备达到2.35 ms，分别越过旧4.900和2.240 ms失败点。共735个接受步，科学失败0；工程故意超限拒绝2次均保留并验证不推进时钟。没有10 ms新轨迹，没有完整同模式OFF/S/O三联。\n",
        "|候选短窗|接受步|最大重叠/nm|最大单次修正/nm|线性／真实重叠差最大/nm|\n|---|---:|---:|---:|---:|",
    ]
    for r in traj:
        report.append(
            f"|{r['case']}|{r['accepted_steps']}|{r['max_overlap_m'] * 1e9:.6f}|{r['max_single_iteration_m'] * 1e9:.6f}|{r['max_frozen_vs_geometric_overlap_m'] * 1e9:.6f}|"
        )
    report += [
        "\n9. **质心与逐泡形态差多少？**只比较旧算法存在保存接受帧的共同前缀，没有旧失败后的真值。\n",
        "|状态|共同前缀/ms|质心最大/nm|固定ID最大/µm|Kabsch最大/µm|Rg相对差/%|长短比相对差/%|旧／新边事件|\n|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for c in morph["comparisons"]:
        v = c["values"]
        report.append(
            f"|{c['reference'].split('/')[-1]}|{c['accepted_common_end_s'] * 1e3:.2f}|{v['centroid_max_m'] * 1e9:.6f}|{v['fixed_ID_max_m'] * 1e6:.6f}|{v['Kabsch_max_m'] * 1e6:.6f}|{v['Rg_max_relative'] * 100:.4f}|{v['aspect_max_relative'] * 100:.4f}|{v['event_count']}|"
        )
    report += [
        "\n新旧接触和迁移率不同；这些差异不能证明哪种更真实，也不等于证明新模型本身不收敛。固定状态求解容差加严100倍通过；新模型dt收敛未运行。PBD修正路径和Newton迭代路径仅作算法工作量描述，不当作同一物理观测。当前保留MORPHOLOGY_NOT_CONVERGED表示尚无完整形态可信度准入。\n",
        "10. **实际成本可接受吗？**CPU小系统验证成本可控，生产速度尚无证据。当前生产同族masked/ordered核的固定状态计时如下；不含受力、壁面、最终接受、输出与批量。\n",
        "|状态|润滑轮数|仅润滑CUDA跨度/µs|润滑＋均值＋预测＋PBD12跨度/µs|\n|---|---:|---:|---:|",
    ]
    for case in CASES:
        for n in [2, 4, 8]:
            a, b = [
                next(v for v in perf if v["case"] == case and v["iterations"] == n and v["scope"] == scope)[
                    "gpu_median_s"
                ]
                * 1e6
                for scope in ["lubrication_only", "lubrication_balance_prediction_PBD12"]
            ]
            report.append(f"|{case}|{n}|{a:.3f}|{b:.3f}|")
    report += [
        "\n每项25次固定输入重放，不推进轨迹。CPU预装配矩阵直接解约52 µs；Python Jacobi两／四／八轮约94／187／372 µs，不能用CPU实现优劣预测GPU算法优劣。原矩阵装配与软接触初始化、接触迭代、GPU输入读回、提交读回均另存。\n"
    ]
    for r in cost["trajectory_costs"]["rows"]:
        t = r["timing"]
        report.append(
            f"{r['case']}：积分＋质量 {t['integration_quality_wall_s']:.6f} s，保存读回 {t['sync_readback_wall_s']:.6f} s，HDF5 {t['hdf5_output_s']:.6f} s，详细诊断序列化 {t['observer_output_s']:.6f} s，驻留主体 {t['resident_body_s']:.6f} s。\n"
        )
    report += [
        f"六个远程冷进程合计 {cost['full_process_wall_s']:.6f} s（包含原始微基准、更正计时、检查与故意拒绝）。它是实际计算进程成本，不是全部研发时长或纯GPU时间；CUDA轨迹跨度包含CPU间隙。未做相同完整时域、相同保存与相同模型检查的无插桩生产对照，不能给出新模式加速倍数。新完整10 ms与独立科学样本均0，合格生产样本单位成本无定义。\n",
        "11. **哪些尚未物理校准？**原近场声学、润滑间隙下限、接触刚度、固定法向与线性化重叠仍是近似。短窗中出现1／16个“真实已分离但线性化仍有排斥”的泡对步，最大2.428／7.047 pN；不是精确几何单边弹簧。真实重叠另行检查，没有偷偷放宽17.8 nm；这些偏差和最终位置下的接触力矩保留于TRAJECTORY_AUDIT。\n",
        "12. **哪些获准入？**离散润滑参考、小系统求解、B2/2kc两静态和两短窗获限定数值PASS；B1跨两状态未通过，B2半刚度／基准刚度未同时通过；新模型dt收敛、完整10 ms、生产GPU批量和AI数据准入均未获得。\n",
        "13. **下一轮是否值得生产GPU模式？**值得继续有限方法研究；现在不宜直接做生产批量集成。先处理真实几何接触与固定法向误差，并验证同一B2模型10/5 µs的逐泡／事件收敛；再用有界矩阵无显式装配或局部预条件GPU方案与CPU参考比较成本。\n",
        f"14. **四组图在哪里？**{BASE / 'figures'}：Fig_F10G_1_lubrication_convergence、Fig_F10G_2_implicit_contact_reference、Fig_F10G_3_dense_cluster_contact、Fig_F10G_4_accuracy_cost_admission。均PNG/PDF；Figure 3另有contact_network补图，使用真实固定ID邻接矩阵。图使用真实半径1.78 µm、等比例3D、英文与实际数据，无虚构轨迹或平滑补帧。数据／源码哈希在FIGURE_PROVENANCE；打开审查在VISUAL_REVIEW。\n",
        "15. **科研进度文档？**已最小插入F10G段落，备份microbubble_aggregation.md.bak_F10G_20261011。before/after、diff、SHA和原标题结构核对在results/fast_pzt/f10g/progress_document；原F10F段落与原文其他字节保留。\n",
        "16. **F10H重点？**先在相同近场声力下研究真实几何单边接触、法向更新和同模型dt收敛，再评估低成本残差受控润滑／接触GPU实现。不要直接增加kc、继续堆Jacobi轮数、扩大种子或将新形态当作实验真值。\n",
        "## 复现、归档与限制\n",
        "WSL使用共享Python，PYTHONPATH=src:.。分析入口 scripts.fast_pzt_f10g.analyze，绘图 figures，验证 verify；不需要重新推进轨迹来重画图。remote从固定Git archive部署唯一Supervisor任务；archive回传并逐项核验实际输入、源码、结果与日志。当前没有GitHub推送请求，只做本地分支归档。\n",
        "数学推导、文献可访问性和兼容记录字段说明见FAST_PZT_F10G_CONTACT_MODEL.md。FIGURE_PROVENANCE与REMOTE_VERIFICATION保留执行提交、SHA与环境。旧接触观测器的部分分配配置字段仍写PBD/graph_chunk=1；真实执行身份由F10G_B2_K2_DT10、OBSERVER和SOFT_STEPS明确，实际不运行PBD或CUDA Graph。不能拿通用记录器的短窗valid_samples_per_gpu_hour当生产吞吐。\n",
    ]
    (ROOT / "docs/FAST_PZT_F10G_DEVELOPMENT.md").write_text("\n".join(report))
    # Read the full current document; insert once before its original epilogue.
    target = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = target.with_name(target.name + ".bak_F10G_20261011")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    insertion = """
**F10G先把润滑方程解清楚。**F10F已经表明，减小步长再减少接触轮数会让这两例更早失败；这次不继续扩大小步搜索，而是把原两轮Jacobi所近似的32泡法向阻力写成96维离线系统。保留原力、CFD、邻接和间隙下限，两轮的块预条件残差仍约为自由速度尺度的18.56%／15.19%。同一首错状态改为四轮后，原12轮PBD的残余由约17.831 nm降到14.842／13.939 nm；离线充分收敛参考则为5.260／5.290 nm。润滑未充分收敛确实影响这些接触预测，但这只是原离散近似的代数检查，不是实际水动力误差测量。

**软接触必须从过阻尼关系重新推导。**沿用 U=U_free+M F_contact，把两泡固定法向的线性排斥写成 delta_next=max(delta_pred,0)/(1+dt kc M_rel)，而非直接照搬[惯性XPBD](https://mmacklin.com/xpbd.pdf)的dt²与质量矩阵。根据保存的法向压紧力／8.9 nm，先登记kc约0.02321 N/m，再比较半倍、原值与两倍；这不是Marmottant径向壳弹性，也不是实测SonoVue接触刚度。90项两泡与18项四泡检查通过独立活跃集参考，成对反力平衡、没有吸引性接触力，但通过小系统求解不等于所有构型都满足原重叠预算。

**有限短窗走过去了，内部排列仍不能当真值。**只把同时通过两个32泡静态门槛的“同一润滑逆矩阵＋2kc隐式接触”接到原GPU声场和CFD上，CPU求解接触，不再叠加PBD。从共同原始初态出发，O/1102完成5 ms，准备1104完成2.35 ms，共735个接受步；最大重叠15.260／7.244 nm，单次修正2.951／0.721 nm，原17.8／178 nm预算保持。与旧算法共同接受前缀比较，质心最大差仅1.566／0.259 nm，但固定ID位置最大差达3.657／2.515 µm，Kabsch约0.917／1.295 µm，边事件也从124／104变为29／72。不同模型不应被要求逐位相同，这些数值更不能证明哪种形态更真实；新模型自身的时间步收敛尚未验证。

**通过预算仍留下明确近似。**当前每步冻结法向并使用线性化重叠，和真实球面重叠最多相差0.137／0.190 nm；有1／16个泡对步在几何分离后仍有小排斥，最大2.428／7.047 pN，已如实记录。因此不是精确的几何单边壳接触模型。两条短窗的积分与质量检查约1.669／0.654 s，详细诊断输出另约4.250／1.854 s；CPU小矩阵和GPU核计时分别保存，不能把插桩或短窗产率冒充生产加速。没有延长10 ms、没有增加新种子、没有改默认，也没有进行AI训练。

> **F10G证据与下一步：**docs/FAST_PZT_F10G_DEVELOPMENT.md、docs/FAST_PZT_F10G_CONTACT_MODEL.md和results/fast_pzt/f10g/保留实际方程、预登记刚度、独立参考、失败状态、两个短窗、完整逐泡与连边数据、实测成本、四组英文PNG/PDF及哈希。结论是限定接触数值检查通过，**MORPHOLOGY_NOT_CONVERGED / PRODUCTION_NOT_ADMITTED / PHYSICS_UNCALIBRATED，training_ready=false**。下一轮先处理法向和真实几何接触，并检查同一模型的时间步收敛，再判断是否值得做低成本GPU集成；不靠任意增大刚度或放宽预算追求成功轨迹。

"""
    anchor = "\n---\n\n### 尾声：故事里还没讲完的部分"
    raw = target.read_bytes()
    text = raw.decode()
    if backup.exists():
        before = backup.read_bytes()
        assert raw == before.replace(anchor.encode(), insertion.encode() + anchor.encode(), 1), (
            "Progress document changed after F10G insertion"
        )
    else:
        assert (
            hashlib.sha256(raw).hexdigest()
            == "f53da1d199086f3204aae2c7d21c4853fca8d67356ca9e377dd7274d360a42ed"
        ), "Review changed document before insertion"
        assert text.count(anchor) == 1 and "**F10G先把润滑方程解清楚。**" not in text
        before = raw
        backup.write_bytes(before)
        target.write_bytes(raw.replace(anchor.encode(), insertion.encode() + anchor.encode(), 1))
    after = target.read_bytes()
    assert after.replace(insertion.encode(), b"", 1) == before

    def titles(b):
        return [s for s in b.decode().splitlines() if re.match(r"^#{1,6} ", s)]

    assert titles(before) == titles(after)
    (out / "before.md").write_bytes(before)
    (out / "after.md").write_bytes(after)
    (out / "insertion.md").write_text(insertion)
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
