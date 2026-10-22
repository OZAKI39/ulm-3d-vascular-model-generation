"""Data-derived reports and minimal, byte-preserving Chinese progress update."""

import re
import difflib
import hashlib
from pathlib import Path
from .common import BASE, ROOT, read, write, sha, identity


def main():
    an = BASE / "analysis"
    runs = read(an / "RUN_ANALYSIS.json")
    comps = read(an / "COMPARISONS.json")
    gate = read(an / "FINAL_GATE.json")
    cost = read(an / "COSTS.json")
    model = read(an / "THREE_MODEL_COMPARISON.json")
    x = runs["NAV_X_DT1"]
    y = runs["NAV_Y_DT1"]
    strong = runs["STRONG_DT1"]
    gx = gate["groups"]["NAV_X"]
    gy = gate["groups"]["NAV_Y"]
    docs = ROOT / "docs"
    lines = [
        "# F10Q three-model common-window comparison",
        "",
        "All controls and initial coordinates match bitwise. Only saved times 0 and 0.1 ms are exactly shared within the 0–0.174 ms legal intersection. No interpolation is used. Different contact and mobility approximations are compared, not reinterpreted as integration error.",
        "",
        "| Pair | COM difference (nm) | Maximum ID difference (nm) | Rg difference (nm) | Relative Q difference |",
        "|---|---:|---:|---:|---:|",
    ]
    for pair, v in model["contrasts"].items():
        lines.append(
            f"| {pair} | {v['COM_m'] * 1e9:.9g} | {v['fixed_id_m'] * 1e9:.9g} | {v['Rg_m'] * 1e9:.9g} | {v['Q_relative']:.9g} |"
        )
    lines += [
        "",
        "The largest cross-model COM contrast is 1.822883789e−5 nm, while Rg differs by up to 15.76497033 nm. This supports limited translation consistency for this initial state, this control and these two saved states only. It does not establish general contact-model independence or any 2/30 ms agreement. Shape is model-dependent.",
        "",
        "Within-model dt and cache differences at exactly the same saved times are separately stored in THREE_MODEL_COMPARISON.json. The original F10N full-window dt COM error is not used as the cross-model difference. Zero-flow OFF is exactly stationary in all inherited records and the new independent three-way OFF repeats; directional acoustic displacement equals displacement relative to that null.",
        "",
        "Model-specific contact quality remains separate: M0 uses its original PBD budget, B2 its original soft-overlap budget, and hard spheres their own complementarity tests. There is no hard-sphere reaction complementarity metric for M0 or B2. All three remain exploratory and uncalibrated.",
    ]
    (docs / "FAST_PZT_F10Q_GROUP_COMPARISON.md").write_text("\n".join(lines) + "\n")
    lines = [
        "# F10Q independent XY navigation",
        "",
        "A failure of the old strong control did not block B. Eight newly optimized controls were registered before simulation and independently completed 0.5 ms. The proposal optimizer reused the original cubic water-field cache, applied four independent complex source weights, checked actual phasor pressure, and used only a frozen-state COM predictor to propose controls. All success claims use actual 32-particle integration.",
        "",
        "Exactly two controls were selected from measured short-window displacement and a documented pressure-risk cost: X_P1200 and Y_P1800. Each was freshly run to 2 ms at dt=1 µs, repeated at 0.5 µs, and repeated with 1.5 µm cache spacing instead of 3 µm. OFF was independently repeated in all three settings.",
        "",
        "| Control | Final COM (µm) | Assigned error (µm) | Basic | COM | Shape | Graph | ID | Target |",
        "|---|---|---:|---|---|---|---|---|---|",
    ]
    for name, r, g in [("OFF", runs["OFF_DT1"], gate["groups"]["OFF"]), ("X", x, gx), ("Y", y, gy)]:
        disp = ", ".join(f"{v * 1e6:.9g}" for v in r["displacement_m"])
        axis = "X" if name == "OFF" else name
        statuses = [
            g[k]
            for k in [
                "BASIC_NUMERICS_PASS",
                "GROUP_POSITION_NUMERICALLY_ADMITTED",
                "GROUP_SHAPE_NUMERICALLY_ADMITTED",
                "GROUP_TOPOLOGY_NUMERICALLY_ADMITTED",
                "PARTICLE_ID_NUMERICALLY_ADMITTED",
            ]
        ]
        lines.append(
            f"| {name} | ({disp}) | {g['target_error_m'][axis] * 1e6:.9g} | "
            + " | ".join(statuses)
            + f" | {g['TARGET_REACHED'][axis]} |"
        )
    lines += [
        "",
        "Original targets remain X=(0.2,0,0) µm and Y=(0,0.2,0) µm, tolerance 0.075 µm. X reaches its original target in this same-model numerical prototype. Y does not. These are no-flow, homogeneous-water numerical controls, not calibrated experimental results.",
        "",
        "## Sensitivity and shape",
        "",
        "| Control / comparison | COM (nm) | ID max (nm) | Rg (nm) | Relative Q | Centered matched RMS/R |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key in ["NAV_X_DT05", "NAV_X_FINE", "NAV_Y_DT05", "NAV_Y_FINE"]:
        v = comps[key]
        lines.append(
            f"| {key} | {v['COM_m'] * 1e9:.7g} | {v['fixed_id_m'] * 1e9:.7g} | {v['Rg_m'] * 1e9:.7g} | {v['Q_relative']:.7g} | {v['matched_RMS_over_R']:.7g} |"
        )
    lines += [
        "",
        f"X has a resolved but tiny Rg change of {(x['rg_initial_final_m'][-1] - x['rg_initial_final_m'][0]) * 1e9:.9g} nm; its X/Y ratio ends at {x['xy_initial_final'][-1]:.12g}, so this is nearly pure translation, not substantial shape control. Its leading in-plane eigenvalues remain degenerate and a unique principal direction is NOT_IDENTIFIABLE.",
        f"Y coarse Rg changes by {(y['rg_initial_final_m'][-1] - y['rg_initial_final_m'][0]) * 1e9:.9g} nm and X/Y ends at {y['xy_initial_final'][-1]:.12g}, but shape differences of 96.4/306.7 nm and ID differences of 395.6/764.6 nm exceed original limits. This is an actual GROUP_POSITION pass with GROUP_SHAPE and PARTICLE_ID failures, not a synthetic metadata fixture. All initial 32-member graphs survive 2 ms; no physical fusion is implied.",
        "",
        "X pressure is 690.39–1113.33 Pa and Y pressure 1742.13–1801.09 Pa across actual sampled steps. Initial proposal caps are not hard runtime amplitude clipping; actual pressure is re-evaluated at every step against the unchanged APECSS range. Source amplitude is not bubble pressure. All source controls and per-step pressure extrema are archived.",
        "",
        "## Exact selected controls",
        "",
        "Order: Y−, Y+, X−, X+. Frequency 2.5 MHz; all sources enabled; duration 2 ms. Full precision is in NAVIGATION_SELECTION.json and each INPUT.json.",
        "",
        "| Control | Source amplitudes (Pa) | Phases (rad) |",
        "|---|---|---|",
    ]
    for row in read(BASE / "NAVIGATION_SELECTION.json")["selected"]:
        e = row["control"]
        lines.append(
            "| "
            + row["axis"]
            + " | "
            + ", ".join(f"{v:.12g}" for v in e["amplitude_Pa"])
            + " | "
            + ", ".join(f"{v:.12g}" for v in e["phase_rad"])
            + " |"
        )
    lines += [
        "",
        "No two-segment/5/10/30 ms experiment was needed: X already reaches the short target, while Y shape/ID convergence should be repaired before claiming dependable shape control. Group COM admission would permit considering a longer experiment, but it is not a requirement or evidence of morphology validity. No budget was spent merely to fill the allowance.",
        "",
        "CONTACT_MODEL_DEPENDENT is established for historical shape contrasts; cross-model dependence of these NEW controls is NOT_VERIFIED because they were not rerun with M0/B2. EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED, TRAINING_READY=false.",
    ]
    (docs / "FAST_PZT_F10Q_XY_NAVIGATION.md").write_text("\n".join(lines) + "\n")
    answers = [
        "开发前有界引用审计未发现新增可证明安全删除的对象，因此没有删除历史文件；清理审计与活跃索引单独提交。后续仅清除新绘图/分析代码中未使用的导入和变量。",
        "F1–F10P报告、M0/PBD、B2、原硬球、R3壁面、F2/APECSS、四PZT缓存、失败状态、历史入口和不可变执行快照全部保留；关键源码SHA复核一致。",
        "实际0.464 ms拒绝提案支持有限切向球面曲率和冻结法向机制；并非仅凭公式猜测。后续0.884 ms的新失败属于活动集循环，不能混为同一原因。",
        "ID 19/21的切向相对滑移9.261203 nm，法向位移−0.011488 nm，法向转动0.149052°；真实曲率间隙0.012046312 nm，二阶近似0.012046333 nm；原反力归一化后复现0.010056962 nm。",
        "每个物理步固定原x_n、F和阻力A，在候选位置重建真实距离与法向、重新线性化同一位移QP并更新反力；每次都从同一x_n形成x_new。独立检查局部KKT、真实法向受力平衡、互补、几何、修正及壁面，无重复推进时间或事后PBD。",
        "先保留三次QP仍不合格的冻结复核；依据残差连续下降，封顶增至五次。原冻结失败在第四次通过；新强轨迹每步1–5次，X完整轨迹1–2次，Y为1–4次。",
        "保持原F10G成对润滑阻力与过阻尼平衡，原声学/F2力仍在步初显式评估，步内只更新接触几何和反力。没有改水动力方程或引入弹簧。",
        "原17.8 nm重叠、178 nm修正、0.01 nm真实互补、局部数学残差、1 nm硬球诊断和原压力/壁面门槛全部不变。",
        "已越过0.464 ms，原初态重新积分的0.5 ms强控制通过；进一步的合法前缀为0.883 ms。",
        "没有完成强控制2 ms。0.884 ms提案出现active_set_cycle并拒绝，按登记顺序不再做强完整半步/细缓存；路线B照常运行。",
        f"新强轨迹已接受前缀：最大重叠{strong['accepted_max_overlap_m'] * 1e9:.8g} nm，真实互补{strong['accepted_max_true_comp_m'] * 1e9:.8g} nm，真实法向迁移率残差{strong['accepted_max_residual_m']:.8g} m，局部QP残差{strong['accepted_max_local_residual_m']:.8g} m，最大力残差{strong['accepted_max_force_residual_N']:.8g} N。拒绝提案单独保留，不能计入合法轨迹。",
        "强控制未做完整2 ms步长/缓存对照，逐泡、质心及形态完整一致性均NOT_RUN/未准入；不能用局部一步质量或旧前缀比较替代。",
        "精确共同保存时刻为0与0.1 ms，无插值。M0/B2、M0/硬球、B2/硬球最大COM差分别约0.000012427、0.000005803、0.000018229 nm。",
        "形态不一致：对应Rg差约13.255825、2.509145、15.764970 nm，最大相对Q差约0.5453%；同ID跨模型差可达18.445099 nm。跨模型差异与各自dt/cache差分开记录。",
        "只能支持该初态、该强控制和两个共同保存时刻下团级平移很接近；不能证明一般模型不敏感，更不能外推2或30 ms。",
        "实际独立运行了8个新的中等X/Y候选，均完成0.5 ms；选择X_P1200和Y_P1800各做完整2 ms、半步与细缓存三条独立轨迹。",
        f"OFF严格零位移；X质心位移({x['displacement_m'][0] * 1e6:.9g},{x['displacement_m'][1] * 1e6:.9g},≈0) µm，目标误差{gx['target_error_m']['X'] * 1e6:.9g} µm；Y质心位移(≈0,{y['displacement_m'][1] * 1e6:.9g},≈0) µm，目标误差{gy['target_error_m']['Y'] * 1e6:.9g} µm。目标均未改变。",
        "X真正达到原0.2 µm目标：完整窗口、基础质量和COM准入通过，误差≤0.075 µm。Y未达到。成功范围仅为当前未校准模型中的短距XY数值原型，不是实验导航校准。",
        "有：Y是真实的团级位置通过而逐泡失败案例；dt/cache最大同ID差395.608440/764.640172 nm，超过原50 nm。",
        "有：同一Y案例位置通过而形状失败；Rg差96.417259/306.729687 nm，相对Q差4.0428%/12.1646%，匹配分布也超限。X各层数值检查均通过，X/Y各自初始团始终连通。",
        "未尝试5 ms或两段控制。X短距目标已达到，Y的形态/ID误差尚未解决；先解决这些问题比延长轨迹更有价值。5/10/30 ms均NOT_RUN。",
        f"核心科学数值执行共{cost['total_numerical_s']:.6f} s：路线A {cost['routes_wall_s']['A']:.6f} s，路线B {cost['routes_wall_s']['B']:.6f} s；保留额度未耗用。远端四进程{cost['remote_process_s']:.6f} s，峰值RSS {cost['peak_RSS_MiB']:.3f} MiB。额外只读最终核验时间独立见VERIFICATION.json；不把研发/绘图耗时或CPU传输余量称为纯GPU时间。",
        "四组PNG/PDF在本独立工作区results/straight_tube_exp1/f10q/figures/，附绘图源码、FIGURE_PROVENANCE.json、逐图SHA、PDF渲染及VISUAL_REVIEW.json。",
        "进度文档备份为/home/lzy/projects/microbubble_aggregation.md.bak_F10Q_20261014；仅在原F10P段后追加中文F10Q段，保留全部历史、标题、公式和引用，before/after/diff/SHA及结构审计在results/straight_tube_exp1/f10q/progress/。",
        "下一轮先用可重放冻结状态解决活动集循环，并针对Y快速内部重排检查步长与场梯度敏感性；导航优化应加入实际接触/形态风险，再有限扩大X距离或探索较低压力Y控制。不先扩大长轨迹或声称任意形态控制。",
    ]
    text = "# F10Q 开发与科研交付\n\n已在独立分支 `dev/fast-pzt-f10q-sequential-qp-xy-navigation-20261014` 完成两条独立路线。基线 `1c785f4803ca7dc0ea2fe91756b02487c7c75a59`；原F10P工作区未覆盖。\n\n**主要结论：X原短目标达到；Y质心数值通过但目标未达到、形态与ID未收敛；旧强控制仍在0.884 ms拒绝。**\n\n"
    text += (
        "\n\n".join(f"{i}. {s}" for i, s in enumerate(answers, 1))
        + "\n\n## 分层状态\n\n| 项目 | 强控制2 ms | X控制2 ms | Y控制2 ms |\n|---|---|---|---|\n"
    )
    for key in [
        "BASIC_NUMERICS_PASS",
        "GROUP_POSITION_NUMERICALLY_ADMITTED",
        "GROUP_SHAPE_NUMERICALLY_ADMITTED",
        "GROUP_TOPOLOGY_NUMERICALLY_ADMITTED",
        "PARTICLE_ID_NUMERICALLY_ADMITTED",
    ]:
        text += (
            f"| {key} | " + " | ".join(gate["groups"][k][key] for k in ["STRONG", "NAV_X", "NAV_Y"]) + " |\n"
        )
    text += "| TARGET_REACHED | false | true (X) | false |\n\nENGINEERING_PASS由最终独立VERIFICATION.json确认。CONTACT_MODEL_DEPENDENT：历史形态对照成立，新控制跨模型依赖尚NOT_VERIFIED。EXPERIMENTALLY_CALIBRATED=NOT_VERIFIED；TRAINING_READY=false。持续保留EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED。\n\n具体方法与证据见同目录TRUE_GAP_AUDIT、SEQUENTIAL_QP、GROUP_COMPARISON和XY_NAVIGATION文档；原始轨迹、每步QP反力/位置、质量表、源快照、运行日志与哈希均位于results/straight_tube_exp1/f10q/。所有图实际打开检查后才完成最终验收。\n"
    (docs / "FAST_PZT_F10Q_DEVELOPMENT.md").write_text(text)
    # Preserve original bytes, headings, formulas and references by insertion only.
    target = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(target) + ".bak_F10Q_20261014")
    before = target.read_bytes()
    oldsha = "b7bf77bbf1a78d8d89652793768faf006a36ba281b532068328737e2b0123656"
    assert hashlib.sha256(before).hexdigest() == oldsha, "Progress changed; reconcile before writing"
    assert not backup.exists(), "Never overwrite stage backup"
    paragraph = """\n\n**F10Q 把硬球曲面修正与XY团级导航分为互不阻塞的两条路线。**在完整F10P独立工作区上先做有界清理审计，无新增可证明安全删除对象，原模型、失败和历史均保留。实际0.464 ms拒绝状态中，ID 19/21切向相对滑移9.261203 nm、法向转动0.149052°，真实球面曲率间隙0.012046312 nm，经原反力归一化复现0.010056962 nm失败值；二阶近似量级吻合，无配对或单位错误。新Sequential QP在同一物理步固定原声力与成对阻力，只重建候选几何、法向和反力，重复求解不重复推进时间，无弹簧或事后PBD；三次迭代的原失败保留，依据真实法向残差下降将上限封顶为五次，冻结状态第四次通过。原两泡/四泡和新增独立非线性曲面参考通过，所有原门槛不变。新强控制从原初态通过0.5 ms，继续到0.883 ms合法，在0.884 ms因原活动集QP循环拒绝，未完成2 ms，也未冒称强完整步长/缓存收敛。路线B独立复用三模型在0与0.1 ms的精确共同保存状态，无插值；M0/B2/硬球最大质心差约0.00001823 nm，而Rg差可达15.76497 nm，说明该短窗平移接近而形态模型依赖，不能外推长时。八个新的中等XY候选均完成0.5 ms，实际选两条各做2 ms、半步与细缓存，四PZT仍为2.5 MHz均匀水代理、32个固定ID、零流、原F2/APECSS与R3壁面。OFF严格静止；X方案移动约0.200896 µm，距原(0.2,0,0) µm目标约0.000896 µm，基础、团COM、形状、图连通和ID都通过，原0.075 µm容差下真正达到该数值短目标。Y方案移动约0.051776 µm，距原(0,0.2,0) µm目标0.148224 µm，未达到；其COM步长/缓存差仅约0.000100/0.003640 nm，但同ID差395.608/764.640 nm、Rg差96.417/306.730 nm、相对Q差4.04%/12.16%，是实际的团位置通过而形状和逐泡失败案例。两方案的初始32成员图保持2 ms，不是新聚合或融合；X的Rg仅减少0.295292 nm，Y的约127.980 nm下降没有超过完整数值不确定性，不能称可靠形态控制。核心科学数值执行A约13.742 s、B约125.473 s，共139.215 s；额外独立只读核验单列，四个独立Supervisor任务、116个远端文件逐项SHA回传，失败、真实步长/反力/几何及四组英文PNG/PDF归档并逐图目视检查。未运行分段或5/10/30 ms：X短目标已到，下一步先解决活动集循环及Y内部重排敏感性，再以接触风险约束优化Y或有限扩大X距离。数值可重复仍不等于实验物理校准，当前继续 **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED，training_ready=false**；详见独立F10Q工作区docs/FAST_PZT_F10Q_DEVELOPMENT.md及results/straight_tube_exp1/f10q/。\n"""
    marker = b"\n---\n\n### " + "尾声".encode()
    assert before.count(marker) == 1
    insert = paragraph.encode()
    position = before.index(marker)
    after = before[:position] + insert + before[position:]
    backup.write_bytes(before)
    target.write_bytes(after)
    out = BASE / "progress"
    out.mkdir(exist_ok=True)
    (out / "before.md").write_bytes(before)
    (out / "after.md").write_bytes(after)
    (out / "change.diff").write_text(
        "".join(
            difflib.unified_diff(
                before.decode().splitlines(True),
                after.decode().splitlines(True),
                fromfile="before.md",
                tofile="after.md",
            )
        )
    )
    headings = lambda b: re.findall(rb"^#{1,6} .*$", b, re.M)
    checks = dict(
        original_is_exact_subsequence=after[:position] + after[position + len(insert) :] == before,
        headings_unchanged=headings(before) == headings(after),
        backup_exact=backup.read_bytes() == before,
        one_insert_only=True,
    )
    write(
        out / "UPDATE_AUDIT.json",
        dict(
            all_checks_pass=all(checks.values()),
            checks=checks,
            before_sha256=oldsha,
            after_sha256=sha(target),
            backup=str(backup),
            heading_count=len(headings(before)),
            inserted_sha256=hashlib.sha256(insert).hexdigest(),
            identity=identity(),
        ),
    )
    print("Wrote three reports and preserved progress history:", checks)


if __name__ == "__main__":
    main()
