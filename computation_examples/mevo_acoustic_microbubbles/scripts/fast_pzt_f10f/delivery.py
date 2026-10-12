"""Actual-data final report and byte-preserving append to the progress document."""

import difflib
import hashlib
import re
from pathlib import Path
from .common import BASE, ROOT, js, sha, write


def main():
    contact = js(BASE / "analysis/CONTACT.json")
    morph = js(BASE / "analysis/MORPHOLOGY.json")["comparisons"]
    cost = js(BASE / "analysis/COST_AND_YIELD.json")
    text = """# F10F 开发与最终证据

工程实现和有界研究交付完成。Small Steps 科学准入失败，不能作为默认生产模式。
继续 EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false。

## 身份与工程实现

基准741ad7d8683bfbb2fc9a270d9e92a218dead5e54，独立分支
dev/fast-pzt-f10f-smallsteps-contact-20261011。
科学与首轮软件检查提交9615a5a1a4e165493e7a4cee28f446a74c983734；
补充Graph首错和连边检查提交4f0e178e690d7ff7518abb0b64ebcefaac356fc7。
服务器固定源码位于/workspace/mevo_acoustic_microbubbles/f10f/.execution/。
Supervisor标签checks_v1、short_v1、failure_checks_v1，未中断其他项目；
整卡检查启动前空闲、运行采样至多一个compute PID、结束为空。
宿主PID在容器内可能没有进程名，不据此伪造进程归属。
输入13项、远程源码快照和回传文件逐项核验SHA-256。

新增f10f_smallsteps.py复用已有物理引擎；F10C recorder仅增加默认关闭的显式接口。
每个物理子步完整更新场、响应、一次力、二次力、两轮局部润滑和实际CFD，
保留原PBD、壁面及边界处理。没有使用多速率缓存、软接受、冻结团或新物理力。
每次Graph只推进一个子步；178nm单轮预算提前到GPU接受前检查，仍用接触错误位2048。
失败不推进接受时钟、不提交位置和活跃状态；正确保存奇数轮tmp2最终提案。
最后合法状态、准确时间、所有子步标量账本与末32步完整逐轮向量分别保存。

控制必须落在固定子步网格，非对齐时刻明确拒绝，不能悄悄延迟。
单泡/双泡独立NumPy参照、Graph开关逐位位置一致、子步内部控制切换、
真实连边跨越时刻、真实超修正拒绝和首错后不继续均通过。
SS2/SS4的Graph首错提案与诊断直接执行逐位相同。
B0重放两例的原始接受前缀逐位等于F10E，且复现原首错时刻。
B>1在新模块明确禁止，F10D仍是默认批量路径；没有声称已验证新模式的批量成员隔离。

## 接触数值结果

| 案例 | 模式 | 首次拒绝结束/ms | 本步预测侵入/nm | 最后残余/nm | 本步最大单轮修正/nm |
|---|---|---:|---:|---:|---:|
"""
    for r in contact["rows"]:
        text += f"| {r['case']} | {r['mode']} | {r['first_rejection_s'] * 1e3:g} | {r['final_predicted_nm']:.6f} | {r['final_residual_nm']:.6f} | {max(r['final_passes_nm']):.6f} |\n"
    text += """
全部被接受子步满足原17.8nm/178nm预算，表中是首次拒绝提案，不是接受状态。
六条重放墙面修正均为零。首错表的预测值对应不同历史和时刻，不是同状态对照。

完全同初态的O第一步：B0预测85.544019nm，SS2为50.535469nm，
SS4为32.978429nm，分别减小40.9246%、61.4486%。但6/3轮修正后仍超限。
同一个B0拒绝预测做只读Jacobi重放：一轮后O为66.191495nm、
准备为55.275530nm，十二轮后都约17.831nm；NumPy逐轮坐标与GPU差<1e-15m。
局部重放没有推进时间，也不是SS12实验。残余可能中途增大，不能假设逐轮单调。

## 逐泡形态：只比较共同接受前缀

| 准备1104比较 | 共同前缀/ms | 质心最大/nm | 固定ID最大/nm | Kabsch最大/nm | Rg相对最大/% | 长短比相对最大/% |
|---|---:|---:|---:|---:|---:|---:|
"""
    for r in morph[:3]:
        v = r["values"]
        name = r["reference"].split("_")[-1] + " / " + r["candidate"].split("_")[-1]
        text += f"| {name} | {r['accepted_common_end_s'] * 1e3:.2f} | {v['centroid_max_m'] * 1e9:.6f} | {v['fixed_ID_max_m'] * 1e9:.6f} | {v['Kabsch_max_m'] * 1e9:.6f} | {v['Rg_max_relative'] * 100:.6f} | {v['aspect_max_relative'] * 100:.6f} |\n"
    text += """
SS2对B0超过178nm固定ID和89nm Kabsch门槛；SS4的较短前缀虽通过位置门槛，
接触路径、连通成员和边事件仍失败。SS2/SS4互比也未过接触路径和边事件。
B0/SS2在2.10ms内边变化91/94；B0/SS4在1.70ms内61/57。
最大团比例相同不能替代成员相同。完整门槛和每帧误差见analysis/MORPHOLOGY.json。

准备过程观察到合并；2.10ms内B0/SS2各一件持续合并，另25/26件右删失；
1.70ms不足以验证2ms持续性，未观察到分裂。持续性准入明确NOT_EVALUATED_RIGHT_CENSORED，
不是零事件就自动PASS。这里是几何连通分裂/合并，不是物理破裂/融合。

O/1102的SS只有初态，因此对B0、C1的演化形态均NOT_EVALUATED。
没有生成B0最后合法4.89ms（最后保存4.85ms）之后的参考。
历史C1仅为数值参考：共同4.85ms前缀中质心差0.882nm，
固定ID差2.046µm、对齐RMSD0.832µm。不能据此降低本轮标准。
缺少同模式匹配OFF控制，声学配对增量NOT_EVALUATED。
总体MORPHOLOGY_NOT_CONVERGED；软件检查不等于微泡壳接触/近场力实验校准。

## 实际成本和有效产出

"""
    text += f"""六条科学诊断主体合计{cost["diagnostic_body_s"]:.6f}s，共用静态准备{cost["static_prepare_s"]:.6f}s。
三个冷进程（含软件检查与失败）合计{cost["cold_all_process_s"]:.6f}s。
科学尝试1836个物理步、1830个接受步、6个首错步；
软件检查另123个尝试步（97+26），不作为独立科学样本。

| 诊断积分每10µs墙钟成本/ms | B0 | SS2 | SS4 |
|---|---:|---:|---:|
| O/1102 |3.109566|5.478180|8.210752|
| 准备1104 |2.983564|4.309574|7.069930|

上表含CPU观察且停止时长不同，不是同窗口无插桩速度比较。
每10µs完整计算中声场/成对力调用1/2/4，CFD调用2/4/8，
润滑2/4/8，PBD均12轮；实际全历程数量见CONTACT.json。
CUDA跨度包含主机间隙，不是纯kernel时间。没有候选通过短窗门槛，
因此没有开展正式子kernel计时、无插桩完整轨迹、SS12或批量测试。

旧C1（5µs/12轮）有一条10ms完成、历史主体2.212250s；
SS2/SS4无完整结果，只能说C1有数值可完成证据，不能给出匹配速度倍数，
也不能认定C1形态或物理通过。

新完整10ms轨迹0、新独立seed0、复用独立seed家族2、完整同模式OFF/S/O三联0。
新合格轨迹与三联每预留GPU小时均0；单个新合格样本成本因零产出而无定义。
不同数值模式保留原sample_id和split_group，不能混成同一种训练标签。

仅作历史累计：F10D/E/F仍只有原三条完整控制轨迹，进程成本
{cost["cumulative_F10D_E_F_process_cost_interval_s"][0]:.6f}–
{cost["cumulative_F10D_E_F_process_cost_interval_s"][1]:.6f}s，
约{cost["cumulative_F10D_E_F_existing_trajectories_per_hour_interval"][0]:.2f}–
{cost["cumulative_F10D_E_F_existing_trajectories_per_hour_interval"][1]:.2f}条/预留GPU小时。
这是有限样本历史账本，不是F10F新产量、纯GPU时间或稳定生产率。

## 交付、复现与F10G

四组English PNG/PDF位于results/fast_pzt/f10f/figures/：
Fig_F10F_1_smallsteps_contact、Fig_F10F_2_morphology_convergence、
Fig_F10F_3_cost_and_yield、Fig_F10F_4_method_admission。
使用实际半径1.78µm、统一3D尺度和固定ID颜色；无虚构平滑帧。
已实际打开PNG及PDF页面渲染；初版底部图注重叠已修复，最终哈希记入VISUAL_REVIEW.json。

/home/lzy/projects/microbubble_aggregation.md只做最小插入，保留原字节、标题、公式和引用；
备份为microbubble_aggregation.md.bak_F10F_20261011。
before/after、插入、diff、SHA及标题检查位于results/fast_pzt/f10f/progress_document。

复现分析与绘图：PYTHONPATH=src:. python -m scripts.fast_pzt_f10f.analyze；
同样运行figures、verify模块。科学计算通过remote模块提交新的唯一Supervisor标签，
使用固定Git快照；不要为重建图片重复科学轨迹。
重要源码、数据、失败、日志和图已回传WSL，并在本地分支归档；本任务不包含GitHub推送。

F10G优先诊断密集接触中局部润滑的迭代收敛，区分迁移率误差、接触投影与近场吸引。
接触预判可随后有界研究，保持原接受限制和逐泡验收。当前没有生产化任何SS候选的证据。
"""
    (ROOT / "docs/FAST_PZT_F10F_DEVELOPMENT.md").write_text(text)
    audit = ROOT / "docs/FAST_PZT_F10F_SMALL_STEPS_AUDIT.md"
    old = audit.read_text().split("\n## Experimental outcome")[0]
    audit.write_text(
        old
        + """
## Experimental outcome

All four SS2/SS4 attempts failed earlier than the corresponding B0 failure.
Original contact budgets, correct odd-round proposal capture, exact B0 prefixes,
CPU/GPU oracles and graph first-error checks passed engineering validation.
The contact and morphology gates did not admit a candidate for production.
This two-case negative result is not a refutation of inertial XPBD or an
experimentally calibrated claim about microbubble shells. Actual tables and
cost limits are in FAST_PZT_F10F_DEVELOPMENT.md and analysis JSON files.
"""
    )
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = path.with_name(path.name + ".bak_F10F_20261011")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    insertion = f"""
**F10F把一个大步拆小，是否就能少做几轮接触？**这次参考[Macklin等2019年的Small Steps](https://mmacklin.com/smallsteps.pdf)，但论文主要讨论惯性XPBD，当前泡泡用的是过阻尼Stokes平移和同步Jacobi位置修正，不能直接搬用论文的性能倍数。保留原声场、近场力、两轮局部润滑、CFD和壁面，每个真实子步都重新计算，再分别尝试10 µs／12轮B0、5 µs／6轮SS2和2.5 µs／3轮SS4；三者每10 µs都是12轮接触，却需要不同次数的完整受力和流场查询。

**侵入变小了，少做接触迭代却让失败更早。**准备1104原在2.240 ms拒绝，SS2和SS4分别提前到2.110和1.745 ms；O/1102原在4.900 ms拒绝，小步候选却在第一个5或2.5 µs子步就失败。从完全相同的O初态看，预测侵入由85.544 nm减到50.535和32.978 nm，约减少40.9%和61.4%，但6轮和3轮后仍余20.920和20.361 nm，超过原17.8 nm。单轮修正都没超过178 nm，墙面修正为零。每个首错立即停止，最后合法状态与拒绝提案分开保存，没有把非法子步推进到大步末尾再称通过。

**质心接近，仍然不等于每颗泡泡排列稳定。**准备1104的SS2与B0共同接受到2.10 ms，最大质心差仅0.00718 nm，Rg和长短比相对差约0.1575%和0.3338%，固定ID位置最大差却达255.470 nm、对齐RMSD达94.579 nm，超过事先登记的178和89 nm门槛。SS4的1.70 ms短前缀位置差较小，但接触路径、成员及连边仍不一致；SS2与SS4互相比也没有通过完整形态准入。O的小步没有新的接受状态，不能拿共同初态当作形态通过。准备过程虽观察到合并，多数事件还没观察满2 ms，持续性保留为右删失，不能宣称没有变化。

**这条小步路线本轮不进入生产。**软件检查、CPU/GPU单泡和双泡对照、控制切换、实际连边时间、Graph首错停止和B0历史回归通过，但四条小步科学尝试都没过原失败窗口。因此按门槛停止，没有继续SS12、10 ms或批量扩展。本轮六条诊断重放共尝试1836个物理步，含检查的三个冷进程总计{cost["cold_all_process_s"]:.3f} s，新增完整10 ms合格样本和完整OFF/S/O三联均为0。诊断每10 µs的耗时增大，但包含逐步读回且终止时长不同，不能写成正式速度倍数；旧5 µs／12轮C1仍只有一个完整案例，也没有因此通过逐泡形态校准。

> **F10F证据与下一步：**独立分支、docs/FAST_PZT_F10F_DEVELOPMENT.md、docs/FAST_PZT_F10F_SMALL_STEPS_AUDIT.md和results/fast_pzt/f10f/保存全部接受前缀、首错、逐轮数据、真实成本、四组英文PNG/PDF及源码和回传哈希。结果是 **工程检查通过，Small Steps接触准入失败／MORPHOLOGY_NOT_CONVERGED**，继续 **EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**。F10G优先检查密集接触下局部润滑的迭代收敛与近场吸引，再有界研究接触预判；不要靠无界增加子步数或放宽重叠预算掩盖问题。

"""
    anchor = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    if backup.exists():
        before = backup.read_bytes()
        assert path.read_bytes() == before.replace(anchor, insertion.encode() + anchor), (
            "Progress changed; inspect"
        )
    else:
        before = path.read_bytes()
        assert (
            hashlib.sha256(before).hexdigest()
            == "ceff73a68b497c1e4e50378341d09a4c695d40f785ea45ff99c010c88735a887"
        )
        assert before.count(anchor) == 1
        backup.write_bytes(before)
        path.write_bytes(before.replace(anchor, insertion.encode() + anchor))
    after = path.read_bytes()
    assert after.replace(insertion.encode(), b"", 1) == before

    def headings(s):
        return re.findall(r"^#{1,6}.*$", s.decode(), re.M)

    assert headings(before) == headings(after)
    (out / "before.md").write_bytes(before)
    (out / "after.md").write_bytes(after)
    (out / "F10F_INSERTION.md").write_text(insertion)
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
    write(
        out / "VERIFICATION.json",
        dict(
            status="PASS",
            original_bytes_preserved=True,
            headings_preserved=True,
            headings=headings(before),
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(path),
            inserted_sha256=sha(out / "F10F_INSERTION.md"),
        ),
    )
    print("Report and minimal progress update complete")


if __name__ == "__main__":
    main()
