"""Minimal append preserving every original byte and all existing headings."""

from pathlib import Path
import difflib
import re
from .common import BASE, sha, write

PATH = Path("/home/lzy/projects/microbubble_aggregation.md")
BACKUP = PATH.with_name(PATH.name + ".bak_F10I_20261012")
BEFORE_SHA = "14b35ac8c531092e391a44a68a5d86a6a202b91181413639fee10eb463ca0015"
ADDITION = """

**F10I把130 µs差异拆开看。**沿用F10H几何B2的O/1102、5 ms、10/5 µs接受轨迹，独立重建全部泡对：29／29条事件身份与同泡对开关次序一致，但12条超出原10 µs门槛，另有一处跨泡对顺序反转。ID7–24最后断边前后，10 µs轨迹在1.650到1.660 ms由177.999758变为178.047093 nm；5 µs在1.785到1.790 ms由177.992241变为178.070134 nm。两者只以约4.73／15.58 µm/s离开178 nm图阈值，nm级间隙差就会放大成事件时间差。这个阈值仍有正表面间隙，不是球面接触或真实黏附。

**边开关敏感，整团在这个窗口没有分开。**原2.100R几何图保持首位，离线只改成2.095R／2.105R，事件数由29变为28／27，最大时差由130变为30 µs；这只是定义敏感性，两个扰动仍没有通过原严格门槛。三个阈值下32个固定ID始终是同一连通分量。原位置与形状的稳定证据保留：固定ID最大差10.412 nm、Kabsch 4.638 nm、Rg／长短比／张量相对差0.00636%／0.03765%／0.05983%。2 ms持续性下没有观察到分裂或合并，右删失0；因为没有事件可测，不能称持续分裂预测已经通过。

**不是把接触方程再解紧一点就能消除时差。**对14个真实保存状态收紧到1e−14 m全部收敛，1e−15 m只有13／14，另一个停在1.00258e−15 m并如实保留失败。紧解位置最多改变约0.0000975 nm，活动接触和边集合不变。为检查累计误差，只从原初态与完整控制补跑一条同10 µs、1e−14 m容差的2 ms轨迹，200步全部接受，逐泡最大位置变化0.00559 nm，23条前缀事件时间全不变，和5 µs比较仍有130 µs差。证据支持时间离散路径差与慢速过阈值的放大，不能据此宣称全部误差来源或时间阶数已解决；本轮没有运行2.5 µs。

**GPU真的试了，但没有进入生产。**矩阵自由阻力、块Jacobi-PCG和真实几何单边罚力原型先在两泡／四泡匹配CPU，两个原本超重叠预算的四泡负例也保持拒绝。真实32泡状态的位置与CPU仅差约2.33e−19 m，实际方程残差2.61e−15 m，但内部保守收敛量未达预登记值，线搜索耗尽，基础PCG也未达自身残差目标。单次GPU尝试35.46 ms，其中标量同步22.24 ms；CPU接触调用12.80 ms但未含阻力装配，因此没有报公平加速倍数。工程启动错误和数值失败均保留，按条件停止，没有追加长轨迹或批量。

> **F10I证据与下一步：**docs/FAST_PZT_F10I_DEVELOPMENT.md、docs/FAST_PZT_F10I_EVENT_DIAGNOSTICS.md、docs/FAST_PZT_F10I_GPU_FEASIBILITY.md及results/fast_pzt/f10i/保存预登记、29事件、全部接受状态、容差传播、GPU失败、成本、四组英文PNG/PDF和输入／执行／回传哈希。原 **RAW_EDGE_TIMING_FAIL / MORPHOLOGY_NOT_CONVERGED** 不改；团级独立记录 **GROUP_LEVEL_NUMERICAL_EVIDENCE**，持续事件预测为 **PERSISTENT_TOPOLOGY_NOT_EVALUATED**，GPU为 **GPU_B2_NOT_ADMITTED**，继续 **EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**。F10J先在同一固定状态审计GPU真实残差、停止语义和同步成本，数值与公平计时通过后再讨论集成；另行预登记任务级持久拓扑验证，不为消除130 µs无限提高精度。
"""


def main():
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if BACKUP.exists():
        assert sha(BACKUP) == BEFORE_SHA, "Unexpected existing backup"
        before = BACKUP.read_bytes()
    else:
        assert sha(PATH) == BEFORE_SHA, "Progress document changed since related section was read"
        before = PATH.read_bytes()
        BACKUP.write_bytes(before)
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    addition = ADDITION.encode()
    after = before.replace(marker, addition + marker)
    assert PATH.read_bytes() in [before, after], "Unexpected concurrent edit"
    assert after.replace(addition, b"", 1) == before

    def headings(b):
        return re.findall(r"^#{1,6} .+$", b.decode(), re.M)

    assert headings(before) == headings(after)
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
    PATH.write_bytes(after)
    write(
        out / "VERIFICATION.json",
        dict(
            status="PASS",
            path=str(PATH),
            backup=str(BACKUP),
            before_sha256=sha(BACKUP),
            after_sha256=sha(PATH),
            all_original_bytes_preserved=True,
            heading_structure_unchanged=True,
            heading_count=len(headings(before)),
            only_minimal_F10I_append=True,
        ),
    )
    print("Progress document appended; backup, before/after, diff and heading checks saved.")


if __name__ == "__main__":
    main()
