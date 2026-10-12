"""Minimal insertion with a verified original-byte backup and unchanged headings."""

from pathlib import Path
import difflib
import re
from .common import BASE, sha, write, js

PATH = Path("/home/lzy/projects/microbubble_aggregation.md")
BACKUP = PATH.with_name(PATH.name + ".bak_F10J_20261012")
BEFORE_SHA = "b24cce46e52627b8039e283c3050ff6bcfd3cc2c299223ead6505af719714918"


def main():
    gate = js(BASE / "analysis/FINAL_GATE.json")
    cost = js(BASE / "analysis/PERFORMANCE_SUMMARY.json")
    assert gate["mathematical_convergence"] == "GPU_FIXED_STATE_NUMERICAL_PASS"
    addition = f"""

**F10J找到了GPU停住的具体原因。**原32泡状态逐位复现F10I失败。基础PCG和第一次接触方向都在128轮后未达到自身残差目标；基础递推／真实残差约1.255e−12／1.254e−12，目标5.190e−13，两者并没有严重脱离。原阻力和方向矩阵保持正定，块预条件器合法。最后一次方向本身合格，但相近总能量计算的舍入让Armijo把实际下降步误判为上升，20次回溯全被拒绝。原1.850e−14 m保守量不是实际位置误差，也不需要把1e−14 m门槛调宽。

**只修求解数值，不改变接触物理。**沿用原求解器，把能量比较写成数学等价的稳定增量，PCG保留原绝对／相对目标、增加有界迭代余量，并同时核验重新计算的真实残差。原18项小系统检查通过，其中两项超预算四泡负例仍拒绝；真实32泡也通过，保守界约7.886e−16 m，实际位置方程残差约1.107e−16 m，相对独立参考的位置／接触力差约1.156e−16 m／5.852e−20 N。重叠17.49757 nm、最大单次修正1.93169 nm，原17.8／178 nm预算及2kc不变。故意不足的PCG、非有限输入、重合球心与严重重叠都保留拒绝能力；这些检查仍只是同一探索性离散模型的数值证据。

**同步减少了，GPU仍没有更快。**复用已接受候选的残差，并使用实际Warp 1.18支持的有界条件Graph后，所有19个参考状态的位置与修正版逐位一致，中途标量读回35→17→12。可是读回等待包含PCG计算，不能都当作可删除的复制开销。三组相同初态、准备、原精度、质量检查和输出的完整成本中位数：旧失败GPU约{cost["medians"]["original"]["median_wall_ms"]:.2f} ms、修正版{cost["medians"]["corrected"]["median_wall_ms"]:.2f} ms、优化版{cost["medians"]["optimized"]["median_wall_ms"]:.2f} ms；包含阻力装配且满足同一精度目标的CPU约{cost["medians"]["cpu"]["median_wall_ms"]:.2f} ms。优化GPU仍约贵{cost["new_gpu_over_cpu_cost"]:.2f}倍，因此本轮停止GPU B2生产集成，保留独立模块和CPU研究参考，没有追加长轨迹或批量。

**持续拓扑有候选，还没有B2验证。**对历史保存数据沿用2.1R图和2 ms持续性重建，找到F5关声后417.5 ms的真实初态，旧M0轨迹在418.54 ms发生持续分裂，可供下一轮登记约4 ms短窗。32个ID、原坐标、控制、CFD、壁面与哈希已保存；但该状态贴壁，M0也不是几何B2，必须先核验合法重启、壁面耦合与原预算，不能拿旧分裂时刻当B2真值。原几何B2的三条历史窗口仍始终是一团，本轮没有新增拓扑轨迹，持续分裂／合并预测未验证。

> **F10J证据与下一步：**docs/FAST_PZT_F10J_DEVELOPMENT.md、docs/FAST_PZT_F10J_GPU_SOLVER_AUDIT.md、docs/FAST_PZT_F10J_PERFORMANCE.md与results/fast_pzt/f10j/保存旧失败逐步残差、稳定能量修正、实际GPU参考与负例、三组公平计时、独立进程成本、历史候选、四组英文PNG/PDF和执行／回传哈希。分层记录 **GPU_FIXED_STATE_NUMERICAL_PASS / GPU_B2_PRODUCTION_NOT_ADMITTED**；原团级数值证据保留，严格连边时间仍FAIL，**PERSISTENT TOPOLOGY VALIDATION NOT RUN / EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**。F10K优先用CPU参考准备有界的持续拓扑与近壁耦合验证，不继续追逐130 µs瞬时边时差或缺乏成本优势的GPU生产化。
""".encode()
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if BACKUP.exists():
        assert sha(BACKUP) == BEFORE_SHA
        before = BACKUP.read_bytes()
    else:
        assert sha(PATH) == BEFORE_SHA, "Unexpected concurrent edit; preserve it"
        before = PATH.read_bytes()
        BACKUP.write_bytes(before)
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    after = before.replace(marker, addition + marker)
    assert PATH.read_bytes() in [before, after]
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
            only_minimal_F10J_append=True,
        ),
    )


if __name__ == "__main__":
    main()
