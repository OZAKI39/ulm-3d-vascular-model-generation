"""Backup and add exactly one paragraph, preserving all historical bytes."""

from pathlib import Path
import re
import difflib
from .common import BASE, read, write, sha


def main():
    assert read(BASE / "analysis/FINAL_GATE.json")["engineering_status"] == "PASS"
    assert read(BASE / "figures/VISUAL_REVIEW.json")["status"] == "PASS"
    p = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(p) + ".bak_F10N_20261013")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "CHECK.json").exists():
        rec = read(out / "CHECK.json")
        assert sha(p) == rec["after_sha256"] and sha(backup) == rec["before_sha256"]
        return
    assert sha(p) == "326f81f32f86668599a4457bdd01294b5f3225b8abbb14b06f4a8a210d7b2880"
    assert not backup.exists()
    before = p.read_bytes()
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    text = "\n**F10N 把研究重心转回四片PZT对微泡位置和形态的有限控制。**按本轮明确的简化，运动腔体视为充水、背景流严格为零，保留F10M-R3真实尺寸几何、双顶部孔、封闭X端及GPU壁面查询，暂不增加PDMS反射、液面、液滴、注射和声流。四枚7 mm圆片采用2.5 MHz均匀水近场积分代理，独立幅相、启停与持续时间可调；同一复压力的一致梯度经过独立面积积分和CPU/GPU核对，3及1.5 µm声场缓存复用于全部控制，没有重跑k-Wave或外推APECSS。32泡分散初态与具有接触网络的4×4×2合法演示初态分别保存，后者不是声学聚集产物；位置目标在运行前固定为(3,2,1) µm，窗口2 ms，形态目标为X/Y方向长短比从1增至1.25。OFF严格静止，BASELINE和OPTIMIZED质心位移模分别约0.0597及0.2282 µm，优化将目标误差从3.7417降至3.5414 µm，未达到导航目标，也没有证明独立Z平移。演示团保持几何连通，但优化位置控制的减半步长同ID差达6.416 µm、加密场差达3.041 µm；质心差很小不能证明内部排列可靠，该形态明确标为未收敛。弱形态控制使Rg增约11.107 nm、X/Y比增约0.4245%，超过本次步长/缓存差异，却远未达到25%目标；这些变化只属未校准接触模型中的数值响应。分散初态在同一优化控制下由32变为8个分量，不能称已聚成完整团。21条有限试验中20条完成、1条在0.4 ms有效前缀后触发原接触预算，失败提案保留，预算未放宽。21项CPU和11项GPU接口检查通过，CPU预检及缓存约7.54 s、远端整组56.42 s，其中轨迹积分26.98 s，整卡采样峰值480 MiB；145个服务器文件回传SHA一致，四组英文PNG/PDF及逐泡敏感性证据已归档。下一步先在XY可控子空间、较低驱动与逐泡误差门槛下研究短目标和少量分段反馈，再用真实轨迹优化去质心形状张量，不继续无必要的材料细化。当前保持 **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED，training_ready=false**；证据见独立F10N工作区的docs/FAST_PZT_F10N_DEVELOPMENT.md和results/straight_tube_exp1/f10n/。\n"
    added = text.encode()
    i = before.index(marker)
    after = before[:i] + added + before[i:]
    assert after[:i] + after[i + len(added) :] == before
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    backup.write_bytes(before)
    (out / "before.md").write_bytes(before)
    p.write_bytes(after)
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
    write(
        out / "CHECK.json",
        dict(
            status="PASS",
            path=str(p),
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(p),
            one_paragraph=True,
            headings_unchanged=True,
            all_original_bytes_preserved=True,
            diff_sha256=sha(out / "change.diff"),
        ),
    )


if __name__ == "__main__":
    main()
