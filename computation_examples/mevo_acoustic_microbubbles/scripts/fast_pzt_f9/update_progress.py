"""One-time minimal F9 insertion; preserve historical narrative and heading order."""

from pathlib import Path
import difflib
from .common import BASE, js, sha, write


def main():
    assert js(BASE / "FINAL_GATE.json")["ledger_closure"] == "PASS"
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = path.with_name(path.name + ".bak_F9_20261011")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "UPDATE.json").exists():
        assert sha(path) == js(out / "UPDATE.json")["after_sha256"]
        return
    before = path.read_text()
    assert "**F9 把每一步的位移记成一笔账。**" not in before
    if backup.exists():
        assert backup.read_bytes() == path.read_bytes(), "Existing backup must not be overwritten"
    else:
        backup.write_bytes(path.read_bytes())
    block = """**F9 把每一步的位移记成一笔账。**这次不改声音、血流或推开泡泡的方法，只在 GPU 里记录：血流带了多少、一次和二次声力各推了多少、润滑与平均平衡改了多少、泡泡接触推开多少、墙面拉回多少，最后与真正接受的位移相加核对。自由运动、平面墙和真实 Y 管检查都通过；打开审计后，固定控制的两种步长和反馈的 10 µs 轨迹都逐帧、逐位重现 F8。这里是同一条数值轨迹的账本，不是把墙面那一项减掉，就得到了“没有墙时会怎么走”。

**固定声音下，为什么明明向下推，团却向上走？**10 µs 的 Y1 中，300 ms 内血流的平均 Z 位移是 −0.00547 µm，一次声力是 −0.15688 µm，墙面带符号的净修正是 +0.31723 µm，合起来正好得到实际的 +0.15489 µm；二次力、润滑平衡和中央接触修正在整团平均上接近零。后者仍会改变内部排列，不能说它们没作用。墙面累计修正长度为 0.31757 µm，是另一个量。现在确认的是几何约束怎样组成数值位移，还不是测到了真实近壁推进力。

**步长减半，墙面的这笔账变大了。**Y1 的净 Z 墙面修正变为 +0.52646 µm，比粗步长多 +0.20922 µm，实际末端 Z 位移多 +0.20918 µm。差异在约 90–110、190–210 ms 等多个阶段增长；原网格的两个面没有改，主要接触泡与面切换时刻却不同，细步长还多出 ID 23 的接触。接触推开与墙面拉回会反复交替。这个证据说明结果对当前接触和墙面处理敏感，不能只凭修正大就宣布实现错误，也不能凭两种步长声称已经收敛。

**关声对照补齐了，导航结论仍不改。**Y0 的 10／5 µs 过关位置只差约 0.00105 nm；Y1、Y2 分别差约 218.69、13.56 nm。反馈相对关声多出的 +Z 进展约 0.1065 µm，超过这组配对的数值门槛；但相对固定声音没有额外进展，自己在两种步长下也都没进入 0.15 µm 的容差。Y2 自身的账本还显示一次声力确实向 +Z 推了约 0.09731 µm，墙面另加约 0.05463 µm，所以也不能说所有上移都是墙面造成的。当前结论是 **WALL_SENSITIVE / NAVIGATION_NOT_RESOLVED**，可以确认开声组保持连通，仍为 **COHESION_ONLY / PHYSICS_UNCALIBRATED**。

> **F9 证据：**`docs/FAST_PZT_F9_DEVELOPMENT.md` 与 `results/fast_pzt/f9/`，包含精确 GPU 位移账本、墙面前后位置与三角面、必要的 Y0 半步对照、审计重放及四组英文 PNG/PDF。F8 的历史标签保留；重放不是新的独立样本，当前不能据此批准 AI 训练。

"""
    anchor = "\n\n\n---\n\n### 尾声：故事里还没讲完的部分"
    assert before.count(anchor) == 1
    after = before.replace(anchor, "\n\n" + block + anchor, 1)
    anchor2 = "\n---\n\n### 一句话总结整个故事"
    assert after.count(anchor2) == 1
    cost = """
9. **F9 的审计成本与下一步：**六个串行科学进程累计约 4.19 分钟；其中只新增一条 Y0 半步对照，另外三条完整轨迹是审计复现。同一段 30 ms 近壁运动，关闭审计积分约 2.46 秒，开启后约 4.25 秒，增加约 73%，所以生产中默认关闭。GPU 只在原保存点批量传回账本，没有每步把全部坐标搬到 CPU；关键科学文件约 68.8 MB 已归档回 WSL。下一步先建立完整状态重启，再用这些实际接触窗口分开检查步长、接触与墙面处理顺序、原三角面近似的影响。RPY 接触仍未准入，暂不做昂贵长轨迹；也没有偷偷增大声力、改变目标或放宽预算。
"""
    after = after.replace(anchor2, cost + anchor2, 1)

    def headings(s):
        return [line for line in s.splitlines() if line.startswith("#")]

    assert headings(before) == headings(after)
    remaining = iter(after.splitlines())
    assert all(any(line == candidate for candidate in remaining) for line in before.splitlines())
    (out / "before.md").write_text(before)
    (out / "after.md").write_text(after)
    (out / "UPDATE.diff").write_text(
        "".join(
            difflib.unified_diff(
                before.splitlines(keepends=True),
                after.splitlines(keepends=True),
                fromfile=str(backup),
                tofile=str(path),
            )
        )
    )
    path.write_text(after)
    write(
        out / "UPDATE.json",
        dict(
            status="APPLIED",
            path=str(path),
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(path),
            headings_unchanged=True,
            heading_structure=headings(after),
            historical_lines_preserved_in_order=True,
            changes="Two F9 insertions only; original F8 statements remain explicitly historical. No formula, reference or heading rewritten.",
        ),
    )


if __name__ == "__main__":
    main()
