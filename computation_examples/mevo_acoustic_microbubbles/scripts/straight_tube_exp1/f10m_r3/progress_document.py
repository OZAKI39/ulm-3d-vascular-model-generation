"""One additive paragraph, original historical bytes and headings preserved."""

from pathlib import Path
import re
import difflib
from .common import BASE, read, write, sha


def main():
    assert read(BASE / "analysis/FINAL_GATE.json")["engineering_status"] == "PASS"
    p = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(p) + ".bak_F10M_R3_20261013")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "CHECK.json").exists():
        rec = read(out / "CHECK.json")
        assert sha(p) == rec["after_sha256"] and sha(backup) == rec["before_sha256"]
        return
    assert sha(p) == "a82e7d451e9025f5073c8621b57816fab8a670a974537c356558c4cc4587c5e1", (
        "Progress document changed since inspection"
    )
    assert not backup.exists()
    before = p.read_bytes()
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    text = "\n**F10M-R3 按新确认的四片PZT、双顶部孔和两端封闭继续修正实验装置。**32.5×9×9 mm的PDMS块、居中内径1.6 mm主通道和两侧3.7 mm厚度保留；现在四枚7 mm发声面分别位于Y负、Y正、X负、X正外表面中心，法向均朝内。两个顶部孔直径均1.5 mm、沿Z向下与主通道形成同一个空腔，关于X=0镜像，注入后保持几何开放，X两端则始终机械封闭；孔轴暂按端部发声面到孔中心5.4 mm放在x=±10.85 mm，孔中心/边缘测距基准仍为PROVISIONAL_MEASUREMENT_REFERENCE。六邻接连通性检查通过，孔没有穿透底部，端PZT面到相邻孔体积最小间隙为4.65 mm。新的CPU/GPU壁面查询排除交汇处被打通的旧圆柱壁，检查真实暴露侧壁、两端封闭面及交线距离，并保留有限泡半径与原0.1R修正预算；原开放X端逻辑只留给历史算例。24项CPU检查、四源权重输入检查和15项GPU检查通过，16泡0.1 ms无声零流轨迹逐泡静止、重叠及墙修正为零；一次通过的整组远端验证约7.91 s，OFF积分输出约1.55 s，采样整卡峰值440 MiB。孔内填充、顶部自由液面和外部介质仍未知，X端密封材料、厚度和胶层也未知，不能据机械封闭自动设PDMS或理想刚性声学壁；未知材料图保留NaN并拒绝求解。内部z+R<1.5 mm仅是数值验证范围，绝非指定液位，越界停止且不反射，孔内测试采用显式局部有液体演示假设。没有新增声波传播求解，四路正式声场、ON及导航均NOT_RUN，旧三路场未冒充四路；失败预检、Warp编译问题、修订图件和实际结果均保留，四组英文PNG/PDF、代码与回传SHA已归档。下一步先核对双孔液位与测距、端部封装、介质和单频声压标定，更新材料与资源准入、匹配APECSS及四路新场，再设计指定目标对照。当前仍为 **EXPLORATORY / MATERIAL_ACOUSTICS_UNCALIBRATED / SOURCE_PRESSURE_UNCALIBRATED / PHYSICS_UNCALIBRATED，training_ready=false**；证据见独立F10M-R3工作区的docs/FAST_PZT_F10M_R3_DEVELOPMENT.md与results/straight_tube_exp1/f10m_r3/。\n"
    i = before.index(marker)
    added = text.encode()
    after = before[:i] + added + before[i:]
    assert re.findall(rb"^#{1,6} .*$", before, re.M) == re.findall(rb"^#{1,6} .*$", after, re.M)
    assert after[:i] + after[i + len(added) :] == before
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
