"""Append one Chinese paragraph, preserving all previous bytes and headings."""

import difflib
import re
from pathlib import Path
from .common import BASE, read, write, sha


def main():
    gate = read(BASE / "analysis/FINAL_GATE.json")
    assert gate["engineering_status"] == "PASS"
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(path) + ".bak_F10M_20261013")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "CHECK.json").exists():
        rec = read(out / "CHECK.json")
        assert sha(path) == rec["after_sha256"] and sha(backup) == rec["before_sha256"]
        return
    before = path.read_bytes()
    assert sha(path) == "e825c9b81cad96241ea8326e120a15cc0bf28dc1245c62097f2fea87bb04d0cb", (
        "Document changed since inspection"
    )
    assert not backup.exists()
    text = "\n**F10M 把暂定薄管改成了实验中的PDMS通道。**F10L的20 mm薄壁管、9 mm圆片和均匀水场只保留为历史演示；现在按用户确认的32.5×9×9 mm长方体PDMS块建模，直径1.6 mm的针倒模通道位于横截面正中央并沿长度贯通，两侧PDMS厚度为3.7 mm，三枚SMD07T05R412WL陶瓷发声面改为7 mm、贴合左右及X正端，上下没有PZT。本轮又确认顶部有直径1.5 mm、向下连通主通道的注入口，注入后开放，而X正端始终封闭；5.4 mm暂按端部发声面到孔中心解释，孔轴x=10.85 mm，横向居中及孔底盲深仍需核对。材料图已把顶部孔与主通道并为液体区域，端片覆盖主通道开口的关系没有被透明压力源掩盖。多介质k-Wave接口支持空间声速、密度及吸收图，并拒绝混用不同吸收幂指数；参考Xu等2020年论文中10:1、85°C固化SYLGARD184的3–7 MHz测量数据，内外水介质仍为演示，未包含完整固体剪切波或真实PZT电压响应。六个小型三维传播对照及一个微型三圆片接口检查实际通过，吸收幅比误差约0.00347、PDMS相位误差约0.0122 rad；这些不是实际装置声场。完整3 MHz装置按最低波长约6.19格估计需64.8 GiB，超出当前4090预算，因此三路正式底稿和ON轨迹均标NOT_RUN，没有粗化网格硬算。更新后的主通道32泡1 ms OFF轨迹已在GPU完成，逐泡位置完全不变、流速和梯度严格零、接触及墙面修正为零；孔交汇区运动尚未准入，程序在进入该邻域前停止，避免把真实开口当成圆柱壁。四组英文PNG/PDF、修正前后记录、失败预检、源码与回传哈希已归档。下一步先核对孔口参考、X负端条件、PDMS品牌固化、内外介质和单频声压标定，准备匹配APECSS及可承受成本的真实装置场，再预登记远离孔口的指定目标对照，随后研究形态指标，不保证任意导航或任意形状。当前仍为 **EXPLORATORY / MATERIAL_ACOUSTICS_UNCALIBRATED / SOURCE_PRESSURE_UNCALIBRATED / PHYSICS_UNCALIBRATED，training_ready=false**；证据见独立F10M工作区的docs/FAST_PZT_F10M_DEVELOPMENT.md与results/straight_tube_exp1/f10m/。\n"
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    index = before.index(marker)
    added = text.encode()
    after = before[:index] + added + before[index:]

    def heading(d):
        return re.findall(rb"^#{1,6} .*$", d, re.M)

    assert heading(before) == heading(after) and after[:index] + after[index + len(added) :] == before
    backup.write_bytes(before)
    (out / "before.md").write_bytes(before)
    path.write_bytes(after)
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
            one_paragraph=True,
            all_original_bytes_preserved=True,
            headings_unchanged=True,
            before_sha256=sha(backup),
            after_sha256=sha(path),
            path=str(path),
            backup=str(backup),
            diff_sha256=sha(out / "change.diff"),
        ),
    )


if __name__ == "__main__":
    main()
