"""Minimal addition to the existing Chinese research narrative, with byte audit."""

from pathlib import Path
import difflib
import re
import shutil
from .common import BASE, write, sha

ADDITION = """**F10C 回到每一颗泡泡：快一些，也要保留内部运动。**G1只能预测团中心，不能产生真实的重排和分散。这次继续让32个固定ID泡泡各自受力、移动、接触，沿用原声场、响应表、CFD和M0算法，没有把团锁成刚体。整体位置加相对位置、平均速度加相对速度，只是每步的诊断分解；完整10 ms轨迹与原逐泡路径逐位相同，分解本身没有带来加速。

**一张显卡可以同时算几群彼此无关的泡泡。**这次沿用已有的批量坐标结构和CUDA Graph，让同一控制家族共用声场及控制时刻，各群保留自己的32泡坐标、力、接触和连边。1、4、16、64条的批量结果，都与相同输入逐条计算逐位一致。不同控制日程先分组再执行，还不是同一批中完全异步的控制器。某个成员触发设备错误会停止整批，其他成员不会被悄悄标成成功。

**慢在哪里，需要实际量。**有界Warp活动计时中，十二轮接触约占47%，局部润滑约占25%，声场和CFD查询占比较小。这是插入CUDA Event后的活动跨度，不是纯kernel时间。旧CUDA Graph确实已开启；短对照有图约59.63 ms、无图约60.21 ms，结果完全相同，没有把已有功能算作新发明。批量变大以后，CPU形态统计和文件输出也占了明显成本。

**批量有可见的收益，但独占速度仍要补测。**保持完整逐泡输出，16条逐条计算约15.58 s，一次批量约2.91 s；64条分别约59.25 s和7.62 s。换算为实际经过时间，后者约30238条/小时。冷进程64条约11.62 s，约19821条/小时。后续显存审计发现，另一项目在测试期间启动并占用同一GPU，因此这些只是共享负载下的观测，不能批准独占显卡的产量或“稳定7.78倍”结论。没有干预别人的求解；独占性能验收保留为待验证。

**少查几次外部场，团中心仍准，内部排列却会变。**尝试每两步更新声音或血流时，距离、二次力、润滑、接触、管壁和最终CFD合法性仍每10 µs检查，控制变化和单元/连边变化会刷新缓存。仅缓存声音的团级最大质心差约2.37 nm，但原相位组逐泡最大差约3.31 µm、对齐RMSD约1.52 µm，没过事先登记的形态门槛。仅缓存血流的原相位组在9.16 ms出现17.865 nm提案重叠，超过原17.8 nm预算，被拒绝；接受前缀和首次失败提案都保留。没有放宽阈值，也没有继续做每四步更新或把缓存升为默认。

**当前能更方便地产生逐泡数值样本，还不能称AI真值。**批量保留每泡身份、完整坐标、团大小和张量、质量、实际控制、源码与输入哈希；本轮场景没有观察到连通团分裂或合并。平移和复制初态只是性能检查，同初态、不同控制和数值模式仍归在同一个训练/测试家族，不能用重复样本制造数据独立性。当前继续是 **EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**。

> **F10C证据与下一步：**`docs/FAST_PZT_F10C_DEVELOPMENT.md`、`docs/FAST_PZT_F10C_PERFORMANCE.md`与`results/fast_pzt/f10c/`保留零容差批量对照、缓存拒绝、真实冷进程、共享资源审计和四组英文PNG/PDF。下一轮先取得空闲GPU窗口做有限独占配对复测，再优化批量形态统计和读回；逐泡接触、集体迁移率、独立初态与实验校准仍要分别验证。不能只凭质心准确或动画看起来聚成团，就把内部排列当成真实物理。

"""


def main():
    source = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = source.with_name(source.name + ".bak_F10C_20261011")
    assert not backup.exists(), "Existing backup must not be overwritten"
    before = source.read_text()
    assert "**F10B 尝试" in before
    marker = "---\n\n### 尾声：故事里还没讲完的部分"
    assert before.count(marker) == 1
    shutil.copy2(source, backup)
    after = before.replace(marker, ADDITION + marker)
    assert after.replace(ADDITION, "", 1) == before
    headings = lambda text: re.findall(r"^#{1,6} .*$", text, re.M)
    assert headings(before) == headings(after)
    source.write_text(after)
    out = BASE / "progress_document"
    out.mkdir(exist_ok=False)
    (out / "before.md").write_text(before)
    (out / "after.md").write_text(after)
    (out / "change.diff").write_text(
        "".join(
            difflib.unified_diff(
                before.splitlines(True), after.splitlines(True), fromfile="before.md", tofile="after.md"
            )
        )
    )
    write(
        out / "VERIFICATION.json",
        dict(
            status="PASS",
            original_path=str(source),
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(source),
            headings_before=headings(before),
            headings_after=headings(after),
            original_bytes_preserved_after_removing_insertion=True,
            scope="Minimal evidence-supported insertion; exclusive performance explicitly pending.",
        ),
    )


if __name__ == "__main__":
    main()
