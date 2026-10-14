"""One paragraph only; preserve all original document bytes and headings."""

from pathlib import Path
import re
import difflib
from .common import BASE, read, write, sha


def main():
    assert read(BASE / "analysis/ENGINEERING_VERIFICATION.json")["status"] == "PASS"
    assert read(BASE / "figures/VISUAL_REVIEW.json")["status"] == "PASS"
    p = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(p) + ".bak_F10O_20261013")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "CHECK.json").exists():
        r = read(out / "CHECK.json")
        assert sha(p) == r["after_sha256"] and sha(backup) == r["before_sha256"]
        return
    assert sha(p) == "a4a8db73d25242ddafb1f970da84db6f01fe0cc8b265ef25feaba473957f6f9f"
    assert not backup.exists()
    before = p.read_bytes()
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    text = """
**F10O 先清理真正无效的代码，再复用已有B2检查四PZT强控制。**从完整F10N独立建立工作区，实际删除一个没有调用的Warp类型转换函数（5行），引用、动态入口及原F10N只读验证通过后单独提交；CPU B2、R3暴露壁面、四路水声场、历史Y管/三PZT复现依赖和所有科学数据均保留，并建立精简活跃入口索引。没有另写接触求解器：直接调用F10H/F10K原几何单边软接触及成对润滑阻力，刚度0.04641319518378658 N/m、数学残差1e-13 m及接触/墙预算不变；声场、F2近场、APECSS、固定32个ID、2.5 MHz、严格零流和封闭X端与F10N一致。14项实际CPU/GPU输入和拒绝检查通过，每组先测150步，OFF/弱/强约3.34/3.64/5.84 ms每步，再开展有界短窗对照。强控制在dt=1及0.5 µs、3及1.5 µm缓存下都于0.175 ms提案超过17.8 nm重叠预算，数学方程已收敛而软接触几何不准入；合法前缀仅0.174–0.1745 ms，失败提案保存，未放宽预算或改变刚度。共同0.174 ms前缀同ID步长差约0.0244 nm、缓存差0.000744 nm；在同一已保存0.1 ms时刻，原M0与B2步长差分别约1.542及0.0201 nm，说明短前缀改善，仍不能解决F10N完整2 ms的6.416 µm问题。OFF和弱控制完成2 ms，弱控制同ID步长/缓存差约0.004545/0.000178 nm；B2的Rg减少0.1842 nm、X/Y比仅增0.005779%，而原M0为Rg增加11.1069 nm、X/Y增0.4245%。两者接触及迁移率近似不同，不能将跨模型差异叫作纯数值误差；小幅B2变化超过所测差异仍未达到1.25形态目标，连通网络也不代表融合。新X=(0.2,0,0)及Y=(0,0.2,0) µm、容差0.075 µm在试跑前登记，旧(3,2,1) µm目标不动；因强控制质量门槛失败，新XY优化、分段反馈及5/10/30 ms均未运行，不能据短前缀+X移动声称独立XY导航或外推长时间成功。两个数值进程约72.91 s，峰值RSS约529 MiB；12条短计算中9条完成各自窗口、3条强控制失败，62个远端文件逐项SHA回传，四组英文PNG/PDF和失败图证据归档。下一步应在较低驱动或相消控制的可准入范围内先完成2 ms逐泡复核，再研究独立X/Y短目标与小段反馈，形态收益必须改用B2实际张量评估；不为获得PASS重写求解器或调整刚度。均匀水声场、近场响应、接触及墙水动力仍未实验校准，保持 **EXPLORATORY / PHYSICS_UNCALIBRATED，training_ready=false**；详见独立F10O工作区docs/FAST_PZT_F10O_DEVELOPMENT.md和results/straight_tube_exp1/f10o/。
"""
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
