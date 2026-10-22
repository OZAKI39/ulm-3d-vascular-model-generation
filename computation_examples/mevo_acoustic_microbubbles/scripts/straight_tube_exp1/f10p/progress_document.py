"""Minimal byte-preserving Chinese addition after scientific and visual checks."""

from pathlib import Path
import difflib
import re
from .common import BASE, read, write, sha


def main():
    assert read(BASE / "analysis/ENGINEERING_VERIFICATION.json")["status"] == "PASS"
    assert read(BASE / "figures/VISUAL_REVIEW.json")["status"] == "PASS"
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = path.with_name(path.name + ".bak_F10P_20261014")
    before = path.read_bytes()
    assert sha(path) == "de8c62cf98e1f42df192fed007fe87f9c1fbe4cca9093ddd2a02bbec8f9fdd7f"
    assert not backup.exists(), "Never overwrite an existing stage backup"
    addition = """
**F10P 把硬球非穿透与团级、逐泡数值准入分开检查。**从完整F10O独立建立分支，核查活跃引用后确认本轮无新增安全删除对象，索引整理先单独提交；M0、B2、历史失败、场与响应表全部保留。新方法按固定球形、无摩擦、非粘附假设，在原成对润滑阻力的过阻尼受力平衡中加入单边速度约束，没有接触弹簧或事后PBD；当前法向局部QP的对称正定前提逐步检查，不宣称完整多球最终位置唯一。两泡/四泡六个独立参考及错误退出检查通过；首次OSQP因3000次迭代上限在第4步拒绝，保留失败后以同一QP的有界活动集直接解复核，未改变质量阈值。32泡、四枚7 mm PZT、2.5 MHz、均匀水缓存、F2/APECSS、R3封闭X壁面及严格零流原样复用。150步和0.25 ms强控制通过，已越过F10O在0.175 ms的17.8 nm软重叠失败点；硬球直接求解的最大实际重叠仅约8.47e−13 nm浮点舍入量、壁修正0。但强dt1向2 ms推进时在0.464 ms因真实球面反力—间隙互补等效量0.010057 nm超过本轮运行前登记的0.010000 nm而拒绝，合法前缀0.463 ms；原重叠/178 nm预算未放宽，局部平面QP残差很小也不能覆盖这一失败。减半步长的有限诊断通过0.5 ms，细缓存仍拒绝于0.464 ms；共同0.463 ms同ID步长差21.647 nm、缓存差0.05214 nm，COM差约4.24e−7/0.001764 nm，Rg步长差2.586 nm，仅属前缀证据。OFF和弱控制三组均完成2 ms：OFF严格静止，弱控制质心位移约1.619 nm，Rg减少0.146854 nm、X/Y比到1.00004925，形态效应超过所测数值差异而远未达到1.25目标。B2相同弱控制的Rg减少0.184224 nm、M0增加11.1069 nm，接触模型依赖仍明确存在。新分层规则分别验收基础质量、团位置、去质心形状/匹配分布、持续网络和固定ID；OFF/弱控制各自通过完整窗口相应数值检查，本轮没有真实的团级通过但ID失败硬球案例，代码独立性测试不冒充科研样本。初始32成员网络保持2 ms并非新聚合或融合；强短窗不足2 ms，不能验证持续事件。X/Y原0.2 µm目标及0.075 µm容差保留，因强完整窗口基础质量失败，新XY优化和5/10/30 ms均未运行；强dt1合法终点X/Y误差约0.14833/0.19751 µm，不能称目标达成。全部科学数值进程约44.856 s、峰值RSS约537 MiB，12条有限运行中的3次失败保留，77个远端文件逐项SHA回传，四组英文PNG/PDF和原始诊断已归档并实际目视检查。下一轮先解决有限步反力与真实球面几何一致性及快速内部重排，再做完整2 ms团级准入和XY导航，不增加力、放宽容差或外推长轨迹来制造成功。当前仍为 **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED，training_ready=false**；详见独立F10P工作区docs/FAST_PZT_F10P_DEVELOPMENT.md及results/straight_tube_exp1/f10p/。
""".encode()
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    after = before.replace(marker, addition + marker)
    assert after.replace(addition, b"", 1) == before

    def headings(v):
        return re.findall(rb"^#{1,6} .*$", v, re.M)

    assert headings(before) == headings(after)
    backup.write_bytes(before)
    out = BASE / "progress_document"
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
    path.write_bytes(after)
    write(
        out / "CHECK.json",
        dict(
            status="PASS",
            path=str(path),
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(path),
            heading_count_before=len(headings(before)),
            heading_count_after=len(headings(after)),
            all_original_bytes_preserved=True,
        ),
    )
    print(read(out / "CHECK.json"))


if __name__ == "__main__":
    main()
