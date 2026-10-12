"""Insert the requested concise F10K update, preserving every original byte."""

from pathlib import Path
import difflib
import hashlib
import re
from .common import BASE, ROOT, js, sha, write


def main():
    assert js(BASE / "ENGINEERING_PREFLIGHT.json")["status"] == "PASS"
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(path) + ".bak_F10K_20261012")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    before = path.read_bytes()
    if (out / "CHECK.json").exists():
        record = js(out / "CHECK.json")
        assert sha(path) == record["file_sha256"][str(path)]
        assert sha(backup) == record["file_sha256"][str(backup)]
        print("Existing F10K insertion and backup verified; no rewrite.")
        return
    assert (
        hashlib.sha256(before).hexdigest()
        == "d01919e64ccbaeb1b8048b3ce764ec793d9f8c38d6000bbd15785aeeb211839b"
    ), "Progress document changed since inspection; reread before editing."
    assert not backup.exists(), "Do not overwrite an existing backup."
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    insertion = """

**F10K 把问题从求解器速度重新拉回科学适用范围。**F10J 的“GPU 比 CPU 慢”只针对 32 泡、单个固定接触问题，不代表所有大泡团。本轮保留 M0 的快速团级响应证据，也保留 B2 的有限位置与形状一致性；数值通过仍不等于真实接触、近场声力和壁面水动力已经校准。

**原 F5 的 417.5 ms 分裂候选没有被硬接到 B2。**独立原壁面查询显示泡面几乎贴壁，没有可分辨穿壁，但当前内部管腔 B2 不具备该近壁适用性，因此拒绝。有限检查已有状态后，选中 F10D prepare_1101 的 2.1 ms 真实远壁状态；没有平移或制造新泡团。新缓存 M0 重放原下一个 50 µs 保存间隔，坐标逐位一致，但新 B2 仍明确是跨模型初始化的探索性反事实，而非原历史的精确 B2 续算。

**终于有了一个持续图聚合的同模型步长对照，仍没有持续分裂正例。**新的 3 ms B2 窗口分别用 10/5 µs，31 泡团与 ID 5 都在局部 0.290 ms 连通，保持 2.710 ms；成员一致，最大固定 ID 差约 27.57 nm，Kabsch RMSD 12.15 nm，接触预算通过。但瞬时边为 46/47 条、配对时间差可达 60 µs，仍 FAIL；旧 12/29 条不达标也保留。图聚合不是泡融合，不能据此准入任意形态控制或 AI 科学数据。

**实际经过目标平面前，声学作用时间比 10 ms 短得多。**真实 F10A OFF/S/O 保存轨迹约在 2.024 ms 经过下游 20 µm 平面；O−OFF 额外 Z 位移仅约 20.32/20.53 nm，低于原 50 nm 作用门限，虽有正确方向仍未建立目标导航。声场在已观察的 10 ms 内仍有效，2.024 ms 是过面前的机会时间，不是场域总寿命。下一步预登记内部母支 X+20 µm、Z+1 µm 的固定目标，保留原历史目标，运行后不能迁就结果移动目标。

**有限规模测试支持继续研究大型 GPU B2，但还不能外推到真实千泡团。**N=32 CPU/GPU 约 13.47/24.17 ms；单一连通合成 N=128 为 151.31/17.33 ms，N=512 为 2028.32/97.68 ms，两端都通过原质量门限。后两者是规则格点，512 的对称性显著减少迭代；不是独立科学样本或 B64 并行轨迹。1024 因原 CPU 全泡对稠密存储估计约 24.5 GiB、超过本轮 4 GiB 数组预算而未运行。求解加预热共 7.44 s、完整规模阶段 10.04 s；N512 进程约 3.49 GiB，GPU 进程采样峰值 418 MiB。

**下一轮先处理真实科学缺口。**优先完成上述固定目标的控制能力与到达误差基线，再寻找合法远壁持续分裂正例、研究严格连边步长敏感；大型不规则团 GPU 与更节省内存的 CPU 参考随后有限推进。没有启动长轨迹、重求声场或 CFD，也没有扩大 AI 准入。两条短 B2 的审计输出成本分别约 7.57/12.46 s；含缺包启动失败和元数据恢复在内，六个服务器进程共约 41.25 s，失败和原始文件均已保留。

> F10K 原始报告、范围矩阵和规模报告：`docs/FAST_PZT_F10K_DEVELOPMENT.md`、`docs/FAST_PZT_F10K_SCIENTIFIC_SCOPE.md`、`docs/FAST_PZT_F10K_SCALING.md`；独立工作区 `computation_examples/_archives/mevo_fast_pzt_f10k_20261012/computation_examples/mevo_acoustic_microbubbles`。图件及分层验收位于 `results/fast_pzt/f10k/`。持续标记 **EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false**。
""".encode()
    index = before.index(marker)
    after = before[:index] + insertion + before[index:]

    def headings(b):
        return re.findall(rb"^#{1,6} .*$", b, re.M)

    assert headings(before) == headings(after)
    assert after[:index] + after[index + len(insertion) :] == before
    # Recheck immediately before mutation to avoid overwriting concurrent work.
    assert path.read_bytes() == before
    backup.write_bytes(before)
    (out / "before.md").write_bytes(before)
    path.write_bytes(after)
    (out / "after.md").write_bytes(after)
    diff = "".join(
        difflib.unified_diff(
            before.decode().splitlines(True),
            after.decode().splitlines(True),
            fromfile="before.md",
            tofile="after.md",
        )
    )
    (out / "change.diff").write_text(diff)
    files = [path, backup, out / "before.md", out / "after.md", out / "change.diff"]
    write(
        out / "CHECK.json",
        dict(
            status="PASS",
            original_bytes_preserved=True,
            original_heading_count=len(headings(before)),
            unchanged_headings=[h.decode() for h in headings(before)],
            insertion_bytes=len(insertion),
            original_bytes=len(before),
            after_bytes=len(after),
            file_sha256={str(p): sha(p) for p in files},
            baseline_review="Read original F10I/F10J narrative and preserved all earlier paragraphs and epilogue.",
            worktree=str(ROOT),
        ),
    )
    print("Progress update backed up; exact original bytes and heading structure preserved.")


if __name__ == "__main__":
    main()
