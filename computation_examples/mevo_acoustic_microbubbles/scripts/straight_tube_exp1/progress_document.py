"""One Chinese paragraph, byte-preserving backup and structure audit."""

import difflib
import re
from pathlib import Path
from .common import BASE, ROOT, read, write, sha


def main():
    gate = read(BASE / "analysis/FINAL_GATE.json")
    assert gate["engineering_status"] == "PASS"
    path = Path("/home/lzy/projects/microbubble_aggregation.md")
    backup = Path(str(path) + ".bak_F10L_20261013")
    out = BASE / "progress_document"
    out.mkdir(exist_ok=True)
    if (out / "CHECK.json").exists():
        record = read(out / "CHECK.json")
        assert sha(path) == record["after_sha256"] and sha(backup) == record["before_sha256"]
        return
    before = path.read_bytes()
    assert sha(path) == "d66d3cdb88ba7d28d9058ab6c6ef47d2a933591b958a0c4424cafe39628835d5", (
        "Document changed after inspection; reread before mutation"
    )
    assert not backup.exists()
    dg = read(BASE / "dynamics/DYNAMICS_GATE.json")["cases"]
    ac = read(BASE / "acoustic/ACOUSTIC_GATE.json")
    on = dg["THREE_PZT_ON"]
    magnitude = sum(v * v for v in on["centroid_displacement_m"]) ** 0.5 * 1e9
    insertion = (
        "\n**F10L 先把电脑实验室对准正在搭建的无流直管装置。**此前的Y管主要研究输运、分叉和近壁数值问题；新的实验把微泡注入直管后不加主动血流，因此单独建立暂定内径1.6 mm的有限直管，保留三枚9 mm圆形PZT：左右各一枚、端部一枚，上下不加PZT。约3 mm暂按发声面到管外壁或端面计算；20 mm管长、0.1 mm壁厚及2.5 MHz均只是演示值，真实壁厚和材料仍记为未知。三路均匀水k-Wave复声场已经实际重算，管壁只约束泡泡位置，没有冒充真实材料反射；背景速度与梯度严格为零，也不读取旧Y管CFD。保留原M0完整逐泡模型、有限泡半径壁面约束和原接触预算，32个半径1.78 µm的演示泡完成OFF与三路同相开启的各1 ms短轨迹，步长1 µs；OFF逐泡位置完全不变，ON质心移动约"
        + f"{magnitude:.2f}"
        + " nm、连通分量由"
        + str(on["components_start_end"][0])
        + "变为"
        + str(on["components_start_end"][1])
        + "，未据此宣称成团或导航成功。每路输入300 Pa只是数值源幅值，实际泡处峰值声压约"
        + f"{on['pressure_peak_range_Pa'][0]:.1f}–{on['pressure_peak_range_Pa'][1]:.1f}"
        + " Pa；三路声源准备及求解约"
        + f"{ac['whole_acoustic_s']:.1f}"
        + " s，另保留首次依赖下载超时。四组英文PNG/PDF、实际轨迹、质量检查和回传哈希均已归档。仍需实验室确认内外径、管长壁厚、内外介质、有效发声面、间距基准、端口、频率电压相位、声压校准及泡径浓度；下一轮先据实预登记可分辨的指定目标与对照，再研究形态指标和控制能力，不保证任意三维导航或任意形状。当前继续标为 **EXPERIMENTAL_GEOMETRY_PENDING / MATERIAL_ACOUSTICS_UNCALIBRATED / SOURCE_PRESSURE_UNCALIBRATED / PHYSICS_UNCALIBRATED，training_ready=false**。证据见独立F10L工作区的docs/FAST_PZT_F10L_STRAIGHT_TUBE_DEVELOPMENT.md与results/straight_tube_exp1/f10l/。\n"
    )
    marker = "\n---\n\n### 尾声：故事里还没讲完的部分".encode()
    assert before.count(marker) == 1
    idx = before.index(marker)
    added = insertion.encode()
    after = before[:idx] + added + before[idx:]

    def headings(data):
        return re.findall(rb"^#{1,6} .*$", data, re.M)

    assert headings(before) == headings(after)
    assert after[:idx] + after[idx + len(added) :] == before
    assert path.read_bytes() == before
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
            one_paragraph_inserted=True,
            all_original_bytes_preserved=True,
            before_sha256=sha(backup),
            after_sha256=sha(path),
            backup=str(backup),
            path=str(path),
            worktree=str(ROOT),
            headings_unchanged=True,
            headings=[h.decode() for h in headings(before)],
            files={p.name: sha(p) for p in out.iterdir() if p.is_file()},
        ),
    )


if __name__ == "__main__":
    main()
