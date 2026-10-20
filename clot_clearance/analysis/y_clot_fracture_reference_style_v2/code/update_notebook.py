"""Append a verified visualization note, preserving existing scientific history."""
from pathlib import Path
import datetime
import difflib
import hashlib
import json
import os

R = Path(__file__).resolve().parents[3]
TASK, RUN = "y_clot_fracture_reference_style_v2", "20261020T094240Z"
O = R / "verification" / TASK / RUN
F = R / "results" / TASK / RUN


def main():
    verification = json.loads((O / "VALIDATION.json").read_text())
    review = json.loads((O / "VISUAL_REVIEW.json").read_text())
    assert verification["status"] == review["status"] == "PASS"
    path = Path("/home/lzy/projects/clot_rupture.md")
    before = path.read_bytes()
    original = before.decode()
    title = "## Y-Clot V4 Fixed-Camera Fracture Replay"
    assert title not in original
    prefix = "temp_storage/clot_github_sync_20261015/clot_clearance"
    section = f"""{title}

**按参考图样式重绘真实碎裂（{RUN}）**：用户明确要求碎裂过程而非绕Z轴展示。本轮只读取最新Y形血块V4已接受的0–3.6ms状态；相机位置、方向、缩放全程固定，世界Z轴竖直。白底六面网格、浅/深红强弱相材料点、蓝灰固定锚点、真实存活键；不使用连续外壳。所有节点位置和材料相区保持原数据，失效键不再绘制，图形标记尺寸不是节点积分体积。

提供原血块区域固定特写与全程全部节点完整视野两个版本。特写的世界坐标窗口不随碎块移动，移出窗口的材料点数明确标注；该视野裁切不改变全图连通分量统计。完整视野覆盖全部接受态节点，便于检查离开特写区域后的实际运动。两版均有2560×2048高清MP4及1600×1280演示版。视频按保存态慢放，右上角显示真实物理时间，没有坐标插值、位移放大或新增计算。当前边界和锚点仍规则，只有412内部点初始扰动；重绘不改变离散范围。

本轮输入SHA保持，视频完整解码、固定相机、Z竖直与全视野节点覆盖检查通过。未新增材料、声流、CFD或输运验证。见[视频及关键帧]({prefix}/results/{TASK}/{RUN}/OPEN_RESULTS.html)、[重放说明]({prefix}/analysis/{TASK}/REVIEW.md)及[验收记录]({prefix}/verification/{TASK}/{RUN}/VALIDATION.json)。

"""
    marker = "# ⑥ 外加剪切"
    position = original.index(marker)
    after = (original[:position] + section + original[position:]).encode()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    audit = O / "notebook_update" / stamp
    audit.mkdir(parents=True, exist_ok=False)
    (audit / "before.txt").write_bytes(before)
    (audit / "after.txt").write_bytes(after)
    (audit / "update.diff").write_text("".join(difflib.unified_diff(
        original.splitlines(True), after.decode().splitlines(True),
        fromfile=str(path) + " before", tofile=str(path) + " after")))
    before_sha = hashlib.sha256(before).hexdigest()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before_sha, "Concurrent document edit"
    temporary = path.with_name(path.name + ".fracture_replay.tmp")
    temporary.write_bytes(after)
    os.replace(temporary, path)
    metadata = dict(status="UPDATED", before_sha256=before_sha,
                    after_sha256=hashlib.sha256(after).hexdigest(), backup=str(audit))
    (audit / "HASHES.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (O / "NOTEBOOK_UPDATE.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("NOTEBOOK_UPDATED", audit)


if __name__ == "__main__":
    main()
