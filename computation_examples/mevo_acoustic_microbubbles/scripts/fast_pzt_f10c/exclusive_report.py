"""Update performance evidence and the existing narrative after the bounded replay."""

from datetime import datetime
from pathlib import Path
import difflib
import re
import shutil
from .common import BASE, ROOT, js, sha, write


def diff(before, after):
    return "".join(
        difflib.unified_diff(
            before.splitlines(True), after.splitlines(True), fromfile="before.md", tofile="after.md"
        )
    )


def main():
    evidence = js(BASE / "analysis/RESULTS.json")["exclusive_performance"]
    assert evidence["completion"]["status"] == evidence["resource_audit"]["status"] == "PASS"
    rows, cold = evidence["batch"], evidence["cold"]
    assert evidence["completion"]["engineering_2x_pass"]
    b16, b64, c64 = rows[2], rows[3], cold[3]
    audit = evidence["resource_audit"]
    samples = [jsline for jsline in (BASE / "exclusive_v1/RESOURCE_SAMPLES.jsonl").read_text().splitlines()]
    import json

    times = [datetime.fromisoformat(json.loads(line)["time_UTC"]) for line in samples]
    gap = max((b - a).total_seconds() for a, b in zip(times, times[1:]))
    commit = js(BASE / "exclusive_v1/IDENTITY.json")["execution_commit"]
    summary = f"""## 独占复测已完成（2026-10-11）

用户选择等待另一项求解结束。Y-clot全部退出后执行 `exclusive_v1`，科学源码 `{commit}`，未暂停或终止其他项目。预先登记 `exclusive_performance.json`，原物理协议、时间步长、完整逐泡输出和数值预算不变。

全程资源记录覆盖 {audit["sample_count"]} 次采样，采样间最大间隔 {gap:.3f} s；起止CUDA进程为空，运行期间最多1个CUDA计算进程，没有其他MEVO/Y-clot科学Supervisor任务。证据为 `exclusive_v1/RESOURCE_AUDIT.json` 和 `RESOURCE_SAMPLES.jsonl`。这是有采样间隔的独占证据，不能排除短于间隔的外部活动；监测开销未从计时扣除。

同口径串行与批量的85/85成员逐位相同，新冷进程85/85成员也逐位相同。B16与B64均超过原≥2倍工程门槛，标记 `THROUGHPUT_IMPROVED`。各规模只做一组有限配对，没有持续一小时试验或性能置信区间；每小时数值为本组墙钟外推，限32泡、10 ms、201完整帧、同控制家族。重复/平移初态仍为 `BENCHMARK_ONLY`，不是独立AI训练样本。

| B | 串行总时长/s | 批量总时长/s | 等效每条/s | 实测加速比 | 合格数 | 合格轨迹/GPU小时 |
|---:|---:|---:|---:|---:|---:|---:|
"""
    for r in rows:
        summary += f"|{r['B']}|{r['serial_total_s']:.5f}|{r['batch_total_s']:.5f}|{r['batch_total_s'] / r['B']:.5f}|{r['resident_speedup']:.3f}|{r['B']}|{r['valid_samples_per_gpu_hour']:.0f}|\n"
    summary += """
静态资源驻留，计时包含模型/图构建、逐泡积分和质量、读回、形态统计及完整文件输出。全批比较和测试驱动额外开销记录在 `exclusive_v1/COMPLETE.json`，不当作生产轨迹开销。单条B1只有约5%的单次观测差，不能据此宣称稳定单条加速。

| B | 新进程完整时长/s | 合格轨迹/GPU小时 | Warp高水位/MiB | 整卡采样峰值/MiB |
|---:|---:|---:|---:|---:|
"""
    for r in cold:
        summary += f"|{r['B']}|{r['full_subprocess_wall_s']:.5f}|{r['valid_samples_per_gpu_hour']:.0f}|{r['warp_allocator_high_water_MiB']:.3f}|{r['peak_observed_gpu_MiB']:.3f}|\n"
    t = c64["timing"]
    summary += f"""
冷启动是实际新Python进程和CUDA上下文；保留磁盘编译缓存。整卡NVML值包含本项目context/驱动，不是Warp分配；10 ms显存采样也不是无遗漏瞬时峰值。B64冷进程积分/质量 {t["integration_quality_wall_s"]:.3f} s、读回 {t["sync_readback_wall_s"]:.3f} s、形态/诊断准备 {t["morphology_quality_output_prepare_s"]:.3f} s、HDF5 {t["hdf5_output_s"]:.3f} s。完整输出的CPU分析仍是下一步优先优化项。

本次补测只关闭性能证据缺口；多速率形态拒绝、原逐泡物理未校准和 `training_ready=false` 均不改变。下面完整保留首次共享负载测量，不能与本次独占值混合。

"""
    perf = ROOT / "docs/FAST_PZT_F10C_PERFORMANCE.md"
    old = perf.read_text()
    assert "## 独占复测已完成" not in old
    old = old.replace("## 结论与资源限制", "## 首次测量：共享负载与资源限制", 1)
    old = old.replace(
        "**本轮不能批准独占 GPU 性能门槛。**",
        "**首次测量不能批准独占 GPU 性能门槛；上面的后续独占复测另行验收。**",
        1,
    )
    old = old.replace("未来复测应记录全程占用", "后续独占复测已补全全程采样记录", 1)
    old = old.replace(
        "下一轮性能任务必须先取得不重叠的 GPU 窗口，再少量同口径串行/批量配对，区分首次资源准备、驻留边际、异步读回与完整输出。",
        "本次独占补测已取得不重叠的GPU窗口并完成有限同口径配对；下一轮优先优化CPU形态统计、读回和完整输出。",
        1,
    )
    heading, rest = old.split("\n", 1)
    perf.write_text(heading + "\n\n" + summary + rest.lstrip())

    dev = ROOT / "docs/FAST_PZT_F10C_DEVELOPMENT.md"
    old = dev.read_text()
    old = old.replace(
        "独占性能验收尚未完成。",
        f"等待其他项目结束后，独占复测已完成：B16/B64实测{b16['resident_speedup']:.2f}/{b64['resident_speedup']:.2f}倍，通过≥2倍工程门槛。",
        1,
    )
    old = old.replace(
        "不再用其批准“独占 GPU ≥2 倍”门槛。",
        "未用其批准独占门槛；后续准入依据独立的 `exclusive_v1` 全程资源记录与同口径复测。",
        1,
    )
    old = old.replace(
        "是共享 GPU 上实际完成数/实际墙钟时间，不是独占吞吐保证。",
        "早期记录是共享负载观测；`exclusive_v1_batch/cold` 为独占复测，区别见清单的 performance_scope。两者均不构成稳定产量保证。",
        1,
    )
    old = old.replace(
        "实际 B、共享负载串行/批量/冷进程、有效数、进程分配/整卡显存及限制。",
        "实际B、独占串行/批量/冷进程、有效数、进程分配/整卡显存及采样范围。旧共享负载图保留在归档提交f49ecba。",
        1,
    )
    old = old.replace(
        "先在没有其他 GPU 作业时做有限配对性能复测，再批准独占吞吐。共享负载下 B64 的形态分析/输出已比积分更贵，",
        "独占复测和≥2倍工程验收已完成；独占B64的形态分析/输出仍比积分更贵，",
        1,
    )
    old += f"\n## 独占性能补充证据\n\n`exclusive_v1`及其batch/cold子目录保留85条真实串行、四个完整批量和四个冷进程，以及{audit['sample_count']}次资源监测。数值比较全部逐位一致，科学提交`{commit}`。B64驻留批量{b64['batch_total_s']:.5f} s，折算{b64['valid_samples_per_gpu_hour']:.0f}条/GPU小时；冷进程{c64['full_subprocess_wall_s']:.5f} s，折算{c64['valid_samples_per_gpu_hour']:.0f}条/GPU小时。限定范围和全部规模见性能报告。新进度文档另存`.bak_F10C_exclusive_20261011`，原备份未覆盖。\n"
    dev.write_text(old)

    readme = ROOT / "scripts/fast_pzt_f10c/README.md"
    old = readme.read_text()
    old = old.replace(
        "All performance rows overlapped another project after the initial idle check.\nRead `docs/FAST_PZT_F10C_PERFORMANCE.md` before using them. Exclusive performance is\n**NOT_VERIFIED**, despite bitwise numerical batch PASS. Do not interrupt unrelated jobs.",
        "The original batch/cold performance rows overlapped another project. They remain\nshared-load observations. `exclusive_v1` waited for that project to exit, then replayed\nserial/batch and fresh-process measurements with continuous sampled resource evidence.\nIts B16/B64 results pass the 2x engineering gate; see the performance report for\nmeasured values and scope. Do not interrupt unrelated jobs.",
    )
    old += "\nExclusive replay (new label, idle device required):\n\n```bash\npython scripts/fast_pzt_f10c/remote.py exclusive --label exclusive_new --timeout 240\n```\n\nThe driver stops its own child and rejects timing if sampled contention is found.\nIt runs at most one CUDA process at a time; the external monitor observes compute\nPID cardinality and science Supervisor tasks, not pure GPU kernel duration.\n"
    readme.write_text(old)

    source = Path("/home/lzy/projects/microbubble_aggregation.md")
    out = BASE / "progress_document"
    verification = js(out / "VERIFICATION.json")
    assert sha(source) == verification["after_sha256"], "Narrative changed since the prior audit"
    before = source.read_text()
    backup = source.with_name(source.name + ".bak_F10C_exclusive_20261011")
    assert not backup.exists()
    new_paragraph = f"**批量的独占收益已补测。**等另一项目结束后，保持完整逐泡输出，16条逐条计算约{b16['serial_total_s']:.2f} s，一次批量约{b16['batch_total_s']:.2f} s，加速{b16['resident_speedup']:.2f}倍；64条分别约{b64['serial_total_s']:.2f} s和{b64['batch_total_s']:.2f} s，加速{b64['resident_speedup']:.2f}倍。后者折算约{b64['valid_samples_per_gpu_hour']:.0f}条/GPU小时，冷进程约{c64['full_subprocess_wall_s']:.2f} s、{c64['valid_samples_per_gpu_hour']:.0f}条/GPU小时。85条串行/批量和85条冷启动对照全部逐位一致，GPU进程和其他科学任务持续采样未见冲突。每个规模只有一组测量，每条仍是10 ms的32泡基准，不代表长期产量或独立训练数据。此前另一项目并行时的观测另行保留，不混作本次独占结论，也没有中断别人的求解。"
    old_paragraph = next(
        p for p in before.split("\n\n") if p.startswith("**批量有可见的收益，但独占速度仍要补测。**")
    )
    after = before.replace(old_paragraph, new_paragraph, 1)
    after = after.replace(
        "下一轮先取得空闲GPU窗口做有限独占配对复测，再优化批量形态统计和读回；",
        "独占有限配对复测已通过，下一轮优先优化批量形态统计和读回；",
        1,
    )
    headings = re.findall(r"^#{1,6} .*$", before, re.M)
    assert headings == re.findall(r"^#{1,6} .*$", after, re.M)
    shutil.copy2(source, backup)
    extra = out / "exclusive_update"
    extra.mkdir(exist_ok=False)
    (extra / "before.md").write_text(before)
    (extra / "after.md").write_text(after)
    (extra / "change.diff").write_text(diff(before, after))
    source.write_text(after)
    (out / "after.md").write_text(after)
    original = (out / "before.md").read_text()
    (out / "change.diff").write_text(diff(original, after))
    write(
        extra / "VERIFICATION.json",
        dict(
            status="PASS",
            backup=str(backup),
            before_sha256=sha(backup),
            after_sha256=sha(source),
            headings_before=headings,
            headings_after=headings,
        ),
    )
    verification.update(
        after_sha256=sha(source),
        scope="Original F10C insertion plus two minimal performance/next-step paragraph updates after exclusive replay; original backup preserved",
    )
    write(out / "VERIFICATION.json", verification)


if __name__ == "__main__":
    main()
