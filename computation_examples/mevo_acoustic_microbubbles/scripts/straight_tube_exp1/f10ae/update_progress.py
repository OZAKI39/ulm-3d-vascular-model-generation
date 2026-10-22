"""Minimal local progress update, with exact before/after and structure evidence."""
import difflib
import re
from pathlib import Path
from .common import BASE, ROOT, sha, write, identity

def main():
    path = Path('/home/lzy/projects/microbubble_aggregation.md')
    backup = path.with_name(path.name + '.bak_F10AE_20261019')
    assert not backup.exists(), 'Existing backup is immutable; inspect before resuming.'
    before = path.read_text()
    assert '### 2.7 ' not in before
    backup.write_bytes(path.read_bytes())
    out = BASE/'evidence/progress_doc'
    out.mkdir(parents=True, exist_ok=False)
    (out/'before.md').write_bytes(path.read_bytes())
    after = before.replace('更新：2026-10-18。', '更新：2026-10-19。', 1)
    after = after.replace('[F10AC + F10AD 工程][sim-root]；[多频率与形变探索报告][acad-report]', '[F10AE 工程][sim-root]；[指定形态与恢复报告][ae-report]', 1)
    oldroot = '/home/lzy/projects/computation_examples/_archives/mevo_qd_f10ac_ad_20261018/computation_examples/mevo_acoustic_microbubbles'
    after = after.replace(oldroot, str(ROOT), 1)
    after = after.replace('当前分支为 `dev/f10ac-qd-f10ad-multifrequency-20261018`，基于已核查的 F10AB `9d304c22817599f85079ada5a42b2e143f3fe055`', '当前分支为 `dev/f10ae-targeted-shape-recovery-20261019`，基于已核查的 F10AC/F10AD `baa47dbe95c78e1800aebe50bc3efc5faf07d5fa`', 1)
    after = re.sub(r'^\[sim-root\]: .*$', '[sim-root]: '+str(ROOT), after, flags=re.M)
    after = after.replace('均匀水圆片声源代理，2.5 MHz；复用四路复声压缓存', '均匀水圆片声源代理，2.5 MHz 主研究；2.4/2.6 MHz 对照各用匹配缓存与响应', 1)
    after = after.replace('主要为 1–2 ms；32 泡部分场景已验证至 5 ms', '主要为 1–2 ms；本轮 64 泡两段动作共 4 ms；32 泡部分场景曾验证至 5 ms', 1)
    section = '''### 2.7 指定形态、恢复与鲁棒性（F10AE）

继续复用 F10AA 逐泡引擎，限定 64 泡、无流、均匀水代理，先登记目标再搜索幅相。官方 APECSS 支持将历史 2.5 MHz、原壳参数和半径 1.78 µm 的低压范围有限扩展至 **0.001 Pa**；0.0001 Pa 的精细求解超时，未准入。历史表与轨迹不变，搜索比较全程保留原压力域，新的低压扩展单独验证。

Ax/BoTorch 与 Sobol 各有 32 次观察，包含 8 个共同起点；本次 Ax 首次达到发现集目标需 10 次观察，Sobol 为 26 次。这是单次有限预算比较。选定 T1 动作使发现初态长宽比增加约 **2.997%**、质心移动 **2.132 nm**，但独立验证及最终测试均只命中 **2/3**；较小形变目标 T2 的独立验证为 **2/3**、最终测试为 **3/3**。两者都未达到预登记的跨初态完整确认条件。

不能将长宽比增加直接称为拉伸。两个验证初态出现超过数值差异的全局主轴增长 **1.186 / 2.218 nm**，但已分成 3 / 2 个邻近分量，其最大 62 泡子群仍收缩。另一个最终测试初态出现连通态约 **7.808 nm** 主轴增长，未对该留出初态追加数值复核，保留为未确认观察。第二段主动控制也没有实现恢复：选定案例的 Rg 误差由约 **0.694%** 增至 **0.746%**，固定 ID 内部结构误差同样增加，不能称为可逆形状控制。

本轮 **230 次逐泡尝试，227 条完整、3 条中止**；630,000 个完整轨迹时间步完成独立审计，66 组半步长/缓存比较通过。核心科学计算与审计约 **18.24 分钟**，图和视频生成另计约 1 分钟，不含安装、开发与回传。12 组 PNG/PDF、3 段真实逐泡 MP4、失败记录和输入哈希均已保存。频率和壳参数变化会使目标值偏离；目前仍是短窗口、未校准的均匀水代理，`training_ready=false`。

![实际主轴变化与初态敏感性][ae-figure]

[25 项中文结论与完整限制][ae-report] · [两段动作与未恢复的真实视频][ae-video] · [开发和复现入口][ae-code]。服务器独立目录为 `/workspace/mevo_acoustic_microbubbles/f10ae_target_20261019`，大型逐步原始数据保留在该目录。下一步应在多个发现初态上优化，再用独立测试判断能否可靠控制；恢复研究应先寻找可重复的反向形变机制，并补做切换瞬态和真实四源声场标定。

'''
    after = after.replace('## 3. 人工标注', section+'## 3. 人工标注', 1)
    result = str(BASE)
    after += '\n[ae-report]: '+result+'/FINAL_CHINESE_REPORT.md\n[ae-figure]: '+result+'/figures/Fig_F10AE_06_absolute_principal_axis_extension.png\n[ae-video]: '+result+'/videos/Video_F10AE_2_deformation_and_recovery.mp4\n[ae-code]: '+result+'/EXECUTION_AND_REPRODUCTION.md\n'
    old_titles = re.findall(r'^#{1,6} .+$', before, re.M)
    new_titles = re.findall(r'^#{1,6} .+$', after, re.M)
    assert [x for x in new_titles if not x.startswith('### 2.7 ')] == old_titles
    old_vision = before.split('## 3. 人工标注', 1)[1].split('[sim-root]:', 1)[0]
    new_vision = after.split('## 3. 人工标注', 1)[1].split('[sim-root]:', 1)[0]
    assert old_vision == new_vision
    path.write_text(after)
    (out/'after.md').write_bytes(path.read_bytes())
    (out/'changes.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True), fromfile='before.md', tofile='after.md')))
    write(out/'CHECK.json', dict(original=str(path), backup=str(backup), before_sha256=sha(backup), after_sha256=sha(path), old_titles=old_titles, new_titles=new_titles, title_structure_pass=True, vision_content_unchanged=True, historical_sections_preserved=True, identity=identity()))
    print('PROGRESS_DOCUMENT_UPDATED', path, flush=True)

if __name__ == '__main__':
    main()
