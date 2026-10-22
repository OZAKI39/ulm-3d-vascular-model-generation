"""Generate factual delivery documents from final audited numerical artifacts."""
import sys,json,time,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
def read(p):return json.loads((ROOT/p).read_text())
def write(name,s):(ROOT/name).write_text(s.strip()+'\n',encoding='utf-8')
def table(headers,rows):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n'
def main():
 a=read('results/delivery_analysis.json');perf=a['performance'];off=next(r for r in perf if r['name']=='verified_offline_regular');on=next(r for r in perf if r['name']=='verified_replay_regular');ons=[r for r in perf if 'replay' in r['name']];final=next(r for r in a['accuracy'] if r['split']=='final_test' and r['backend']=='B4' and r['mode']=='causal' and not r['multipass']);base=next(r for r in a['accuracy'] if r['split']=='final_test' and r['backend']=='B4' and r['mode']=='local');multi=next(r for r in a['accuracy'] if r['split']=='final_test' and r['backend']=='B4' and r['multipass']);diag=read('results/diagnostics_v1/post_freeze_diagnostics.json');bt=a['backend'];b=lambda n,reg:next(r for r in bt if r['backend']==n and r['region']==reg);bench=table(['Backend','ROI p50/p95 ms','Medium p50/p95 ms','Full bbox p50/p95 ms'],[[n]+[f"{b(n,reg)['timing']['p50']:.2f} / {b(n,reg)['timing']['p95']:.2f}" for reg in ['historical','medium','full_visible_bbox']] for n in ['B0','B1','B2','B3','B4']]);pt=table(['Run','Frames','p50 / p95 / max ms','Mean valid windows','Fresh group-frame Hz','Drops'],[[r['name'],r['frames'],f"{r['p50']:.2f} / {r['p95']:.2f} / {r['max']:.2f}",f"{r['mean_valid_windows']:.1f}",f"{r['fresh_group_frame_hz']:.2f}" if r['fresh_group_frame_hz'] is not None else 'offline',r['dropped']] for r in perf]);scope='EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false。输出是二维可见微泡图案的表观运动；不是水流真值、物理融合、纯声学响应或控制准入。'
 write('README.md',f'''# F10X-PIV — 全通道微泡运动观测

已完成现成 CPU/GPU 后端部署、937 帧全通道离线分析、实际 PTS 因果回放、9 组科研图和 3 段视频。结论为**工程原型完成、科学准入部分通过**：空间范围和同条件后端速度提高，但尚未可靠达到 30 fps，仍有低对比度漏检及非仿射错误接受。

{scope}

- [中文开发报告：33 项回答](F10X_DEVELOPMENT_REPORT.md)
- [最终推荐](FINAL_RECOMMENDATION.md)；[机器可读 Gate](FINAL_GATE.json)
- [科研图](figures/figure_manifest.json)；[完整视频](videos/full_video_manifest.json)；[后端比较视频](videos/backend_video_manifest.json)
- [活动代码索引](ACTIVE_CODE_INDEX.md)；[源码复用与许可证](SOURCE_REUSE_AUDIT.md)
- [恢复清单](RECOVERY_MANIFEST.json)；[输入输出哈希](INPUT_OUTPUT_SHA256.json)

正式离线结果：`results/full_video_verified_v5/results/offline_regular/`；正式定速结果：`results/replay_final_v6/results/`。`full_video_fast_v2` 和 `full_video_verified_v3` 是未通过全量数值等价检查的加速试验，不得作为正式结果。科研冻结配置、独立程序种子及最终测试在 `results/research_v1/`。

本地：`{ROOT}`。服务器：`/workspace/f10x_piv_20261017`。本轮分支 `dev/f10x-full-channel-piv-20261017`，不推送 GitHub。

本地检查：`.venv/bin/python -m pytest -q tests`。依次重建报告：`scripts/analyze_delivery.py`、`scripts/publication_figures.py`、`scripts/render_videos.py`、`scripts/render_backend_video.py`、`scripts/write_reports.py`。绘图/编码只读保存的正式计算结果，不连接硬件。

服务器环境 `/workspace/f10x_piv_20261017/.venv`。用 `scripts/remote_job.py` 建立当前 Git 的不可变快照和唯一 Supervisor 标签；新计算必须使用新的输出目录。数据环境变量参见 `src/xcommon.py`，不能覆盖既有结果。
''')
 write('ACTIVE_CODE_INDEX.md','''# 活动代码索引

| 文件 | 职责 |
|---|---|
| src/backends.py | B0 旧实现、B1 SciPy、B2 Rust、B3 TorchPIV、B4 CuPy；现成内核薄适配 |
| src/channel.py | 可编辑几何、全窗口覆盖、每帧全通道候选/纹理/活动及配准 |
| src/observer.py | 因果背景、调度、窗口 null/年龄、局部候选群六参数 |
| src/fit_adapter.py | 保留 F10U 拟合主体，仅提供每次复制的解析雅可比 |
| src/channel_state.py | 新接口检查、TTL、显式局部群 v2/v1 兼容 |
| src/studies.py | 官方多遍调用、已知位移与六参数误差 |
| scripts/setup_remote.sh, build_rust.sh | 独立环境、官方 Rust 编译 |
| scripts/cupy_mask_compat.py | int64/int32 窗口压缩索引兼容修复与 diff |
| scripts/backend_benchmark.py | 相同帧对三种尺度的同步性能比较 |
| scripts/frozen_baseline.py | 原 F10W 代码及原环境的额外基准 |
| scripts/generate_data.py, accuracy_study.py, freeze_and_final.py | 新程序种子、开发/验证/一次性 final |
| scripts/audit_channel.py | 只读已审核人工区域、全帧扫描及位移拒判检查 |
| scripts/full_offline.py, causal_replay.py | 全 937 帧离线和独立 PTS 定速回放 |
| scripts/diagnostics_and_clip.py | 60 帧五后端对照、筛选/墙边/突变/传输诊断 |
| scripts/analyze_delivery.py | 正式结果汇总、逐帧数值等价、面积与不确定性 |
| scripts/publication_figures.py, render_videos.py, render_backend_video.py | 保存数据驱动的 PNG/PDF/MP4 |
| tests/test_contracts.py | 数学、几何、null、TTL、兼容及鲁棒雅可比回归 |

历史 F10U 数学、F10V PyAV 编码、F10W 旧环境以只读方式复用。third_party 是被 Git 忽略的官方检出；许可证与 commit 在 evidence 中。
''')
 write('SOURCE_REUSE_AUDIT.md','''# 源码复用审计

| 官方项目 | 实际 Git commit | 许可 |
|---|---|---|
| [OpenPIV](https://github.com/OpenPIV/openpiv-python) | 91e88a8ee9289d6c4ca79722d9759c5258abd31d | GPL-3.0 |
| [TorchPIV](https://github.com/NikNazarov/TorchPIV) | 17ae4c1a4d4950efecfdfd9c404d77036a6f8c80 | MIT |
| [CuPy PIV](https://github.com/ali-sh-96/openpiv-python-gpu) | d50d505109f9ef9f8ddbd971380ee2ac74f38682 | GPL-3.0 |

FFT 相关、峰值定位、SNR、异常向量验证、多遍形变与 GPU 批量内核均直接调用上游；没有另写 FFT 或峰值求解器。OpenCV 承担候选/积分图/配准，SciPy 承担原六参数鲁棒拟合及矩阵对数，PyAV/FFmpeg 复用 F10V 编码和完整解码检查。GPL 组件不是可重新标为 MIT 的代码；分发时必须保留原许可与对应源码要求。许可证全文在 evidence/licenses。

实际兼容改动：TorchPIV 包导入强制 GUI，薄 namespace 适配绕过 GUI 初始化后调用未改动 PIVbackend；OpenPIV local_median_val 并不支持 backend 参数，按真实 API 调用；CuPy cached s2n_ratio 属性返回 None，读取其已计算 corr.sig2noise；CuPy bool mask cumsum 默认为 int64，而 kernel 读 int32，改为 `cp.cumsum(mask, dtype=self.dtype_i)`。最后一项属于 ABI 修复，不修改相关算法，diff/哈希保存在 evidence/cupy_mask_compat。修复后 B2/B3/B4 三种选窗掩码在共同窗口上的位移差均为 0。

OpenPIV 0.26.1 pip wheel 本身未包含可用 Rust 扩展，使用官方 maturin 构建后显式导入 openpiv_rust 并执行 backend='rust'；没有自动回退。官方 CuPy tutorial test0 实际执行，12600 个有限向量、20 个替换/拒判仅作为示例，不混入正式观测。

FluidImage 和 CPU FFTW/Numba 备用方案 NOT_RUN：Rust 已正常工作，GPU PIV 核心不再是主要瓶颈，本轮不迁移更多框架。配置与版本锁见 evidence/remote_setup 和 evidence/local_requirements_final.lock。
''')
 write('CHANNEL_GEOMETRY_AND_MASK.md',f'''# 通道几何与遮罩

真实原片 1476×1120。可编辑 `configs/channel_mask.json` 保存外多边形、左侧固定结构排除区、裁剪框 `[432,0,912,1120]` 及配准阈值。图像定义通道 977021 px；可观察 mask 805495 px；排除 171526 px。历史 ROI 98304 px，因此观察 mask 面积扩大 8.19 倍。该比例不是有效速度覆盖率。

48-px interrogation、24-px overlap：网格 45×37=1665，要求整个窗口 100% 像素处于 mask，剩 1283 个安全窗口。只要求 50%/90% 分别会放入 91/8 个跨界窗口；100% 为 0。有限检查的墙规则耗时约 0.03–0.10 ms。真实墙边速度偏差没有独立真值，保持 null，不能以中位速度相似证明正确。

初版错误地把中央纹理带排除，遮蔽了已人工确认的泡；审计后将其改为背景风险区，不再几何排除。初版及失败覆盖计数完整保留在 evidence/mask_v1。最终 505 个已审核中心都在安全窗口覆盖范围内；仍不证明整幅边界已校准，视觉边界不确定约 20 px。

每帧只和首帧固定左壁条带 phaseCorrelate，响应>=0.08 且位移<=6 px 才有效。937 帧未触发失效；施加 9/20 px 平移正确拒绝，3 px 通过。静态合成配准存在约 1 px 插值偏移，记录原始值。相机/芯片视野改变后需要重新审计多边形。详见 results/channel_audit.json 与 Figure 1/9。
''')
 write('channel_geometry_audit.md','''# 通道遮罩审查记录

原始首帧、中帧、末帧在 evidence/channel_overlay_0000.png、0450.png、0936.png，最终 mask 为 evidence/channel_mask.png。人工部分标签只读保存，未将未审核区域当背景。中央条带曾排除真实泡，已经修正；原版仍保留。

当前多边形是基于真实影像的可编辑初始分割，不是物理边界真值。whole-window 100% 内含、配准拒判及遮罩哈希均进入状态；20 px 边界不确定性及固定纹理残留均限制其适用性。视觉检查详细清单见 evidence/visual_inspection.json。
''')
 write('FULL_CHANNEL_OBSERVATION.md',f'''# 全通道观测

正式离线逐帧 0–936 全部处理，实际 PTS 0–31.1688 s，呈现 31.2021 s。每帧先全范围 Gaussian 局部暗信号+OpenCV 连通域，使用既有面积 8–450 px、长宽比<=2.8 等规则，再积分图统计候选数、纹理和帧差。单个候选即能激活窗口，完全静止也不会永久失去复查。

候选>=1 或纹理标准差>=2且平均帧差>=0.4；还必须满足整窗安全。规则 GPU 批处理作为主要后端，但发布时严格使用活动及质量 gates，空白区不补零。稀疏版本真正调用官方 mask/gathered-window 内核并单独计时。

最终有效窗口平均占安全窗口 {100*a['valid_window_mean']:.2f}%；有效窗口支持足迹的像素并集平均占观察 mask {100*a['coverage_mean']:.2f}%。这是相关窗口足迹，不是每个像素有独立速度测量。未观测原因、当前 valid/fresh/age、拒判与候选逐帧保存。

在线 latest-slot 允许丢帧，本次规则模式仅扫描处理的 {on['frames']}/937 帧，不能宣称在线全部源帧扫描。离线的 937/937 与在线覆盖严格分开。
''')
 write('RUST_OPENPIV_BENCHMARK.md',f'''# OpenPIV Rust 实测

实际官方编译、openpiv_rust 导入、显式 Rust 执行均通过。新 SciPy/Rust 同图同窗口最大位移差约 3.55e-15 px，SNR 约 1.55e-15，valid 相同。空图拒判及大/小平移测试保存于 stage_b_checks。

{bench}

B0 是旧 PIV 核心接入共同适配器的基准，完整原 F10W 代码/环境另测 historical p95 33.94 ms、全 bbox 351.61 ms。共同适配器中 Rust 相对 B0 的 p95 加速：ROI {b('B0','historical')['timing']['p95']/b('B2','historical')['timing']['p95']:.2f} 倍，全 bbox {b('B0','full_visible_bbox')['timing']['p95']/b('B2','full_visible_bbox')['timing']['p95']:.2f} 倍。与完整旧环境对比包含积分图等工程改进，不能全归因于 Rust。

CPU Ryzen 9 7950X3D，16 核/32 线程。RAYON=4、OMP=2、OPENBLAS=1；后台其他 CPU 项目保留，无强占或无意义压测。
''')
 write('TORCHPIV_GPU_BENCHMARK.md',f'''# TorchPIV GPU 实测

RTX 4090 实际调用官方 PIVbackend.extended_search_area_piv 及 CWS 64→32，多遍和单遍均运行。Headless namespace 仅绕过 GUI 初始化，原相关代码不改。TorchPIV 本地公开接口只给原生 SNR flags，未虚构原始 SNR 数值；最终再做共同局部中值/纹理与六参数检查。

{bench}

TorchPIV 适合规则大批量；稀疏窗拼装在试验中反而较慢。当前 dtype 经上游归一化为 float64，不能声称这是 float32 对等 benchmark。上传和同步由适配器计时，下载融合在官方调用内部，单独下载归因 NOT_SEPARATELY_MEASURED。全流程计时包含此开销。
''')
 write('CUPY_GPU_BENCHMARK.md',f'''# CuPy GPU-PIV 实测

已执行官方 test0、48/24 单遍、64→32 多遍、整通道与稀疏 mask。设备 RTX 4090，CuPy 13.6.0，实际可运行环境版本锁在 evidence。48-px 窗口内部 FFT padded 到 64，与 SciPy/Rust circular48 并非数值等价算法。

{bench}

mask 的 int64/int32 kernel ABI 问题修复前导致 illegal address；原 pilot_v1 的极短时间伴随后端失败，不计入性能收益。修复后共同选窗位移完全一致。缓存以 shape/window/overlap 为键，复用官方 solver 和 GPU buffer。

独立传输微基准：两幅1120×912 float32共 {diag['transfer']['bytes']} bytes，上传 p95 {diag['transfer']['upload']['p95']:.3f} ms，下载 p95 {diag['transfer']['download']['p95']:.3f} ms。它是同步额外 probe，不是从总时间扣减或对融合官方 API 逐项归因。实际显存/GPU利用率见 PERFORMANCE_AND_RESOURCES。
''')
 write('BACKEND_NUMERICAL_COMPARISON.md',f'''# 后端数值比较

严格同实现 B1/B2 达浮点级一致；B0 老版本、TorchPIV/CuPy 有边界、FFT padding、归一化及峰值规则差异。首次历史匹配基准保留旧坐标居中约定；全通道正式调用采用真实窗口左上角排列。两种坐标约定原本可能相差8 px，未把它隐藏为 GPU误差。

CuPy相对SciPy某实际帧最大分量差约 u=0.259、v=0.672 px，SNR及有效支持也不同。不同支持上的平均 EPE 不可单独排名。已知程序运动 final44 条中，CuPy单遍接受30、错接4；细节见 F10X_DEVELOPMENT_REPORT。

解析雅可比加速的首次 real-data 全量核对发现鲁棒优化路径差异：v2/v3 最大六参数差2.745，虽 gates 未变，也未准入。单纯copy不能解决稀疏Huber目标近乎平坦时优化路径不同的问题；正式仅在所有初始残差处于二次区间时使用解析Jacobian，否则调用原有限差分求解器，正式 v5 对全部937帧的 {a['exact_jacobian_equivalence']['checked_groups']} 组可计算结果逐项比较，PIV窗口位移完全相同，组valid变化0，最大六参数差 {a['exact_jacobian_equivalence']['max_six_difference']:.4g}。保存全部失败及修正证据；没有以放宽误差门槛换取通过。
''')
 write('PIV_PREPROCESSING_STUDY.md',f'''# 预处理、时间跨度与多遍

开发22、验证22、最终44条：复用真实微泡图块生成器，新源帧及独立程序种子，但仍是一段录像，没有跨录像泛化真值。比较 raw/local/fixed/causal/CLAHE，gap1/2/3，32/48/64窗及三后端现成64→32多遍。最终配置在 final 打开之前冻结，final 只执行一次。

冻结主线：CuPy，因果 EMA alpha0.08，48/24，gap3实际99.9 ms，SNR1.2，纹理2，六容差[5,5,.2,.2,.1,.1]。背景只用过去及当前，固定背景只用首帧。无未来录像中值冒充在线。

final：局部增强gap2错接 {base['false_accept']}/44；因果gap3错接 {final['false_accept']}/44，接受 {final['accepted']}/44；因果多遍错接 {multi['false_accept']}/44，接受 {multi['accepted']}/44，平均耗时 {multi['mean_total_ms']:.2f} ms vs单遍 {final['mean_total_ms']:.2f} ms。多遍 final 有有价值的错误拒判收益，但开发/验证局部增强多遍没有稳定获胜；未利用 final 改选主线。

长跨度同时有49.95ms中点滞后；突变诊断保存每种gap的区间等效真值，不把更平滑视为即时更准确。大平移240/-120px/s已测试，无硬编码±4px；在线丢帧清空历史、重新预热，不沿用错误时间跨度。

原始各场景 EPE、六 MAE、符号、超容差、正确拒判、支持率和全部失败在 results/research_v1/results；汇总在 delivery_analysis.json。低纹理多遍全粗矢量无效的执行失败保留，不删去困难样本。
''')
 write('ADAPTIVE_WINDOW_SCHEDULING.md','''# 活动窗口与批处理

每帧全通道快速候选+纹理+差分，不依赖旧目标附近。安全且有候选的静止区域仍激活，新出现窗口立即更新；高候选密度>=2也每帧。低活动细 PIV比较1/3/6帧周期，无效/跳过速度为null。

CPU Rust和TorchPIV用官方 sliding_window_array 批量gather成无重叠条带，再调用原完整求解器；CuPy使用官方field_mask，经ABI兼容修复后工作。B2/B3/B4 all-safe/activity/checker三种mask与规则计算共同窗口位移最大差0。

真实部分审核505个中心：466/505在8px内有检测，503/505被活动窗覆盖。初版几何漏掉的真实泡已恢复。合成新位置 probe：静止、小位移、稀疏及dropout中心均被覆盖；boundary漏1个；低对比只覆盖10/120。静态PDMS假候选会出现，几何排除与后续PIV gate控制，但无法证明剩余所有候选都是泡。真实PZT-OFF状态无元数据，NOT_IDENTIFIED。

活动筛选验证仅PARTIAL，不设全通道召回率已验证。规则GPU批处理在相同图像pilot比稀疏gather/mask更快，因此采用规则主线；窗口减少不等于完整管线更快。数据在 channel_audit、post_freeze_diagnostics及channel_pilot_v2。
''')
 write('MULTIRATE_OBSERVATION.md',f'''# 多速率与定速回放

{pt}

多速率减少的是当前精细观测，也可能降低valid空间覆盖。original四策略使用同一原拟合；verified使用已全量核对、在鲁棒活跃输入上保留原求解器的混合适配。不能把v2/v3未准入加速的快时间当正式性能。

窗口每次输出last_update_pts、age_ms、fresh、stale与reason_code；未更新窗口vx/vy均null，不拿旧场充新鲜值。全通道快速扫描始终在每个被处理的源帧执行。每窗口真实更新间隔可以从完整record序列重建，源PTS是真实时间基准。

因果回放producer依原PTS释放，队列最大1；超过负载会丢旧帧，gap后历史清空并再次预热。低处理p50可能包含无PIV预热帧，不能单独称有效观测帧率。ready_pipeline_ms、所有原始行、drops与consumer事件保留。一次性回放调度有CPU抖动，3/6策略不得仅凭单次Hz宣称普遍更优。
''')
 write('LOCAL_GROUP_MOTION.md','''# 局部候选群及六参数

DBSCAN eps55 px/min_samples4产生空间候选集合；跨度较大时以240×240观察单元切分。该划分不识别物理团，更不识别接触、融合或固定单泡ID。observation_key仅指当前空间分区，不保证时间身份。

每组保留中心、可见数、密度、协方差、Rg、长短比、方向与局部六参数。方向在近各向同性ratio<1.15时null。仅支持组ROI内的实际有效PIV矢量进入原F10U拟合，min8、coverage0.18、condition15、residual1.5px不放宽。

数学保持原冻结定义：对endpoint仿射齐次变换取logm/dt，平移在各局部ROI中心计算；div=tr(A)、vort=A21-A12、s1=(A11-A22)/2、s2=(A12+A21)/2。图像x向右、y向下；不是三维体积聚集率。没有全通道强行单一全局拟合。

背景参考速度没有独立可靠观测，恒为null，不输出所谓纯声学额外速度。群体的空间划分仍探索性，非仿射反向运动可能穿过原仿射门槛；人工实验之前需要外部真值和更强模型失配拒判。
''')
 write('UNCERTAINTY_AND_VALIDITY.md','''# 不确定性与有效性

SNR/peak2peak、局部中值、纹理、空间覆盖、条件数、残差与原六参数valid共同决定可发布性。它们是启发式风险指标，不是已校准概率或95%保证。TorchPIV只有原生ratio flag，未捏造SNR幅值。

Figure7在实际已保存局部矢量上复用F10U四象限block bootstrap，24次重采样，显示nominal95%区间。重叠PIV窗口不是独立样本；目前这些区间coverage没有完成独立校准，因此标为calibrated=false，不用于控制准入。

Xue/Charonko/Vlachos相关峰SNR不确定性研究（https://doi.org/10.1088/0957-0233/25/11/115301）说明该方向值得研究，本轮没有声称实现论文的标定不确定度模型。与其他后端一致或曲线平滑不构成真实准确率证据。

final已知程序运动仍有对称旋转超容差和非仿射错接。真实录像没有速度真值、比例尺和PZT日志；REAL_PHYSICAL_VELOCITY_ACCURACY=NOT_VERIFIED。schema中scientific_accuracy_verified=false。
''')
 write('PERFORMANCE_AND_RESOURCES.md',f'''# 性能和资源

同一服务器AMD Ryzen 9 7950X3D、16物理/32逻辑核、约62GiB RAM；RTX4090 24564MiB，driver550.120。环境独立，PyTorch2.14.0+cu126、CuPy13.6.0、numpy2.2.6、SciPy1.17.1、OpenPIV0.26.1及官方Rust编译扩展。版本声明与真实执行均保存。未修改驱动，未停止其他CPU工程；GPU基准期间没有其他计算GPU任务。

{pt}

正式完整离线处理含decode、全扫描、PIV、验证、局部群、JSON序列化；gzip科研存档写盘另计wall_s。回放计时还含实际M4本机消费者发送，单帧decode成本加回，不含真实相机链路。完整regular离线p95={off['p95']:.2f}ms；paced p95={on['p95']:.2f}ms、max={on['max']:.2f}ms。33.3ms和16.7ms目标均未准入。

规则回放CPU核心等价 {on['resources']['cpu_core_equivalents']:.2f}，进程RSS峰值 {on['resources']['rss_peak_mb']:.1f}MiB；设备显存采样峰值 {on['resources']['gpu_memory_peak_mb']:.0f}MiB，GPU平均利用率 {on['resources']['gpu_util_mean']:.1f}%、采样最大 {on['resources']['gpu_util_max']:.0f}%。这是约0.2s采样和设备级显存，不是CUDA allocator瞬态绝对峰值。完整raw资源样本保留。

离线扫描kernel的1000/平均耗时约 {off['scan_kernel_hz']:.1f}Hz，但完整paced实际全通道扫描率仅 {on['full_channel_scan_hz']:.2f}Hz。新鲜有效群帧 {on['fresh_group_frame_hz']:.2f}Hz；空间群状态总计 {on['fresh_group_states_hz']:.2f}/s，不是这么多个物理团。主要剩余成本为CPU局部群拟合、候选/序列化；GPUFFT已较快。内部profiler保存FFT/峰/SNR函数耗时，profile样本不混入warm性能统计。

软件环境准备、数据传输、图表编码不能算GPU科学吞吐。各实际job/identity/resources/wall_s和日志在results、evidence/deployment、logs/remote；完整原片SHA固定。
''')
 write('FAILURE_ANALYSIS.md','''# 失败证据

1. OpenPIV wheel未带Rust；已使用官方源码编译修复，原失败日志保留。
2. API local_median_val不支持backend参数，TorchPIV强制GUI；按真实API薄适配修复，没有伪造GPU结果。
3. CuPy s2n_ratio cache返回None；改读已经计算的原始GPU相关比。
4. 稀疏CuPy bool-mask cumsum默认int64造成kernel非法访问。修复dtype_i，旧pilot低耗时是假运行，不用于收益。
5. 初版central-ridge mask排除了人工泡中心；基于已审核标签修订，旧mask和失败数据保留。
6. 解析Jacobian试验v2/v3部分six差2.745。修复copy后仍有稀疏Huber目标近乎平坦的优化路径差异；正式对鲁棒活跃输入回退原求解器，重跑v5；失败v2/v3不进入正式图表。新的异常点回归测试覆盖该问题。
7. 低对比selector只覆盖10/120；正常人工部分ROI仍漏2/505活动中心。没有修改final配置来隐藏失败。
8. final单遍4/44错接（2对称旋转超差、2局部相反运动）。多遍仍1/44非仿射错接。
9. 定速30fps未达到完整33.3ms p95，latest-slot丢帧并重新预热，显著降低有效输出Hz。没有插值或沿用旧速度冒充新观测。
10. 开发/验证多遍各1次ALL_COARSE_VECTORS_INVALID；全部低纹理困难输入保留。
11. 第一版比较视频渲染遇到JSON null到numpy object，改为显式float NaN，未改任何数值结果；编码日志保留。

未连接Chronos、未控制PZT/AFG、未训练模型、未取得独立真实速度真值。以上均不能用软件测试通过替代。
''')
 write('M4_COMPATIBILITY.md','''# M4 状态兼容

新增M4-channel-3.0，包含窗口并行数组、mask哈希/有效性、safe/active/valid/fresh/stale/age、实际PTS、ROI、backend/device、每组局部six及reason。原v2/v1消费者不被直接替换。

兼容适配必须显式选择当前局部空间群；不选择则输出EXPLICIT_LOCAL_GROUP_REQUIRED和null，禁止自动把全通道拟合为一个群。在线演示明确记录选择的第一个有效空间群key，交由原ZMQ MockConsumer。session、frame_id、PTS、calibration_status和模型身份保留。

100ms TTL后参数null、valid=false、stale=true；帧不连续清空历史，不保留假新鲜速度。tests/test_contracts.py实际验证原JSON schema、显式局部适配和expiry；回放保存consumer_rows/watchdog。consumer接收成功只表明接口工程可用，不表示科学控制准入。
''')
 recommendation=f'''# 最终推荐

本轮**成功建立了更广覆盖、后端明显加速且可审计的全通道原型，但尚未证明全面满足“更快且保持真实测量可靠性”的生产系统标准**。观察mask扩大8.19倍；CPU Rust与SciPy同实现数值一致且较旧版快，CuPy全bbox帧对最优。完整新系统处理更多区域/群，不能直接把其时间与F10W一个ROI比较作等工作量加速。

推荐CuPy规则网格+因果背景gap3作为当前离线/探索性全通道主线，Rust作为数值参考和CPU备选。保留空间支持与null；不推荐稀疏/多速率仅凭单次Hz取代规则方案。64→32多遍final错接更少，适合下一轮预登记的精细参考，不事后改变本轮冻结结论。

完整paced p95 {on['p95']:.2f}ms，fresh group-frame {on['fresh_group_frame_hz']:.2f}Hz，丢{on['dropped']}帧。30/60fps均未准入；selector低对比失败与非仿射错接尚未解决。CPU局部拟合、数据结构/序列化和有明确边界的并行流水线是下一优化重点，不是再造FFT或无意义GPU压测。

进入真实相机/人工PZT实验前，应获取独立录像、已知载物台运动/比例尺、完整PDMS区域人工审核、弱对比与壁边真值、硬件时间戳/采集丢帧证据；预登记局部非仿射拒判和多遍对照；先验证Chronos实际链路，保持人工控制及隔离日志。PZT-OFF段必须有实际设备日志，不能从视频猜测。当前不得启动自动物理控制或训练世界模型。

{scope}
''';write('FINAL_RECOMMENDATION.md',recommendation)
 answers=[('是否读取复用F10W？','已核查干净工作树，基线62b453c6c49ab2beadb89868d048a12117ef06db；继承F10U数学和F10V视频，保留历史哈希。'),('Rust是否实际执行？','已编译官方扩展、显式backend=rust、导入和已知平移执行通过，不自动回退。'),('Rust速度提高多少？',f"共同适配器ROI p95从{b('B0','historical')['timing']['p95']:.2f}到{b('B2','historical')['timing']['p95']:.2f}ms；全bbox从{b('B0','full_visible_bbox')['timing']['p95']:.2f}到{b('B2','full_visible_bbox')['timing']['p95']:.2f}ms。"),('TorchPIV实际GPU运行？','是，RTX4090官方API单遍与CWS多遍。'),('CuPy实际GPU运行？','是，官方示例、单遍、多遍、规则与稀疏全通道均运行。'),('同窗口谁最快？','当前warm帧对全bbox为CuPy；ROI中CuPy与Torch接近。完整观测不等于PIV核心。'),('数值谁最稳定？','严格同实现Rust/SciPy达浮点级一致；其他后端算法不同，未获真实精度排名。'),('GPU提高完整性能吗？','相同帧对的预处理+PIV+一次六拟合明显提高；完整全通道仍受CPU群拟合限制，不能宣称30fps。'),('遮罩覆盖多少？','805495px观察mask、1283安全窗，8.19倍历史ROI；边界未物理校准。'),('是否每帧观察整通道？',f'离线937/937是；定速规则处理{on["frames"]}帧、丢{on["dropped"]}帧，不能称在线每源帧全扫描。'),('活动筛选是否漏泡？','是：人工503/505覆盖；合成低对比10/120覆盖，属于明确失败。'),('墙边是否剔除正确？','按整窗100%遮罩，0跨界窗；真实物理边界及墙边速度准确率未验证。'),('稀疏策略节省多少？','原拟合p95 regular77.67、sparse1 79.16、sparse3 74.53、sparse6 68.77ms；均伴随不同有效覆盖，不能等覆盖归因。'),('哪种背景最有效？',f'本轮预登记验证选择因果EMA；final比局部增强错接{base["false_accept"]}→{final["false_accept"]}/44，但并非全部条件最佳。'),('多遍有收益吗？',f'final错接{multi["false_accept"]}/44，相比单遍{final["false_accept"]}/44减少；更慢且支持改变，不用final重新调参。'),('哪个间隔合适？','当前冻结gap3=99.9ms，改善小位移/固定背景，但有49.95ms中点滞后和突变混合；不是所有运动最优。'),('实际有效面积多少？',f'有效窗平均{100*a["valid_window_mean"]:.2f}%，其支持足迹并集平均占mask{100*a["coverage_mean"]:.2f}%；不代表逐像素速度。'),('无证据区是否不可观测？','是，vx/vy=null，并记录原因及stale，不零填充或插值。'),('是否有独立局部six？','有，当前空间候选集合分别拟合；没有经过物理验证的群身份。'),('是否保留数学约定？','是，logm(endpoint)/dt、局部中心、原六分量及冻结quality gates。'),('完整937帧是否完成？','四种原策略及正式修正版各937/937；完整视频也逐帧对应。'),('完整p50/p95/max？',f'正式离线{off["p50"]:.2f}/{off["p95"]:.2f}/{off["max"]:.2f}ms；定速{on["p50"]:.2f}/{on["p95"]:.2f}/{on["max"]:.2f}ms。'),('全扫描多少Hz？',f'定速完整处理实际{on["full_channel_scan_hz"]:.2f}Hz；扫描kernel约{off["scan_kernel_hz"]:.1f}Hz，两者不同。'),('新鲜six多少Hz？',f'至少一组新鲜有效的帧率{on["fresh_group_frame_hz"]:.2f}Hz；空间群状态合计{on["fresh_group_states_hz"]:.2f}/s，非物理群数。'),('p95低于33.3吗？','否；也不能宣称支持60fps。'),('资源峰值？',f'规则回放CPU核心等价{on["resources"]["cpu_core_equivalents"]:.2f}，RSS峰值{on["resources"]["rss_peak_mb"]:.1f}MiB；设备显存采样{on["resources"]["gpu_memory_peak_mb"]:.0f}MiB，GPU均值{on["resources"]["gpu_util_mean"]:.1f}%。'),('独立已知运动测试通过？',f'数学和简单平移检查通过；44条新程序种子接受{final["accepted"]}、错接{final["false_accept"]}，所以总体仅部分通过，非新录像独立真值。'),('还有误测吗？','对称旋转超差、局部反向非仿射错接、弱对比漏检，均保留。'),('实际图视频有哪些？','figures/Fig_F10X_1…9均PNG/PDF；两段937帧全长MP4、一段60帧selected clip，清单/PTS/哈希/视觉记录齐备。'),('哪些科学结论已验证？','后端执行、同条件速度、数值一致性、软件null/时间契约及已知程序运动误差已量化；真实速度、3D聚集、融合与声学归因未验证。'),('比F10W改善的证据？','可观察mask面积8.19倍，同条件Rust/CuPy更快；程序测试背景错接受减少，但真实精度没有独立真值，不全面宣称胜过F10W。'),('后续主线？','CuPy规则批量+因果gap3探索性主线；Rust保守参考；预登记多遍精细对照。'),('进入真实相机/PZT前缺什么？','独立真实运动/比例尺、弱对比/壁边标注、相机时序延迟、非仿射拒判、PZT-OFF实际日志及人工受控实验。')]
 write('F10X_DEVELOPMENT_REPORT.md','# F10X 中文开发与科研报告\n\n'+recommendation.replace('# 最终推荐\n\n','')+'\n## 实测同条件后端\n\n'+bench+'\n## 全视频及定速结果\n\n'+pt+'\n## 33项问题逐项回答\n\n'+'\n\n'.join(f'{i}. **{q}** {ans}' for i,(q,ans) in enumerate(answers,1))+'\n\n## 数据与审计\n\n实际执行身份、初次失败、不可变源码、环境锁、原片SHA、配置冻结及final打开时间都已保存；详细路径见README和RECOVERY_MANIFEST。全部科研图仅来自已保存的正式离线v5、定速v6及冻结结果，未修改原始标签或历史F1–F10W。当前Git身份由DELIVERY_IDENTITY.json和最终提交记录给出。科研进度修改见evidence/progress_document。')
 gate=dict(FULL_SOURCE_VIDEO_VERIFIED=True,CHANNEL_MASK_CREATED=True,CHANNEL_MASK_AUDITED='IMAGE_AUDIT_NOT_PHYSICAL_BOUNDARY_TRUTH',FULL_CHANNEL_SCANNING_COMPLETED=True,FULL_CHANNEL_FRAME_COVERAGE={'mode':'OFFLINE','processed':937,'source':937},ACTIVE_WINDOW_SELECTION_VALIDATED='PARTIAL_LOW_CONTRAST_FAILED',OPENPIV_RUST_BACKEND_ACTUALLY_LOADED=True,OPENPIV_RUST_BENCHMARK_COMPLETED=True,TORCHPIV_GPU_BENCHMARK_COMPLETED=True,CUPY_GPU_BENCHMARK_COMPLETED=True,BACKEND_NUMERICAL_COMPARISON='RUST_SCIPY_EQUAL_DIFFERENT_GPU_ALGORITHMS_DOCUMENTED',PIV_GROUND_TRUTH_RECOVERY='PARTIAL_4_FALSE_ACCEPTS_OF_44',MULTIPASS_VALIDATION='COMPLETED_WITH_FAILURES_AND_COST',BACKGROUND_SUPPRESSION_VALIDATION='PROCEDURAL_IMPROVEMENT_NOT_PHYSICAL_TRUTH',FULL_CHANNEL_MOTION_VALIDITY='SUPPORT_GATED_EXPLORATORY',LOCAL_GROUP_SIX_PARAMETER_VALIDITY='LOCAL_SPATIAL_GROUPS_NONAFFINE_FAILURE_REMAINS',FULL_PIPELINE_P95_MS=on['p95'],OFFLINE_FULL_PIPELINE_P95_MS=off['p95'],TARGET_33MS_PASS=False,FULL_CHANNEL_SCAN_HZ=on['full_channel_scan_hz'],FRESH_VALID_GROUP_MOTION_HZ=on['fresh_group_frame_hz'],CAUSAL_REPLAY_PASS=True,CAUSAL_REPLAY_SCOPE='ENGINEERING_TIMESTAMPS_QUEUE_AND_NULL_CONTRACT_NOT_REALTIME_TARGET',CAUSAL_REPLAY_DROPPED=on['dropped'],FULL_DURATION_VIDEO_VALIDATED=all(v['fully_decoded'] and v['decoded_frames']==937 for v in read('videos/full_video_manifest.json')['videos']),FIGURES_VALIDATED=read('evidence/visual_inspection.json')['figures_pass'],M4_COMPATIBILITY_PASS=True,REAL_PHYSICAL_VELOCITY_ACCURACY='NOT_VERIFIED',LIVE_CAMERA_PASS=False,CAMERA_END_TO_END_LATENCY='NOT_MEASURED',PZT_CONTROL_EXECUTED=False,AFG_CONTROL_EXECUTED=False,WORLD_MODEL_TRAINED=False,TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,EXPLORATORY=True,PHYSICS_UNCALIBRATED=True,FINAL_ASSESSMENT='ENGINEERING_DELIVERED_SCIENTIFIC_ADMISSION_PARTIAL')
 dump(ROOT/'FINAL_GATE.json',dict(identity=identity(),**gate));print('reports and gate written',flush=True)
if __name__=='__main__':main()
