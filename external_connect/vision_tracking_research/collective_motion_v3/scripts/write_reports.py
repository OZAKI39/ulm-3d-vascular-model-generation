"""Generate the Chinese delivery from actual saved measurements only."""
import sys,json,collections
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
def read(p):return json.loads((OUT/p).read_text())
def write(p,s):
    path=OUT/p;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(s.rstrip()+'\n')
def main():
    audit=read('results/full_video_source_audit.json');full=read('results/full_video_analysis_summary.json');valid=read('results/full_video_validity.json');an=read('results/model_analysis.json');freeze=read('results/model_freeze.json');runtime=read('results/full_video_runtime.json');bench=read('results/MODEL_BENCHMARK_COMPLETE.json');assets=read('evidence/remote_assets.json');official=read('results/UnLite_official_examples.json');ot=read('results/ot_numerical_audit.json');hist=read('results/historical_diagnostics.json');clip=read('results/real_comparison_clip.json');env=read('evidence/server_environment_final.json');video=read('videos/full_video_manifest.json');video3=read('videos/comparison_video_manifest.json')
    replay={m:read('results/causal_replay_'+m+'.json') for m in ['PIV','LRAFT_ST','CoTracker3','POT_unbalanced']}
    header='EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false\n\n'
    table='| Method | Finite global MAE n | MAE (vx, vy, div, vort, s1, s2) | Geometry accepted / wrong | Risk-gated accepted / wrong |\n|---|---:|---|---:|---:|\n'
    for m,r in an['methods'].items():
        table+=f'| {m} | {r["finite_global_mae_n"]} | '+', '.join(f'{x:.4f}' for x in r['finite_global_mae'])+f' | {r["base"]["accepted"]} / {r["base"]["bad_accepted"]} | {r["gate"]["accepted"]} / {r["gate"]["bad_accepted"]} |\n'
    rt='| Method | Processed /937 | Dropped | Full ready p50 / p95 ms | Fresh valid states | Fresh six Hz |\n|---|---:|---:|---:|---:|---:|\n'
    for m,r in replay.items():rt+=f'| {m} | {r["n_processed"]} | {r["n_dropped"]} | {r["full_pipeline_ms"]["p50"]:.2f} / {r["full_pipeline_ms"]["p95"]:.2f} | {r["n_fresh_valid"]} | {r["fresh_six_hz"]:.3f} |\n'
    count_info=f'历史ROI候选数为{full["count_range"][0]}–{full["count_range"][1]}，首尾{full["count_first_last"][0]}→{full["count_first_last"][1]}；候选中心从(151.78,176.14)变为(114.23,145.53) px（ROI局部坐标）。这不是同一物质团的轨迹，也不能把计数变化等同于聚集、融合或通量。'
    fulltext=f'''# 完整录像 PIV 报告

{header}原始文件存在，SHA-256为`{audit['sha256']}`，H.264、1476×1120、27,832,433 bytes。共937帧，PTS 0–31.1688 s，帧间隔33.3 ms，画面呈现时长31.2021 s；容器31.232 s包括容器/音轨时长差异，没有为对齐容器时长虚构视频帧。无异常PTS或完全重复解码图像；没有采集计数器，无法据此排除相机端掉帧。

历史F10U只分析360–479帧，F10V真实回放660–779帧，均为120帧。本轮从0至936逐帧分析，937/937覆盖，无缺失中间时段；有效运动887/937（94.664%），50帧拒判，其中2帧WARMUP、26帧不一致/非仿射向量、22帧支持不足。94.664%是冻结质量门槛通过率，**不是真实准确率**。

固定ROI `[650,300,256,384]`；冻结OpenPIV 48 px窗口、24 px重叠、SNR≥1.2、纹理标准差≥2、至少8向量、覆盖≥0.18、条件数≤15、残差≤1.5 px。三帧只使用首尾两帧及真实PTS差，不是三帧联合神经网络。另存candidate_5、candidate_10、raw_3、local_3共4685条方法/帧记录，不用对照结果修改主方法。

{count_info}

有效片段中vx中位数−0.841 px/s、vy中位数24.394 px/s；vy最大221.777 px/s，末段vx上升、涡量和剪切改变，同时候选和空间支持下降。所有极值保留；缺乏独立运动真值，不能把末段突变认定为真实流场加速。0.067–1.199 s有多段向量不一致；2.165、4.196、4.329 s另有拒判。28.605–30.236 s间支持不足较集中；完整不连续区间列表在`results/full_video_analysis_summary.json`。

全画幅做了每帧候选分布和128 px分块时序证据扫描。图中通道边界是第一帧固定的目视导线，非标定分割；静态条纹可能被检测为候选，候选加时序变化也不等于材料运动。没有宣称已获得全画幅真实速度场。保留固定背景/低纹理/边界候选比例与拒判诊断。

全部离线方法和全画幅扫描耗时{runtime['invocation_wall_s']:.3f} s。两个完整MP4均937帧，1600×1000，原PTS映射误差为0，不插帧；Video3为360–423帧的64帧selected clip，不能代表全视频。

## 完整文件定速回放（与离线覆盖分开）

{rt}
每个方法都实际输入937帧，单槽最新帧，跳帧ID完整保存。断帧清空历史窗口，TTL=100 ms，源结束后过期清空通过。Hz分母为源活动区间31.1688 s；若计入尾帧保持时间，应采用31.2021 s而非改变源帧数。相机未连接，capture timestamp缺失，camera-to-controller延迟NOT_MEASURED。
'''
    write('FULL_VIDEO_PIV_REPORT.md',fulltext);write('results/full_video_summary.md',fulltext);write('FULL_VIDEO_SOURCE_AUDIT.md','# 原片与历史覆盖审计\n\n'+fulltext.split('## 完整文件定速回放')[0]);dump(OUT/'FULL_VIDEO_PIV_SUMMARY.json',dict(source=audit,coverage=valid,analysis=full,scope='QUALITY_VALIDITY_NOT_TRUTH_ACCURACY'))
    write('ZERO_SHOT_MODEL_COMPARISON.md',f'''# 开源零样本模型对照

{header}新开发46、验证92、校准46、最终138条16帧序列，来自原片8个不同帧/4个不相交区域。已知解析运动施加于真实暗泡图块，保留固定纹理、漏检、密集、非仿射和不可观测条件；不是实测六参数真值，也不是微泡动力学仿真。所有模型使用同一13→15帧、66.6 ms真值。CoTracker额外可见0–15历史帧，其他方法仅使用端点；时间上下文不同已明确披露。

只在开发集选择raw/inverted和OT尺度5/8 px。验证集选择经验风险≤5%、至少10接受的质量阈值；校准集只统计误差分位数。`FINAL_OPEN_RECORD.json`证明配置、权重和阈值冻结后才开启最终集。最终集未用来调整阈值。误差容差保持[5,5,0.2,0.2,0.1,0.1]。独立分块不等于跨录像独立，复制样本不能当成独立实验。

{table}
MAE只在有有限原始估计且有全局真值的样本上计算，分母不同；缺失和拒判均另报，不能凭Partial OT的18个简单样本低MAE认定最优。主四法共同89个样本上，PIV/Unbalanced OT/LRAFT-ST/CoTracker3平移MAE分别约(1.554,1.361)/(0.322,0.442)/(2.639,0.971)/(0.360,0.429) px/s，CoTracker在同一支持上的四项形变误差最低。Unbalanced OT与CoTracker有继续研究价值，但不构成真实场景精度证明。

25%诊断覆盖取35/138：PIV、CPD、Unbalanced OT、CoTracker均0/35错误；LRAFT 1/35；UnLite32/35；Partial OT只有18个基础接受，不能达到该覆盖。风险曲线是事后描述，不作为新部署阈值。有限同录像样本的0错误不代表总体风险为0。

UnLite官方示例8对均已实际运行，EPE为{min(r['epe_px'] for r in official['rows']):.3f}–{max(r['epe_px'] for r in official['rows']):.3f} px；同一实现迁移到微泡后，基础114个接受中111个超过容差或错误接受，验证找不到合法阈值，最终正式输出0。官方样例可运行不等于微泡适用。

F10U既有18类seed7201/source840、F10V旧final指定16例共34例×7方法=238次额外推理，均实际运行，全部只标记历史诊断。没有把已暴露数据称为新盲测。
''')
    provenance='# 官方来源、权重与复用\n\n'+header
    for name,r in assets['third_party'].items():
        provenance+=f'## {name}\n\n官方仓库：{r["url"]}\n\n固定提交：`{r["commit"]}`；许可：`{r["license_status"]}`。\n\n'
    provenance+='CoTracker LICENSE.md为CC-BY-NC 4.0；本地科研验证，禁止据此假定商用许可。LRAFT和UnLite未发现LICENSE，授权不明；第三方源码/权重没有纳入本项目Git交付。POT为MIT，直接调用官方求解器。OpenPIV、pycpd、F10U几何转换器均复用。\n\n'
    for r in assets['files']:
        if '/models/' in r['path']:provenance+=f'- `{r["path"]}`：{r["bytes"]} bytes，SHA-256 `{r["sha256"]}`。\n'
    provenance+='\nLRAFT权重来自官方README链接的Google Drive学生Dataset1（ID 18tTvU8WjquJLY84GHrxrgeHQXpfZ5UQh）；CoTracker来自官方Hugging Face facebook/cotracker3/scaled_online.pth；UnLite使用官方仓库自带权重。下载脚本和日志保留。UnLite仅将CuPy旧memoize/compile_with_cache接口适配为memoize/RawModule(backend=nvcc)，CUDA核内容、网络和权重未改变；原失败、diff、before/after SHA保存在evidence/unlite_compat。\n'
    write('PRETRAINED_MODEL_PROVENANCE.md',provenance);dump(OUT/'models/INDEX.json',dict(official=assets['third_party'],weights=[r for r in assets['files'] if '/models/' in r['path']],redistribute_weights=False))
    numerical=[r for r in ot['rows'] if 'numerical error' in r['warning'].lower()];changes=[r for r in ot['rows'] if r['base_valid']!=r['previous_base_valid']]
    write('OPTIMAL_TRANSPORT_VALIDATION.md',f'''# POT匹配与数值复核

{header}POT 0.9.7官方Partial/Unbalanced熵正则求解器；空间平方距离、面积比、峰值差异组成代价，最大位移300 px/s×dt+2 px，超限代价100；reg=0.05、reg_m=1、partial mass=0.85，最多500步。开发选择distance_scale=8。行质量、匹配峰≥0.5、归一化熵≤0.6和真实有限性检查后交给原F10U拟合。

最终138例×2变体复核中，明确numerical errors共{len(numerical)}例（Partial 96、Unbalanced 22），全部拒判。18个接受的Partial例存在中间divide警告，但最终矩阵有限、质量约0.85、求解日志误差接近机器精度；保留警告，不把所有warning混为有效收敛证据。数值错误拦截修正前后接受集合变化{len(changes)}例。冻结测试结果未被重新挑选。

Unbalanced最终87/138接受、0个已知错误；CPD风险门槛后58/138、1错误，PIV98/138、3错误。真实整片Unbalanced仅66个新鲜有效状态，不能把合成成功外推到真实录像。端点数量变化/孤立假点能使低正则核出现下溢，本轮没有在最终集上改用更容易通过的超参数。应在新开发数据上研究稳定/log域求解器及显式unmatched项，再用新测试集评价。

软对应不是长期材料ID；对称结构、漏检和固定纹理仍有歧义。tests/test_contracts.py验证已知平移和异常点下安全拒判。
''')
    write('SELF_SUPERVISED_ADAPTATION.md','''# 自监督适配决策

状态：NOT_RUN。本轮没有训练或更新任何模型权重，没有自监督收益可报告。

先完成了真实零样本实验。LRAFT和CoTracker有运动敏感性，但分别存在固定背景假接受和当前短窗口在线吞吐不足；UnLite虽能运行官方示例，在微泡迁移上幅值偏差及错误接受严重。已读取官方UnLite multiscaleUnsupervisorError（双向光度、流一致性、二阶平滑），机制存在；没有重新发明损失或网络。

在单录像、无真实材料运动真值的当前条件下，减小光度损失也可能只是追踪PDMS背景。选择保留零样本基线，将新增独立录像/已知载物台运动验证和CoTracker持续流式缓存列为优先事项；不为完成训练次数而改动已冻结模型。本轮任务把适配设为有明确依据时的可选阶段，因此如实记NOT_RUN，图表不制造“适配后”结果。
''')
    write('SIX_PARAMETER_VALIDATION.md','''# 六参数数学与时间定义

所有方法复用collective_motion_v1/src/affine.py。有限端点变换H的生成元为logm(H)/dt；v0在固定ROI几何中心，图像x右、y下。div=A11+A22，vort=A21−A12，s1=(A11−A22)/2，s2=(A12+A21)/2。平移px/s，其余1/s。没有用协方差变化替代速度场。

共同误差限[5,5,0.2,0.2,0.1,0.1]未变。密集流先在真实检测点取样，经前后向一致性≤1.5 px后拟合。CoTracker用官方online predictor的可见性/置信门槛获取短对应。OT采用低歧义软对应。局部流无可信仿射/空间支持不足则逐字段null。

8项本地与服务器合同测试通过：非方形ROI有限logm解析真值、流场xy顺序和越界NaN、M4 v2与v1适配、无效值必须null、断帧窗口重置、TTL、真实本机ZMQ消费者以及OT已知平移/数值失败。OpenCV多通道border NaN行为已显式修复；旧版本任一通道NaN同样被finite过滤，该修复不改变已冻结有效支持。
''')
    write('M4_COMPATIBILITY.md','''# M4观测兼容

复用F10V的M4-collective-2.0 Schema、VersionedGate及to_v1；未修改原1.0工程。全部正式定速回放逐消息schema.validate，并实际发送到本机MockConsumer。使用单槽最新帧、host monotonic clock、源PTS，断帧清空窗口，重复/逆序/过期拒判，断流后100 ms将运动和分布过期清空；不能复用旧状态冒充新观测。

available_capture_timestamp=null，科学准确率未验证，hardware_control_executed=false，training_ready=false。相机到控制器端到端延迟NOT_MEASURED。最初回放因缺calibration_version被Schema阻止，未被当成完整回放；修复后重新完整执行。权重SHA最初每帧重复计算的慢版本保留在results/replay_before_hash_cache_fix，正式版模型加载时计算一次。
''')
    resources=f'''# 真实运行时间与资源

{header}独立服务器根`/workspace/f10w_vision_20261017`、独立.venv，不覆盖F10V或物理环境。RTX4090 24 GiB，CPU 16物理/32逻辑核；完整环境及有限32 MiB读写探针见evidence/server_environment_final.json，读缓存不等于磁盘峰值性能。其他CPU科研作业未停止；这些是实测该服务器负载下的性能，未宣称独占CPU实验。

全视频离线PIV+5配置+全帧扫描：{runtime['invocation_wall_s']:.3f} s；新322样本生成{read('data/DATASET_COMPLETE.json')['wall_s']:.3f} s；完整零样本分阶段基准{bench['wall_s']:.3f} s；历史诊断{hist['wall_s']:.3f} s；真实比较片段{clip['wall_s']:.3f} s；UnLite官方8例{official['wall_s']:.3f} s。四个修正后完整定速回放合计{sum(r['wall_s'] for r in replay.values()):.3f} s。两个全长视频渲染+完整解码{video['wall_s']:.3f} s，比较视频{video3['wall_s']:.3f} s。各阶段分开报告，不把安装、网络传输、人工分析和不同阶段重叠时间混成“GPU训练时间”；初次失败与慢回放日志保留。

基准GPU框架峰值分配{bench['peak_cuda_allocated_mb']:.1f} MiB，设备采样最高2324 MiB / 93%利用率；采样是设备范围，不是每个模型独占计数。所有必要Supervisor任务已终止，无隐含后台训练。

{rt}
PIV回放保留历史质量门槛；AI/OT增加验证风险门槛，基准表另给所有方法统一风险门槛对照。完整管线包括候选检测、分布、模型/求解、几何、Schema和发布；all-frame与ready窗口分布分别保存。推理延迟不等于端到端相机延迟。CoTracker采用官方16帧online窗口、step8，但当前适配器每个短窗重新初始化；断帧后重建16帧，首次可输出PTS=0.4995 s。它是因果短时原型，不是优化后的持续在线吞吐上限。OT首个调用包含库冷启动，最大约547 ms；稳态快不代表高有效更新率。
'''
    write('RUNTIME_AND_RESOURCES.md',resources)
    real_counts={m:sum(all(r['published_estimate']['valid'].values()) for r in rows) for m,rows in clip['methods'].items()}
    write('REAL_CHRONOS_OBSERVATION.md',f'''# 真实Chronos录像观测

{header}原始937帧都已完成PIV，详细见FULL_VIDEO_PIV_REPORT。整段因果回放均以源PTS定速，实际处理/丢弃可追溯。未连接真实相机，LIVE_CAMERA_PASS=false。

固定预登记selected clip为360–423帧，CoTracker只看当前及之前16帧，不使用未来；PIV/OT/LRAFT使用过去第2帧至当前。64个展示帧中各方法质量/风险门槛接受数：{real_counts}。CoTracker局部轨迹与可见泡群运动相符是观察性描述，不是已核实单泡材料身份。LRAFT或OT拒判时视频显示null，不把其低可信原始预测绘成可信速度。

{rt}
REAL_SIX_PARAMETER_ACCURACY=NOT_VERIFIED，跨实验泛化未验证。候选中心、密度和2D散度均不代表3D物理团收缩率。
''')
    write('GENERALIZATION_AND_FAILURES.md','''# 泛化、失败和保留限制

SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED。八个素材帧、四个空间区、独立随机种子仍共享一次录像的照明、背景、光学系统；同一源帧上的复本不是独立实验。全录像PIV及可视化已经使用所有源帧，因此这是已知解析运动目标的留出测试，不是原始场景从未接触的完全盲测。新最终测试只开一次用于科学比较，后续OT重复只做数值审计，历史F10U/F10V仅诊断。

保留：PyArrow26/NumPy1.26不兼容（隔离环境降到PyArrow19）；数据脚本语法错误；初始非仿射场不满足预设残差条件（在任何模型选型前把生成强度45→90 px/s，旧不完整数据保留）；UnLite旧CuPy API失败和最小diff；Schema缺少必需字段；OpenCV双通道边界NaN；逐帧重复权重哈希拖慢；OT核下溢/数值错误；全部科学误判和未通过模型。

UnLite在官方粒子样例有效但微泡域幅值不足；LRAFT固定背景条件6个中基础误接受5个，风险门槛后只保留1个；CoTracker合成成绩改善但16帧重建造成新鲜输出很低；Unbalanced OT合成低误差却在真实录像频繁数值/支持拒判。没有为好看改预测、裁掉失败帧或放宽容差。实时参数不确定度没有在独立真实真值上校准。
''')
    write('SOURCE_AND_REUSE_AUDIT.md','''# 代码与历史审计

实际从F10V完整交付80440a8fc5f15716a2a5835b8eb2d509e20a30c1创建隔离worktree和dev/f10w-pretrained-full-video-20261017，未用较旧F10U同步覆盖最新工作。原F10V工作树启动时干净；服务器物理仓库66d4ceb未改。适用目录未发现额外AGENTS；服务器运行前已读/etc/vast-agents-guide.md并使用独立Supervisor标签。

F10U CandidateDetector/distribution/OpenPIV/pycpd/affine/state直接只读复用；F10V只复用图块splat、M4 v2、资源和视频编解码辅助。新代码只提供官方模型适配、数据协议、运行编排、结果审核与可视化；没有重写RAFT、CoTracker、Sinkhorn或视频编解码器。third_party被Git排除。冻结执行快照及代码SHA见各COMPLETE/EXECUTION文件。

原片、人工标注、F10U/F10V、相机/AFG和物理历史共5214个受保护文件逐一复核，结果见evidence/final_integrity.json。科研总进度文档仅在全部实测后按用户要求备份并添加F10W段落，未删历史结论。没有GitHub push。
''')
    qa=[
    ('指定MP4是否存在，基本信息是什么？',f'存在。SHA-256 `{audit["sha256"]}`；1476×1120、H.264、937帧；画面呈现31.2021 s，容器31.232 s，实际PTS 0–31.1688 s。'),
    ('此前是否已全片处理？','没有。F10U为360–479、F10V回放为660–779，各120帧。'),
    ('本轮是否完整PIV？','是。固定历史ROI冻结单遍PIV，0–936每帧均有处理记录。'),
    ('解码、分析、不可观测各多少？','937解码、937分析、887质量有效、50拒判（2 warmup、26不一致向量、22支持不足）。真实准确率未知。'),
    ('是否有未分析中间区间？','没有。缺帧ID列表为空，覆盖率100%；定速回放的丢帧另报。'),
    ('群中心、密度、数量如何变化？',count_info),
    ('六参数显著变化？','vy中位24.394 px/s、vx中位−0.841 px/s。初段有大波动，末段vx、涡量/剪切明显改变且支持下降。不能把未标定的极值解释为真实流场变化。'),
    ('PIV何时失效？','0.067–1.199 s有多段不一致；2.165、4.196、4.329 s亦有失败；28.605–30.236 s有多段空间支持不足。全部离散区间已保存。'),
    ('PDMS背景是否造成异常？','存在风险，原灰度/局部增强/候选增强结果差异保留；真实片没有材料真值，不能把每个异常归因于背景。已知运动的固定背景测试中LRAFT基础误接受5/6、UnLite6/6，证实这些方法会受背景干扰；本轮该条件PIV和CoTracker为0/6，但样本小。'),
    ('两个完整MP4是否合格？','均937帧，1600×1000，完整解码、时间戳一致、首尾及失败帧视觉核查通过，没有插帧。'),
    ('LRAFT官方权重是否运行？','已严格加载官方学生Dataset1并完成开发/验证/校准/最终/历史/真实测试。最终风险门槛后55/138接受，其中1错误。'),
    ('CoTracker3是否能追踪真实短时运动？','官方online模型已在真实候选上给出因果短时对应；比较视频能看到与泡群一起移动的点。缺少独立材料ID真值，不能称已验证真实泡追踪准确率。'),
    ('OT是否更合理？','已知运动下Unbalanced OT在共同89例的平移误差约0.322/0.442 px/s，低于PIV；最终87/138接受0已知错误。但真实片新鲜有效只有66次，不能全面替代PIV。'),
    ('UnLite运行与表现？','官方8例EPE 0.078–0.240 px；微泡最终基础114接受中111错误，风险门槛后0接受。环境成功，迁移科学准入失败。'),
    ('哪些仍有近零或错误运动？','UnLite明显幅值不足，部分近零；LRAFT有固定背景误匹配；CPD对称点集/错误对应、PIV支持不足/非仿射、CoTracker漏检/变形失败均保留。'),
    ('哪些不需长期ID？','全部当前方法都不要求长期单泡ID：PIV/光流局部位移、CoTracker短时对应、CPD/OT端点集合配准。'),
    ('已知真值下谁最准确？','不能只按不同支持上的MAE排名。在主四法共同89例上，Unbalanced OT的vx最低，CoTracker的vy和四个形变分量最低；CoTracker总体覆盖较高。'),
    ('相同有效率谁错误接受最少？','25%（35/138）质量排序诊断覆盖下PIV、CPD、Unbalanced OT、CoTracker均0错误；LRAFT1、UnLite32；Partial OT达不到该覆盖。不能声称唯一胜者。'),
    ('完整延迟谁最低？','当前Unbalanced OT ready管线p50/p95约7.34/12.14 ms，但包括快速拒判，不能据此称有效观测最好；PIV36.26/45.88 ms，LRAFT46.35/55.63 ms，CoTracker79.66/85.56 ms。'),
    ('是否自监督及有无收益？','NOT_RUN，无收益可报告。已读官方UnLite自监督机制；当前优先处理可观测性、真实真值与持续在线缓存，不在同录像上为训练而训练。'),
    ('新AI真正优于PIV吗？','CoTracker在受控合成运动的共同支持上有精度改善，但没有真实真值，而且在线新鲜输出1.60 Hz低于PIV21.40 Hz。尚不能声称整体优于PIV。'),
    ('是否缺独立录像与真实真值？','仍缺。SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED；真实六参数准确率NOT_VERIFIED。'),
    ('资源和计算时间？',f'RTX4090实际运行三个模型，基准峰值CUDA分配1516.8 MiB、设备采样最高2324 MiB。全片离线146.744 s、数据生成7.109 s、正式模型基准63.626 s、修正后四路定速回放{sum(r["wall_s"] for r in replay.values()):.3f} s；其他阶段、冷启动和失败分别保留于资源报告。'),
    ('已验证与探索性分别是什么？','已验证文件完整性、逐帧覆盖、解析运动数学、官方推理可运行、当前样本比较、Schema/TTL、视频编码与哈希；真实物理准确、跨录像泛化、材料ID和3D聚集解释均未验证。'),
    ('下一轮优先AI还是标定？','优先新增独立显微录像、像素/时间/载物台已知平移旋转真值；保留PIV主线。并行做CoTracker官方持续缓存/查询更新、OT稳定求解和固定背景遮罩，再用新留出数据评价；有可验证收益再考虑有限自监督。'),
    ('是否允许自动PZT控制？','没有依据。相机直播未测、端到端延迟未测、真实六参数未验证；PZT/AFG均未操作、世界模型未训练、scientific_control_admitted=false。')]
    final='# F10W最终中文报告（26项）\n\n'+header+'**结论：没有模型已被证明能够可靠替代PIV作为当前真实微泡群观测主线。继续保留PIV；CoTracker3和Unbalanced OT作为有潜力的研究对照。**\n\n'
    for i,(q,a) in enumerate(qa,1):final+=f'{i}. **{q}** {a}\n\n'
    final+='## 已知运动测试\n\n'+table+'\n## 完整定速回放\n\n'+rt
    write('FINAL_RECOMMENDATION.md',final);write('F10W_DEVELOPMENT_REPORT.md',final)
    write('README.md',f'''# F10W-Vision

{header}已完成全937帧冻结PIV、三种官方预训练模型+POT/CPD的实际对照、四路全文件定速回放、11组PNG/PDF、3段MP4与审查。结论：暂不替换PIV。

- [26项最终中文报告](FINAL_RECOMMENDATION.md)
- [完整录像报告](FULL_VIDEO_PIV_REPORT.md)
- [模型对照](ZERO_SHOT_MODEL_COMPARISON.md)
- [代码与复现入口](ACTIVE_CODE_INDEX.md)
- [机器Gate](FINAL_GATE.json)
- 本地视频位于`videos/`；两段full_duration均937帧，比较视频明确为64帧selected clip。

本地目录：`{ROOT}`；分支`dev/f10w-pretrained-full-video-20261017`。服务器独立目录`/workspace/f10w_vision_20261017`。没有GitHub push，第三方完整仓库、权重、大缓存和原片不纳入Git；恢复索引见REMOTE_ARTIFACTS.json和models/INDEX.json。

隔离环境复现：本地`.venv/bin/python -m unittest discover -s tests -v`；已下载结果下运行`figures_full.py`、`analyze_models.py`、`full_videos.py`、`comparison_video.py`重绘。所有脚本位于scripts。使用新的F10W_OUTPUT_ROOT运行新实验，不能覆盖原结果。模型推理需服务器CUDA环境；最终锁定依赖在evidence/server_bundle/setup/requirements_final.lock。NumPy1.26使用PyArrow19，不能直接升到26。
''')
    write('ACTIVE_CODE_INDEX.md','''# 活跃代码与复现入口

| 范围 | 入口 | 说明 |
|---|---|---|
| 路径/身份/只读复用 | src/wcommon.py | 环境变量F10W_OUTPUT_ROOT/VIDEO/THIRD_PARTY/MODELS |
| 原片审计 | scripts/audit_source.py | PyAV完整解码、PTS、逐帧SHA |
| 完整离线PIV | scripts/full_video.py | 冻结历史参数，5配置，每50帧保存可验证chunk |
| 读取大记录 | src/full_records.py | 分块SHA校验后流式读取 |
| 官方模型 | src/motion_models.py | LRAFT/UnLite/CoTracker thin adapters；POT官方函数 |
| 新素材协议 | configs/model_protocol.json | 46/92/46/138，选型与最终测试隔离 |
| 数据 | scripts/generate_model_data.py | 解析运动+真实图块，禁止称物理真值 |
| 统一评价 | src/benchmark.py / scripts/model_benchmark.py | 先开发、验证、校准冻结，再final；已有freeze拒绝重调 |
| 实际定速回放 | scripts/causal_replay.py | 全原片、latest slot1、断帧重置、M4 v2→v1 |
| 历史/真片/官方例 | scripts/extra_diagnostics.py | real/historical/ot/official；历史仅诊断 |
| 图表 | scripts/figures_full.py / analyze_models.py | 保存数据重绘，无预测修饰 |
| 视频 | scripts/full_videos.py / comparison_video.py | 复用PyAV/既有Writer，源PTS、完整解码 |
| 合同测试 | tests/test_contracts.py | 数学、越界、OT失败、Schema、TTL、ZMQ |
| 环境和资产 | scripts/remote_asset_audit.py / audit_server_environment.py | 官方SHA、许可证、锁定依赖、有限I/O探针 |
| 发布快照 | scripts/deploy_snapshot.py | Git archive不可变目录、独立Supervisor标签，无GitHub push |
| 最终核查 | scripts/verify_delivery.py | 全帧、manifest、分割、历史保护、视频/图表 |
| 报告/进度 | scripts/write_reports.py / update_progress.py | 原始结果出表，保留中文历史及diff |

恢复大结果：`python scripts/restore_assets.py fullvideo_delivery.tar.gz`及`python scripts/restore_assets.py models_complete.tar.gz`；脚本先校验SHA，只补缺失文件，不覆盖已有不同内容。

大数据可读取入口：pyarrow.parquet.read_table('results/full_video_piv_per_frame.parquet')；局部向量用full_video_vectors.parquet（dx/dy为区间像素位移，速度需除实际PTS差）。完整候选/图像质量/各方法对象在full_video_chunks，远程映射与SHA见REMOTE_ARTIFACTS.json。原始录像直接使用指定本地文件，模型权重用models/INDEX.json恢复。F10V历史依赖仍在相邻collective_motion_v2，不需要拷贝或改写。
''')
    dump(OUT/'results/report_generation.json',dict(identity=identity(),input_sha256={p:sha(OUT/p) for p in ['results/model_analysis.json','results/full_video_analysis_summary.json','results/model_freeze.json']},scientific_recommendation='KEEP_PIV',selfsupervised='NOT_RUN'))
if __name__=='__main__':main()
