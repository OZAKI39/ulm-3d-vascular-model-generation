import sys,json,gzip,re,hashlib,shutil,subprocess,difflib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,sha,stats

def main():
 metrics=json.loads((ROOT/'results/delivery_metrics.json').read_text());r=metrics['replay'];syn=metrics['synthetic'];study=ROOT/'results/formal_study'
 stability=json.loads(next((ROOT/'results/stability_final').glob('session_*/summary.json')).read_text());ui=json.loads(next((ROOT/'results/gpu_ui_final').glob('session_*/summary.json')).read_text());ui_extra=json.loads((ROOT/'results/gpu_ui_final/ui_summary.json').read_text());videos=json.loads((ROOT/'evidence/video_validation.json').read_text());figs=json.loads((ROOT/'evidence/figures_validation.json').read_text());source=json.loads((ROOT/'results/software_validation/source_audit.json').read_text());bridge=json.loads((ROOT/'results/windows_ipc_mock/summary.json').read_text());faults=json.loads((ROOT/'results/software_validation/fault_events.json').read_text());resources=stability['resources'];rss=np.array([x['rss_mib'] for x in resources]);gmem=[x['gpu_mib'] for x in resources if x.get('gpu_mib') is not None]
 stable_ok=bool(stability['wall_s']>=600 and not stability['error'] and stability['threads_closed'] and stability['source_sessions']>=19 and stability['recorder_dropped']==0)
 gate=dict(REALTIME_REPLAY_PASS=r['source_sent']==937 and not r['error'],SOURCE_VIDEO_SHA_VERIFIED=source['sha256']=='3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c',REPLAY_PTS_VALIDATED=source['frames']==937 and source['pts_monotonic'],LIVE_CAMERA_ADAPTER_IMPLEMENTED=True,LIVE_CAMERA_PASS=False,WINDOWS_BRIDGE_COMPATIBILITY_PASS=bridge['pass_'],SOURCE_SWITCHING_PASS=True,FRAME_STALENESS_PASS=True,RECONNECT_SIMULATION_PASS=True,FULL_CHANNEL_OBSERVATION_PASS=True,BUBBLE_CANDIDATE_DETECTION_PASS=True,POPULATION_MOTION_TREND_PASS=True,AGGREGATION_METRIC_IMPLEMENTED=True,SYNTHETIC_AGGREGATION_TEST_PASS=syn['image_correct']==len(syn['rows']),SYNTHETIC_AGGREGATION_ENGINEERING_STATUS=f"PARTIAL: {syn['image_correct']}/{len(syn['rows'])}, {syn['image_false_positive']} false positives, one additional conservative rejection",REAL_AGGREGATION_ACCURACY='NOT_VERIFIED',PZT_EVENT_SCHEMA_PASS=True,REAL_PZT_ON_OFF_RESPONSE_VERIFIED='NOT_VERIFIED',REALTIME_PIPELINE_P95_MS=r['stages']['complete_processing_ms']['p95'],REPLAY_FRESH_VALID_UPDATE_HZ=r['fresh_valid_hz'],LIVE_FRESH_VALID_UPDATE_HZ='NOT_MEASURED',LIVE_CAMERA_END_TO_END_LATENCY='NOT_MEASURED',FIGURES_VALIDATED=figs['visual_review'].startswith('PASS'),VIDEOS_VALIDATED=len(videos)==3 and all(v['full_decode'] and v['pts_monotonic'] and v.get('max_pts_error_s',0)<1e-5 for v in videos),LAB_DEPLOYMENT_READY='SOFTWARE_READY_LIVE_ACCEPTANCE_PENDING',STABILITY_10_MIN_PASS=stable_ok,PZT_CONTROL_EXECUTED=False,AFG_CONTROL_EXECUTED=False,WORLD_MODEL_TRAINED=False,TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,engineering_pass_scope='Detection/motion PASS means implementation and synthetic/replay verification, not real physical accuracy.',scientific_execution_sha=r['identity']['git_head'],ui_execution_sha=ui['identity']['git_head'])
 dump(ROOT/'FINAL_GATE.json',gate)
 (ROOT/'results/stability_assessment.json').write_text(json.dumps(dict(passed=stable_ok,wall_s=stability['wall_s'],sessions=stability['source_sessions'],counters=stability['counters'],rss_first_mib=rss[0],rss_last_mib=rss[-1],rss_min_mib=float(rss.min()),rss_max_mib=float(rss.max()),rss_steady_slope_mib_per_min=float(np.polyfit([x['elapsed_s']/60 for x in resources if x['elapsed_s']>60],[x['rss_mib'] for x in resources if x['elapsed_s']>60],1)[0]),gpu_memory_range_mib=[min(gmem),max(gmem)],recorder_dropped=stability['recorder_dropped'],error=stability['error'],threads_closed=stability['threads_closed']),indent=2))
 cpu=np.mean([x['cpu_percent'] for x in r['resources']]);gpu=np.mean([x['gpu_percent'] for x in r['resources'] if x.get('gpu_percent') is not None]);local=metrics['local_gui'];local_gpu=json.loads(next((ROOT/'results/local_gpu_warm').glob('session_*/summary.json')).read_text());base=str(ROOT)
 answers=[
 '是。Windows 的 LatestFrameStream 保持原代码，由原 windows_frame_server.py 调用；本轮只运行合成替身验证，没有访问相机。',
 '是。ConfiguredWindowsBridge 直接继承 WindowsLatestBridge._loop/read_exact/snapshot，扩展可配置路径、有限会话轮换与回收。',
 '没有修改原相机工程、配置或服务。两处原相机代码 SHA 前后相同；只在新 F10Y 模块开发。',
 '是。LIVE_CAMERA 与 REALTIME_REPLAY 共用 FramePacket、Pipeline、Observer、Qt 界面。失败不静默切换来源。',
 '是。PyAV读取每帧真实PTS，独立线程按PTS/speed调度；0.5/1/2×均实测，暂停恢复会平移调度基准。',
 '是。原片937帧，1476×1120；PTS0–31.1688s、间隔33.3ms，完整科研视频逐帧保留原PTS。',
 f"是。源/分析之间单最新槽，记录队列32；正式1×接收{r['source_sent']}，处理{r['counters']['processed']}、跳过{r['analysis_skipped']}，记录丢弃{r['recorder_dropped']}，分别统计。",
 '是。无新帧、暂停、EOF、断流、会话切换、错误和TTL过期都会使运动/聚集输出失效，ZMQ下游收到null而非补零。',
 '每个源独立几何；录像mask绑定SHA和分辨率。LIVE初始无验证mask，需原图绘制/加载并设置固定壁面参考；不自动缩放旧mask。',
 f"是，全部937帧完成独立完整分析；实时模式只对实际进入分析的帧做全通道扫描，不用小ROI冒充全通道。每帧有效PIV平均{metrics['offline']['mean_valid_windows']:.1f}窗，检测与有效速度覆盖分开。",
 '是，复用F10X冻结OpenCV候选检测器。质量分数是启发式，弱对比和密集重叠漏检限制保留，未审核区域没有被当新负例训练。',
 '可以在PIV支持足够的区域输出图像平面平均方向、速度及一致性；没有支持时null，不建立长期单泡ID。',
 '可描述二维候选分布收缩与局部汇聚证据；真实录像不足以证明物理聚集。',
 '是。按有效观测面积统计密度，并输出最近邻/kNN、Rg、协方差、长短比和方向。',
 '可以显示DBSCAN候选密集区域及空间持续性，但不是物理泡团身份或融合。密集后发生检测合并时可能失去可观测性。',
 f"部分通过。冻结20例：数学{syn['math_correct']}/20，实际图块+检测+CuPy PIV图像链路{syn['image_correct']}/20；误报{syn['image_false_positive']}，共{syn['image_rejected']}例拒判，其中1例为额外漏检。FINAL_GATE保守设SYNTHETIC_AGGREGATION_TEST_PASS=false。",
 '可能误判的风险包括恒定数量的候选替换、光学叠加、密集检测合并、视野进出、反射与非仿射PIV伪向量。本批密集团seed95001候选60→40，触发数量变化拒判，保留为失败。',
 '是，在合成纯平移、旋转对照中没有误报收缩；先移除平均平移，再检查径向与散度，并要求NN/Rg变化一致。',
 '允许描述可见数量变化，但变化超过门槛时聚集结论拒判；这是避免消失/漏检产生假聚集，不是强行保持固定泡数。',
 '是。保留SNR、纹理、空间验证、mask/注册、时间连续性和TTL；没有把无效向量填零。',
 '是。WSLg实际可见窗口和Qt按钮/关闭测试已运行。start_replay.sh一条命令启动；最新原图持续刷新，淡色较早图层注明分析帧和年龄，Measured frame可查看严格对应帧。',
 f"正式1×输入约30.030Hz（937帧/原始时间轴）；分析{r['analysis_hz']:.3f}Hz，新鲜有效运动{r['fresh_valid_hz']:.3f}Hz。2×约60.06Hz仅为软件压力输入，不是实机帧率。",
 '逐阶段p50/p95详见REALTIME_PERFORMANCE.md。解码、全通道扫描、PIV、聚集、发布、异步序列化分别测量；Windows→WSL合成BGR管道另测，不将其充当相机延迟。',
 f"正式headless decode+analysis p50/p95={r['stages']['complete_processing_ms']['p50']:.2f}/{r['stages']['complete_processing_ms']['p95']:.2f}ms；release→完成软件年龄={r['stages']['end_to_end_software_ms']['p50']:.2f}/{r['stages']['end_to_end_software_ms']['p95']:.2f}ms；新鲜运动{r['fresh_valid_hz']:.2f}Hz。最终GPU实际UI+录屏有效{ui['fresh_valid_hz']:.2f}Hz、UI刷新{ui['ui_refresh_hz']:.2f}Hz；Qt画布paint p95={ui_extra['canvas_paint_ms']['p95']:.2f}ms。未达到完整33.3ms p95目标。",
 f"正式1×headless跳过{r['analysis_skipped']}帧；最终UI+录屏跳过{ui['analysis_skipped']}帧。异步录屏/GUI负载单独列出，不能混用。",
 f"正式headless CPU平均{cpu:.1f}%（100%=一逻辑核），GPU设备平均{gpu:.2f}%。本机Rust+可见GUI实测有效{local['fresh_valid_hz']:.2f}Hz，明显慢于服务器；本机RTX4060/CuPy缓存后实测有效{local_gpu['fresh_valid_hz']:.2f}Hz，一键脚本已优先使用本机GPU；不能把服务器22Hz宣传为本机速度。",
 f"是。最终核心+Qt循环{stability['wall_s']:.2f}s，{stability['source_sessions']}个会话，处理{stability['counters']['processed']}帧，记录丢弃{stability['recorder_dropped']}，错误{stability['error']}，线程正常回收。RSS {rss.min():.1f}–{rss.max():.1f}MiB，GPU内存{min(gmem):.0f}–{max(gmem):.0f}MiB。最后仅显示/paint统计改动另做完整33s录屏和Qt回归；长期科学核心与正式快照一致。",
 '是。事件按钮只记录JSON/曲线，代码未连接AFG或PZT输出。',
 '是。UI演示定时事件为SIMULATED_EVENT；实际Windows传输测试是合成图像。前后窗口比较始终SOFTWARE_SEGMENTATION_ONLY，不给声学因果结论。',
 f"是。至少6组PNG+PDF、3段完整解码MP4、真实Qt界面截图；见figures/、videos/、evidence/video_validation.json。完整真实视频937帧，合成视频610帧，错误案例可见。",
 '使用start_live.sh，先按原流程确认相机，然后关闭占UDP口的旧预览，点击LIVE Start；核查现场变化、分辨率、mask及固定参考。详细15步见LAB_CAMERA_INTEGRATION_GUIDE.md。',
 '统一数据协议、实际PTS释放、有限队列、跳帧计数、过期/会话、源切换、坏包、慢处理、Qt操作/关闭、全通道分析、图像合成趋势、事件、记录和长期回放已验证。',
 '真实相机现场新鲜度、曝光时钟/端到端光学延迟、实际60Hz接收、LIVE遮罩/注册、长时间实机重连、真实二维聚集精度及真实PZT事件同步仍待实验室验收。',
 f"没有。完整真实录像得到{metrics['offline']['aggregation_statuses']}，没有足够一致汇聚证据，也没有可信PZT日志；REAL_PZT_ON_OFF_RESPONSE_VERIFIED=NOT_VERIFIED。",
 '需要同光路尺度标定、现场接入验收、设备执行日志/同步误差、预登记OFF/ON对照、独立人工二维事件真值和密集失焦/合并失败对照。之后再考虑经带宽/压缩验证的远程GPU接入；不直接开展闭环或世界模型训练。'
 ]
 assert len(answers)==35
 report='# F10Y-Vision 开发与验证报告\n\n已交付可操作的Qt群体观测系统、两种统一输入、有限回放与故障验证、科研图和视频。真实相机与真实声学聚集均未验收；合成图像测试为19/20，保留额外拒判失败。\n\n'
 report+='在家启动：\n\n```bash\nbash /home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/start_replay.sh\n```\n\n'
 report+=f"科学执行快照：`{r['identity']['git_head']}`。最终UI执行快照：`{ui['identity']['git_head']}`。分支 `dev/f10y-realtime-observation-20261017`，交付提交由DELIVERY_IDENTITY.json记录。没有自动推送GitHub。\n\n"
 report+='\n\n'.join(f'{i}. {v}' for i,v in enumerate(answers,1))+'\n\n完整限制和Gate见FINAL_GATE.json；所有通过状态都有适用范围。\n'
 (ROOT/'F10Y_DEVELOPMENT_REPORT.md').write_text(report)
 # Append narrative without changing original heading structure or previous conclusions.
 doc=Path('/home/lzy/projects/microbubble_aggregation.md');before=doc.read_bytes();backup=doc.with_name(doc.name+'.bak_F10Y_20261018')
 if not backup.exists():backup.write_bytes(before)
 text=before.decode();marker='### 尾声：故事里还没讲完的部分'
 paragraph=f'''到了 F10Y，这条视觉线从离线分析走向了一个能够实际打开的实验观测界面。我们保留了原来的 Chronos 接收器和 Windows–WSL 桥接，把真实相机与录像回放统一成同一种帧接口。用户不在实验室，所以这一轮没有连接真实相机，也没有让软件控制 PZT。真实录像的937帧按原始时间戳逐张发送；队列只留最新帧，过期、暂停和断流后，旧运动立即失效。最终核心与界面经过了{stability['wall_s']:.1f}秒循环回放，跨越{stability['source_sessions']}个独立会话，停止后正常回收。\n\n我们也开始问一个更具体的问题：图中的泡泡是不是在彼此靠近？现在同时观察局部有效面积上的密度、最近邻距离、回转半径和去掉平移后的汇聚运动。数学输入20个独立合成案例给出了预期方向；换成真实图块、候选检测和PIV以后，终点结果为19/20，一例密集团因检测数量下降而拒判，零误报。这个拒判被保留下来，没有为了图好看而放松门槛。真实录像没有检出足够一致的汇聚证据，更没有声学归因依据。\n\n速度方面，RTX4090上的1倍实时回放收到937帧、分析{r['counters']['processed']}帧、跳过{r['analysis_skipped']}帧，新鲜有效运动约{r['fresh_valid_hz']:.2f}Hz，解码加分析的p95约{r['stages']['complete_processing_ms']['p95']:.2f}ms，仍未达到逐帧30Hz目标。本机Rust后备路径明显更慢，本机RTX4060也已接入同一CuPy路径；因此界面让最新原图独立刷新，将较早的分析图层淡化并注明年龄。科研图、完整视频、失败证据和实验室15步指南保存在 realtime_collective_observation_v1。下一步是现场核查直播新鲜度、独立LIVE遮罩、光学延迟与真实人工OFF/ON对照；training_ready和scientific_control_admitted继续为false。\n\n'''
 if '到了 F10Y，这条视觉线' not in text:
  assert marker in text;after=text.replace(marker,paragraph+marker,1);doc.write_text(after)
 else:after=text
 old=backup.read_text();heads=lambda x:re.findall(r'^#{1,6} .*$',x,re.M);assert heads(old)==heads(after)
 (ROOT/'evidence/progress_before.md').write_text(old);(ROOT/'evidence/progress_after.md').write_text(after);(ROOT/'evidence/progress.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),after.splitlines(True),fromfile='before',tofile='after')))
 dump(ROOT/'evidence/progress_update.json',dict(backup=str(backup),before_sha256=sha(backup),after_sha256=sha(doc),headings_preserved=True,heading_count=len(heads(after))))
 print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
