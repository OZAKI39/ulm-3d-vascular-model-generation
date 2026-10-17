"""Write reproducible reports from actual saved results; pending metrics remain null."""
from pathlib import Path
import json,sys,hashlib,datetime,subprocess
ROOT=Path(__file__).resolve().parents[1]

def load(p):return json.loads((ROOT/p).read_text())
def write(name,text):(ROOT/name).write_text(text.strip()+'\n')
def main():
 A=load('results/live_opencv_20s/summary.json');B=load('results/live_cellpose_20s_resumed/summary.json')
 P=load('results/paired_comparison/summary.json');audit=load('results/annotation_audit_current.json')
 gui=load('evidence/annotation_repair_sam_smoke.json');videos=load('evidence/video_validation.json')
 source=load('data/candidate_manifest.json')['rows'][0]
 gate={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'LIVE_CAMERA_PASS':True,
       'ANNOTATION_TOOL_PASS':gui['PASS'],'EXPERT_LABELS_READY':False,'EXPERT_GROUND_TRUTH_PENDING':True,
       'annotation_snapshot_utc':audit['utc'],'expert_confirmed_frames':audit['expert_confirmed_frames'],'draft_frames':audit['draft_frames'],
       'candidate_frames':audit['candidates'],'BUBBLE_DETECTION_VALIDATED':False,'SHORT_TRACKING_VALIDATED':False,
       'GROUP_SEGMENTATION_VALIDATED':False,'GROUP_SHAPE_NOT_ADMITTED':True,
       'REALTIME_PROCESSING_PASS':{'opencv':'ROI_SOFTWARE_THROUGHPUT_PASS_ONLY','cellpose':'FAIL_60FPS_REFERENCE'},
       'strict_every_frame_16_7ms_pass':False,'validated_bubble_update_hz':0,'validated_group_update_hz':0,
       'CROSS_SESSION_GENERALIZATION_NOT_VERIFIED':True,'END_TO_END_LATENCY_NOT_MEASURED':True,
       'EXPOSURE_TIMESTAMP_UNKNOWN':True,'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED',
       'HARDWARE_CONTROL_DISABLED':True,'CLOSED_LOOP_CONTROL_NOT_RUN':True,'LED_TEST_NOT_RUN':True,'SAM2_NOT_RUN':True,
       'finetuning':'NOT_RUN_EXPERT_LABELS_PENDING','training_ready':False,
       'accuracy_metrics':None,'status':'ENGINEERING_DELIVERY_COMPLETE; EXPERT_EVALUATION_PENDING',
       'executions':['results/live_opencv_20s','results/live_cellpose_20s_resumed','results/paired_comparison'],
       'failed_or_interrupted':['logs/install_plan.log','logs/annotation_gui_smoke.log','results/live_cellpose_20s/INTERRUPTED.json']}
 write('FINAL_GATE.json',json.dumps(gate,indent=2))
 write('SMART_TRACKING_V2_DEVELOPMENT.md',f'''# F10T-Vision-2 开发与修复记录

本轮基于上一轮视觉调查继续开发，Git 基线 `9c22d207be5318f9d45b014caaf52bd04d345328`，分支 `dev/f10t-vision2-20261016`。原 Windows Chronos、WSL Chronos、Tektronix 和历史研究保留；无明确安全清理对象，未删除历史。65 个原设备核心文件 before/after SHA 检查见 evidence。两个环境均通过 pip check。

## 标注故障及修复

初版官方 micro-sam Annotator2D 在显示窗口前同步初始化模型，换帧又关闭、重建 Qt 窗口并重复初始化；提示点和类别层需要手动反复选择，退出缺少自动保存。用户反馈实际卡住，无法按原说明完成。用户所述整次 WSL 异常退出原因没有系统级证据，标为未证实，不把它全部归为某一个库。

修复采用任务允许的 napari 原生 Labels 后备路径，保留五类数据格式。默认手工描画不加载 SAM；同一窗口切换图片，中文分步按钮，自动新编号，草稿定时 / 换帧 / 关闭保存，保存错误阻止退出，恢复上次序号，单实例写入锁。原版 SAM 通过 micro-sam API 在可取消 QProcess 中运行，90 秒超时，结果校验图片 ID 和 token，接受后清空提示点，不覆盖已有实例。没有重写绘图、分割或关联算法。

实际 Qt 窗口检查通过：两个实例的多边形描画、保存重开、跨分辨率切换与返回、自动保存、SAM 接受及取消。换帧约 {gui['checks']['frame_switch_s']:.3f} s；SAM 任务约 {gui['checks']['sam_wall_s']:.2f} s，事件循环最大采样间隔 {gui['checks']['event_loop_max_gap_s']:.3f} s。测试数据隔离于 evidence，修复当时两张草稿的六个文件逐项保持原 SHA；用户随后新增的正式修订不受此检查限制。操作说明见 [ANNOTATION_READY.md](ANNOTATION_READY.md)。

## 实际视觉工作

仅运行 OpenCV + Norfair 和 Cellpose 3.1.1.3 cyto3 + Norfair。原 SAM 只辅助人工标注，不作为第三条在线检测路线；没有 SAM 2、LED、AFG 或 PZT 控制。早期合并安装因 NumPy 1 / 2 依赖冲突失败，改为两个隔离环境；WSLg OpenGL 通过 xcb + glx 修复。所有失败日志保留。

真实输入复用原 Windows `LatestFrameStream`。仅增加有界 OS 管道 GET/response 桥，原始 BGR 或明确 ROI，无缩放。Windows 与 Linux 单调时钟不直接相减，以发送端解码后年龄 + WSL GET 往返 + WSL 后续年龄构成保守软件年龄上界。它不包含未知曝光、相机内部缓存或网络前段延迟。

每端仅保留最新图像；落后跳帧。Norfair 无固定方向、流速或 F10S 先验，最近 2–5 个观测点用于短时速度；丢失返回 null 坐标，时间跳变 / session 改变 / 大间隙重置 ID。零件分割的面积、Rg、Q 只存为 UNREVIEWED_IMAGE_REGION，正式 group_shape=null。稀疏 LK 仅为图像运动，不是液体速度。

## 有限执行与复现

真实相机各测 20 s；另有约 12 s 的稀疏实时候选采集。Cellpose 首次正式试验被宿主重启中断，无总结，保留 INTERRUPTED；恢复使用新目录，不覆盖失败。成对录像测试 12 帧 × 全图 / ROI × 两种方法，另有 90 帧、约 3 s 的严格时间顺序离线回放，共耗时 {P['wall_s']:.3f} s（含模型初始化 / 输出）。没有通过视频播放速度宣称 Cellpose 实时。

执行身份在每个 results 子目录的 execution_identity.json。测试时 HEAD 与所有实际源码 / protocol SHA 均保存；有少量脚本在提交前运行，因此必须同时查 Git 状态和文件 SHA，不能声称每次执行工作区都干净。代码已分别提交本地分支，未自动推送。

后续还补上过期帧入口和流结束时的显式空状态发送；性能表来自此前不可变记录，没有伪称新版本重新跑过同一测量。软件故障契约检查覆盖这些输出规则。图件、视频及其脚本和哈希均保留。

两段视频在首张人工标签正式确认前生成，其中“尚无专家参考”的文字反映生成当时的状态；随后收到的训练帧确认已纳入固定审计快照、图 1 和图 3。视频中的检测和跟踪候选仍未通过独立验证，不据此修改原始输出。
''')
 write('DATASET_SPLIT_AND_PROVENANCE.md',f'''# 数据划分与来源

原始录像：`{source['source_file']}`。SHA-256：`{source['source_sha256']}`。原始图像 1476×1120，录像实际 PTS 约 30.03 fps；不是 60 fps 实时流，更不代表传感器采集帧率。没有物镜倍率、适配器、焦平面、像素到 µm 标定，保持 null。

72 张录像候选：train 帧 0–270、validation 360–570、test 660–936，各 24 张；每组含一个连续 8 帧短片段（120–127、450–457、780–787）。组间保留时间间隙，不随机拆分相邻帧。完整候选 manifest 保存原始 PTS / time_base、图像和视频 SHA、尺寸、ROI、类别审核状态。

另 12 张原始 1280×1024 实时帧约每秒采样一次，共 84 张。直播接收器没有暴露传感器 PTS，因此 source_pts=null；另存 Windows host decode time / UTC、相机配置、session、原始 PNG SHA。source_sha 在直播帧中指 PNG 文件，不冒称原始流文件哈希。直播数据全部 QUARANTINE_UNASSIGNED，未经新分组和审核不得训练或评价。现有样本是否包含真实团、突然转向、不同焦平面等需要专家确认；不能把 12 张相近直播画面当作 12 种独立条件。

当前 {audit['expert_confirmed_frames']} 张正式确认、{audit['draft_frames']} 张草稿。原 20/21 正例匹配没有混入真值。保存整个源图，任何 ROI 推理另记 xywh，未改变检测分辨率；演示视频只在显示时放大 2 倍。

同视频内部初步划分不证明跨实验泛化：CROSS_SESSION_GENERALIZATION_NOT_VERIFIED。下一次采集应补独立注射轮次 / 实验日、清晰单泡与团、焦平面变化、亮度变化、边缘进出、急转弯和聚集 / 分离片段，并记录真实光学配置。

## 人工完成后继续

先运行 `.venv/bin/python -B scripts/review_labels.py` 核查实际文件、SHA、完整审核和类别冲突。用户完成后可用 `--freeze-name <唯一名字> --export-train` 固定所确认修订并导出 train / validation 的标准 Cellpose PNG + masks TIFF。测试集不会导出到训练目录；含 UNCERTAIN / OUT_OF_FOCUS 的帧默认排除训练，防止标准 Cellpose loss 把未知区域误当背景。未确认帧一律不导出。

评估规则已在 protocol.json 预登记：单泡 6 px 中心匹配、团实例 IoU 0.5，Hungarian 一对一；不确定 / 失焦区域显式忽略，预测与忽略区重叠超过 0.5 时拒判。阈值仅可在训练 / 验证组选择。当前只有部分训练帧确认，验证/测试真值仍不足，尚未运行检测或短轨迹准确率评价，也未据测试结果调参。
''')
 rows=[]
 for scope in ['FULL_FRAME','ROI']:
  for route in ['opencv','cellpose']:
   r=P['scope_results'][f'STATIC_PAIRED/{scope}/{route}'];rows.append(f'| {scope} | {route} | {r["n"]} | {r["detect_ms"]["mean"]:.2f} | {r["detect_ms"]["p95"]:.2f} | {r["candidates"]["mean"]:.2f} |')
 write('OPENCV_VS_CELLPOSE.md','''# 两条识别路线的实际结果

OpenCV 复用原局部背景 / 对比度及连通域实现，Norfair 负责关联。额外保留较大连通域作为未审核区域；面积分成“小 / 大”只是候选分组，不是学到的单泡 / 团语义。原阈值保留，本轮没有在无真值测试集上调参。

Cellpose 3.1.1.3 的 cyto3 是本轮唯一 AI 检测模型，作为成本较低的 Cellpose 系替代。灰度 channels=[0,0]，diameter=30（内部尺度比为 1），无增强，flow_threshold=0.4、cellprob_threshold=0、min_size=15。预测 score 是 cell probability 的映射，不是已校准的微泡置信度。没有训练或微调。

以下为同一录像、相同原始像素、相同帧、相同 RTX4060 Laptop/CPU、同一计时方法的成对检测；每次 CUDA 同步后结束计时。全图为 1476×1120，ROI 为 x650 y350 w180 h350。这里只测 detection，不把全图 detection 时间误称完整在线时间。

| 范围 | 路线 | 帧数 | 平均检测 ms | p95 ms | 平均候选数 |
|---|---|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

90 帧因果离线回放的完整引擎 p95：OpenCV {P['scope_results']['CAUSAL_REPLAY_OFFLINE/ROI/opencv']['engine_ms']['p95']:.2f} ms，Cellpose {P['scope_results']['CAUSAL_REPLAY_OFFLINE/ROI/cellpose']['engine_ms']['p95']:.2f} ms。它是按顺序实际计算的离线对照，不是实时输入计时。

目视实际叠加：OpenCV 选出很多小对比区域；cyto3 候选更少，常偏向较大的模糊结构。这个差异**不能证明**后者背景误检少，也不能证明前者微泡召回高。当前没有完整人工真值，Precision、Recall、F1、定位 / 分割和 ID 误差均 NOT_EVALUATED。直播 Cellpose 在当前 ROI 没有可用候选，并频繁超出 100 ms 软件年龄门槛；失败保留。

当前工程选择：保留 OpenCV + Norfair 作为实时候选管线；Cellpose 暂不承担 60 fps 每帧输出。待人工标签确定实际漏检和误检后，才判断是否进行一次有限微调，不能仅凭“新模型”或候选数量宣布优胜。

许可证及权重身份见 [MODEL_AND_ENVIRONMENT.md](MODEL_AND_ENVIRONMENT.md)。
''')
 lines=['# 真实相机在线性能','', '两次 LIVE_CAMERA 试验各约 20 s，真实 Windows 原接收器、原始 256×384 ROI，未缩放。两次是先后运行，场景及后台负载可能不同；不作为公平准确率对比。成对录像基准另列。OpenCV 试验时人工标注 GUI 仍在运行，Cellpose 恢复试验在 WSL 重启后运行；GPU 指标为整卡共享负载。','', '| 路线 | 帧数 | 处理 Hz | 新鲜候选 Hz | 验证微泡 Hz | 处理 p95 ms | 处理最大 ms | 跳过源帧 | 处理后过期 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for d in [A,B]:
  lines.append(f'| {d["route"]} | {d["processed_frames"]} | {d["processed_hz"]:.3f} | {d["fresh_candidate_hz"]:.3f} | 0 | {d["stage_ms"]["software_ms"]["p95"]:.3f} | {d["stage_ms"]["software_ms"]["maximum"]:.3f} | {d["unprocessed_source_frames_in_span"]} | {d["stale_result_rejected"]} |')
 for d in [A,B]:
  lines+=['',f'## {d["route"]} 分项毫秒','', '| 项目 | mean | p50 | p95 | maximum |','|---|---:|---:|---:|---:|']
  fields={**d['stage_ms'],'GET_roundtrip_ms':d['get_roundtrip_ms'],'pipe_write_ms':d['windows_pipe_write_ms'],
          'mock_transport_ms':d['mock_transport_ms'],'mock_consume_ms':d['mock_consume_ms']}
  for k,v in fields.items():
   if v:lines.append('| '+k+' | '+' | '.join(f'{v[t]:.3f}' for t in ['mean','p50','p95','maximum'])+' |')
  lines+=['','资源（CPU 为一核 100% 基准；整卡 GPU 含其他进程）：','', '| 项目 | mean | p50 | p95 | maximum |','|---|---:|---:|---:|---:|']
  for k,v in d['resources'].items():
   if v:lines.append('| '+k+' | '+' | '.join(f'{v[t]:.3f}' for t in ['mean','p50','p95','maximum'])+' |')
  lines+=['',f'Windows接收 {d["windows_receiver"]["received_frames"]} 帧，失败 {d["windows_receiver"]["failures"]} 次，nominal 60 fps；处理区间 host decode 频率 {d["input_host_decode_hz"]:.3f} Hz，可能包含启动缓冲追赶，不等同于曝光帧率。管道 / 接收跳过 {d["bridge_or_receiver_skipped"]}、视觉最新槽覆盖 {d["vision_latest_overwrites_in_processing_span"]}。录制丢帧 {d["recorder_drops"]}。Mock 收到 {d["mock_received"]} 条，均 NO_VALID_TARGET；其中过期 / 无图拒绝 {d["mock_stale_rejected"]} 条，包含启动心跳，不与处理后过期数混用。']
 lines+=['','## 解释边界','', 'OpenCV 的 ROI 平均吞吐和 p95 达到 60 fps 参考，但最大耗时超过 16.7 ms，存在跳帧，不能称每一帧硬实时。Cellpose 不通过 60 fps 参考。新鲜候选只是按软件时钟仍可用的图像区域，专家已验证微泡 / 团更新率仍均为 0。','', 'decoder 只暴露完成解码的时间戳，单独每帧解码耗时 NOT_INSTRUMENTED；没有伪造曝光到解码延迟。软件新鲜度使用 sender_age + WSL GET 往返 + WSL 后续耗时的保守上界，不直接相减异机单调时钟。GET 包含等待下一帧，不是纯复制开销。','', 'EXPOSURE_TIMESTAMP_UNKNOWN；CAMERA_TO_CONTROLLER_END_TO_END_LATENCY=NOT_MEASURED。没有 LED 测试，也没有实际 PZT / AFG 控制。']
 write('LIVE_INFERENCE_BENCHMARK.md','\n'.join(lines))
 write('GROUP_AND_BUBBLE_ADMISSION.md','''# 单泡与团级准入

| 层次 | 当前状态 | 依据与限制 |
|---|---|---|
| 实时相机入口 | LIVE_CAMERA_PASS | 原 Windows 接收器与 WSL 有界桥，实际两条20秒试验 |
| 人工工具 | ANNOTATION_TOOL_PASS | 修复后实际 Qt 多对象、保存恢复、换帧、SAM 后台与取消 |
| 专家标签 | EXPERT_GROUND_TRUTH_PENDING | 详见FINAL_GATE.json的带时间快照；验证/测试标签尚缺 |
| 单泡识别 | NOT_EVALUATED | score 未校准；面积分类不是语义真值 |
| 短轨迹 | NOT_EVALUATED | Norfair软件关联可用；无人工连续ID真值，不声称改善断裂 |
| 团轮廓 / 形状 | GROUP_SHAPE_NOT_ADMITTED | group_shape=null；原始图像区域形态单独保存，不进正式接口 |
| 团附近单泡 | NOT_ADMITTED | 有像素邻域几何辅助函数；没有人工准入团，因此正式输出null |
| 并入 / 离开事件 | NOT_EVALUATED | 当前不产生已验证事件；未来只能称APPARENT_INCORPORATION_CANDIDATE，不能推断物理融合 |
| 软件实时 | OpenCV ROI throughput PASS | 并非所有帧低于16.7ms，也非身份准确率PASS |
| 曝光端到端延迟 | NOT_MEASURED | 无LED测量；曝光时间未知 |
| 硬件控制 | DISABLED | 仅本地ZeroMQ mock；无设备控制API |

单泡位置采用当前观测；丢失不输出预测旧坐标。速度为最近观测的像素/源时间秒，直播源时间为host decode时间，不能解释为标定物理速度。LK为局部图像运动而非液体流速。跨session、逆向时间、过大帧间隔重置轨迹。原始mask的质心、面积、Rg、Q和长短比通过独立矩形软件检查，但形状公式正确不证明其输入是真实泡团。

当前训练和科学控制状态：training_ready=false。后续需冻结完整人工标签、用未调参测试组验证单泡和团分割、用短序列检查真实关联，才考虑解除相应语义门控。本轮不会自动将模型预测转成专家标签。
''')
 answers=[
 ('是否复用原接收器','是，F:\\CHRONOS_control 的 LatestFrameStream，未重写RTSP。'),
 ('相机是否正常输入','是，两条正式20秒LIVE_CAMERA试验均完成；接收失败数均0。'),
 ('标注界面能否打开','修复后已实际打开并通过多边形、保存重开、换帧、后台SAM和取消检查。'),
 ('怎样标注','按ANNOTATION_READY.md第2节：类别→新建轮廓→沿边界点击→完成轮廓→保存；整帧检查后才确认。'),
 ('候选多少帧','84：72录像+12实时候选。'),('正式确认多少帧',f'{audit["expert_confirmed_frames"]}，另有{audit["draft_frames"]}张草稿；软件测试不算。'),
 ('是否避免相邻帧泄漏','是，录像按独立时间块划分并留间隙；直播隔离，跨实验泛化未验证。'),
 ('OpenCV表现','ROI速度较快，但选出很多小对比区域；真实微泡Precision/Recall尚未测得。'),
 ('Cellpose表现','cyto3实际运行；候选少且常偏向较大模糊结构；直播ROI无可用候选，速度不满足60fps。'),
 ('是否微调','没有，缺少正式专家标签，不能用自动输出冒充训练真值。'),
 ('谁漏检少','NOT_EVALUATED；候选数量不能当召回率。'),('谁背景误检少','NOT_EVALUATED；需要完整负例和忽略区标签。'),
 ('谁适合60fps','当前保留OpenCV作实时候选管线；仅通过此ROI的吞吐参考，不是识别准入。'),
 ('Norfair断裂是否改善','尚无法评价；短ID真值未完成，软件ID数量不是微泡数。'),
 ('突然转向是否可跟踪','算法未固定运动方向，但真实转向准确率尚未验证。'),
 ('团中心与形状有人工验证吗','没有，正式group_shape=null。'),
 ('团附近单泡能可靠输出吗','尚未准入；像素几何函数可用，但缺少可靠单泡和团掩码。'),
 ('是否有并入事件','未评价；不把目标消失或面积增大当作融合，事件保持不确定。'),
 ('处理p95多少',f'真实在线完整软件：OpenCV {A["stage_ms"]["software_ms"]["p95"]:.3f} ms，Cellpose {B["stage_ms"]["software_ms"]["p95"]:.3f} ms。'),
 ('有效视觉更新率',f'处理率{A["processed_hz"]:.3f}/{B["processed_hz"]:.3f} Hz；新鲜候选{A["fresh_candidate_hz"]:.3f}/{B["fresh_candidate_hz"]:.3f} Hz；已验证微泡均0 Hz。'),
 ('是否未做LED','是，LED_TEST_NOT_RUN。'),('端到端延迟','继续NOT_MEASURED；曝光时间未知。'),
 ('是否确保无PZT自动控制','是，仅mock，无真实AFG/PZT命令。'),
 ('四组科研图','figures/Fig_VisionV2_1_detection_and_labels、2_online_tracking、3_accuracy、4_live_latency，均PNG/PDF。'),
 ('视频及播放验证','videos/Video_VisionV2_live_tracking.mp4：1202帧、约20.118秒；method_comparison.mp4：90帧、约2.997秒。均1280×960 H264，逐帧解码与单调时间戳检查通过，并抽帧目视检查。'),
 ('哪些指标待标注','PR/F1、FP/FN、定位、团Dice/IoU和形态误差、真实短ID关联、转向/失焦表现、团邻域和并入候选。'),
 ('推荐方案及理由','先保留OpenCV+Norfair实时候选管线；人工标注+原SAM辅助；Cellpose作为待验证离线备选，不强行逐帧闭环。'),
 ('下一轮最小工作','先确认少量完整train/validation/test帧与8帧短轨迹，冻结修订，定量定位漏检/误检；如有证据再进行一次有界微调。补清晰且独立实验日样本。')]
 write('FINAL_RECOMMENDATION.md','# F10T-Vision-2 中文交付报告\n\n标注工作台已修复，后续工程验证已继续完成。正式科学识别评价等待人工标签；没有宣称闭环控制完成。\n\n'+'\n\n'.join(f'{i}. **{q}：**{answer}' for i,(q,answer) in enumerate(answers,1))+f'''\n\n两条完成的实时测量窗口合计约 {A['wall_measured_s']+B['wall_measured_s']:.2f} s；成对录像计算 {P['wall_s']:.2f} s。候选采集、冷启动、GUI检查、安装下载、开发绘图另计；未保存的中断进程耗时未知，不虚构“全部研发耗时”。

分层机器记录：[FINAL_GATE.json](FINAL_GATE.json)。详细操作：[ANNOTATION_READY.md](ANNOTATION_READY.md)。性能：[LIVE_INFERENCE_BENCHMARK.md](LIVE_INFERENCE_BENCHMARK.md)。结果和脚本 SHA见 evidence/delivery_manifest.json。
''')
 write('MODEL_AND_ENVIRONMENT.md','''# 模型、来源与环境

- Vision：Python3.12独立venv安装层，只读引用上级研究venv；OpenCV-headless4.10.0.84、Norfair2.3.0、Cellpose3.1.1.3、Torch2.11.0+cu128。有效锁定见evidence/requirements-vision-effective.lock。不修改Windows Python3.14环境。
- Annotation：独立Python3.12 venv，micro-sam1.8.14、napari0.9.2、PyQt6，NumPy2.2.6。独立锁定见requirements-annotation.lock，避免与Norfair的NumPy1.x要求冲突。两环境pip check通过。
- 原SAM ViT-B仅辅助标注。官方权重URL：https://dl.fbaipublicfiles.com/segment_anything/sam_vit_b_01ec64.pth ，SHA-256 `ec2df62732614e57411cdcf32a23ffdf28910380d03139ee0f4fcbe91eb8c912`。原SAM官方代码与权重Apache2.0；SA-1B数据有单独研究条款，本轮未下载训练数据。实际采用micro-sam加载，不安装/运行SAM2。
- Cellpose唯一模型cyto3，官方下载 https://www.cellpose.org/models/cyto3 ，SHA-256 `2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c`。3.1.1.3代码LICENSE为BSD3条件；该版官方README明确cyto3等模型训练涉及CC-BY-NC数据。不能从代码BSD许可推断所有训练数据或商业用途无额外条件。未下载Cellpose训练数据，未将用户视频上传第三方。许可原文快照和URL/SHA保留在evidence/licenses。
- micro-sam MIT，napari BSD3，Norfair BSD3，OpenCV Apache2.0；具体依赖还有各自许可。官方项目：https://github.com/computational-cell-analytics/micro-sam ，https://github.com/napari/napari ，https://github.com/tryolabs/norfair ，https://github.com/opencv/opencv 。micro-sam源克隆只用于核查，实际运行来自锁定的PyPI包，不能混淆两者版本。
- cyto3选择依据：[Cellpose3官方模型文档](https://cellpose.readthedocs.io/en/v3.1.1.1/models.html)。仅运行一个轻量替代模型，没有竞赛式扫描。

权重、环境和外部仓库未加入Git。中文界面使用用户Windows系统已有微软雅黑字体，通过用户字体目录符号链接供WSLg显示，没有改动设备工程。
''')
 write('ACTIVE_CODE_INDEX.md','''# F10T-Vision-2 活跃入口

| 入口 | 用途 |
|---|---|
| start_annotation.sh / scripts/annotate.py | 已修复的持续窗口、napari原生绘图、类别与版本化保存 |
| scripts/sam_worker.py | 可取消原SAM后台辅助；不含SAM2 |
| src/annotation_store.py | 五类掩码、确认关口、SHA、Cellpose导出 |
| scripts/prepare_candidates.py / acquire_live_candidates.py | 预登记录像候选、独立直播采集；拒绝覆盖 |
| scripts/review_labels.py | 核查人工标签；完成后冻结具体修订和导出 |
| scripts/windows_frame_server.py / src/live_bridge.py | 原Windows接收器的有界原始图像IPC |
| src/vision_engine.py | 两条候选检测路线、Norfair、图像矩、稀疏LK |
| src/runtime.py / scripts/run_live.py | 最新帧、过期拒绝、mock、资源与真实直播日志 |
| scripts/compare_routes.py | 固定同帧同ROI成对测试、严格因果离线短回放 |
| scripts/render_delivery.py | 四组真实PNG/PDF和两段MP4；--figures-only只更新图 |
| tests/test_contracts.py | 标注I/O、冲突、导出、形态、时间/session、mock契约 |
| scripts/write_delivery_reports.py | 从真实日志生成交付文档与分层关口 |
| scripts/verify_delivery.py | 最终文件身份、视频、标签保护及交付检查 |

本地分支dev/f10t-vision2-20261016；基线9c22d207be5318f9d45b014caaf52bd04d345328。原相机 / AFG工程不修改，历史未清理。旧版GUI仅保留evidence/annotation_repair_20261016作故障证据，不再是启动入口。
''')
 print('Reports written from actual logs; expert labels remain pending.')
if __name__=='__main__':main()
