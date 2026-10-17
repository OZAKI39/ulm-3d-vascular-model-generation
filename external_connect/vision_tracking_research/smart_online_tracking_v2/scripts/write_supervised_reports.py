"""Publish the finite, frozen-ROI study without erasing earlier deliveries."""
from pathlib import Path
from datetime import datetime, timezone
import json, sys, shutil, difflib, re

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from annotation_store import sha

TAG='supervised_20261016_v1'
OUT=ROOT/'results'/TAG


def main():
    summary=json.loads((OUT/'summary.json').read_text())
    paced=json.loads((OUT/'paced_replay/summary.json').read_text())
    replay=json.loads((OUT/'replay/summary.json').read_text())
    training=json.loads((OUT/'training.json').read_text())
    before=OUT/'documents_before';before.mkdir(exist_ok=False)
    names=['SMART_TRACKING_V2_DEVELOPMENT.md','ANNOTATION_READY.md','MICROBUBBLE_ANNOTATION_PROTOCOL.md',
           'DATASET_SPLIT_AND_PROVENANCE.md','OPENCV_VS_CELLPOSE.md','LIVE_INFERENCE_BENCHMARK.md',
           'GROUP_AND_BUBBLE_ADMISSION.md','MODEL_AND_ENVIRONMENT.md','FINAL_RECOMMENDATION.md','FINAL_GATE.json']
    for name in names:shutil.copy2(ROOT/name,before/name)
    models=['opencv','cellpose_pretrained','cellpose_tuned']
    table='| 方法 | TP / FP / FN | Precision | Recall | 中心 F1 | 掩码 F1 (IoU≥0.5) | 中心误差均值 / p95 (px) |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for model,label in zip(models,['OpenCV','cyto3 预训练','cyto3 微调（拒绝采用）']):
        m=summary['metrics'][f'{model}/validation/ROI'];c=m['center']
        table+=f'| {label} | {c["tp"]} / {c["fp"]} / {c["fn"]} | {c["precision"]:.4f} | {c["recall"]:.4f} | {c["f1"]:.4f} | {m["mask"]["f1"]:.4f} | {m["center_error_px_mean"]:.3f} / {m["center_error_px_p95"]:.3f} |\n'
    performance='| 因果回放方法 | 处理帧 / 输入帧 | 处理 Hz | 新鲜候选 Hz | 软件 p50 / p95 / 最大 (ms) | 最新帧覆盖 | 结果过期 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for route,label in [('opencv','OpenCV'),('cellpose','cyto3 预训练')]:
        v=paced[route];t=v['stage_ms']['software_ms']
        performance+=f'| {label} | {v["processed_frames"]} / {v["source_frames"]} | {v["processed_hz"]:.3f} | {v["fresh_candidate_hz"]:.3f} | {t["p50"]:.3f} / {t["p95"]:.3f} / {t["maximum"]:.3f} | {v["latest_slot_overwrites"]} | {v["stale_result_rejected"]} |\n'
    gate=dict(utc=datetime.now(timezone.utc).isoformat(),status='OFFLINE_STUDY_COMPLETE; LIVE_CAMERA_AND_SCIENTIFIC_ADMISSION_PENDING',
        result_directory=str(OUT),LIVE_CAMERA_PASS=False,live_reason='USER_CONFIRMED_CAMERA_DISCONNECTED; 0 frames in both attempted runs',
        historical_live_pass_preserved=True,ANNOTATION_TOOL_PASS=True,
        EXPERT_LABELS_READY=True,expert_labels_scope='REVIEWED_ROIS_ONLY: 2 train frames / 2 validation frames; no test or track/group truth',
        expert_reviewed_source_frames=4,expert_confirmed_full_frames=0,candidate_frames=84,
        train_reviewed_rois=5,train_single_masks=210,validation_reviewed_rois=5,validation_single_masks=133,
        validation_eligible_single_masks=129,EXPERT_GROUND_TRUTH_PENDING=True,
        BUBBLE_DETECTION_VALIDATED=False,SHORT_TRACKING_VALIDATED=False,GROUP_SEGMENTATION_VALIDATED=False,
        GROUP_SHAPE_NOT_ADMITTED=True,group_shape=None,
        REALTIME_PROCESSING_PASS={'opencv':'PILOT_NATIVE_30FPS_CAUSAL_REPLAY_SMALL_ROI_ONLY; LIVE_60FPS_NOT_RETESTED',
                                  'cellpose':'FAIL_100MS_FRESHNESS_IN_PACED_REPLAY'},
        strict_every_frame_16_7ms_pass=False,validated_bubble_update_hz=0,validated_group_update_hz=0,
        CROSS_SESSION_GENERALIZATION_NOT_VERIFIED=True,END_TO_END_LATENCY_NOT_MEASURED=True,
        EXPOSURE_TIMESTAMP_UNKNOWN=True,CAMERA_TO_CONTROLLER_END_TO_END_LATENCY='NOT_MEASURED',
        HARDWARE_CONTROL_DISABLED=True,CLOSED_LOOP_CONTROL_NOT_RUN=True,LED_TEST_NOT_RUN=True,SAM2_NOT_RUN=True,
        finetuning='ONE_CYTO3_100_EPOCH_RUN_COMPLETED; REJECTED_BY_PREREGISTERED_VALIDATION_RULE',
        selected_cellpose='PRETRAINED_CYTO3',training_ready=False,
        accuracy_metrics={k:v for k,v in summary['metrics'].items() if '/validation/' in k},
        failed_or_limited=['Fine-tuned model loses recall','Full-frame Cellpose inference loses recall',
                           'No camera connected','Cellpose results stale in paced replay','No human group/track/test labels'])
    (OUT/'FINAL_GATE.json').write_text(json.dumps(gate,ensure_ascii=False,indent=2))
    (ROOT/'FINAL_GATE.json').write_text(json.dumps(gate,ensure_ascii=False,indent=2))
    report=f'''# F10T-Vision-2：基于现有人工区域标注的离线续研

本轮完成了冻结真值、两条检测路线的实际评价、一次有限 cyto3 微调、同帧因果回放、最新帧软件链路测量、四组科研图和两段真实 MP4。相机未连接，用户要求先完成离线部分；未继续尝试接流。**离线交付完成不等于单泡识别、团形态或控制准入通过。**

## 1. 真实数据与评价边界

固定快照：`results/expert_freezes/vision2_supervised_20261016_v1/`。训练为第6、7张（原视频帧90、108），5个已审核区域、210个单泡；验证为用户新完成的第25、26张（原视频帧360、374），5个区域、133个单泡。没有整帧完整审核，也没有测试、团、短轨迹ID的人工真值。原第1–5张的部分标签仍是草稿，全部保留，不参与本轮训练或计分。

验证133个标签中4个因6 px边界中心范围或裁切规则排除，最终129个可评价单泡。先求最大有效匹配数，再在有效匹配内优化距离／IoU；中心距离≤6 px与掩码IoU≥0.5分别计分。边界、裁断和忽略区的排除逐ROI记录；没有不确定／失焦标签时拒判准确率为NOT_EVALUATED，不写成0%错误。5个验证ROI来自**2张源图**，不能算5个独立样本。未审核的框外不视作负例。

分组仍为train源帧0–270、validation360–570、test660–936，保持时间间隙和连续片段分组。本轮没有用test调参；验证参与最终候选选择，所以不能把验证成绩称独立测试结果。同一录像、少数区域仍有场景偏倚，不能证明跨实验日泛化。人工标签包含模型辅助后人工修正的来源记录，不是盲法、双人独立标注；需保留这一可能偏倚。

## 2. 检测比较和失败结果

{table}
OpenCV在这些区域漏检10个，cyto3预训练漏检52个；FP分别35与33，相差很小且没有独立统计检验。FP表示已审核区域内未匹配候选，不能进一步断言全是某类背景伪影。微调模型FP为0，但漏检114个、召回率仅11.6%，不能称性能更好。OpenCV中心F1较高，预训练cyto3掩码F1较高；定位和边界恢复是不同任务。

额外使用完整1476×1120原图推理，仍**只在同一已审核ROI评分**：OpenCV TP119/FP35/FN10，预训练cyto3 TP1/FP0/FN128，微调cyto3 TP0/FP0/FN129。全图归一化与局部裁剪的输入分布可能不同，这是需要验证的解释，尚未完成归一化消融，不能直接认定唯一原因。更不能把局部表现推广成整帧准确率。

## 3. 有界训练与模型选择

OpenCV仅在训练ROI比较对比阈值6/8/10，选择10。cyto3直径参考由训练单泡面积等效直径中位数得到11.113245963 px，min_size=5，channels=[0,0]；未逐泡强制同径。仅使用既有cyto3 3.1.1.3，不引入第三条自动检测路线。

训练基线中心F1<0.90或掩码F1<0.80触发一次训练。复用官方train_seg，5个训练裁片、100 epochs、batch2、SGD lr0.1、weight_decay1e-4、bsize224、种子20261016、180秒预算；没有给优化器验证或测试数据。只比较最终epoch，未挑最好epoch。预登记采用条件为验证中心F1至少增加0.02且FP不多于基线+2；实际不满足，继续选用原预训练模型。

实际训练 {training['wall_s']:.3f} s，完整静态试验 {summary['total_wall_s']:.3f} s（训练包含其中）。损失、失败指标、检查点保留在`training.json`、`all_rows.json`和`models/`；失败检查点SHA-256为`{training['weights_sha256']}`。禁止把该检查点悄悄切换为在线模型。

后台原SAM提示分割是另一个用途：既有自动服务在2张验证图、64个有提示单泡上的候选平均IoU约0.536，原版约0.432，按自身门槛允许辅助掩码。它仍需要逐泡提示及人工确认，不能当作自动检测性能，也没有改变本轮Cellpose拒绝采用结论。没有运行SAM2。

## 4. 实际回放、软件性能和相机状态

原Windows Chronos接收器保持不变。本次两条接流尝试均取得0帧，保留原日志；用户确认相机未连接。旧20秒直播成功记录是历史，不作为本次LIVE_CAMERA_PASS依据。

第一组严格顺序回放：源帧360–419，550×560原像素ROI，每路60帧，共处理 {replay['wall_s']:.3f} s。OpenCV engine p95 {replay['routes']['opencv']['timing_ms']['engine_ms']['p95']:.3f} ms，cyto3 p95 {replay['routes']['cellpose']['timing_ms']['engine_ms']['p95']:.3f} ms。只用于离线检测/短时关联展示，不能把按源PTS播放的2秒视频当作实时吞吐证明。

第二组按源PTS约30.030 fps释放帧，源帧360–659，原有180×350像素ROI，单槽最新帧、100 ms新鲜度上限，逐路运行约10秒。包含推理、Norfair、形态、局部图像运动和序列化发送；mock消费与传输另外测量，**不把分段p95相加作为端到端p95**。

{performance}
OpenCV最大软件时间28.028 ms，因此即使p95低于16.7 ms，也不是逐帧60 Hz硬实时。Cellpose全部54个结果超过100 ms，mock中没有新鲜候选；246个源帧被最新帧机制覆盖，不积压旧画面。两路经过验证的生物单泡和团更新率均为0。约30.079 Hz的计数除墙钟估计略高于30.030源帧率，是首末帧/短窗口端点口径差异，不代表生成额外帧。

OpenCV mock传输p95 {paced['opencv']['mock_transport_ms']['p95']:.3f} ms、消费p95 {paced['opencv']['mock_consume_ms']['p95']:.3f} ms；cyto3对应 {paced['cellpose']['mock_transport_ms']['p95']:.3f}/{paced['cellpose']['mock_consume_ms']['p95']:.3f} ms。CPU/RAM和GPU日志均保存。GPU为设备整体采样，标注工作台仍运行，并非独占性能测试。模型warmup另计；没有更改相机固件、系统网络或原设备工程。

Norfair仍使用短时观测关联，失去观测时坐标为空，session或大时间间隙重置。没有使用未来帧或控制仿真先验。没有人工轨迹真值，软件ID数量和LOST观测次数不能转化成真实ID断裂率，也不能证明急转弯跟踪改善。

## 5. 科研图、视频与复现

- 四组PNG/PDF：`figures/{TAG}/Fig_VisionV2_1_detection_and_labels.*`、`Fig_VisionV2_2_online_tracking.*`、`Fig_VisionV2_3_accuracy.*`、`Fig_VisionV2_4_live_latency.*`。后者明确展示当前直播不可用与离线计时。
- 两段MP4：`videos/{TAG}/Video_VisionV2_CAUSAL_REPLAY_tracking.mp4`、`Video_VisionV2_CAUSAL_REPLAY_method_comparison.mp4`。1280×960、H.264、每段60帧约1.998秒，保留原PTS间距、不插帧、不利用未来帧。人工轮廓仅在真实源帧360和374的已审核区域出现，其余帧明确未审核。
- 绘图脚本：`scripts/render_supervised_delivery.py`。原始输出：`predictions/`、`replay/records.jsonl`、`paced_replay/*`。输入、源码、输出SHA与视频逐帧解码检查见`render_identity.json`、`video_validation.json`及最终`delivery_manifest.json`。
- `scripts/run_supervised_study.py`、`replay_supervised.py`、`benchmark_supervised_replay.py`保留有限实验入口；正式输出目录禁止覆盖，再跑需选择新TAG和新冻结数据。当前GUI与标注原文件没有被重启或改写。
- 检查：评价匹配/忽略区7项、原消息/跟踪契约8项通过。ffprobe本地命令不可用，采用已安装PyAV/FFmpeg库读取元数据并逐帧完整解码；记录该检查工具替代，不伪称运行ffprobe。

## 6. 准入与最小下一步

当前推荐先完善OpenCV＋Norfair小ROI候选路线，cyto3保留为离线轮廓和人工标注参考。不要仅因训练损失下降就替换模型；下一次开发优先固定归一化/ROI策略、保留未知区域、扩大少量独立人工审核样本，再做有界比较。

用户无需重标第25、26张。最小补充是第49、50张（test源帧660、678）各1–2个小区域，框内逐一补漏、纠错并“确认当前区域”，框外可继续不标；测试图不进入后台训练。随后若要验证短轨迹，选择第56–63张连续8帧（780–787）的同一小区域，标注可确认的短ID，模糊/遮挡不猜。团级评价还需要真正的团掩码与不确定区域，不把粘连暗斑或图连通直接称物理融合。

本次无需等待这些新标签才交付离线工作。相机接回后再复测固定ROI的实际输入、fresh-candidate更新率和软件p95。LED延迟、SAM2、AFG/PZT自动控制均未运行，CAMERA_TO_CONTROLLER_END_TO_END_LATENCY=NOT_MEASURED，group_shape=null，training_ready=false。

## 7. 原任务28项回答

| 序号 | 回答 |
|---:|---|
| 1 | 复用原Windows Chronos接收器，没有重写RTSP。 |
| 2 | 本次无实时输入；用户确认相机未连接，后续按要求只做离线。 |
| 3 | 既有标注界面可用且本轮保持运行；未强制关闭用户窗口。 |
| 4 | 框选小区域→批量候选→删错补漏→确认当前区域；N新增泡，P画笔，E标签橡皮，R候选橡皮。详见ANNOTATION_READY.md。 |
| 5 | 84张候选：72张录像、12张历史实时隔离帧。 |
| 6 | 本次冻结4张源图中的10个已审核区域；整帧审核数0。训练210、验证133个单泡。 |
| 7 | 保留时间块/连续片段分组；没有随机拆相邻帧，没有使用test训练。 |
| 8 | OpenCV验证中心P0.773/R0.922/F1 0.841，TP119/FP35/FN10。 |
| 9 | cyto3预训练P0.700/R0.597/F1 0.644，TP77/FP33/FN52；掩码F1 0.469优于OpenCV 0.254。 |
| 10 | 完成一次100 epoch微调；训练基线未达预登记触发阈值。微调验证退步，拒绝采用。 |
| 11 | 当前已审核ROI内OpenCV漏检最少；不能推广到全部图或其他实验。 |
| 12 | 预训练cyto3 FP33略少于OpenCV35，缺少背景分类真值；微调FP0伴随FN114，不能称更好。 |
| 13 | OpenCV更有在线使用潜力；本次30 fps离线回放不构成60 fps直播通过。 |
| 14 | 没有人工短ID，ID断裂改善NOT_EVALUATED。 |
| 15 | 没有人工急转弯片段真值，不宣称已验证；不以固定方向约束轨迹。 |
| 16 | 团中心和形状尚无人工验证，group_shape=null。 |
| 17 | 不能可靠输出科学团邻域分布；候选几何不能替代准入。 |
| 18 | 本轮没有确认的并入事件；事件输出为空，不把消失/遮挡解释为融合。 |
| 19 | 小ROI完整软件p95：OpenCV9.958 ms、cyto3 206.005 ms；大ROI仅engine p95：46.736/1058.136 ms。 |
| 20 | 小ROI离线fresh候选约30.079/0 Hz；已验证单泡/团均0 Hz；本次实时更新率未测得。 |
| 21 | 没有LED延迟测量。 |
| 22 | 相机端到端延迟继续NOT_MEASURED。 |
| 23 | 硬件控制禁用，仅本机mock消费者，没有PZT自动控制。 |
| 24 | 四组英文PNG/PDF在figures/{TAG}/。 |
| 25 | 两段真实回放MP4在videos/{TAG}/；逐帧解码、时间戳与抽帧检查留档。 |
| 26 | 独立test、跨实验日、团、短ID、急转弯及不确定/失焦拒判准确率仍待人工真值。 |
| 27 | 优先OpenCV＋Norfair小ROI：当前漏检较少、处理更快；中心与轮廓精度分开验收。 |
| 28 | 先在第49/50张各审核1–2小块；按目标补8帧短ID或团掩码，再接回相机复测。 |

## 参考与来源

- [Cellpose 3官方训练说明](https://cellpose.readthedocs.io/en/v3.1.1.1/train.html)：本轮复用train_seg及推荐SGD设置；保留少样本失败，不据模型名称推断准确性。
- [Norfair官方代码](https://github.com/tryolabs/norfair)：沿用成熟短时关联。
- 模型来源、BSD代码与cyto3训练数据NC条款区别、环境锁定见MODEL_AND_ENVIRONMENT.md。本轮未将录像或标注上传第三方。
'''
    (OUT/'REPORT.md').write_text(report)
    sections={
      'SMART_TRACKING_V2_DEVELOPMENT.md': '已基于第25、26张的人工审核区域完成冻结评价、一次有限微调、两类因果回放及新图表/视频。相机未连接，按用户要求完成离线部分。新增评价模块与测试；原相机控制、标注原数据和历史结果保留。',
      'ANNOTATION_READY.md': '第25、26张已用于本次冻结评价，无需为了继续开发而重新标注。当前训练为第6、7张，验证为第25、26张，各5个完整审核的小区域；没有把整张图当背景。\n\n下一步最省工作：在“更多 / 训练”跳到**第49张、第50张**，各框1–2个小区域，删错、补漏后点“确认当前区域”；框外不必标。它们是test，不会用于后台训练。若要验证轨迹，再审核第56–63张同一小区域的8帧短ID。没有可靠团掩码前，团形态仍不准入。\n\n本轮独立试跑的Cellpose微调效果变差，没有替换批量标注模型；后台SAM提示分割的微调状态是另一回事。',
      'MICROBUBBLE_ANNOTATION_PROTOCOL.md': '本轮只评价冻结的10个已审核ROI。验证133个单泡中4个边界/裁断排除，129个计分；中心6 px匹配与mask IoU≥0.5分开。一对一先最大化有效匹配数；排除项逐ROI保存。不确定类样本不足，不报拒判准确率。未经审核的外部区域不作为背景。',
      'DATASET_SPLIT_AND_PROVENANCE.md': '新冻结集vision2_supervised_20261016_v1含4张源图：train90/108（210泡）、validation360/374（133泡），均为区域审核。test未使用。5个验证ROI只算2张源图；模型辅助后人工审核，非盲法双人真值。原第1–5张仅部分标签，不作为全图真值。',
      'OPENCV_VS_CELLPOSE.md': table+'\n只在5个已审核验证ROI计分，129个有效真值。一次cyto3微调被拒绝采用；全图推理时cyto3在这些ROI仅检出1/129，失败保留。没有整帧准确率或独立测试通过。',
      'LIVE_INFERENCE_BENCHMARK.md': '本次相机未连接，两个接流尝试0帧，LIVE_CAMERA_PASS=false。以下均为约30.030 fps录像因果回放，不是历史60 fps实时流。\n\n'+performance+'\n小ROI为180×350；软件含发送但不把独立mock传输p95相加。大ROI550×560的engine p95为46.736/1058.136 ms。无曝光时间戳，无LED，本次端到端延迟NOT_MEASURED。',
      'GROUP_AND_BUBBLE_ADMISSION.md': 'EXPERT_LABELS_READY仅适用于2 train+2 validation源图的局部单泡掩码。BUBBLE_DETECTION_VALIDATED、SHORT_TRACKING_VALIDATED、GROUP_SEGMENTATION_VALIDATED仍false；没有test、团或短ID真值。group_shape=null，control_permission=false，training_ready=false。',
      'MODEL_AND_ENVIRONMENT.md': f'本轮复用既有cyto3 3.1.1.3及环境，未安装新模型家族。一次100 epoch SGD微调耗时{training["wall_s"]:.3f} s，验证失败，仍选预训练。拒绝检查点SHA `{training["weights_sha256"]}`；权重只在本地results内，未入Git。后台原SAM提示掩码与自动检测分开评价。',
      'FINAL_RECOMMENDATION.md': '优先推进OpenCV＋Norfair的小ROI候选路线：本次验证漏检10/129，小ROI软件p95约9.96 ms；但有35个FP、轮廓F1仅0.254，尚不科学准入。预训练cyto3轮廓较好、漏检较多且太慢；微调候选拒绝。下一步补两张test的小区域，检查归一化/裁剪一致性，再复测真实相机；不进行更多盲目训练或硬件控制。'
    }
    for name,text in sections.items():
        original=(before/name).read_text()
        first,rest=original.split('\n',1)
        notice=f'\n\n> **最新续研：2026-10-16，已完成离线人工ROI评价。** 以下较早记录按历史保留；当前状态以文末续研和[完整报告](results/{TAG}/REPORT.md)为准。\n'
        appendix=f'\n\n## 2026-10-16 已有标注的离线续研（{TAG}）\n\n{text}\n\n完整参数、失败证据、28项回答和输出路径见[本轮报告](results/{TAG}/REPORT.md)。\n'
        (ROOT/name).write_text(first+notice+rest+appendix)
    # Minimal global progress insert, with heading and old-paragraph preservation checks.
    global_path=Path('/home/lzy/projects/microbubble_aggregation.md')
    old=global_path.read_text()
    backup=global_path.with_name(global_path.name+'.bak_F10T_Vision2_20261016_supervised')
    if backup.exists():raise FileExistsError(backup)
    shutil.copy2(global_path,backup)
    block='''**2026年10月16日的离线续研已用上新完成的验证标注。**第1–5张早期图片只有部分泡、没有检查全图，已保留标签并撤销全图审核。本次冻结第6、7张的5个已审核区域（210个单泡）作训练，第25、26张的5个区域（133个单泡）作验证；框外不当背景。4个验证对象按边界/裁断规则排除，实际129个计分，5个区域只算2张源图。OpenCV中心匹配TP/FP/FN为119/35/10，P/R/F1为0.773/0.922/0.841；预训练cyto3为77/33/52、0.700/0.597/0.644，但mask IoU≥0.5的F1较高（0.469对0.254）。这只是同录像局部试验，没有独立test或跨实验日通过。

**一次有限Cellpose微调失败，结果没有被替换或美化。**仅用5个训练裁片完成100 epoch、约15.61秒，整个静态比较约73.08秒；微调召回降到0.116，漏检114个，按预登记规则拒绝采用，保留权重与损失。完整原图推理时预训练cyto3在同一验证ROI只检出1/129，说明局部结果不能推广到全图；归一化差异只是待验证解释。原SAM提示掩码的自动微调是人工辅助用途，不能与自动检测混为一谈。

**本次相机未连接，因此继续的是实际录像因果回放。**两个接流尝试0帧，用户确认设备未连接后不再重试；旧实时成功记录仍保留。550×560像素ROI的60帧顺序回放，OpenCV/cyto3 engine p95约46.74/1058.14 ms。原180×350小ROI约30.03 fps、300帧最新帧回放中，OpenCV处理300帧，软件p95约9.96 ms、最大28.03 ms、新鲜候选约30.08 Hz；cyto3处理54帧、覆盖246帧，p95约206.01 ms且全部结果过期，新鲜候选0 Hz。不能把录像播放速度或p95当逐帧60 Hz实时通过。两路已验证单泡/团更新率仍0，短ID断裂、急转弯、团掩码和团邻域尚无人工真值，group_shape=null。

新增四组英文PNG/PDF和两段标记CAUSAL_REPLAY的H.264视频，使用原PTS、不插帧，保留输入/输出SHA、源代码身份、完整逐帧解码与视觉检查；详见smart_online_tracking_v2/results/supervised_20261016_v1/REPORT.md。科研进展是完成了可复现的局部精度评价及失败训练，尚未通过科学或控制准入。下一步只需在第49、50张各审核1–2个小区域形成独立test，再按目标补连续8帧短ID或团掩码，接回相机后复测。未做LED、SAM2或PZT/AFG控制，端到端延迟NOT_MEASURED，training_ready=false。

'''
    marker='---\n\n### 尾声：故事里还没讲完的部分'
    assert old.count(marker)==1
    new=old.replace(marker,block+marker)
    headings=lambda s:re.findall(r'^#{1,6} .+$',s,re.M)
    assert headings(new)==headings(old) and new.replace(block,'')==old
    global_path.write_text(new)
    evidence=OUT/'progress_document';evidence.mkdir()
    (evidence/'before.md').write_text(old);(evidence/'after.md').write_text(new)
    (evidence/'diff.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before',tofile='after')))
    (evidence/'verification.json').write_text(json.dumps(dict(backup=str(backup),before_sha256=sha(backup),after_sha256=sha(global_path),
         heading_structure_unchanged=True,all_old_text_preserved=True,heading_count=len(headings(old))),indent=2))
    (OUT/'document_diff.patch').write_text(''.join(''.join(difflib.unified_diff((before/n).read_text().splitlines(True),(ROOT/n).read_text().splitlines(True),fromfile='before/'+n,tofile=n)) for n in names))
    print('Reports, current gate and minimal progress update saved; historical documents backed up.')


if __name__=='__main__':main()
