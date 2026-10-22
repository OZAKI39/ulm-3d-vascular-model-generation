# s1-8 显式支撑：使用与输出

> 本文以下内容是 **legacy_explicit_support_v1 历史说明**，V1 已退出生产选择。当前入口默认运行 V2，打印姿态可改变，最新输出在 `vascular_support_planner/v2/`。请使用 [FDM V2 使用说明](FDM_SUPPORT_V2使用说明.md)，不要按本文旧的“三候选/27 格试片/固定打印姿态”解释新结果。

继续使用原 s1-8 入口，不存在新的 s1-9 阶段。正式方向为 `CUSTOM_LITERATURE_AWARE_SUPPORT`：检测必须支撑的区间，筛选下表面接触位置，建立与修复显式支撑图，生成可折断底脚，再关闭 Bambu 自动支撑进行真实切片。血管、半径、端口、2 mm 盒体和姿态保持不变。

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-8_vascular_support_planner.py --stage all
```

程序从最新成功的 `support_aware_roi_60x50` manifest 自动发现输入。默认结果目录为：

`outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/`

先读其中的 `vascular_support_summary.json` 与 `vascular_support_report.md`。只有全部生产门槛通过，才会有 `production_selected` 和 `READY_FOR_HUMAN_REVIEW_AND_COUPON_TEST`。`NEEDS_ADJUSTMENT` 下导出的几何属于诊断候选，不是完整模型的打印推荐。

## 分阶段运行

|参数|作用|
|---|---|
|`--stage plan`|表面、真实层分析、FreeFloating、SWC 融合与弧长筛选|
|`--stage geometry`|生成三个由实际线宽派生的支撑候选，检查间隙、底脚、自身可打印性与拆除|
|`--stage removal`|复用已有显式几何，重新检查整组件折颈后的取出路径|
|`--stage coupon`|生成 27 格试片，并实际关闭自动支撑切片|
|`--stage slice`|对已有三个候选做完整装配与独立支撑的 Bambu 真实切片|
|`--stage report`|从已保存结果重建报告与 18 张主要 QC 图|
|`--mode benchmark --stage all`|单独重跑旧 Bambu/Orca/Guided 对照；结果不会被选为生产方法|

首次运行用 `all`；分阶段命令依赖此前阶段的文件。输入、工艺或配置不匹配时，缓存检查会要求重新运行 `plan`。`--source` 指定已存在的上游 ROI 目录，`--output` 指定独立输出位置，`--config` 指定 YAML。不要把输出指向上游模型目录。

## 主要文件

|位置|意义|
|---|---|
|`geometry/custom_support.stl`|与原生产模型同一打印坐标的显式支撑；只有批准脚垫接触盒内底板|
|`geometry/final_mold_with_breakaway_support.stl`|冻结原模型加显式支撑的装配体；未通过门槛时仅供诊断|
|`geometry/frozen_mold.stl`|原生产模型副本，便于分离检查|
|`geometry/support_tips.vtp`、`support_graph.vtp`|尖端与支撑图数据|
|`custom_support_auto_support_off.3mf`|当前诊断/入选装配的真实 Bambu 切片，自动支撑关闭|
|`custom_support_only_auto_support_off.3mf`|显式支撑与其合法底板基础的独立切片，不借助血管或侧壁支撑|
|`coupon/support_breakaway_coupon.stl`|三种颈宽 × 三种尖端宽度 × 三种顶部间隙的 27 格试片|
|`coupon/bambu/sliced.3mf`|试片真实 Bambu 切片|
|`coupon/coupon_test_plan.md`|试片观察、折断与拆除记录方法|
|`custom_candidates/custom_00…02/`|各参数候选、真实切片、间隙、自身悬空和整组件移除记录|
|`tables/`|原始检测、融合区间、筛选尖端、边、根部、组件移除与对比表|
|`QC/`|18 张显式路线主要图；清单见 visual_qc_manifest.json；旧图保留但不混作新结果|
|`implementation_manifest.json`、`artifact_manifest.json`|本轮代码、依赖版本、测试统计及交付物大小、SHA256|
|`verification.json`|最终 STL、3MF CRC、实际工艺与输入不变性的复核记录|
|`benchmark/guided_v1/`|方向修正前的代码结果、报告、测试与 3MF，全部 BENCHMARK_ONLY|

不要把历史 `BG001_RMCA_60x50_guided_support*.3mf` 误认为当前显式路线的生产推荐。不要单独移动、缩放或自动落床 `custom_support.stl`；装配 STL 已包含共同坐标。

## 代码职责

`support_custom_detection.py` 调用原四源检测器，将证据映射为命名分支区间。`support_contact_planner.py` 完成确定性的弧长覆盖筛选。`support_custom_graph.py` 建立向下合并的支撑图、校正真实尖端间隙，并使用竖直立柱修复失败合并。偏移底脚用于避让下方血管，支柱截面保持与打印层平行，避免倾斜端面制造悬空薄片。`support_breakaway.py` 生成受控脚垫、窄颈、渐扩树干及试片。`support_removal.py` 在删除窄颈后检查整组件的有限确定性移动，使用 FCL 距离与运动界避免只检测离散端点。`support_custom_qc.py` 检查实际 STL 精度下的闭合体、实体相交、自身悬空和真实刀路特征。`support_custom_pipeline.py` 调度生产路线，原 `vascular_support_pipeline.py` 保留旧基准实现。

盒外端口与“根部只能在盒内底板”的冲突必须作为明确失败保留；不能在盒外偷偷落床、以侧壁为根、删掉危险区间或改原模型。探针能接近某个树干不代表整个支撑能取出，所以新的生产判定使用整组件扫掠路径。

## 工艺、依赖与复核

优先读取本机的 X1 Carbon 0.2 mm、ABS Vessel、0.10mm Standard Vessel。当前线宽为 0.22 mm，因此三种尖端/窄颈宽为 0.44/0.66/0.88 mm；下次运行仍以实际工艺解析值为准。支撑参数见 `config/vascular_support_planner.yaml` 的 `custom_support`。

切片器负责 G-code；程序的命令白名单只允许本地操作。原生/Guided Orca 基准的外部工具安装入口仍为 `./run_python.sh tools/prepare_support_engines.py`。开源来源、原作者 ZIP 实际内容及许可边界见根目录 `open_source_support_audit.md`。

```bash
./run_python.sh -m pytest -q tests/test_s1_8_vascular_support_planner.py
./run_python.sh tools/verify_s1_8_support_outputs.py
```

所有永久夹具位于 `tests/data/support/`。输入模型、SWC、端口、盒体与 s1-2 至 s1-7 的 SHA256 结束时全部复核。完整打印前应先执行 coupon；软件不验证真实折断力，也不自动发送打印任务。

文件复核中的 `artifact_integrity_passed` 只表示交付物完整、工艺与源文件对应，不表示完整模型达到打印条件。生产是否通过仍看 `production_ready`、`production_status` 和 `hard_gate_failures`。
