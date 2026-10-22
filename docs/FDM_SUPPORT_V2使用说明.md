# FDM 血管支撑 V2：入口与结果阅读

本轮继续 s1-8，没有新建 s1-9。它为已经冻结的 BG001 RMCA 血管、四个端口及 2 mm 浇注盒寻找打印姿态和可分段移除的显式支撑。打印姿态可以改变，浇注姿态不变；不缩放、不改半径、不移动局部血管、不改盒壁。

## 运行

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-8_vascular_support_planner.py --stage all
```

默认模式为 `v2`。旧命令中的 `--mode custom` 也映射到 V2，不会重新生成 V1 锥柱。`--mode benchmark` 仅运行原生对照。V1 原代码、测试及已生成结果保留，旧生产入口会明确拒绝执行。

完整流程包含昂贵的碰撞、布尔、分层与本机 Bambu 切片，不能按秒级交互程序预期。后续单独检查已有结果：

```bash
./run_python.sh s1-8_vascular_support_planner.py --stage removal
./run_python.sh s1-8_vascular_support_planner.py --stage slice
./run_python.sh s1-8_vascular_support_planner.py --stage coupon
./run_python.sh s1-8_vascular_support_planner.py --stage report
./run_python.sh -m pytest tests/test_s1_8_vascular_support_planner.py -q
```

完整 `all` 已包含主干分段细化。分阶段开发可在 `plan` 后执行 `--stage refine`，给保留的四个候选补充按层高量化的主干断颈，再运行 `removal` 和 `slice`。这一步保留 `before_trunk_segmentation/` 并重新检查几何和自身可打印性，不把“增加断颈”直接当作可拆除证明。规划被中断时可用 `--stage resume`；它核对上次输入、工艺和配置，复用已完成候选或保存的构件，继续未完成姿态。

`orientation` 只作粗筛；`plan` 完整重建 12 个规划候选。后续阶段读取输入/配置指纹，变化时要求重新规划。真实切片适配器校验几何、工艺及切片器指纹，命中时复用已验证归档。不会上传、发送或启动打印机。

`config/vascular_support_planner.yaml` 中的 `custom_support`、固定 5 mm 和固定根部限制属于保留的 V1 配置，不能解释为 V2 的生产约束。V2 的实际策略另写入 `v2_policy_runtime.json`：支撑数量不设上限、分叉只增加惩罚、根部允许打印板与合法盒面；主门槛是逐层完整多边形支撑性。

## 输入、处理与输出

输出根目录为：

`outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/v2/`

输入由 `support_inputs.discover()` 定位：当前 `support_aware_roi_60x50` 下 2 mm 盒体、完整倒模件、血管及端口实体、优化 SWC、语义映射和浇注变换。本机保存的 X1 Carbon / 0.2 mm / ABS Vessel / 0.10mm Standard Vessel 配置经既有适配器读取，不从旧 P1S 配置推断。

|文件或目录|用途|
|---|---|
|`vascular_support_v2_report.md`|先读：15 个关键问题、真实结果、未通过的条件|
|`vascular_support_v2_summary.json`|机器可读状态；`production_selected` 为空时不得当作生产通过|
|`orientation_candidates.csv`|128 个刚体姿态的粗筛指标及拒绝原因|
|`planning_results.json`|12 个完整规划结果，保留失败项|
|`top4_candidate_ids.json`|进入真实切片的四个候选；后续阶段固定这四个，避免验证后悄悄换人|
|`candidates/Oxxx/`|每个候选的几何、检测、支撑图、失败原因和切片归档|
|`print_transform_v2.json`|冻结浇注坐标 → 诊断候选打印坐标，比例 1|
|`casting_restore_transform_v2.json`|严格逆变换，用于恢复浇注姿态|
|`support_v1_postmortem.md`|V1 零合并的来源追踪；明确区分历史记录与冻结源码重放|
|`V1_vs_V2_topology.csv`、`v1_vs_v2.csv`|拓扑与负担比较；V1 弧长覆盖与 V2 层面积覆盖不可直接混为一项|
|`QC/`|20 个主检查图，另含正、侧、顶和剖切视图|
|`coupon/`|4 个间隙、4 个桥跨、轨道、Y 合并、窄颈及细圆杆试片和真实切片|
|`test_results.xml`、`test_run.log`|永久测试证据|
|`source_protection.json`|冻结上游 SHA256 检查|

每个候选包含 `custom_support.stl`、`mold_without_support.stl`、`mold_with_custom_support.stl`、`support_roofs.stl`、`support_trunks.stl`、`support_breakaway_parts.stl`。装配切片成功后保存 `candidate_Oxxx_X1C.3mf`，可在 Bambu 打开检查；不要把无支撑整件与有支撑整件同时叠放打印。

`bambu_validation.json` 同时记录切片成功、Auto Support OFF、真实顶垫覆盖、构件保留、根部起始路径和血管刀路；不能只看 `slicer_success`。`bambu/floating_warning_localization.csv` 是几何与 G-code 的独立关联，不冒充原生警告坐标。

`support_removal_plan.md/json` 给出模拟断颈后的拆除顺序；路径按浇注开口定义，未拆碎片仍作为障碍。`TRAPPED_REMOVAL_FRAGMENT` 表示给定路径集内未找到安全路线，不能通过改变颜色或隐藏几何消除。

`sensitivity/roof_gap_*.stl` 只是不同间隙的接触顶垫候选，不是完整打印件。当前 0.20 mm 是回退参数，标记 `PHYSICAL_GAP_UNVALIDATED`。桥跨 2/4/6/8 mm 的软件比较不等于实物校准；默认 2 mm 标记 `BRIDGE_LIMIT_UNCALIBRATED`。只有实物通过且 `config/fdm_bridge_calibration.json` 的工艺指纹吻合才读取校准上限。

校准文件的 `profile` 必须与 `v2_policy_runtime.json` 中的 `bridge_calibration_signature` 完全一致，包括打印机、工艺名称、喷嘴、层高、线宽、材料、桥参数、速度和温度。另需真实填写 `physical_test_passed: true` 和实测通过的 `maximum_bridge_span_mm`；仅有软件切片结果不能填写为通过。

## 后续开发位置

`vascular_processing/support_v2_orientation.py` 处理姿态；`support_v2_layers.py` 复用旧分层实现并做完整多边形覆盖；`support_v2_planner.py` 负责顶垫、共享分支、合法根部和拒绝记录；`support_v2_primitives.py` 生成 FDM 量化构件；`support_v2_removal.py` 负责断颈及带顺序的扫掠；`support_v2_validation.py` 读取实际切片刀路；`support_v2_coupon.py` 和 `support_v2_review.py` 生成试片与报告。

已有三角面、FreeFloating、SWC 和层检测代码及 V1 合成测试继续保留。新合成几何位于 `tests/data/support/v2/`，测试继续写入 `tests/test_s1_8_vascular_support_planner.py`。不要修改 s1-2、s1-3 的交互视觉来掩盖支撑问题。
