# BCF V2：内墙支撑、完整侧墙与固化后切割

当前生产入口仍是 **s1-9**，没有新建 s1-10。完整 ABS 血管芯及四根延长端口都放在四面封闭 PLA 墙内。PLA 没有底板，培养皿直接形成 PDMS 底面。端口由内墙长出的竖直薄片承托，另用独立插片防止明显上浮；支撑和插片留下的印痕随后随外围 PDMS 一起切掉。

上一版穿墙端口、U 槽、外部锁销和密封套结构已退出默认生产路线。V1 的代码、测试和输出原位保留作为历史，不要混用其 3MF。

## 如何使用

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-9_vascular_casting_fixture.py
```

默认完整运行包括：只读输入核对、9 组小幅刚体摆放比较、三种托座间隙候选、2/3 线宽真实切片比较、装配与切割检查、STL/STEP、21 张 QC 和中文报告。不会重新计算 SWC、语义迁移、ROI、VascularMD，也不会更改现有可视化或控制打印机。

配置：[internal_wall_supported_fixture_v2.yaml](../config/internal_wall_supported_fixture_v2.yaml)。临时试验建议指定 `--output /绝对路径/实验目录`。`--geometry-only`、`--no-step` 只用于诊断，结果不能通过正式交付门槛；`--no-render` 可省略图像，但最终审阅应使用完整输出。

## 输入与输出

输入由 `vascular_processing/support_inputs.py` 发现，是 `support_aware_roi_60x50/candidates/candidate_0008/manufacturing/core_with_ports.stl` 及其端口 manifest、SWC 和映射。原始 STL 的 SHA256 必须为当前受保护值；本轮输出不回写上游目录。

新结果集中在：

```text
outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/
  vascular_casting_fixture/
    ... V1 历史结果，不再用于生产 ...
    internal_wall_supported_fixture_v2/
```

以下路径相对于新结果目录：

|路径|作用|
|---|---|
|vascular_core/frozen_ABS_core_with_extensions.stl|原坐标、原字节的独立 ABS 源副本|
|vascular_core/assembly_placed_ABS_core.stl|只施加整体刚体变换的装配参考副本|
|assembly/casting_core_transform.json|源坐标与夹具坐标转换；夹具 Z=0 为皿面|
|frame/closed_wall_bottomless_frame_with_support_blades.stl、.step|四壁无底 PLA 框架、下托座、内侧限位和切割肋|
|frame/frame_with_wall_supports.3mf|当前实际 PLA 框架切片，Auto Support OFF|
|retainers/port_*_anti_float_retainer.stl、.step|四个独立防浮插片，文件本身在装配坐标中|
|retainers/retainers_print_layout.stl|已将各插片的平整侧面朝下，分开放在打印板上|
|retainers/anti_float_retainers.3mf|真实独立插片切片，Auto Support OFF|
|assembly/assembled_casting_fixture.step|独立零件装配；血管保留原三角面，不与 PLA 合并|
|assembly/assembled_review_only.stl|仅供审阅的多零件文件，不能直接作同材料一体打印|
|diagnostic/PDMS_before_trim.stl、PDMS_after_trim.stl|真实布尔计算的切割前后 PDMS，仅作诊断|
|trim/trim_planes.vtp|有限切面，包括端口主切面和局部切口侧面|
|trim/trim_alignment_ribs.stl|导向肋；其内侧面与主切面重合|
|trim/trim_strategy.json|每面墙的局部切口/整侧条带比较、实际切面清单和移除检查|
|trim/trimmed_port_sections.csv|真实截面积、等面积直径和暴露检查|
|tables/|摆放、端口直段、接触、装配和切割明细|
|config/clearance_study/|0.11、0.22、0.33 mm 下托座候选 CAD 与真实间隙比较|
|bambu/frame_2lines、frame_3lines|0.44、0.66 mm 薄片真实切片、日志、逐层刀路保留检查|
|bambu/retainers_v2|实际插片切片、G-code、归档配置与日志|
|QC/|21 张模型、装配、切割前后及实际切片检查图|
|reports/、protection/|测试、实现文件哈希与源文件/V1 历史保护检查|
|internal_wall_supported_fixture_summary.json|机器可读状态与各项门槛|
|internal_wall_supported_fixture_report.md|先回答 20 个制造问题的中文报告|

## 当前几何与装配逻辑

内腔约 **99.17 × 76.68 × 41.01 mm**；墙厚 2.00 mm，向外法兰宽 5.00 mm。完整芯仅平移，端口末端距对应内墙约 1.50 mm，中央血管最低处距皿面 10.00 mm。含法兰所需平坦圆形底面直径至少 **145.01 mm**；这是平底有效直径，不是培养皿标称外径。实际培养皿尚未测量。

采用本机 **X1 Carbon / 0.2 mm 喷嘴 / PLA Basic / 0.10 mm 层高**，实际线宽 0.22 mm。2 线宽和 3 线宽都通过真实刀路保留检查，因此采用 0.44 mm 薄片，切割肋保持同厚。当前框架估时约 5 小时 6 分，四个插片合计约 12 分；这些是切片估计，并非实物计时。

O2 为主基准：较窄下托弧和带正间隙的末端轴向挡块；其余端口使用较宽双托弧，并保留 ±1.00 mm 轴向浮动。完整芯先从上方竖直放入，插片再从上方放下并向内墙滑动 1.60 mm，进入内部斜面限位。只竖直落下而不完成这一步滑动，不能认为已经防浮。

四根已冻结端口各有约 9 mm 保守直段。主切面在其中点，所以最终 PDMS 外表面至保守曲率过渡还保留 4.50 mm 直段。支撑、插片和肋都留在主切面的外侧牺牲区。每个局部切口通常有三刀，O2 靠转角仅需两刀；实际有效切面共 11 个，以 JSON/VTP 为准。

固化后应先沿肋的内侧面和切口侧面切开，随后移走外围 PDMS 与 PLA，最后取出中央 PDMS。不要试图先从未切开的 PDMS 中强拔薄片。几何证明不包含真实黏附、材料剥离力或手工切割误差。

## 代码分工与测试

|代码|职责|
|---|---|
|s1-9_vascular_casting_fixture.py|当前 V2 命令入口|
|casting_fixture_v2_pipeline.py|串联输入、CAD、检查、切片、图像及报告|
|casting_fixture_v2_geometry.py|整体摆放、四壁无底 CAD、托座、插片、竖直肋及公差候选|
|casting_fixture_v2_audit.py|连续装配路径、接触位置、抗浮扰动、PDMS 布尔差、有限切面和倒扣检查|
|casting_fixture_v2_slicing.py|解析真实 Bambu 刀路，逐层检查薄片与肋是否保留|
|casting_fixture_v2_review.py|独立 QC 和中文报告，不修改 s1-2/s1-3 的 UI|
|casting_fixture_v2_synthetic.py|十类永久合成数据，包括不足直段、碰撞和错误切割的负例|
|casting_fixture_profiles.py|共享的实际预设读取与仅本地切片工具|

模块均位于 `vascular_processing/`。无 `v2` 后缀的旧 geometry/audit/pipeline 保留供 V1 历史测试使用；新生产入口不调用其旧结构生成器。

```bash
./run_python.sh -m pytest -q tests/test_s1_9_vascular_casting_fixture.py
./run_python.sh -m pytest -q tests
```

底部密封、防浮强度、打印收缩、手工切割及 PLA/PDMS 脱离尚未物理验证。当前软件状态是 `READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_SEAL_PENDING`，不是制造认证。实测培养皿后可填写 `config/casting_dish_measurements.json` 的 `dish_inner_flat_diameter`、`dish_wall_height`、`dish_floor_thickness`（mm）再核对，勿填猜测值。
