# 全 ABS 牺牲浇注系统：入口、边界与输出

当前入口仍为 **s1-10**，使用已选中的 **90° COMPACT_C_RETAINER_V2** 插片。2026-10-04 的尺寸调整将侧壁和底板统一为 **2 mm**，完整血管表面到内底面、顶部开口平面的净距均设为 **5 mm**。盒体外尺寸为 **103.17 × 80.68 × 31.01 mm**；顶部仍开放。

本轮适度缩短盒体和托座竖直部分。血管芯、半径、端口路线、端口长度、C 插片几何及其配合保持不变，装配整体沿 Z 方向下移 6 mm；不改变既有血管或 ROI 交互界面。旧模型及检查记录原位保留。

## 六组件同板打印

当前六件同板交付文件为 [BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf](../outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/all_abs_sacrificial_casting_system/compact_c_retainer_v2/mold_2mm_clearance_5mm/print_plate_all_components/BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf)。在 Bambu Studio 中以 **“打开项目”** 方式载入，保留项目预设和逐对象支撑设置。

盒体位于打印板前部，血管芯位于后部，I1/O1/O2/O3 四个插片独立排在右侧。所有组件保持原有打印姿态和尺寸，分别落床，使用 **X1 Carbon / 0.2 mm / ABS / 0.10 mm 层高 / 按层打印**。只有血管芯开启树状支撑；盒体和插片关闭自动支撑。它们是可分别选择的六个对象，打印后再装配。

原生 Bambu Studio 02.07.01.57 已实际切片成功，无关键警告；归档、几何回读、逐对象支撑与参考挤出路径覆盖均通过。估计 **12 小时 48 分钟、32.93 g ABS**，不等同于实测打印结果。10 mm brim 预留范围及血管额外支撑预留范围互不相交，最小分隔约 8.57 mm；实际模型和支撑挤出均在可用打印板内。

只重新排列并切片现有已验收组件，可运行：

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-10_all_abs_sacrificial_casting_system.py --print-plate
```

实现模块为 `vascular_processing/all_abs_print_plate.py`，测试为 `tests/test_s1_10_all_abs_print_plate.py`。该模式读取下述 `D/` 中的最新外模和插片打印姿态，并复用形状不变的冻结血管打印参考；不重新建模。输出位于 `D/print_plate_all_components/`，其中 `single_plate_layout.png` 为布局预览，`single_plate_summary.json` 为总检查记录，`bambu_all_abs_plate/` 保留实际切片命令、配置、日志和 G-code。使用 `--output <新目录>` 可另存；不与 `--config` 或跳过检查的选项组合。不会上传或启动打印机。

## 运行

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-10_all_abs_sacrificial_casting_system.py
```

默认配置为 `config/compact_abs_mold_2mm.yaml`，重新生成新高度的外模、装配 STL/STEP 和 PDMS 诊断，检查连续装入、配合余量、有限抗浮姿态、壁侧溶剂通路、双 SWC 可达性及原裁切面，实际切片新外模和四个插片，并输出 5 张检查图和中文调整说明。

打印配置继续为 X1 Carbon / 0.2 mm / ABS，外模和插片 Auto Support OFF。切片输入、配置及二进制签名相同时复用校验后的结果。不会上传或启动打印机。

`--config <yaml>` 指定实验配置，`--output <新目录>` 指定新输出。`--geometry-only`、`--no-step`、`--no-render` 会因缺少必要证据保留 `NEEDS_ADJUSTMENT`。

若明确需要复现上一轮仅调整 C 插片的历史流程，使用：

```bash
./run_python.sh s1-10_all_abs_sacrificial_casting_system.py --config config/compact_c_retainer_v2.yaml
```

历史流程保留原 3 mm 底板、10/12 mm 上下净距，不是当前默认的盒体尺寸。

## 输入与修改边界

以下路径均相对于项目根目录：

```text
B = outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/all_abs_sacrificial_casting_system
C = B/compact_c_retainer_v2
D = C/mold_2mm_clearance_5mm
```

`B/` 是原全 ABS 阶段，`C/` 是上一轮紧凑 C 插片阶段，**`D/` 是当前调整后的模型**。

读取父级布局、配置、CAD 和清单，冻结上游 `support_aware_roi_60x50/candidates/candidate_0008/` 的输入：

- `manufacturing/core_with_ports.stl`：完整血管与四端口；不缩放、不重新拟合。
- `vascular_core.stl`：中央血管树，用于禁止辅助接触和裁切损伤检查。
- `manufacturing/geometry.json`：端口半径、直段、轴向、SWC 锚点与已接受连接。
- `compensated.swc`、`VascularMD/compensated_vmd_smooth.swc`：双 SWC 分支可达性。
- `compensated_edges.csv`、`roi_casting_transform.json`：分支和坐标关系。
- `C/retainers/port_*_C_retainer.step`：实际已接受 C 插片，只改变统一装配位置。

O2 主基准、其余次基准、0.22 mm 上方间隙、0.12 mm 导向间隙以及 ±1 mm 次基准相对配合余量均继承原值。新盒体根据完整血管表面高度计算；若特定输入的插片或导向过高，会扩大顶部余量以保留至少 1 mm 净空，不通过缩短插片来迁就盒体。本例最终上下净距恰为 5 mm。

## 当前输出

|位置|用途|
|---|---|
|`D/outer_mold/ABS_sacrificial_outer_mold_2mm.stl/.step`|2 mm 侧壁及底板的新外模，含适配高度后的原壁生托座|
|`D/outer_mold/ABS_sacrificial_outer_mold_2mm_X1C_ABS.3mf`|新外模的真实无自动支撑打印文件|
|`D/retainers/port_<I1/O1/O2/O3>_C_retainer.stl/.step`|四个 C 插片的新装配坐标；几何未变|
|`D/retainers/C_retainers_X1C_0p2_ABS.3mf`|四插片独立打印布局与实际切片|
|`D/assembly/all_ABS_casting_system_C_retainer_2mm.step`|六个独立组件的新总装|
|`D/assembly/all_ABS_C_retainer_mold_2mm_review.stl`|完整审阅网格，不用于联合打印|
|`D/vascular_core/assembly_placed_ABS_core.stl`|新装配位置的血管芯|
|`D/vascular_core/BG001_RMCA_frozen_ABS_core.stl`|源字节完全一致的原始坐标血管副本|
|`D/assembly/casting_core_transform.json`|新旧变换及共同 Z 平移量|
|`D/diagnostic_PDMS/PDMS_after_all_ABS_removed.stl`|新高度下溶解后的 PDMS 诊断|
|`D/diagnostic_PDMS/PDMS_after_trim.stl`|外围裁切后的最终 PDMS|
|`D/QC/`|同尺度新旧对比、厚度与净距截面、新总装、裁切 PDMS、实际刀路共 5 张图|
|`D/tables/`|尺寸变化、装配、间隙、轴向余量及壁侧通路|
|`D/solvent_access/`|实际接触证据、去除插片节点后的端口图、双 SWC 可达性|
|`D/trim/`|有限裁切面、端口暴露、保留直段与辅助孔腔检查|
|`D/bambu/outer_mold/`、`D/bambu/retainers/`|真实切片命令、配置、归档、G-code 和特征保留检查|
|`D/reports/run.log`、`D/reports/pytest_s1_10.log`|本次生成和回归日志|
|`D/source_protection.json`|源数据和历史文件开始/结束 SHA256|
|`D/mold_revision_summary.json`、`D/盒体调整说明.md`|完整数值结果及中文说明|

**使用新盒体时应采用 `D/` 的装配坐标。** 原 `B/` 或 `C/` 的血管和插片坐标未下移，直接与新盒体叠加会造成装配高度错误。冻结血管本体的独立打印参考文件仍可使用 `B/vascular_core/ABS_core_reference.3mf`，其形状与单独落床后的姿态不因总装平移而改变。

`C/` 保留上一轮三种 C 包覆角比较、18 张检查图、34 项检查结果及完整中文研究报告。旧长杆最小归档仍为 `B/legacy_retainer_v1/`。这些历史文件不应误认为当前新高度盒体。

## 代码职责

下表模块均位于 `vascular_processing/`。

|入口或模块|职责|
|---|---|
|`s1-10_all_abs_sacrificial_casting_system.py`|默认进入 2 mm 盒体调整；显式旧配置可复现 C 插片阶段|
|`compact_mold_revision.py`|保护历史；保持 XY、旋转和血管形状；调整 Z 放置、内底和盒沿；重建外模、装配、孔腔并切片|
|`compact_mold_review.py`|新旧尺寸对比、截面和总装等 5 张 QC，中文调整说明|
|`compact_c_geometry.py`|只读加载原布局和 CAD；C 弧、短颈及短舌参数|
|`compact_c_audit.py`|C 插入、延长直段最近点、轴向范围、有限运动、壁侧及插片自身可达性|
|`compact_c_slicing.py`|C 轴向端面落床摆放；实际 C 壁、导向和通道刀路检查|
|`all_abs_geometry.py`|既有闭底外模、壁生托座、真实下表面接触和浅导向几何内核|
|`all_abs_audit.py`|既有连续装入、双 SWC、PDMS 布尔和裁切检查|
|`all_abs_slicing.py`|真实 Bambu 配置、归档 CRC、几何及自动支撑核查|
|`compact_c_pipeline.py`、`compact_c_review.py`|上一轮 C 插片比较及 18 张 QC，保留为独立可复现流程|
|`tests/test_s1_10_all_abs_sacrificial_casting_system.py`|原 V1 和 C 插片的永久回归|
|`tests/test_s1_10_mold_2mm_revision.py`|实际 STL 射线测厚、上下表面净距、刚体变换、新切片及既有功能回归|

## 工艺边界

上下 5 mm 为当前合理收紧空余区域的设计值，尚不是经过浇注实验优化的材料参数。四条壁侧通路继续由“外部浴—外墙—窄片—接触区—延长端口”建立，插片后端接触只服务其自身去除。名义接触若被 PDMS 薄膜隔开，实物通路可能失效。

裁切面保持相同几何平面，其有效高度随新 PDMS 域缩短；切后辅助空腔为数值零，四端口各保留 4.50 mm 已确认直段。血管芯和插片禁止蒸汽平滑；外模可选处理仍未校准。当前结论仅支持软件审阅，不代表打印强度、真实抗浮力、溶解和 PDMS 行为已验证。

## 测试

```bash
./run_python.sh -m pytest -q tests/test_s1_10_all_abs_sacrificial_casting_system.py tests/test_s1_10_mold_2mm_revision.py
```

先查看 `D/盒体调整说明.md` 与 `D/QC/`，再进行实物干装配和工艺确认。
