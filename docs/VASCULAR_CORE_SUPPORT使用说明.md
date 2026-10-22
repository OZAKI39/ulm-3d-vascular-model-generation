# s1-11：独立血管芯的显式打印支撑

本阶段只读取 s1-10 当前冻结的 BG001 RMCA 血管芯及其四条延长端口，为其生成独立 ABS 支撑。外模、夹具、插片、PDMS、溶剂通路和裁切设计均不参与支撑计算，也不会重新生成。

> 本文记录 s1-11 初始建模流程。当前默认入口已进入后续局部闭合与安全拆除检查，详见 [本轮说明](S1_11覆盖闭合与安全拆除.md)。旧报告可使用 `--stage legacy-review` 重建。

## 初始阶段入口

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-11_vascular_core_support.py
```

支持 `--stage prepare`、`geometry`、`slice`、`review`；默认 `all`。`--output <目录>` 可另存一轮结果。已有几何和实际切片用于续跑；修改支撑参数或生成算法后，应指定新的输出目录，以避免把旧候选误作新结果。

程序执行 48 个确定性刚体姿态的粗筛，前 6 个姿态各尝试 0.20/0.25 mm gap 与 2/3 线宽组合，前 3 个姿态进入实际 Bambu 切片。未通过时继续执行逐层补撑和根据实际缺口进行的局部修正。所有修正只作用于支撑。失败候选、数值问题和未覆盖区域均留在记录中，不会自动改成通过。

## 输入

通过 `all_abs_sacrificial_casting_summary.json` 和 `reports/source_discovery.json` 追溯冻结 STL、补偿 SWC、VascularMD 平滑 SWC、端口曲线和坐标变换。使用最新的 s1-10 冻结血管副本，核对其 SHA256 与原始制造源一致。程序没有 s1-8 候选输入路径。

打印参数从本机最新保存的 `0.10mm Standard Vessel`、`Bambu ABS Vessel` 和 X1 Carbon 0.2 mm 配置读取。只在本次导出的配置副本中关闭 Auto Support；用户保存的预设文件保持不变。引擎版本、二进制及配置哈希归档于输出的 `bambu/profiles/`。

## 输出

根目录：

```text
outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_core_support_production/
```

|位置|用途|
|---|---|
|`source/`|冻结源副本与数据来源，源 STL 原始字节保留|
|`orientations/ranking.json`|48 个姿态、刚体变换和粗筛指标|
|`candidates/<ID>/`|独立支撑 STL、血管打印副本、下表面目标和逐候选检查|
|`candidates/<ID>/bambu/`|原生切片、3MF、G-code、命令、配置及日志|
|`selected/frozen_vascular_core.stl`|仅经过刚体放置的血管打印副本|
|`selected/custom_vascular_support.stl`|独立的显式 ABS 支撑|
|`selected/vascular_core_plus_support_review.stl`|仅拼接的审阅网格；未布尔合并血管与支撑|
|`selected/BG001_RMCA_vascular_core_supported_X1C.3mf`|原生切片工程；血管与支撑为两个 normal part，Auto Support OFF|
|`selected/print_transform.json`、`inverse_transform.json`|正向和逆向刚体变换，缩放为 1|
|`tables/candidate_comparison.json`|实际切片候选的覆盖、区域、支撑自身和体积比较|
|`tables/removal_scorecard.json`|拆除各项原始指标及硬门槛，不使用黑箱总分|
|`QC/`|22 张检查图；第 16、17 张保留实际未覆盖位置|
|`vascular_core_support_summary.json`|状态、硬门槛、覆盖、拆除、切片与保护结果|
|`vascular_core_support_report.md`|中文运行报告，逐项回答制造与验证问题|
|`vascular_support_removal_guide.md`|根据几何生成的剪口及片段移除顺序|
|`validation/`、`reports/`|真实水平管回归、源文件哈希、测试和运行日志|

**先读取状态，再使用模型。** `NEEDS_ADJUSTMENT` 时，`selected/` 仅是诊断候选，不能作为已验收生产文件。中心线达到 100% 并不代表逐层区域、支撑自身和拆除检查全部通过。当前实物去除力、表面质量及分支损伤风险均未校准。

在 Bambu Studio 中以“打开项目”载入 3MF，保留 part 关系和关闭自动支撑的设置。本工具只在本地生成文件，不上传或启动打印机。

## 代码职责

|模块（`vascular_processing/`）|职责|
|---|---|
|`core_support_inputs.py`|只读发现、哈希保护、当前预设、SWC/端口分析、刚体姿态粗筛|
|`core_support_geometry.py`|实际 STL 下表面射线、连续顶面、分段幕墙、投影相交处理、逐层补撑及底脚|
|`core_support_slicing.py`|两独立 part 的 3MF、原生切片、实际线宽刀路、完整下表面及逐层区域检查|
|`core_support_repair.py`|根据实际缺口加局部支撑，拒绝与血管相交的补块，不放宽判定|
|`core_support_removal.py`|3 mm 探针、有限虚拟剪切、包覆角和保留其余障碍物的连续移除扫掠|
|`core_support_travel.py`|只检查实际非挤出移动，并仅与此前已打印材料比较|
|`core_support_fixtures.py`|九类合成输入、支撑及真实 30 mm 水平管回归|
|`core_support_pipeline.py`、`core_support_review.py`|运行编排、硬门槛、候选记录、检查图与中文报告|

逐层检查将真实向下表面三角形裁切到对应 gap 高度区间，再与实际挤出覆盖及完整需求多边形求交；不把切片层顶高度直接当成血管下表面。同一位置更高且间隙过小的挤出会导致检查失败；不能用较低幕墙层冒充合格顶面。逐层检查也只接受暴露的实际顶面。G-code XY 小数精度和 STL 浮点误差单独记录，不用扩大容差追求 100%。

## 回归

```bash
./run_python.sh -m pytest -q tests/test_s1_11_vascular_core_support.py
```

永久合成数据位于 `tests/data/vascular_core_support/`。30 mm 水平管的实际切片回归使用位于管外的独立层连续性见证柱，避免理想水平模型的全局空层被 Bambu 拒绝；该柱不计入任何下表面承托指标，也不加入真实 BG001 模型。

当前阶段采用严格门槛：即使全部程序回归测试通过，只要真实模型覆盖不足 100%、存在自身悬空或拆除失败，仍输出 `NEEDS_ADJUSTMENT`，命令返回退出码 2。这是制造验收未通过，不等同于程序崩溃。失败的 STL 导出尝试记录在 `reports/generation_errors.json`。历史运行日志中的中间指标不替代最终 JSON 和报告。
