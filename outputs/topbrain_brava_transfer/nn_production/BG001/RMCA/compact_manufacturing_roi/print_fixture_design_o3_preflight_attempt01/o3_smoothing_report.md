# O3 局部表面连续性检查

## 1. 原来哪里有问题

圆口虽然对齐了，但原血管侧壁与人工直管的表面方向并不一致，所以连接处仍有一圈折肩。已有封口是模型的终端闭合面，不是必须保留的解剖结构。当前 5 mm 中心线平均方向与旧 cap normal 相差 None°。旧实现保留为 BASELINE_G0_ALIGNMENT 对照，本轮没有用它生成新的生产几何。

## 2. 为什么旧测试没有发现

旧测试检查位置、圆口轮廓、中心、轴向、闭合与体积，即“接没接上”。这些是位置连续（G0）的证据，不能证明侧壁切向连续（G1），也就没有回答“接得顺不顺”。近乎零的轮廓误差与轴向误差仍可伴随明显折肩。

## 3. 这次用了什么方法

本轮首先实测既有 VMTK 环境，并完成跨环境 VTP 文件往返；通过状态为 True。主环境没有导入 VMTK。指定环境实际为 VMTK 1.5.0 / VTK 9.2.6，没有 SetInterpolationModeToRamp，也没有 SetPreserveCrossSectionShape。现有 Slicer 的已安装扩展为空，不能提供替代 VMTK 实现。官方新版源码存在这些 API，但它不是本机已安装的运行时。全部版本、方法列表、源码提交与文件哈希保存在 api_audit 和 VMTK/bridge_smoke 中。

按照用户第 60 条“仍没有：停止并报告”，当前停止代码为 **O3_PREFLIGHT_ERROR**。没有安装新的 VMTK、没有自行实现 ramp、没有把 thinplatespline 或 linear 当成 ramp。尚未进入正式剪切、平滑过渡、CadQuery 接续或 Boolean 合并。

[官方 VMTK 实现](https://github.com/vmtk/vmtk/blob/master/vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.h)与 [Slicer ClipVessel](https://github.com/vmtk/SlicerExtension-VMTK/blob/master/ClipVessel/ClipVessel.py)的固定提交副本和实际本机版本分开记录；不能因为网页有新方法就认定本机可调用。

## 4. 改了多少原血管

本轮实际剪切 0 mm、局部修改 0 mm。三个预定回退距离从实际直径计算，见 O3_transition_candidates.csv，均标为未执行。没有生成 VMTK 平滑段或新 handoff；对应长度、位置和几何指标为未产生，而非零误差通过。

## 5. 修复前后有多大改善

尚无修复后模型，所以不能报告改善百分比。以下均直接读取现有 aligned STL 得到，未重建基线。截面使用旧 cap-normal 坐标系的 -1 至 +1 mm，步长 0.25 mm；它仅用于基线诊断，不代替任务要求的未来新中心线截面验收。

| 指标 | 当前基线 | 修复后 |
| --- | --- | --- |
| 相邻面法向中位数 | None° | 未生成 |
| 相邻面法向 p95 | None° | 未生成 |
| 相邻面法向最大值 | None° | 未生成 |
| 相邻截面最大面积变化 | None | 未生成 |
| O3 对其他血管间距 | None mm（既有记录） | 未生成 |

特征边由 vtkFeatureEdges 提取；闭环是否绕管一周由周向绕数检查。法向跳变使用成熟网格库的相邻三角面，局部窗口的中心、轴向和尺寸保存在 JSON。曲率只作为可视化，不用单个噪声三角面判失败。局部窗口和百分位可能遗漏窄小但明显的肩部，仍需结合线框图检查；它们不是解剖学标准。

当前图片仅有真实“修复前”图，未生成虚假的前后对比或最终 STL。没有用角度变化制造改善效果。

## 6. 是否影响其它端口

I1/O1/O2 的既有数据和原模型未修改。保护文件 1269 个，哈希不一致数 0。单端口实体此前没有独立文件，因此仅从已选参数确定性重放并保存单端口审计副本，未进行路线搜索；副本哈希和参数哈希均在 JSON 中明确标明来源。不能把这些副本冒称为原先已有的单端口 STL。

## 7. 是否可以进入下一阶段

**NEEDS_ADJUSTMENT**。本轮不是 O3 修复成功。必须先获得能实际调用官方 ramp 的运行时，再执行局部剪切、候选比较与真正的连续性验收。由于用户要求当前条件不满足时停止，本轮没有新增生产 STL、没有 box/core 一体式合并，也没有 Bambu 最终切片。既有模型与可视化保持原样。
