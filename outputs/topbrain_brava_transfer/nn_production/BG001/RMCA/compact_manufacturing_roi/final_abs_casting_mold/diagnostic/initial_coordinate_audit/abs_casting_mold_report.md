# BG001 RMCA 一体式 ABS 牺牲倒模件评估

## 1. 这次做了什么

将已验收的血管芯与新生成的完整盒体作实体合并。盒体四侧及底部封闭，顶部完整开放；没有装配间隙孔、顶盖或顶部横梁。旧开放框架作为历史探索保留，本轮生产目标改为 OPEN_TOP_FIVE_WALL_CASTING_BOX。没有重新生成血管、调整半径、重做端口路径或接口光顺。

## 2. 最终是不是一个 ABS 整体

是。连通分量 1；封闭检查 True，流形检查 True，退化三角形 0。合并后的血管材料损失 0.0 mm³，盒体材料损失 0.0 mm³。原始血管 SHA256 为 `b9fb2d8a57b58b877f940f6c2486f923d93ac71b42793593f4092180dcb7e35e`。

## 3. 盒子结构

BOX_DIMENSIONS_REUSED：沿用已审核尺寸，本轮没有重新优化或扩大盒体。

| 项目 | 数值 |
|---|---|
| 内部尺寸 | 78.29 × 81.85 × 56.39 mm |
| 外部尺寸 | 84.29 × 87.85 × 59.39 mm |
| 侧壁厚度 | 3.0 mm |
| 底厚 | 3.0 mm |
| 顶部开口 | 78.29 × 81.85 mm |
| 开口面积 | 6407.78 mm² |

上行射线与内腔实体交集检查确认顶部无封板。单独盒体 STL 和 STEP 位于 reference/，最终制造件位于 final/。

## 4. 四个端口

| 端口 | 原固定侧壁 | 重叠体积 (mm³) | 轴向交叠 (mm) | 实体连接 |
|---|---|---:|---:|---|
| I1 | -Y | 3.383817 | 3.000 | PASS |
| O1 | +Y | 3.190049 | 3.000 | PASS |
| O2 | -X | 3.654329 | 3.000 | PASS |
| O3 | +Y | 3.711491 | 3.000 | PASS |

要求交叠长度至少达到壁厚的 90%，并且交叠体积为正。端口在完整侧壁中形成实体连接，没有悬在装配孔里。四个血管衔接口的局部几何另存于 tables/frozen_interface_qc.csv。

## 5. PDMS 最薄区域

全局最薄非连接区域位于 140_141 与 O3 之间，距离 1.687713 mm，几何分级为 THIN_LOCAL。O3 间距 1.687713 mm，相对验收值变化 0.000002 mm。距离由 FCL 对最终交付三角面计算，排除同分支及真实分叉、端口衔接的局部连接区；没有把连接处的零距离当作 PDMS 薄区。

对内壁厚度，四个合法穿墙区域先排除 2.0 mm 的局部轴向带。带外最小值 2.0 mm 位于排除带边缘，属于该定义的边界值；中央原生血管到壁面的最小厚度为 10.0 mm。未来 PDMS 占据内腔减去血管芯的空间，诊断体积 361.04 mL。诊断 VTP 只用于理解外形和管腔关系。

## 6. 内部打印支撑

需要内部支撑。竖直基准先进行了关闭支撑、开启支撑的实际切片；其顶部通路代理识别到 28 个阻挡风险邻域，因此继续搜索姿态。选定姿态中，30 个潜在支撑邻域有 30 个找到顶部通路，其中 30 个通过直径 3.0 mm 探针检查；剩余阻挡风险 0 个。顶部通路始终在倒模坐标下计算，侧壁从未被当作移除出口。

实际支撑类型取自当前 Bambu 工艺，为 未确认。根据真实挤出路径在内腔中的分布，盒内支撑估计使用耗材 None g、对应挤出体积 None mm³。这是按 G-code E 值和实际耗材直径、密度计算的估计，不是独立支撑实体体积。

Probe accessibility is a geometry proxy, not a guarantee that real support can be removed. 当前没有独立支撑网格，状态为 BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE；实际拆除仍需要人工审核和实物试验。

## 7. Bambu 实际切片与打印姿态

实际状态：BAMBU_VALIDATION_PENDING。版本 02.07.01.57，打印机 Bambu Lab P1S，喷嘴 0.4 mm，ABS 配置 ['Bambu ABS @BBL P1S 0.4 nozzle']，工艺 0.20mm Standard @BBL X1C。选定切片预计耗时 None s、总耗材 None g。

共检查 52 个确定性姿态，采用显式字典序比较阻挡风险、内部支撑风险面积、顶部可达性、外部支撑、总支撑量、高度和占地；真实切片后用实际支撑挤出量细化同等几何候选的排序。最终姿态 SIDE_y90_yaw0，整体包围盒 59.39 × 107.85 × 94.29 mm，高度 94.29 mm。打印后按照 casting_restore_transform.json 恢复顶部朝上再灌注。

四个端口、主要血管分支、四壁与底板均使用真实非支撑挤出路径作采样保留检查，结果见 tables/slicer_feature_preservation.csv；3MF 内部网格也单独核对刚性变换。真实 Bambu 原生预览保存为 QC/15_bambu_slice_preview.png，补充支撑路径图来自实际 G-code。没有发送、上传或启动打印。工艺阈值来源及未输出文本的真实 --help 探测均保留在 bambu/。

## 8. 当前最重要的风险

1. 顶部可达性仅是几何代理；实际树状支撑可能需要分段切除，必须确认不会牵拉或折断细血管。
2. O3 附近约 1.69 mm 的 PDMS 薄区虽然通过既定回归门槛，仍需实物验证。
3. 侧放可能显著增加盒内壁支撑，需结合真实切片预览判断清理工作量。
4. 尚未验证 ABS 实际打印、PDMS 灌注、丙酮去除和流动实验。

## 9. 最终状态

**NEEDS_ADJUSTMENT**。这是供人工审阅的结果，不是制造验证完成。

受保护文件 1701 个，哈希不匹配 0 个。失败项：['SLICER_DROPPED_VASCULAR_FEATURE', "TypeError: render.<locals>.save() got multiple values for argument 'lines'"]。限制与警告：['BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE: actual toolpaths inspected, no independent support solid; top removal remains a geometry proxy.', 'BAMBU_METADATA_WARNINGS: anomalous first-layer time metadata retained but not used.']。永久测试在 tests/test_s1_6_abs_casting_mold.py，测试结果保存在 test_results.xml。生成入口为 s1-6_abs_casting_mold.py；复运行时使用 --output-root 指向新目录，以保护所有历史结果。
