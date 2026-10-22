# 支撑顺序拆除检查指南

状态：**NEEDS_ADJUSTMENT**。

当前未完整通过时，以下表格只用于核对已计算动作，不能按顺序强行完成剩余拆卸。整棵血管不得作为把手，禁止扭动大组件、沿血管轴线拖拽或横拉分叉。先剪远离血管的连接，再取下 bulk，最后只处理短接触残片；未验证的 bulk/roof 分离不会被写成成功步骤。

坐标为打印坐标系，单位 mm。接近方向指从目标向工具所在空间的向量；取出方向为组件移动向量。检查图 13–15 显示前两步及最后可完成状态。

|步骤|动作及真实 ID|夹持/剪切位置 (mm)|接近或取出方向|附近血管|禁止动作|
|---|---|---|---|---|---|
|1|CUT_STRATEGIC_CONNECTOR SC37|[166.128, 164.12, 16.229]|[-0.0, 1.0, 0.0]|PORT_O2|向血管内推、轴向拉、整体扭动、横拉分叉|
|2|CUT_STRATEGIC_CONNECTOR SC28|[102.521, 142.84, 15.12]|[-0.0, 1.0, -0.0]|PORT_O1|向血管内推、轴向拉、整体扭动、横拉分叉|
|3|CUT_STRATEGIC_CONNECTOR SC25|[109.63, 104.588, 17.02]|[1.0, -0.0, 0.0]|PORT_I1|向血管内推、轴向拉、整体扭动、横拉分叉|
|4|CUT_STRATEGIC_CONNECTOR SC39|[163.403, 153.877, 11.32]|[0.0, -1.0, -0.0]|PORT_O3|向血管内推、轴向拉、整体扭动、横拉分叉|
|5|CUT_STRATEGIC_CONNECTOR SC31|[90.175, 141.513, 15.12]|[-0.0, -1.0, 0.0]|PORT_O1|向血管内推、轴向拉、整体扭动、横拉分叉|

## 必须停止的位置

- **M00**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：['SC12', 'SC03']。禁止把此组件整体当作杠杆拉下。
- **M01**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：['SC22', 'SC14']。禁止把此组件整体当作杠杆拉下。
- **M02**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：['SC03']。禁止把此组件整体当作杠杆拉下。
- **M03**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：[]。禁止把此组件整体当作杠杆拉下。
- **M04**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：[]。禁止把此组件整体当作杠杆拉下。
- **M05**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：[]。禁止把此组件整体当作杠杆拉下。
- **M06**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：[]。禁止把此组件整体当作杠杆拉下。
- **M07**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：[]。禁止把此组件整体当作杠杆拉下。
- **M08**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：['SC22']。禁止把此组件整体当作杠杆拉下。
- **M09**：LONG_CONTACT_COMPONENT_REQUIRES_VERIFIED_BULK_AND_SHORT_ROOF_CUTS；尚未可达的相关剪口：['SC12', 'SC14']。禁止把此组件整体当作杠杆拉下。

尚未完成的剪口：SC12, SC03, SC22, SC14。

本指南没有假设 ABS 的实际剪断力、剥离力或细管承载能力已验证。
