# BG001 RMCA：P1S 0.4 mm / ABS 血管芯与支撑

状态：**NEEDS_ADJUSTMENT**。只生成血管芯与显式支撑，不包含外模或 C 插片。

真实切片使用 **Bambu Lab P1S / 0.4 mm / ABS**，当前保存工艺为 `0.20mm Standard @BBL X1C`。
这个官方共享工艺名称包含 X1C，但其兼容列表明确包含 P1S；机器启动/结束 G-code 使用 P1S 机器配置。
层高 0.2 mm、常规及支撑线宽 0.42 mm，自动支撑关闭。
血管原始 STL 逐字节保留，打印姿态沿用已有 O002。支撑根据新工艺重新生成，未放大或增粗血管。

## 文件

`selected/` 包含原始血管、打印坐标血管 STL、独立支撑 STL、可编辑工程 `_project.3mf`、
原生真实切片 `_sliced.3mf`、`plate_1.gcode` 和坐标变换。
在 Bambu Studio 中选择“打开项目”，打开 `_sliced.3mf` 可看真实刀路，打开 `_project.3mf` 可编辑。
不要把两个 STL 分别自动落板；二者共用同一打印坐标，血管通过显式支撑承托。

`QC/` 保存组合、独立血管、独立支撑图和原生切片预览；`candidates/` 保存每次生成与切片日志。
`regeneration_summary.json` 保存配置、门槛、源保护和交付文件 SHA256。

## 本次实际检查

- 真正缺承托：9 个；最大缺失弧长 0.626754 mm。
- 支撑宏观失败：0；浮空血管起始区域：0。
- 实际血管刀路特征保留：True。
- 几何生成阶段遗留路由错误：3 个，逐项记录于摘要。
- 薄连接验证：True；完整安全拆除路径：False。
- 未通过门槛：TRUE_MISSING_SUPPORT_ZERO, NO_DANGEROUS_SPAN, COMPLETE_SAFE_REMOVAL, SHORT_FINAL_CONTACT_PIECES, CRITICAL_WARNING_ZERO, GEOMETRY_ROUTING_ERRORS_ZERO。

切片成功不代表实物制造验收。只有全部软件门槛通过时才标为可首次试印；未通过时这些文件是诊断候选。
ABS 实物表面、拆卸力与血管损伤尚未实测。本次没有向打印机发送文件或启动打印。
历史 X1 Carbon 输出、外模、C 插片及用户 GUI 配置均未改写。
