# FDM 血管支撑 V3：开发与结果使用

当前 `s1-8_vascular_support_planner.py` 默认执行 **FDM_SEGMENTED_RIB_SUPPORT_V3**。没有新增 s1-9。V1、V2、原生 Bambu Organic 基准和历史测试均保留；V2 入口只读取历史报告。s1-2、s1-3 的可视化未改动。

## 运行入口

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-8_vascular_support_planner.py --stage all
```

它读取当前已冻结的 BG001 RMCA 完整血管、端口及 2 mm 浇注盒，读取本机保存的 X1 Carbon / 0.2 mm / ABS Vessel / 0.10mm Standard Vessel 工艺。只进行刚体变换，比例始终为 1。流程包括 8 个姿态的区间薄肋规划、前 4 个候选的拆除与真实 Bambu 验证，以及中文报告和 20 张 QC 图。首次运行较慢，检测和真实切片均有缓存。不会发送、上传或启动打印机。

分阶段命令：

```bash
./run_python.sh s1-8_vascular_support_planner.py --stage orientation
./run_python.sh s1-8_vascular_support_planner.py --stage resume
./run_python.sh s1-8_vascular_support_planner.py --stage removal
./run_python.sh s1-8_vascular_support_planner.py --stage slice
./run_python.sh s1-8_vascular_support_planner.py --stage report
./run_python.sh -m pytest tests/test_s1_8_vascular_support_planner.py -q
```

`plan` 重新规划几何；`resume` 复用指纹一致的几何并完成未规划候选。修改几何策略后应重新规划，再运行拆除和切片；不要只更新报告来宣称旧刀路对应新几何。`slice` 检查真实挤出路径，单线窄颈丢失时保留之前的候选并自动尝试两线窄颈、三线肋身。不存在物理强度通过的含义。

`--mode v2` 只读取 V2 历史报告。旧 V2 使用说明中的默认命令属于历史版本，以本文为准。`--mode custom` 现在是 V3 的别名。`--mode benchmark` 的新运行输出放在 V3 的 `native_benchmark`，不会覆盖原生历史基准。

## 输入、逻辑与输出

输入仍由 `support_inputs.discover()` 定位，不复制或重建原血管、端口及盒体。检测继续使用既有 `support_detection.py`、`support_custom_detection.py`、`support_freefloating.py` 和 `support_v2_layers.py`，作为冻结检测基础。

规划从连续危险区间出发，不从散点支柱出发。区间关联 SWC 弧长、半径、曲率、切向、端口语义以及真实欠支撑层多边形。薄肋厚度是配置线宽的整数倍；顶部窄屋脊逐层扩宽，工具窗用拱形闭合避免依赖长桥。几何上可切的窄颈与物理上可折断的结构是两回事。

输出根目录：

```text
outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/v3/
```

| 路径 | 内容与用途 |
|---|---|
| `vascular_support_v3_report.md` | 首先阅读；17 个问题、实际结果和未通过的门槛 |
| `vascular_support_v3_summary.json` | 机器可读状态；查看 `production_selected`，不要只看切片成功 |
| `candidates/Oxxx/`、`candidates/Nxxx/` | 各姿态的实际几何、冻结检测证据、切片与拆除记录 |
| `geometry/` | 当前最佳诊断候选的几何副本；不代表生产通过 |
| `tables/` | 姿态、区间、薄肋、屋脊、切割点、碎片、拆除路径、真实覆盖与警告 |
| `QC/` | 20 张检查图；16 为实际覆盖，17 为未覆盖区域，11 为真实肋中面纵剖 |
| `sensitivity/line_widths/` | 1、2、3 条线结构的实际 Bambu 保留性试验；无实物试验 |
| `sensitivity/full_rib_variants/` | 全薄肋的线宽和 0.15/0.20/0.25 mm 间隙敏感性模型，仅供诊断 |
| `validation/` | 历史/源文件 SHA256、V2 同门槛复算、实际工艺和缓存指纹 |
| `history/` | 修正前的 V3 候选；V1/V2 历史仍在原目录 |
| `production_selected/` | 仅当全部软件门槛通过才创建；不存在时不得将诊断模型当作已验收生产件 |

各候选的重要文件：

- `mold_without_support.stl`：相同刚体姿态下的冻结完整倒模件。
- `custom_support_v3.stl`：显式支撑；`mold_with_custom_support_v3.stl`：倒模件与支撑的合并件。不要在切片中将它们重复叠放。
- `support_ribs.stl`、`support_roofs.stl`、`support_cut_necks.stl`、`support_exterior.stl`：分部检查。没有辅助柱时 `support_pillars.stl` 是空几何，清单注明原因。
- `print_transform_v3.json` / `casting_restore_transform_v3.json`：打印与冻结浇注坐标之间的互逆变换。
- `X1C_project.3mf`：真实切片归档，自动支撑关闭；软件门槛未通过时仅供 Bambu 审阅。
- `bambu_validation.json`：实际覆盖、构件保留、支撑自身实际挤出检查、警告、归档 CRC 和打印估时。
- `support_removal_plan.json`：切断已设计窄颈后的完整碎片拆除顺序。内部碎片必须整体离开实际浇注开口，尚未拆除的碎片继续作为障碍。
- `tool_access.json`：独立的 3 mm 圆柱探针检查。它不等于真实钳子/剪刀已验证可操作。
- `exported_geometry_audit.json`：源模型保留、STL 闭合及盒体接触是否仅发生在批准根部。

## 判断结果

覆盖必须是实际 Bambu 刀路对每个危险层多边形的完整承托，不是设计屋脊的投影覆盖。V3 不以未经校准的桥接补足缺口。STL 层检查与实际挤出层检查分别报告；前者为零而后者非零时仍不通过。

`NEEDS_ADJUSTMENT` 意味着仍有软件失败；`SLICER_VALIDATION_PENDING` 意味着还缺实际切片证据。只有全部软件门槛通过才使用 `READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED`，仍然不代表生产工艺或材料已实验验证。

**NO PHYSICAL COUPON AVAILABLE**。屋脊间隙、桥接能力、断裂/切割力和 ABS 薄肋强度均 **NOT PHYSICALLY VALIDATED**。默认 0.20 mm 是未校准回退间隙；不能依据切片成功判断它是最佳间隙。

## 后续开发位置

| 模块 | 职责 |
|---|---|
| `support_v3_pipeline.py` | 姿态筛选、缓存、阶段调度、真实单线失败后的升级及门槛 |
| `support_v3_planner.py` | 冻结检测调用、层多边形到弧长区间、薄肋规划与拒绝原因 |
| `support_v3_geometry.py` | 薄肋、逐层屋脊、拱形工具窗、根部窄颈、有限切割截面 |
| `support_v3_removal.py` | 3 mm 工具探针、切断后的碎片、按顺序的连续扫掠检查 |
| `support_v3_validation.py` | 真实 G-code 支撑路径隔离、覆盖、自支撑、细构件保留、警告定位 |
| `support_v3_sensitivity.py` | 可复核的全肋线宽/间隙敏感性几何 |
| `support_v3_verify.py` | 导出 STL、批准根部、源模型保留审计 |
| `support_v3_review.py` | 中文报告、V2 同标准复算和 20 张 QC 图 |
| `support_v3_fixtures.py` | 13 类永久合成几何；文件位于 `tests/data/support/v3/` |

测试仍统一放在 `tests/test_s1_8_vascular_support_planner.py`。不要删除旧 V1/V2 测试，也不要修改可视化来隐藏漏支撑、受困碎片或警告。
