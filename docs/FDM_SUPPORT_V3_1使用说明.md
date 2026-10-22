# s1-8 V3.1：开口倒置实验与实际刀路修复

本阶段保留同一个 `s1-8_vascular_support_planner.py` 入口。输入仍为已冻结的 BG001 RMCA 血管、四个端口、2 mm 五面盒体、SWC 分支映射和当前保存的 X1 Carbon / 0.2 mm / ABS Vessel / 0.10mm Standard Vessel 配置。没有修改原有交互可视化，也没有增加 s1-9。

## 执行方式

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-8_vascular_support_planner.py --stage all
```

程序先验证历史 SHA256，再评估 O000、精确倒置及 casting-local X/Y 两轴的 ±15°、±30°、±45°。所有姿态都检查旋转矩阵与名称；倾角绝不是屏幕或相机轴。完整支撑规划与真实 Bambu 本地切片限定为倒置粗筛 Top 6、必须保留的 INV_0，以及正放 O000 对照。

只要倒置候选通过全部软件门槛，就在倒置目录选择生产审阅候选并停止进入 Phase B。若全部失败，则自动用 O000 进行“实际切片—残差测量—局部修复—再次切片”，最多 5 轮修复；连续两轮未覆盖面积改善不足 1% 即停止。停止并不等于通过。

```bash
# 仅生成矩阵与粗筛
./run_python.sh s1-8_vascular_support_planner.py --stage orientation
# 完成 Top 6、INV_0、O000 的支撑几何
./run_python.sh s1-8_vascular_support_planner.py --stage plan
# 根据已完成两阶段结果更新图表与报告
./run_python.sh s1-8_vascular_support_planner.py --stage report
# 从已有检查点继续；已经触发停滞停止时，只更新结果，不无限增加修复轮次
./run_python.sh s1-8_vascular_support_planner.py --stage resume
# 只读历史 V3，不覆盖原始输出
./run_python.sh s1-8_vascular_support_planner.py --mode v3
# 完整永久回归
./run_python.sh -m pytest -q tests/test_s1_8_vascular_support_planner.py
```

## 路径与输出

以下 `S` 代表 `outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner`。

| 路径 | 内容与用途 |
|---|---|
| `S/inverted_study/orientation_matrix.csv` | 矩阵、开口/盒底法向、实际倾角及名称一致性 |
| `S/inverted_study/surrogate_metrics.csv` | 14 姿态的盒底/盒壁负担、血管/端口负担、首层接触及射线机会 |
| `S/inverted_study/candidates/<姿态>/` | 支撑 STL、完整件、实际 3MF/G-code、拆除模拟、逐层分类覆盖和日志 |
| `S/inverted_study/full_candidate_metrics.csv` | 四类实际覆盖、自打印、拆除、切片和失败门槛 |
| `S/inverted_study/QC_INVERTED/` | 14 张倒置实验 QC 图；抬离失败图不代表可实施动作 |
| `S/v3_1/iteration_00/` | 重新切片的 O000 基线 |
| `S/v3_1/iteration_01/` 等 | 每轮局部修复、拒绝原因、实际覆盖、支撑自身残差、实际切片和 QC |
| `S/v3_1/residual_failure_map.json` | 尚未解决的需支撑层区、支撑自身缺口、受困件及失败门槛 |
| `S/v3_1/actual_coverage_before_after.png` | 基线与最终实际 G-code 覆盖，绿色已托住、红色漏托 |
| `S/v3_1/orientation_then_v3_1_report.md` | 先回答倒置实验，再回答闭环修复的中文结果报告 |
| `S/v3_1/support_orientation_and_closure_summary.json` | 两阶段可供后续开发读取的机器汇总 |
| `S/inverted_study/validation/protected_before.json` | 上游及 V1/V2/V3 历史输出的开始时 SHA256 清单 |

如果 Phase A 直接通过，最终汇总与报告放在 `inverted_study/`，不生成未执行的 Phase B 结果。`production_selected/` 仅在全部软件门槛通过后创建；候选目录下存在 STL 或 3MF 本身不代表已经通过。

## 主要实现

| 模块 | 作用 |
|---|---|
| `support_v3_1_orientation.py` | casting-local 刚体矩阵、法向审计、真实开口射线可见性 |
| `support_v3_1_planning.py` | 复用 V3 分段肋，增加有空气间隙的盒底顶棚支撑 |
| `support_v3_1_release.py` | 固定落板支撑，整模 +Z 抬离，必要时 ±5/10/15° 小角度路径与 FCL 连续扫掠 |
| `support_v3_1_validation.py` | 按实际导出实体分根域；血管、端口、盒底、盒壁分别读取真实刀路覆盖 |
| `support_v3_1_residual.py` | GEOS 差集、残差分组、线宽量化、向下传播与停滞判定 |
| `support_v3_1_repair.py` | 局部屋脊/肋修改、历史受困件/切口尝试，接受前检查拆除与工具通道 |
| `support_v3_1_closure.py` | 从第 0 轮开始的有限次实际切片修复闭环 |
| `support_v3_1_pipeline.py` | 严格先 A 后 B、复用检查点、保护历史、通过后才导出 production |
| `support_v3_1_review.py` | 实际数据 QC、四类覆盖表、中文报告及汇总 |

二维操作使用已有 Shapely/GEOS，没有新增手写多边形布尔内核。支撑 STL 浮点精度问题仅清理新支撑；无法可靠导出的局部支撑会记录为拒绝，实际覆盖门槛仍会显示其缺口。实际连通实体包含内部负向空腔表面，这些空腔不作为独立可拆碎片。

所有“覆盖通过”均以实际 Bambu 刀路为依据。0.20 mm 间隙、短桥、断裂力及 ABS 强度仍未经过物理 coupon 校准；3 mm 工具代理与首层接触仅是几何检查。程序只在本地切片，不上传、发送或启动打印。

## 本次真实运行结果

14 个姿态已完成几何粗筛，6 个倒置候选及 O000 对照均完成真实 Bambu 全模型与支撑单独切片。所有倒置候选仍有软件门槛未通过，因此自动执行 O000 的两轮局部修复；连续两轮未覆盖面积改善不足 1%，按规则停止，状态为 `NEEDS_ADJUSTMENT`。

| 实际刀路指标 | iteration_00 基线 | iteration_02 最终诊断件 |
|---|---:|---:|
| 血管与端口合计覆盖率 | 92.060184% | 92.050749% |
| 未完整覆盖层区域 | 337 | 336 |
| 支撑自身欠支撑区域 | 36 | 33 |
| 受困碎片 | 3 | 3 |
| 工具不可达切口 | 8 | 8 |

最终血管覆盖率为 87.884921%，端口为 94.578559%。O000 的盒底和盒壁没有额外需支撑层区域，表中对应的 100% 不能解释为倒置盒体已经解决。实际 X1 Carbon 切片成功，但仍有浮空提示及覆盖、拆除失败项，未创建正式生产候选，不建议据此直接进行完整模型实物试打。

最终 STL、3MF、逐层残差和日志位于 `S/v3_1/iteration_02/`，只是诊断文件。首选阅读 `S/v3_1/orientation_then_v3_1_report.md`，再查看实际覆盖对照图和 `residual_failure_map.json`。`history/export_guard_failed_iteration_01/` 保留的是网格导出修复前的失败尝试，不应当作正式迭代或生产文件使用。

完成结果文件后，可核对压缩包 CRC、实际 G-code、Auto Support OFF、历史源文件 SHA256 以及永久回归测试：

```bash
./run_python.sh tools/verify_support_v3_1.py
```

验证记录保存在 `S/v3_1/validation/artifact_verification.json`、`pytest.log` 和 `pytest.xml`。文件与代码验证通过，只说明本轮结果可复核，不表示模型已经通过制造门槛。

本轮最终验证：209 项永久测试通过，170 项结果文件检查通过；4887 个受保护文件的 SHA256 差异为 0。测试中保留 1 条既有零体积碰撞测试的数值警告，未导致测试失败。
