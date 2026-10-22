# 支撑感知血管 ROI 与自适应盒体

## 当前 2 mm STL 导出

最新厚度调整把已选中 `candidate_0008` 的侧壁和底板改为 **2 mm**。内腔及血管、端口坐标保持，外壁向内退回；盒体外尺寸为 **82.17 × 74.53 × 31.01 mm**。这是对已选中模型的独立 STL 修改，下文完整筛选、支撑比较及原 3MF/G-code 仍记录之前的 3 mm 版本。

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-7_support_aware_roi_60x50.py --stage box-stl --box-thickness-mm 2
```

新 STL 位于默认输出根目录的 `casting_mold_2mm/`：`final_abs_casting_mold_2mm.stl` 是一体模型，`final_abs_casting_mold_2mm_oriented.stl` 已居中落床，`intact_open_top_box_2mm.stl` 是单独盒体。对应厚度参数和检查结果保存在该目录的 `resolved_box_config.json`、`geometry_summary.json` 中。本命令不重选 ROI、重建端口或调用切片器。原切片文件不能作为新厚度版本的切片结果使用。

## 原始支撑感知筛选流程

本阶段入口为项目根目录的 `s1-7_support_aware_roi_60x50.py`。它从冻结的 BG001 RMCA MeVO 中心线重新选择制造子树，生成血管、适度长度的人工延长段和贴合几何范围的五面倒模盒。原有 candidate 0、15 的已验收文件保持只读；本阶段的 `candidate_0001` 等编号表示新的血管子树，不是旧打印姿态编号。

## 输入与处理逻辑

输入由既有 manifest、节点映射和来源哈希定位，包括原始 BraVa BG001 数据、冻结的 `semantic_refined` 与 `d075` MeVO 范围，以及既有铸造坐标变换。程序不重新执行 TopBrain 配准、NN 标签迁移或解剖分类，也不通过整体缩放凑目标尺寸。

处理顺序为：确定性子树枚举 → 中心线支撑风险比较 → 前八个候选的 VascularMD 表面建模 → 新端口和自适应盒体 → 顶部支撑可达性检查 → X1 Carbon 0.2 mm 的实际本地切片 → 中文报告和审核图。端口设计失败会保留原因，并由后续合格候选补充，直至三个完整模型完成切片或候选用尽。

60 × 50 mm 指中央血管在固定铸造 XY 平面的目标范围，不包含人工延长段和盒体。最终以真实表面尺寸核验；达不到允许的 ±3 mm 范围时，报告 `ROI_SIZE_TARGET_NOT_MET` 和 `NEEDS_ADJUSTMENT`，保留实际最佳结果。

## 延长段与盒体如何自适应

参数集中在 `config/support_aware_roi_60x50.yaml`：

| 参数 | 当前含义 |
|---|---|
| `ports.outside_stub_mm` | 墙外接口直段长 6 mm |
| `ports.preferred_maximum_extension_mm` | 优先让每根完整人工延长段不超过 35 mm；仍须通过真实端面、弯曲与碰撞约束 |
| `box.adaptive_unconnected_faces` | 按实际中央血管和人工延长段收紧未连接端口的面 |
| `box.unconnected_face_gap_mm` | 自由侧壁、底面至实际几何包络的目标净距 5 mm |
| `box.top_headroom_mm` | 顶部开放边缘至最高几何的目标高度余量 5 mm |
| `box.wall_thickness_mm` / `bottom_thickness_mm` | 原始筛选与切片比较采用 3 mm；当前独立 STL 导出使用上方 `box-stl` 的 2 mm 参数 |
| `box.margin_*` | 初始端口搜索空间；开启自适应后不再是最终固定盒体余量 |

端口只穿 ±X、±Y 侧壁。每次选路优先减少过长端口，再减少最长段与总长度。自由面收紧时保留端口中心至邻近盒边至少 8 mm 的材料空间；这会使个别面的实际净距大于 5 mm。实际净距、限制原因和 FCL 表面距离均写入 `tables/adaptive_box_faces.csv`。血管与非连接血管、端口之间的绝对间距下限仍为 1.6 mm。

修改制造参数后，`--stage evaluate` 会核对制造输入签名，将受影响的旧模型、姿态和切片移入相应候选的 `history/`，再重新生成。历史 JSON 保留原始绝对路径，重定位关系见 `relocation.json`。冻结的上游血管结果及已归档的比较工艺不随盒体参数变化而重算。

如果人工端口在顶部可达性检查中失败，程序会排除该端口当前的侧壁走向，最多重新选路两次，并在该候选的 `port_route_review.json` 保存具体风险与允许侧壁。这一反馈绑定该候选血管 STL 的 SHA256；换了血管数据就不能沿用。原生血管不会被这项反馈移动或修改。本次三出口候选的 O2 曾在 +Y 走向出现困陷风险，因而改向 -X 重新检查。该记录不会覆盖其他候选的选路，也不放宽间距、弯曲或真实端面连接要求。

## 运行

```bash
cd /home/lzy/projects/vascular_printing

# 使用新输出目录运行完整流程
./run_python.sh s1-7_support_aware_roi_60x50.py \
  --output-root outputs/support_aware_trial_01

# 已有默认目录修改制造参数后，重新生成端口、盒体和切片
./run_python.sh s1-7_support_aware_roi_60x50.py --stage evaluate

# 更新图、表、报告及交付模型副本
./run_python.sh s1-7_support_aware_roi_60x50.py --stage report

# 永久回归测试
./run_python.sh -m pytest -q tests/test_s1_7_support_aware_roi_60x50.py
```

`prepare` 只接受新目录；已有结果可用 `model`、`baseline`、`baseline-upright`、`evaluate`、`report` 分阶段继续。真实切片需要当前本地 Windows Bambu Studio 和已配置工艺。程序仅执行本地切片，不向打印机发送任务。

## 默认输出

默认根目录为：

```text
outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/support_aware_roi_60x50/
```

| 路径 | 内容与用途 |
|---|---|
| `roi_60x50_report.md` | 中文研究报告，先读实际结果与限制 |
| `roi_60x50_summary.json`、`validation_results.json` | 实际尺寸、状态、来源、几何和切片检查 |
| `selected/` | 选中子树 SWC、补偿 SWC、血管 STL/VTP、含新端口的血管芯 |
| `casting_mold/` | 独立盒体、完整一体模型、摆放后的 STL、实际 3MF 和 G-code |
| `tables/port_length_comparison.csv` | 各端口调整前后长度及实际间距 |
| `tables/adaptive_box_faces.csv` | 各面实际净距、收紧幅度、连接或边缘限制 |
| `tables/roi_candidate_summary.csv` | 所有子树，以及进入制造评估者的最终状态 |
| `QC/08_new_ports.png`、`QC/09_port_interfaces.png`、`QC/10_new_box.png` | 延长段、真实端面连接、自适应盒体检查图 |
| `QC/13_baseline_vs_new_support.png` | 相同工艺重新切片后的旧基线与新结果比较 |
| `candidates/` | 每个候选的输入、VascularMD 输出、失败记录、制造结果和姿态切片 |
| `baseline/` | 已验收旧模型的重新切片证据，含当前姿态和正立补充对照 |
| `bambu/evaluation_profile/` | 比较使用的本机打印机、ABS、工艺配置及 SHA256 |
| `history/` | 本阶段调整前结果；不是当前交付入口 |
| `test_results.xml`、`protected_before.json`、`implementation_manifest.json` | 测试、上游保护快照、新实现版本证据 |

盒内支撑材料依据实际 G-code 正挤出量计算。接触面积是支撑路径邻域内的三角面面积估计；顶部探针可达性也是几何检查，不能替代真实打印后的拆支撑、灌注和溶解实验。
