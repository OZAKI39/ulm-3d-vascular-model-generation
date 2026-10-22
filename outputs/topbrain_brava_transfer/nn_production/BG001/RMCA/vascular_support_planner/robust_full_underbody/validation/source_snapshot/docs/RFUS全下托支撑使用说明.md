# RFUS 全下托支撑的开发与使用

RFUS 继续使用 `s1-8_vascular_support_planner.py`，固定 O000（底板朝下、开口朝上）。输入是当前冻结的 BG001 RMCA 血管、SWC、四端口、2 mm 盒体和本机保存的 X1 Carbon 0.2 mm 工艺。既有可视化与输入几何不修改。

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-8_vascular_support_planner.py --stage all
```

程序首先核对历史 SHA256 和真实 V3/V3.1 报告，然后沿全部血管及端口中心线采样。每个下表面点都保留在目标集合中，支撑生成不受旧稀疏风险掩膜控制。曲率、高度和分叉决定分段长度；几何避让依次尝试侧移、双帘与局部竖向结构。每段设计独立底颈和抓取柄。

2/3 条线宽的帘与 1/2 条线宽的底颈形成四组初始候选。若实际刀路丢失薄特征，只对对应部位增加一条线宽；若段不可取出，进行有限的缩短分段、移动分缝与抓取点调整，再重新切片。结果未通过时保留 `NEEDS_ADJUSTMENT`，不会自动改动血管或盒体来制造通过结果。

## 入口与配置

```bash
# 只生成 RFUS 几何
./run_python.sh s1-8_vascular_support_planner.py --mode rfus --stage plan
# 从新模式自己的检查点继续
./run_python.sh s1-8_vascular_support_planner.py --mode rfus --stage resume
# 根据已完成真实结果更新报告和 QC
./run_python.sh s1-8_vascular_support_planner.py --mode rfus --stage report
# 只读历史报告
./run_python.sh s1-8_vascular_support_planner.py --mode v3.1
```

配置位于 `config/robust_underbody.yaml`。线宽和层高来自当前保存的 Bambu 配置，不以提示词里的示例数值覆盖。默认屋脊间隙 0.20 mm，标记为未校准；帘和颈尺寸使用线宽与层高倍数。新的代码模块均以 `support_rfus_` 命名，历史 `support_v*` 模块保持原样。

## 输出在哪里

以下 `S` 指 `outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner/robust_full_underbody/`。

| 路径 | 输入、作用与输出含义 |
|---|---|
| `S/robust_underbody_config.json` | 实际使用的工艺、工程参数和固定姿态 |
| `S/underbody_corridors.csv` | 所有分支及端口的中心线、弧长、半径、切向、曲率、下表面样本 |
| `S/candidates/RFUS_C2_R2_N1/` 等 | C 是帘线数，R 是顶脊线数，N 是底颈线数；保存独立候选的完整证据 |
| 候选后缀 `_P1` / `_RMS` | 局部薄特征保留修复 / 包含分段拆除修复的有限补充候选 |
| `S/geometry/` | 当前选中结果的支撑、顶脊、帘、柄、颈、原模和合并模型 STL；先看其中 README |
| `S/underbody_coverage.csv` | 每个弧长样本在真实 G-code 中是否有下托、是否由已有模型直接承托 |
| 候选的 `actual_required_layer_coverage.csv` | 独立的完整层多边形覆盖检查；与中心线覆盖不是同一百分比 |
| `S/segment_removal.csv` | 工具可达、受困状态、取出顺序与连续扫掠路径 |
| `S/travel_near_miss.csv` | 真实非挤出移动对已打印血管、屋脊、帘和柄的几何风险代理 |
| `S/QC_ROBUST/` | 提示词要求的 20 张人工审阅图；图 14 显示真实覆盖，图 20 显示分段与拆除状态 |
| `S/method_comparison.csv` | Bambu Organic、V3 稀疏肋与 RFUS 在统一完整弧长指标下的比较 |
| `S/robust_full_underbody_report.md` | 开头以中文回答 19 个问题，再列真实成本和方法边界 |
| `S/robust_full_support_summary.json` | 后续代码应读取的状态、选择、指标和文件入口 |
| `S/validation/` | 历史 SHA256、真实输入读取记录、永久测试、真实水平管切片和运行日志 |
| `S/history/` | 本轮开发调试中间件，不能作为正式生产件 |

每个候选的 `bambu/` 保存实际执行命令、stdout/stderr、3MF 和 G-code。`X1C_diagnostic_full_support.3mf` 是本地诊断切片；只有所有软件门槛通过后，才导出 `production_selected/X1C_robust_full_support.3mf`。没有上传、发送或启动打印的步骤。

## 如何判断是否可进入人工审阅

先读 summary 的 `status` 和 `production_selected`。STL 可打开、3MF 能切片，都不代表覆盖或拆除通过。中心线全弧长覆盖、需支撑层覆盖、自打印、实际分段拆除、工具通道、模型特征保留和原生警告是独立检查。

覆盖使用真实 G-code 的逐段 `LINE_WIDTH`，包括变线宽壁；不使用名义 STL 给覆盖率加分。移动近失检查考虑 G-code 顺序，只检查当时已经挤出的材料。3 mm 工具代理与喷嘴圆柱代理都不是实物操作力或完整喷头碰撞证明。

通过全部软件门槛后，允许的状态为 `READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED`。屋脊间隙、底颈行为、帘/柄强度和手工移除力仍需实物校准。

## 回归与真实测试

```bash
./run_python.sh -m pytest -q tests/test_s1_8_vascular_support_planner.py
./run_python.sh tools/verify_support_rfus.py
```

验证器检查历史未变、所需文件、20 张 QC、STL、真实 3MF 的 CRC 和 Auto Support OFF，并运行永久测试；独立水平圆管如已有完整真实切片证据则校验复用，否则执行切片。该 30 mm 圆管试件使用 2 mm 试验底板和一根远离血管投影的竖柱，使屋脊空气间隙所在高度不形成整个试件的空层。底板和竖柱不计入其下托覆盖，稀疏对照也执行真实切片，须从归档 G-code 测出中间约 20 mm 缺口。它只在 RFUS 新目录写日志，不覆盖旧研究结果。Bambu 调用共用本机暂存目录，请顺序运行，避免同时启动多个切片进程。
