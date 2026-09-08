# Mirheo starter

本目录同步自 `/home/lzy/projects/mirheo_starter`，保留现有源码、配置、测试、
可视化入口和环境记录的原始字节。该子工程使用独立的 Python 环境。

## 代码入口

| 内容 | 入口与说明 |
| --- | --- |
| 已验收血管形状、管壁与端口标签迁移 | `py_scripts/import_vessel_geometry.py`、`py_scripts/vessel_geometry/`；[几何核查说明](test_code/README_vessel_geometry.md) |
| 旧 Musubi 工况冻结、SI/DPD 单位与有限端口边界设计 | `py_scripts/prepare_fluid_physics.py`、`py_scripts/fluid_physics/` |
| 有持久预算与缓存的小规模 DPD 标定 | `py_scripts/calibrate_dpd_fluid.py`；[物理与标定说明](test_code/README_fluid_physics.md) |
| 原生 DPD/SDPD 纯液体对照、共享剩余预算与候选独立判定 | `py_scripts/compare_dpd_sdpd.py`、`py_scripts/fluid_comparison/`；[对照测试说明](test_code/README_fluid_model_comparison.md) |
| 离线交互核查页及浏览器检查 | `test_code/review_*.py`、`test_code/check_*_browser.cjs` |
| 原有 RBC 入门示例与结果查看器 | `py_scripts/hello_rbc_flow.py`、`py_scripts/view_rbc_flow.py` |
| 原有 WSL 激活与官方算例启动脚本 | `scripts/` |

源码、配置、测试与环境记录按本地版本同步；仓库另附本说明、CPU 依赖清单、
忽略规则和上游许可证副本。
原始数据、运行结果、GPU 预算账本、截图、生成的 HTML、虚拟环境、Mirheo 第三方
源码和编译库留在本地，由 `.gitignore` 排除。

## CPU 环境与测试

已验证环境为 WSL、Python 3.12.3。以下命令从仓库根目录开始，创建独立环境：

```bash
cd mirheo_starter
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-cpu.txt
mkdir -p test_code/outputs/fluid_physics
```

98 项合成/CPU 测试可在没有真实几何数据或 Mirheo 库的情况下运行：

```bash
.venv/bin/python -B -m unittest \
  test_code.test_vessel_geometry \
  test_code.test_vessel_geometry_migration \
  test_code.test_fluid_physics.UnitsTests \
  test_code.test_fluid_physics.StatisticsTests \
  test_code.test_fluid_physics.ProtectionBudgetTests \
  test_code.test_fluid_physics.EosReportTests \
  test_code.test_fluid_model_comparison.EvidenceTests \
  test_code.test_fluid_model_comparison.SharedBudgetTests -v
```

另有 26 项 `LegacyTests` 和对照模块的 `ModelTests` 读取配置固定的真实来源文件、
已冻结物理工况和单位包。具备这些外部输入后，可以运行全部 124 项测试；
这些测试均不启动 GPU：

```bash
.venv/bin/python -B -m unittest \
  test_code.test_vessel_geometry \
  test_code.test_vessel_geometry_migration \
  test_code.test_fluid_physics \
  test_code.test_fluid_model_comparison -v
```

2026-09-08 同步核查：在独立检出中使用原 WSL 的 Python 环境和配置指定的真实
外部输入，全部 124 项测试通过；Python、JavaScript 和 Shell 语法检查通过。
64 个同步文件及其文件权限与原工程一致。本次增加 13 个对照模块相关文件，
更新 3 个公共模块及本说明；此同步核查没有重新运行 GPU 或浏览器交互。

## 外部数据和运行环境

三个 YAML 配置保留原 WSL 绝对路径、指定 run ID 和来源哈希。几何输入来自
`bloodflow_starter` 的已验收包，旧流体工况来自本地 `ulm_3D_vascular/outputs/`，
第二阶段也引用本地 `mirheo_starter/data/geometry/`。这些输入不包含在此次代码同步中；
缺失或哈希冲突时程序会明确阻断。DPD/SDPD 对照还引用已冻结的
`data/fluid_physics/result_1b008103675eb6ff/` 和旧 DPD 原始任务。
详细来源及许可证说明见各模块的 `SOURCES.md`。

GPU 环境记录在 `metadata/`：Mirheo 提交
`8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`、CUDA 12.6、GCC 12、Open MPI 4.1.6、
单精度构建，以及已有的 CUDA 兼容补丁。`CMakeCache.txt` 和库哈希记录原机器的
实际构建；`python_packages.txt` 是部署初期快照，当前 CPU 依赖见本目录依赖清单。
上游 Mirheo 的 MIT 许可证副本为 [MIRHEO_LICENSE.txt](MIRHEO_LICENSE.txt)。

原有 `scripts/activate_mirheo.sh` 与 `scripts/run_official_case.sh` 固定使用
`/home/lzy/projects/mirheo_starter`。迁移 GPU 运行位置时，需先配置相应的本地环境、
Mirheo 源码及编译库、实际构建记录和数据路径。继续原标定 campaign 时，必须保留其
完整 `runs/fluid_calibration/dpd_round1_20260908/`，包括原预算账本及任务结果。
继续 DPD/SDPD 对照时，还须保留
`runs/fluid_calibration/shared_budget_pool.json` 和
`runs/fluid_model_comparison/dpd_sdpd_remaining_20260908/`，确保旧任务、本轮任务和
未结束预约累计使用同一份 3600 秒授权；缺少账本或未注册的 campaign 会阻断执行。
Windows 浏览器检查脚本也依赖原机器已有的 Chrome/Node 路径。

## 当前科学状态

血管几何及标签迁移已完成。首轮 DPD 标定及本轮 DPD/SDPD 对照的科学结果均为
`PARTIAL`，没有合格液体参数，`selection=null`。SDPD 基准实测黏度点估计更接近
目标，但温度、统计充分性、敏感性及含不确定性的压力范围尚未通过初筛。
短任务执行成本不能用作合格目标液体的效率排名。
边界设计中的原生接口缺口仍标为 `REQUIRES_EXTENSION`。
已有 RBC 文件属于入门示例。此次代码同步与 CPU 测试不增加 GPU 标定任务，
也不表示已完成真实血管流动、SDF 或生产参数验收。
