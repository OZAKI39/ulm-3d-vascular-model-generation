# BraVa 人体脑血管流场与微泡工程

本工程研究 **BraVa BG001 RMCA 四端口人体血管，入口 18 mL/min**。**本目录是唯一日常开发与 Git 根目录，长期分支为 `dev/human-vessel-flow`。** 代码、输入、旧结果和必要诊断证据已独立整理在这里。小鼠血管工程另在 `new_mouse_vessel_flow`，没有把它的血管、流场或粒子队列合并进来。

先看[本次清理报告](docs/CLEANUP_REPORT_ZH.md)、[直接 Git 开发方式](docs/GIT_WORKFLOW_ZH.md)和[当前方法清单](docs/CURRENT_METHODS_ZH.md)，或打开[结果图表与四视角动画](results/reviews/migration_final/OPEN_RESULTS.html)。此前迁移过程见[中文迁移报告](docs/MIGRATION_REPORT_ZH.md)，详细来源在[模块对应表](docs/MODULE_MAPPING.csv)、[函数与类对应表](docs/FUNCTION_CLASS_MAPPING.csv)和[冻结文件索引](docs/CASE_INDEX.md)。本地冻结资料的当前身份清单为 `configs/active_data_manifest.json`。

## 当前可以做什么

读取两套原有流场、核对原始壁面数据、分析旧 1565 个微泡、运行轻量测试、复查已保存的诊断数组、生成核查图表。GPU 和原生工具的本次实际验证范围详见 [验证汇总](validation/DELIVERY_STATUS.json)；不能将软件加载成功等同于新流场精度通过。

本次没有重新计算正式 CFD，也没有运行新固定时间窗的大批量微泡。原来失败、停止或尚未执行的科学检查仍然保留原状态。

## 从哪个文件读起

入口是 [main.py](src/human_vessel_flow/main.py)。按下面顺序阅读即可跟随工作流程；具体数值算法仍调用已经保留的原实现。

| 文件 | 负责什么 | 主要单位或限制 |
|---|---|---|
| [geometry.py](src/human_vessel_flow/geometry.py)、[transforms.py](src/human_vessel_flow/transforms.py) | 读取打印工程提供的已封闭血管，识别四个端口，转换坐标 | 源几何 mm，计算 m；candidate 0/15 是同一物理解的两个姿态 |
| [meshing.py](src/human_vessel_flow/meshing.py)、[mesh_audit.py](src/human_vessel_flow/mesh_audit.py) | 明确调用 SimVascular/TetGen；检查四面体质量、壁面与端口 | 旧权威网格保留；本次没有生成新网格 |
| [boundary.py](src/human_vessel_flow/boundary.py) | 用中心线长度和半径估算阻力，再设计新算例出口压力 | 入口 3e-7 m³/s；压力 Pa；不运行校准 CFD |
| [fem.py](src/human_vessel_flow/fem.py)、[native_runner.py](src/human_vessel_flow/native_runner.py) | 准备 XML，校验软件身份，显式调用原生求解器并检查日志 | 密度 1056 kg/m³；黏度 0.00345312 Pa·s；数值求解仍在 C++/Fortran 中 |
| [flow_field.py](src/human_vessel_flow/flow_field.py)、[wss.py](src/human_vessel_flow/wss.py)、[postprocess.py](src/human_vessel_flow/postprocess.py) | 读取速度、压力，检查网格对应关系，计算流量与壁面剪切 | 速度 m/s；梯度 1/s；原始面片 WSS Pa；展示节点值不能替代原始面片值 |
| [particles.py](src/human_vessel_flow/particles.py)、[population_v1](src/particle_3d/population_v1) | 载入冻结场，准备固定灌注任务，计算微泡接触和运动，保存状态与终态 | dt=0.0005 s；FP64；保留原 CPU 回退、GPU 内核和粒径分布 |
| [validation.py](src/human_vessel_flow/validation.py) | 只读核对几何、流场、WSS、全部旧状态块与原科学判定 | 数据一致性通过不等于局部物理精度通过 |
| [visualization.py](src/human_vessel_flow/visualization.py) | 从旧结果生成核查图表，链接原有动画 | 显示筛选不改变原始轨迹或账本 |

`src/particle_3d` 保留现有包名及真实依赖中的 77 个模块。部分文件名带早期 particle 阶段编号，是当前算法实际使用的函数，不能仅凭名称删除。`vascular_validation` 和 `flow_solver_support` 是项目内部唯一的一套必要辅助实现。

受保护的共享模块还含可选 RBC 和旧 LAMMPS 辅助接口，当前 BraVa 微泡入口没有启用它们。相关旧契约放在 `docs/provenance/inherited_contracts`，不作为 BraVa 的活动物理配置；没有为这些旧接口复制小鼠算例。

## 运行环境与已执行的命令

Python 环境为本目录 `.venv`，无需设置 `PYTHONPATH`，无需激活旧环境。第三方原生软件统一在 `/home/lzy/projects/new_vessel_flow_solvers`，详见 [环境说明](docs/NATIVE_ENVIRONMENT_ZH.md)。

```bash
cd /home/lzy/projects/new_human_vessel_flow
./.venv/bin/human-flow status
./.venv/bin/human-flow native-status
./.venv/bin/python -B -I scripts/check_independence.py validation/新的独立性检查目录
```

本次只读检查曾执行以下命令，结果和日志已保存。`--output` 必须是一个**尚不存在**的目录；复查时请换一个新名字：

```bash
./.venv/bin/human-flow verify --output validation/read_only_v1
./.venv/bin/human-flow plot --output results/reviews/migration_final
bash scripts/run_small_tests.sh validation/small_tests
```

现成的报告无需重跑即可查看。测试使用小型人工场景，不会启动旧队列或正式群体。GPU 核查只读取原 BraVa 场，不产生新轨迹。

## 五类结果请分开理解

| 内容 | 位置 | 应怎样理解 |
|---|---|---|
| 历史正式压力出口 CFD | `cases/balanced_pressure_final` | 旧微泡依赖的冻结场；出口压力经过历史校准，各约分流 1/3 |
| 等流量校准 CFD | `cases/calibration_equal_split` | 出口流量被强制设置；用于校准，不能称自然分流预测 |
| 旧微泡群体 | `results/population/realization_0000_contactfix` | 1565 个接受微泡，458/527/437 出流，141 观察截止，2 数值失败 |
| 新固定灌注协议 | `configs/population/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json` | 灌注 1.240815660 s，再观察 10 s；尚无新正式大批量结果 |
| 研究诊断 | `diagnostics/` | Gmsh V2/P1 局部场检查未通过；P2/P1 停止；V3 仅提案 |

旧群体按出口覆盖情况停止灌注，实际总观察时间是 13.240815660385584 s；新协议总时间是 11.240815660 s。旧结果没有被改名为新协议结果，`scientific_integrity_pass=false` 原样保留。

## 以后建立新算例

配置索引在 `configs/cfd/cases.json`、`configs/cfd/native.json`、`configs/population/environment.json`。几何参数与原始端口信息集中在 `inputs/geometry`，不另建重复几何配置。

`prepare-case --output <新目录>` 只准备中心线阻力压力算例；`prepare-population --output <新目录> --backend cuda` 只准备固定时间窗任务。两者都不会自动启动计算。真正计算需明确调用 `solve --case <新目录>` 或 `run-population --output <新目录>`，本次没有执行这些正式入口。

新计算的输出必须使用独立目录；不能覆盖 `inputs/`、`cases/` 或旧 `results/population/`。若以后更换流场，需同时维护网格、边界、场 SHA 和新协议身份，不能只替换一个文件。原后处理仍保留历史分流与稳态门槛；中心线估算压力尚未通过实际三维分流验证，失败时应保留并审查，不自动转回 CFD 校准。

## 保留的科学问题与旧目录

参考管的 `REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC`、`MIXED_FIELD_ERROR`、`P1_FIELD_TRAJECTORY_IMPACT_LARGE` 均保留。P2/P1 为 `STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT`，曲面压力回填也仍有失败。微泡异常速度、近壁场误差及正式轨迹收敛需要后续科学研究，工程整理没有解决或掩盖这些问题。

2026-10-23 已清除两处经过身份核对的历史 GitHub 同步工作区，停用旧同步入口，并直接在本目录建立 Git。原始 BraVa 工程、共享依赖和服务器工作区保留；本轮未登录服务器。需要人工确认的原工程和备份见[当前待确认清单](docs/cleanup_20261023/PENDING_REVIEW.csv)。GitHub 发布时排除的大文件不能自动视为可删除文件。
