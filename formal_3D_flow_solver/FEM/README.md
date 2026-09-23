# Formal 3D flow solver — FEM

当前完成 Stage 0 审计、Stage 1 体网格管线和 Stage 2 Stokes 核心验证。WSL 的本目录是唯一源码来源。
两个参考工程 `/home/lzy/projects/ulm_3D_vascular` 与
`/home/lzy/projects/ulm_microbubble_traj_gen_2D` 只读，既不导入执行，也不写入。

已实现的核心模型为 SI 单位的 3D steady incompressible Newtonian Stokes：
P2 三分量速度、P1 压力、一个全局 Real 流量约束自由度；PDMS 壁无滑移；
入口仅给总流量，出口大气表压以自然牵引实现。全局乘子必须保存。
Real 采用 block/MixedFunctionSpace 或显式 block assembly，不放入普通 mixed_element。
不含 RBC、微泡、非牛顿黏度或脉动。

Stage 1 的几何、质量与人工审核条件见 [报告](reports/stage01/REPORT.md)。
两个开发网格保存在 `outputs/stage01/{coarse,medium}/mesh/`；`fluid.xdmf` 必须与同目录
`fluid.h5` 一起使用，包含 mesh、cell_tags 和 facet_tags。边界完全冻结，没有表面重网格。
Stage 0 的历史审计结论见 `reports/stage00/REPORT.md`，不重写其输入或结果。
环境报告保留探测前状态；新独立环境的状态另行保存。

```bash
python3 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -e '.[audit]'
.venv/bin/python -B -m pytest -q
# 历史 Stage 0/1 证据已冻结，不要重写其报告或测试 XML。
```

远程配置在 `configs/remote.json`；探测后实际目录在 `remote/connection.local.json`，
每次调用记录该目录。SSH 私钥始终留在 WSL。所有同步均增量且无删除参数；
只同步本项目白名单目录；取回仅日志/输出/报告，绝不取回远程源码覆盖 WSL。
每次执行前验证 WSL 源码哈希。远端环境不建 Git 分支。
重要命令通过 `scripts/record_run.py --stage 1 --label NAME -- COMMAND ...` 保留元数据和独立 stdout/stderr。

Stage 1 管线的顺序为：`prepare_stage01_surface.py` → 单独测试 source/welding →
远端 `mesh_stage01.py --profile development_coarse` → `qc_stage01.py --profile coarse` →
`dolfinx_stage01.py save --profile coarse` → 分别用新 1/2-rank 进程 `reload` →
fetch / `collect_stage01.py --profile coarse` → `render_stage01.py --coarse-only` →
coarse 测试和 `gate_stage01_coarse.py` → medium 同样流程 →
`compare_stage01.py` / 正式图 → 完整 pytest → reference integrity → 报告。
所有远端调用使用现有四个 shell wrapper 的 `--stage 1`，解释器为远端项目内
`remote/.env/bin/python`；远端不装 PyVista，渲染只在 WSL。
`mesh_stage01.py` 和 DOLFINx 保存动作拒绝覆盖已有产物；如需重新生成，应先明确新 run 的保存位置，
保留现有验收证据。完整可复查命令已保存在 `logs/stage01/`，不要重新运行 Stage 0 建立契约的脚本。

Stage 1 的低质量非退化单元及人工审核条件仍然有效；Stage 2 不修改或求解这些真实血管网格。

## Stage 2：只规定入口总流量

验证结果、七张审核图和边界条件的解析兼容性说明见 [Stage 2 报告](reports/stage02/REPORT.md)，
符号与弱形式见 [FORMULATION.md](reports/stage02/FORMULATION.md)。正式入口只有总流量约束和
纯法向牵引，正式出口为零牵引；没有入口速度剖面、压力固定点或压力 nullspace。
完整应力下，有限圆管的这些条件存在端部效应。经用户确认，另设端面切向牵引一致的
Poiseuille 解析验证案例；它是独立验证接口，不改变正式边界条件。

源码按空间、弱形式、边界、求解、几何、诊断、结果 I/O 和派生场分开。
`solve_stokes` 使用 SI、P2/P1/native Real 和 MUMPS；`rho` 仅用于 Reynolds 数诊断。
`configs/stage02_pipe_benchmark.yaml` 是合成数值基准，不是实验配置。

已有结果在 `outputs/stage02/cases/CASE/`：`metadata` 保存配置/来源，`qc` 保存积分与解析比较，
`solution/fields.xdmf` 与 `fields.h5` 供可视化；原始 P2/P1/Real 系数由
`primary_checkpoint.npz` 和对应 `meshes/PROFILE/pipe.xdmf`、`pipe.h5` 共同恢复。
派生场是单元参考中心处采样的 DG0 梯度、应变率、涡量，定义在每个 `export_manifest.json` 中。
全部 12 个案例已在新进程恢复，包含从 2 rank 检查点恢复到 1 rank 的情况。

以下命令只重建 Stage 2 的汇总/图，依赖已取回的计算产物：

```bash
.venv/bin/python -B scripts/summarize_stage02_validation.py
.venv/bin/python -B scripts/render_stage02.py
.venv/bin/python -B -m pytest -q
```

远程使用现有 `remote_sync.sh`、`remote_run.sh`、`remote_fetch.sh` 的 `--stage 2`。
完整每次运行命令及源码哈希保存在 `logs/stage02/` 和
`outputs/stage02/remote_return/logs/stage02/`；求解脚本拒绝覆盖已有 case。
首次执行顺序为：数学文档 → 1/2 rank Real probe → 独立 block 测试 → 三个圆管网格与
分别求解 → 解析收敛与正式条件端部效应 → Q/反向/μ/ρ 回归 → 主场及派生场导出/新进程重载
→ 2 rank medium → WSL 图与完整测试 → 历史完整性审计。具体配置和运行环境随每次运行保存。

**止于 Stage 2。Stage 3 未开始；真实血管求解未验证。**

## Stage 1.5：人工端盖内部重铺试验

[Stage 1.5 报告](reports/stage01_5/REPORT.md) 的结果为 **FAIL**：预先冻结的 A/B/C 均未通过
cap 三维面积相对误差 ≤1e−12 的门槛。原 rim 轻微不共面；逐位保留 rim、在原平面新增内部点后，
端盖标量面积随三角化改变约 1e−9（扩展精度已复核）。端盖三角形改善不能替代几何验收。

没有生成候选体网格、winner 或替代 Stage 1 medium，未执行新网格 DOLFINx round-trip，
未修改 Stage 2 solver，未运行候选 FEM，未开始 Stage 3。八张诊断图明确区分实际表面对比
和因 gate 失败而不存在的体网格对比。Stage 0/1/2 的证据目录及既有核心源码全部冻结。

Stage 1.5 的源码为 `cap_remesh.py`、`cap_selection.py` 及 `scripts/*stage015*.py`；
现有 remote wrappers 接受 `--stage 1.5`，将新证据写入 `stage01_5`。
`acceptance_policy.json` 在候选运行前冻结，原 SHA 保存在 `acceptance_policy_lock.json`；
候选脚本拒绝覆盖已完成记录，不进行 A/B/C 之外的搜索。三个被拒绝的表面保留在
`outputs/stage01_5/candidate_{A,B,C}/surface/`，不得视为通过验收的计算输入。

## Stage 1.6：平面端口契约与稀疏端盖

[Stage 1.6 报告](reports/stage01_6/REPORT.md) 的结果为 **FAIL**。新的 derived planar port
contract 保留原始 3D rim、wall、plane 和 Stage 0 contract SHA，正式端口面积改用固定投影
多边形面积；legacy scalar area 继续公开。投影面积、向量面积、中心、覆盖和 cap 质量均通过。

三个预先冻结的 graded-sizing 候选 sparse_A/B/C 分别产生 **1143 / 975 / 887** 个 cap triangles，
全部超过总预算 800 和各端口预算。没有选出 winner，没有生成新体网格或执行 selected
DOLFINx round-trip；没有重新求解 Stage 2，也没有开始 Stage 3。十张图明确标识失败候选
与未生成的体网格比较，不能将诊断图中的 sparse_C 当作 selected。

规则及 SHA 在 `reports/stage01_6/acceptance_policy.json` 和 `freeze_lock.json`；
新源码在 `src/fem3d/planar_port.py`、`sparse_volume.py`、`scripts/*stage016*.py`。
`mesh_stage016.py` 被 surface gate 阻止而未运行，其 QC 算法通过冻结基线和成本回归测试验证。
现有远端 wrappers 支持 `--stage 1.6`，不建立新的 SSH 路径。候选拒绝覆盖和自动扩展搜索。

已取回的证据在 `outputs/stage01_6/sparse_{A,B,C}/` 和 `remote_return/`；历史 Stage 0/1/1.5/2
目录与全部既有核心源码冻结，Stage 1.5 的 FAIL 不变。仅重跑只读验收可用：

```bash
.venv/bin/python -m pytest tests/test_stage016_*.py -q
```

二十个 Stage 1.6 测试模块永久保留；没有 winner 的两项 round-trip 明确跳过，不计作验收通过。

## Stage 1.7：自适应端盖网格优化

[Stage 1.7 报告](reports/stage01_7/REPORT.md) 的结果为 **CONDITIONAL PASS，等待人工审核**。
四个端口独立使用几何导出的尺寸进行有界搜索；不使用固定 A/B/C 或 cap triangle 数量门槛。
Stage 1.6 planar-port v2 契约及 SHA 保持不变，接受策略在搜索前冻结。

首轮真实体网格 `iteration_00` 即满足全部门槛并以 `FIRST_FEASIBLE_ACCEPTED` 停止：
tetra 从 147,569 增至 147,948，minSICN<0.1 从 153 降至 3，cap-adjacent 从 129 降至 0；
P1/P5/median 均提高，P2 velocity proxy 从 801,405 增至 804,528（+0.390%）。
wall/rim displacement 均为零，投影面积和向量面积身份保持。

独立重跑的几何、连接、质量及决策一致。生产和验证合计生成 2 个体网格，未超过总预算 5。
selected 的新进程 1/2 rank DOLFINx 重载均通过，包括对重载 tetra 重新计算质量。
完整测试 267 passed、5 个历史 skipped；Stage 1.7 的 21 个模块共 50 passed、无跳过。
12 张 WSL 图、决策日志、策略 SHA、完整性审计及终端摘要保存在 `reports/stage01_7/`。

实现为 `src/fem3d/adaptive_{port,surface,volume,qc}.py` 和 `scripts/*stage017.py`。
真实试验保存在 `outputs/stage01_7/surface_trials/`、`iteration_00/`、`selected/`、`determinism/`；
现有 remote wrappers 支持 `--stage 1.7`。已完成 optimizer 与试验拒绝覆盖，原始运行日志永久保留。
依赖已取回 artifact 的只读验收可执行：

```bash
.venv/bin/python -m pytest tests/test_stage017_*.py -q
```

Stage 1.5 和 Stage 1.6 的 FAIL 保留，历史证据和 Stage 2 core 未修改。
尚未执行真实血流求解或 FEM 网格收敛验证；**停止于 Stage 1.7，不自动开始 Stage 3**。

## Stage 1.8：替代网格器评估（已清理）

用户已确认 Stage 1.8 PASS，最终生产结论为 **KEEP_STAGE017**。
[最终报告与审核图](reports/stage01_8/REPORT.md) 保留；按 Stage 3 明确要求，
专用实现、测试、安装包及原始实验数据已删除，说明见 [CLEANUP_NOTICE.md](reports/stage01_8/CLEANUP_NOTICE.md)。
生产网格继续使用 `outputs/stage01_7/selected`，相关历史源码只由已有 Git 历史追溯。

## Stage 3：真实三维血管 FEM 首次求解

**FAIL — SINGULAR_PRESSURE_SUPPORT。** [正式报告](reports/stage03/REPORT.md) 与完整失败证据已保存。
使用原 Stage 1.7 selected 和冻结的 Stage 2 core，实际创建 847,707 个未知量并完成仅组装检查。
唯一一次 4-rank MUMPS 正式分解报数值奇异；独立矩阵诊断确认两个压力零模式，
对应两个低质量 tetra 的全部 P2 速度节点被壁面约束。没有修改网格、添加 pressure pin 或更换 solver。

没有有效流场、checkpoint 或 solution roundtrip；报告图片明确区分真实几何/资源证据与结果不可用状态页。
物性来自只读 3D vascular 当前正式配置；实验 Q 保持 null。
完整测试 288 passed、14 skipped、0 failed，其中 9 项 Stage 3 解相关测试因无解跳过，阶段仍为 FAIL。
43,038 个只读参考条目、2,276 个历史文件和 Stage 2 核心未变。
停止于 Stage 3；不自动开始下一阶段。
