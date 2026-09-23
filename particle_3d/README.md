# Particle-0：冻结 FEM 三维采样

本目录只做 Frozen FEM → 三维位置查询。所有开发文件均位于 `particle_3d/`，与冻结 FEM payload 分离。没有微泡/RBC 实体、时间积分、壁面力、碰撞或 CFD 求解。Particle-1 尚未开始。

## 重复运行

从总 worktree 根目录运行；Python ≥3.11，依赖见 `pyproject.toml`。本机复用已有只读读取环境：

```bash
export PATH=/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin:$PATH
python -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle0
python -B particle_3d/scripts/generate_particle0_report.py --stage all
python -B particle_3d/scripts/finalize_particle0.py --handoff-root /home/lzy/projects/ulm-3d-vascular-model-generation/formal_3D_flow_solver/FEM_SimVascular
```

新环境可先 `python -m pip install -e 'particle_3d[validation]'`；之后所有测试及制图均可离线运行。测试不需要 GPU、SSH、服务器、CFD 或网络。本机没有安装、改动 FEM 环境或调用任何求解器。`requirements-validation-lock.txt` 记录本次读取环境。

每步独立制图：`--stage 0` 到 `--stage 6`。报告/图/CSV/JSON/日志分别保存在 `reports/particle0/`、`figures/`、`data/`、`logs/`。科学采样和几何数据固定 seed；计时不是确定性数据。

不读取真实 FEM 的快速数学测试：

```bash
python -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle0/test_tetra_barycentric.py particle_3d/tests/particle0/test_affine_*.py
```

默认 FEM 根目录由当前 checkout 推导；可通过 `PARTICLE0_FEM_ROOT` 指向另一个**完整且通过 manifest 验证**的冻结 payload。不能指向旧的、没有 handoff 的开发目录。SHA 不符会报告文件、expected 和 actual，并停止载入。

handoff 原测试 `test_main_history_preserved.py` 明确要求运行分支名称仍为冻结分支，而且仅允许当年的 FEM 改动。因此最终脚本在未修改的原 handoff checkout 运行完整官方 18 项测试，在 Particle 开发 worktree 内另行运行两个官方 payload 完整性脚本。没有跳过、修改或放宽该历史测试。

## 查询接口

```python
from pathlib import Path
from particle_3d import FrozenFEMField

field = FrozenFEMField.from_frozen(Path("formal_3D_flow_solver/FEM_SimVascular"))
sample = field.sample([1.0e-4, 6.0e-5, 1.2e-4])
print(sample.inside_lumen, sample.tetra_id)
batch = field.sample_many([[1.0e-4, 6.0e-5, 1.2e-4], [0., 0., 0.]])
```

`FlowSample` 的返回值如下；`FlowBatch` 每项增加一个最前面的查询维度。

| 字段 | 单点形状 | 单位 |
|---|---|---|
| velocity_m_s | (3,) | m/s |
| pressure_pa | scalar | Pa |
| velocity_gradient_s_inv | (3,3) | 1/s |
| vorticity_s_inv | (3,) | 1/s |
| strain_rate_s_inv | (3,3) | 1/s |
| inside_lumen | bool | 无 |
| tetra_id | int | canonical volume 的零起始 cell index |

核心只接受有限的三维坐标并转换为 float64。错误形状、NaN、Inf 抛出 `ValueError`；一个批次含无效坐标时整个批次拒绝，不返回部分结果。空批次 `(0,3)` 合法。

域外返回 `inside_lumen=False`、`tetra_id=-1`，全部物理量为 NaN。精确位于闭合边界上的点允许作为流体场查询；这不代表未来有限半径粒子的中心可以位于墙上。多个单元都包含某点时选择最小 canonical tetra_id。对输入数组作防别名复制，并将源数组及返回数组标为只读。

## 数学约定与数值精度

barycentric weights 是四个顶点对查询位置的贡献权重。令 `A=[x1-x0,x2-x0,x3-x0]`，用 float64 计算 `N[1:]=inv(A)@(x-x0)`，`N0=1-sum(N[1:])`。速度和压力分别按这四个权重求和。没有借用 VTK 的插值结果作为答案。

VTK locator 只枚举包围盒相交的候选单元，最终包含判断由自己的四面体公式完成。体网格决定所有插值顶点和单元编号；flow VTU 只提供按全局节点编号排列的 POINT Velocity/Pressure。flow connectivity 仅用于审计行对应，不用于插值。

gradient 是速度沿各空间方向变化的快慢。定义 `G[i,j]=∂u_i/∂x_j`；实现使用节点速度差与 `inv(A)` 相乘，等价于 `Σ ua⊗∇Na`，减少恒定速度分量的抵消误差。vorticity 是局部旋转趋势，按 curl(u) 的三个分量计算。strain rate 是局部拉伸和剪切的快慢，定义为 `(G+G.T)/2`。

每个线性四面体内 G 为常量；共享面上的速度/压力连续，G 可以有跳变。不对 G 平滑。在面/边/节点上返回所选最小编号单元的 G。

固定几何容差公式：

```text
eps = float64 machine epsilon
h = ||A||inf
kappa = ||A||inf * ||inv(A)||inf
c = max(abs(tetra vertex coordinates))
tau = 64 * eps * kappa * (1 + c/h)
-tau <= every N <= 1 + tau
locator bounding-box padding = max_over_cells(4 * tau * h)
```

`c/h` 覆盖未平移的 SI 坐标相减误差，`kappa` 覆盖单元形状影响，64 为这段固定的减法、3×3 求逆、乘加和权重求和提供浮点运算余量。这个常数从实现开始固定，未依据失败结果放大；并非指定一个物理长度来扩张管腔。包围盒 padding 只扩展候选搜索，不决定 inside。`tau >= sqrt(eps)` 的几何直接拒绝，不尝试修网格。本次 max kappa=41.2181，max tau=1.0611e-10，搜索 padding=3.6917e-16 m。

人工场固定 `atol=rtol=256*eps=5.6843e-14`，各分量约为 1，单元 kappa≈1.94。真实节点误差界为 `4*tau*局部节点值跨度 + 16*eps*局部最大绝对值`；共享面采用 `4*tau*局部最大绝对值`，逐点实际误差和界均保留在 CSV/JSON。以上阈值未因测得结果而放宽。

## 验证边界

本目录验证采样实现，不重新验证 Stage Q 的物理正确性。没有重做 mesh convergence、FEM timestep study 或 CPU/GPU field equivalence。真实场导数没有独立连续真解；报告的最大 gradient error 来自人工已知答案测试。SonoVue sampler 只登记来源，尚未接入正式 particle population。

人工审核状态始终为 `PENDING_USER_REVIEW`。自动测试通过不能替代用户审核。
