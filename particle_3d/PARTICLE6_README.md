# Particle-6：LAMMPS 状态、邻居和二进制重启桥接

P6 只接入存储与查询。P0–P5 的公式、几何、采样、接触及物理时间积分代码保持原提交内容。P5 人工审核确认单独保存在 `reports/particle5/PARTICLE5_MANUAL_ACCEPTANCE.json`，不覆盖原自动验证快照。

## 运行环境

在项目根目录运行；只安装到本项目 `.venv`，不修改系统 Python 或 Frozen FEM venv。

```bash
python3 -m venv .venv
.venv/bin/python -m pip install lammps==2025.7.22.4.0 mpich==5.0.1.post1
```

本次环境复用了已有数值依赖：项目 `.venv/lib/python3.13/site-packages/particle_upstream_readonly.pth` 只读引用 `/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/lib/python3.13/site-packages`。精确 Python、LAMMPS、共享库路径、MPI 能力和包列表在 `reports/particle6/data/00_lammps_environment.json`。LAMMPS wheel 需要 MPICH ABI；桥接构造器仅在本进程预载项目 `.venv/lib/libmpi.so.12`，不建立系统链接、不更改全局环境。版本为 22 Jul 2025 / `20250722`。验证为 WSL CPU、单 rank；多 rank 明确拒绝，尚未验证并行分区与 ghost 通信。

## 接口和责任

- `lammps_state.py`：冻结 V1 编码和 per-atom 属性，`BridgeParticle` 验证有限值、ID、形状和包围半径。`contracts/particle6_state_v1.json` 是机器契约。
- `lammps_bridge.py`：`LammpsParticleBridge` 读写实际原子 ID、位置和 `fix property/atom` 字段。`rebuild()` 只执行 `run 0 post no`，从 Python API 提取邻居并按稳定 ID 排序。读出的状态为物理计算输入。
- `lammps_neighbors.py`：两条查询使用同一 `ValidationNeighborPolicy`。LAMMPS 的查询半径为 center cutoff + skin；standalone 复用原 P4 包围盒查询并施加同一中心距离筛选。其后仍由原 P4/P5 决定真实 gap、接触和近场资格。若查询范围不能覆盖原求解器需要的候选，直接报错。
- `particle6_stepper.py`：调用原 P4/P5 trial 的同一个 Python code object，仅在私有函数 globals 副本中注入邻居查询依赖，不修改历史模块或任何全局绑定。原 P3 接受子步时才写回 LAMMPS；原 P2 `advance_orientation` 推进四元数。两条比较路径独立演化，standalone 直接调用原始工厂。此适配针对锁定依赖提交，升级上游时必须重新验证。
- `particle6_checkpoint.py`：写出实际 LAMMPS `state.restart` 与全局 `particle_sidecar.json`，按 SHA256 检查后在新实例中读入。sidecar 没有逐粒子状态，不重新采样；恢复使用相同 fix ID、属性名称、类型和顺序。Frozen 文件哈希、P0–P5 提交、黏度、模式、时间、步号和查询策略都必须匹配。此确定性验证无运行时 RNG，`rng_states = null` 明确记录这一点。

物理几何来自原 `particle_3d` 类型，不使用 LAMMPS sphere/ellipsoid geometry。为逐位恢复原 P4 椭球旋转矩阵，额外保存 `d2_rotation 9`，避免从四元数重建引入最后一位变化。胶囊轴、半径和柱段长度决定胶囊几何；`RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION`。两种表示按原规则保存，不能把胶囊四元数用作旋转轴。

LAMMPS 使用 `atom_style atomic`。`mass * 1` 是 `LAMMPS_BRIDGE_PLACEHOLDER_ONLY`；原求解器不读取质量、force 或 torque。force 数组在初始化分配后置零，此后每次重建前后只读审计，非零立即失败。atomic 模式没有 torque 数组，诊断记 0 并明确 `UNDEFINED_FOR_ATOM_STYLE_ATOMIC`。所有物理速度来自 P5 阻力求解器或原 P4 验证路径。命令白名单禁止非零 run、积分 fix 和真实 pair potential。

混合形状 100 步与重启使用原 P4/P2 运动学；P5 球体场景单独验证阻力、润滑与同时接触。真实 FEM 两球另行重跑完整原 P5 `simulate_resistance`，逐接受子步核对物理状态、间隙、近场和残差。混合 RBC 验证不构成非球形润滑或真实 RBC 通过血管的证明。

## 复现

```bash
.venv/bin/python -B particle_3d/scripts/run_particle6_validation.py --part synthetic
.venv/bin/python -B particle_3d/scripts/run_particle6_validation.py --part real
.venv/bin/python -B particle_3d/scripts/run_particle6_validation.py --part restart
.venv/bin/python -B particle_3d/scripts/generate_particle6_report.py
.venv/bin/python -B -m pytest -q -p no:cacheprovider particle_3d/tests/particle6
```

`restart` 不覆盖已存在的 bundle；重现前应将旧的 `reports/particle6/checkpoints/{mixed,sphere}` 移到自选归档目录。测试自动使用临时目录，不修改交付 checkpoint。图只从保存的 CSV/JSON 生成；永久测试要求重新生成的 11 张 PNG 与清单 SHA256 完全相同。完整初始/最终测试命令和 JUnit 见 `reports/particle6/logs`。

checkpoint 首次生成时记录工作树源码 SHA256 和当时 HEAD；全部检查通过后先创建源码提交，再用 `finalize_particle6.py --require-pass --bind-source` 核对已测试文件与提交内容并绑定 checkpoint 清单的确切源码提交。此绑定不改变二进制或 sidecar 的任何物理内容。最终验证 JSON 保存源码、CSV/JSON、图、checkpoint 和日志哈希；证据在第二次本地提交保存。

## 范围

未冻结生产 neighbor cutoff、skin、lubrication cutoff 或 timestep。非球形润滑未冻结、真实 RBC 通过未建立、真实 RBC LAMMPS 动力学延期、亚纳米间隙连续介质有效性未建立。固定存储 box 是 V0 验证域，越界报错。无完整悬浮液、CFD 重算、LAMMPS 物理、Particle-7 或 MPI 并行验收。

接口依据：[LAMMPS Python 安装](https://docs.lammps.org/Python_install.html)、[property/atom 与 restart](https://docs.lammps.org/fix_property_atom.html)、[pair zero](https://docs.lammps.org/pair_zero.html)、[Python 邻居访问](https://docs.lammps.org/Python_neighbor.html)。本次以已安装 20250722 版本的实际探针和永久测试确认可用能力，不依赖最新文档中的 full-list 扩展。
