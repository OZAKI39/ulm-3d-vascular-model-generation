# Clot Simulation

当前不规则 Y 血栓 V4 的独立工程。先读 [`clot_sim/main.py`](clot_sim/main.py) 的 `run()`，再按它的显式 import 阅读各个物理模块。模型源码只有 `clot_sim/` 一份；自定义原生计算在 `native/src/`。运行和测试均不读取旧工程源码或数据。

研究基线：**Clot Simulation V4 — Codebase Consolidation**。日常开发规则见 [`AGENTS.md`](AGENTS.md)。

本次仅重构结构，没有重新运行冻结的 3.6 ms 正式轨迹。已通过的等效范围及误差见 [`VALIDATION.md`](VALIDATION.md)。原工程、正式轨迹和原始检查点仍保留。

## 快速使用

本机已创建 `.venv` 并完成原生编译。从项目目录执行：

```bash
cd /home/lzy/projects/new_clot_simulation
source .venv/bin/activate
source native/env.sh
python -m clot_sim.main inspect
python -m pytest -q
```

运行 **3 个短步**；保留原始加载时程，不压缩成快速断裂演示：

```bash
python -m clot_sim.main run --output results/my_smoke --steps 3
```

从该短步检查点继续 1 步，输出仍写在同一个案例中：

```bash
python -m clot_sim.main run --stage continue --output results/my_smoke \
  --checkpoint results/my_smoke/checkpoints/short_validation_step_000003 --steps 1
```

以下命令是将来主动进行完整计算的入口，**本次整理未执行它们**：

```bash
python -m clot_sim.main run --stage pilot --output results/my_case
python -m clot_sim.main run --stage continue --output results/my_case
```

`pilot` 计算至 2 ms，`continue` 从合法检查点计算至 3.6 ms。新计算拒绝覆盖已有 `history.csv`；后验验收失败即停止，失败试算与最后合法状态分开保存。当前原生驱动只允许 **1 个 MPI rank**。

## 在另一台机器安装

需要 Linux、Python ≥ 3.12、C++17 编译器、CMake、OpenMPI 开发库，以及离屏 VTK 所需的系统图形库。版本约束在 `pyproject.toml`；本机完整环境记录在 `configs/requirements.lock.txt`。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
python native/build.py --jobs 4
source native/env.sh
python -m pytest -q
```

`native/build.py` 从项目内的固定版本源码包和自定义 C++ 文件完整编译，不借用旧对象文件。它也安装匹配的 LAMMPS Python 包。不要用普通 PyPI LAMMPS 替代这个自定义求解器。详细构建与来源说明见 [`native/README.md`](native/README.md)。

默认使用 `native/build/liblammps.so`。外部安装可通过 `CLOT_LAMMPS_LIBRARY` 和匹配的 `CLOT_NATIVE_MANIFEST` 明确指定，再 source `native/env.sh`。程序核对加载库与构建记录，拒绝混用其他 `liblammps`。

## 目录与阅读顺序

```text
new_clot_simulation/
├── clot_sim/
│   ├── __init__.py
│   ├── main.py                 # 完整运行流程、预测细分、验收、账本与检查点
│   ├── config.py               # 显式参数和项目相对路径
│   ├── io.py                   # 原子写入、SHA 校验、状态读写
│   ├── geometry.py             # 冻结参考点/键/体积和原生输入导出
│   ├── material.py             # 非均匀体积变分 LPS 能量及内部力
│   ├── damage.py               # 历史最大拉伸与不可逆软化预测
│   ├── loading.py              # +X 切向滑移代理和原始时间包络
│   ├── surface.py              # 当前表面、可加载性、守恒分配、壁面力
│   ├── surface_support.py      # L1 有支撑单纯形几何
│   ├── surface_cells.py        # L2 映射积分胞与 C++ 几何接口
│   ├── vessel.py               # 冻结 Y 网格与严格四面体归属查询
│   ├── wall.py                 # 最近壁面、朝向与开放端口几何
│   ├── contact.py              # 非存活键对的颗粒排斥和接触势能
│   ├── solver.py               # 显式 NativeSolver 和 LAMMPS 回调
│   ├── fragmentation.py        # 连通分量、原始 ID 谱系、体积和形态
│   ├── analysis.py             # 参考统计、能量残差
│   └── visualization.py        # 真实保存状态的静态图和视频
├── native/                     # C++、固定第三方源包、构建脚本及说明
├── configs/                    # 物理参数、环境版本
├── inputs/                     # 唯一参考 NPZ、材料、边界、必要 Y 资产
├── tests/                      # 原版数值基线与独立回归测试
├── results/                    # 新运行、缓存、验证日志及可视化
├── archive/reports/            # 来源/冻结/迁移报告；无运行时依赖
├── README.md
├── MIGRATION.md                # 旧→新函数映射、依赖判断与待清理清单
├── VALIDATION.md
└── pyproject.toml
```

`.venv` 是安装环境；`native/work` 是可重建的编译暂存，不是另一套模型源码。`reference.native` 和 `body.data` 在每次新运行中从唯一的 `inputs/reference.npz` 导出，保留原始 ID；它们是求解器输入格式，不是另一份几何实现。

## 一步计算怎样执行

1. `main.run()` 校验输入 → `geometry.load_reference()` 读取原始 ID、参考位置、键和积分体积。
2. `NativeSolver` 配置质量 `ρ V_i`、原始逐键阈值和原生求解器；`CurrentSurface` 建立 L1/L2 几何与固定血管壁查询。
3. `damage.predict()` 对预计漂移位置只读计算 `κ`、`d`。按原有启动步长和 `Δd ≤ 0.1` 门槛二分预测步长。
4. `NativeSolver.step()` 调用 LAMMPS：NVE 半步速度/位置漂移 → C++ `post_integrate` 中提交一次损伤和断键 → 原生材料力与颗粒接触力 → 当前表面/壁面外力回调 → NVE 末半步速度。
5. 回调使用当前位置、存活拓扑、当前表面法向/面积和收敛后的末半步速度。表面分配满足合力、力矩与虚功配对；被遮挡或不可解析的面不会强行施加载荷。
6. `main.run()` 检查预测/原生一致性、不可逆性、固定锚点、积分、断键、能量账本及壁面穿透门槛。
7. 只有通过验收，才提交 Python 表面历史、累积功和碎块谱系，并保存输出/合法检查点。后验失败不会被静默回滚后重新接受。

Python 的 `material.py` 是独立的能量/力校验实现；真实动力学材料提交由 C++ 执行。二者不是可以互换的求解后端。

## 冻结输入与参数

- 1,927 个血栓点、1,561 个锚点、153,986 条参考键；412 个内部参考点曾扰动。
- 血栓积分体积 0.240875 mm³，质量 0.240875 mg。非均匀 `V_i` 存在参考 NPZ 中，不能替换为节点数乘统一 `Δx³`。
- 名义间距 50 µm，视界 150.5 µm；`K = 100000/3 Pa`，`G = 940 Pa`，`ρ = 1000 kg/m³`。
- 实际逐键 `s0/sc` 在 `inputs/material.npz`；弱相 0.10/0.12，混合 0.175/0.21，强相 0.25/0.30，界面 0.140625/0.16875。
- 名义驱动 1000 Pa；上升至 2 ms，保持至 2.4 ms，降至 2.8 ms，随后保留被动滑移阻力；总观察时间 3.6 ms。
- 时间基步 0.25 µs，最小预测子步为其 1/8。原有启动细分和全部验收条件保留。

这些是当前探索性模型的配置，不是新的材料标定或声流求解。Y 流场文件在本计算中用于固定管腔几何归属；没有重新求解 CFD 或增加双向流固耦合。开放端口到达仍按原逻辑终止并要求显式输运交接。

参考几何在此轮直接读取冻结输入，**不重新随机采样、不重建另一张邻域图**。每条参考键按原 ID 保留，运行中仅依据损伤失效。今后若要改变形状/离散/阈值，应作为新物理案例显式生成并审计输入，而不是让目录整理自动修改它们。

## 可视化

```bash
python -m clot_sim.visualization --states \
  tests/fixtures/initial.npz tests/fixtures/first_failure.npz tests/fixtures/final.npz \
  --output results/my_replay --video
```

也可把自己的 `results/my_case/keys` 目录交给 `--states`。输出真实坐标、当前存活键、红色强弱相、固定相机和竖直 Z 轴。每个提供的保存状态停留 1 秒；不插值、不伪造缺失帧。当前验证视频是 **3 个保存状态的稀疏回放**，不是重新计算的连续断裂轨迹。已生成的图与视频入口：[`results/validation/saved_state_replay/index.html`](results/validation/saved_state_replay/index.html)。

## 日常 Git 工作流

只在本目录修改、测试和提交，不再创建同步副本。当前是独立 Git 仓库，初始提交为 `29202fe`，分支 `main`。原 GitHub 仓库的工程布局和提交根与这里不同；连接远端的方案须先确认，不能直接覆盖原仓库的 `main`。

建议沿用 `https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`，在该仓库只使用一个长期分支 `clot-simulation-v4` 保存本工程，保留现有全部历史分支。方案确认前不配置推送目标、不推送。

平时提交：

```bash
cd /home/lzy/projects/new_clot_simulation
git status --short
git diff
# 按本次实际改动显式选择路径，避免把大型输出一起加入。
git add clot_sim/loading.py tests/test_equivalence.py
git diff --cached
git commit -m "Describe the actual change"
```

上述 `git add` 路径只是示例，不表示需要修改加载方法。远端方案确认后，可在本工作树一次性配置：

```bash
git remote add origin https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git
git branch -m clot-simulation-v4
# 仅当用户明确要求推送时执行首次推送：
git push -u origin clot-simulation-v4
# 此后日常提交仍在相同分支；要求推送时直接：
git push
```

这会在已有远程仓库增加一个长期分支，不合并无关历史、不改写任何旧提交，也不删除旧分支。如果该分支已存在，先检查其历史，不能强制覆盖。

## 科研数据保留与 Git 边界

- Git 保存 `clot_sim/`、`native/src/`、固定源码包、配置、必要输入、回归样本与研究文档。每次检查暂存差异和体积，避免笼统 `git add .`。
- `results/` **不整体忽略**：新案例、账本、状态清单和报告会显示为未跟踪文件，须显式决定哪些纳入 Git。不要为使工作树看起来干净而删除或隐藏正式结果。
- 正式轨迹、原生 restart、完整精度状态和唯一输入原地保留，不自动压缩、降采样或覆盖。大型数据不默认推送；若需另行归档，先由用户决定位置，再在案例的小型清单中记录路径、大小、SHA-256、代码提交、配置及原生库身份。清单不代替数据本身。
- `.gitignore` 只忽略环境、编译目录、Python/测试缓存、可重建的 `results/cache/`，以及已注明可重放的三帧验证视频/图片。小型验证报告和索引仍可提交。不要添加全局 `*.npz`、`*.restart`、`*.log` 或 `results/` 忽略规则。
- `native/build/` 虽不提交，但含当前已验证的本地动态库与构建记录；**被 Git 忽略不等于可删除**。任何删除仍需核查可重建性及是否唯一，并按用户批准清单执行。
- 旧正式结果和服务器 `/workspace/y_clot/` 保持原位。历史路径仅用于出处记录，不是新工程运行依赖。
