# 迁移基线、函数映射与清理清单

## 实际选择的基线

本地原开发根目录：

`/home/lzy/projects/temp_storage/clot_github_sync_20261015/clot_clearance`

基线入口：`analysis/y_clot_irregular_current_method_v4/code/run_case.py`。原 Git HEAD 为 `5036268`，分支 `clot-clearance/y-irregular-v4-20261020`。迁移时记录了 **2,647 项未提交状态**；未提交文件和唯一的已跟踪修改均保留在原处，详细清单和哈希见 `archive/reports/SOURCE_BASELINE.json`。

GitHub 快照是 `clot-clearance/project-sync-20261020`，远端核对 HEAD 为 `032cb75ecf2b99290fceb0f4ef672f91eb3df981`。对原开发目录与本地 GitHub 快照的 analysis/pd_clot/configs/scripts 中 **1,160 个源码及配置文件**逐一比较，内容相同；不是只凭分支名称选择基线。逐文件 SHA 在 `archive/reports/SOURCE_COMPARISON.json`。

服务器实际执行版本也通过独立数值基线导出和模块 SHA 核对。已导入的项目源码中，只有 `pd_clot/portable_paths.py` 与本地不同：本地新增依赖根目录、字体环境变量及字体回退；服务器保留原部署路径映射。它不改变本构或加载，新项目已完全删除这个历史路径搬运层，直接使用 `inputs/` 和系统字体。原生 7 份自定义 C++ 使用当前 V4 文件，不从旧生成脚本重新打补丁。

## 旧代码 → 新代码

旧路径前缀 `V4/` 指 `analysis/y_clot_irregular_current_method_v4/code/`。

| 旧实现 | 新实现 | 处理 |
|---|---|---|
| `V4/run_case.py::main` | `clot_sim/main.py::run` | 保留预测细分、原生推进、所有力学/能量验收及提交顺序；运行目录变成显式参数 |
| `V4/run_case.py::save` | `clot_sim/io.py::save` | 同样先写临时文件再原子替换 |
| `V4/run_case.py::Lineage` | `clot_sim/fragmentation.py::Lineage` | 原始 ID、分裂关系、积分体积、形态与功统计原样迁入 |
| `V4/common.py` 的物理常量/路径 | `config.py` + `configs/model.toml` | 单一参数入口；移除历史根目录和 `V2/V3/Y_OLD` 路径 |
| `common.sha/cv/dump/load_npz` | `io.py` 同名函数 | 显式 import；不再 `from common import *` |
| `common.module/energy_model` | 显式导入 `material.SofteningLPS` | 删除 `importlib` 源码加载 |
| `V4/energy.py::SofteningLPS` 全部方法 | `material.py::SofteningLPS` | 非均匀 endpoint-star 运算及顺序保留；用于独立 Python 校验 |
| `V4/native.py::predict` | `damage.py::predict` | 历史最大正拉伸与渐进软化预测，不能提交原生历史 |
| `V4/native.py::Contact` 全部方法 | `contact.py::Contact` | 活键排除、距离门槛、平均体积平方、势能和力保持不变 |
| `analysis/v3_softening_pilot_v1/native.py::Native` 中实际被 AST 选中的方法 | `solver.py::NativeSolver` | 将 V4 当时进行的文本替换展开为可读源码：非均匀质量、SI 力/能量、field 15 |
| `V4/native.py::native` 中嵌套 `Native.__init__/step/dynamics` | 同一 `NativeSolver` | 合并构造函数和覆写方法；不保留工厂生成类、`exec` 或历史基类副本 |
| `V4/native.py::make_case` | `geometry.write_native_inputs`、显式案例目录 | 从唯一参考 NPZ 导出两种必要原生格式；移除历史符号链接 |
| `V4/shear_field_library.py::ShearField` | `loading.py::ShearField` | 保留 +X 切向投影、三次包络和关断后的被动阻力 |
| `V4/supported_surface.py::CurrentSurface` 的几何、phase、metrics | `surface_support.py::SupportedSurface` | 保留 L1 支撑/开隙/可达性；删除未被当前路线调用的另一套 evaluate/commit/work |
| `V4/core_surface.py::CurrentSurface` | `surface.py::CurrentSurface` 的表面、拓扑、历史方法 | 与 Y 壁面适配合成一个显式类 |
| `V4/current_surface.py::CurrentSurface` | 同一个 `surface.py::CurrentSurface` | 几何阶段转为 `_quadrature_geometry`，表面力阶段转为 `_surface_forces`；Y 归属与壁面计算顺序不变 |
| `V4/mapped_helpers.py` 的 cells/initialize/reference_crossings/mapped_ray | `surface_cells.py` 同名函数 | 直接读取当前参考胞；缓存只写 `results/cache` |
| `mapped_helpers.load_kernel` | `surface_cells.load_kernel` + `native/build.py` | C++ ctypes 接口保留；编译移到显式构建步骤，不在运行时编译 |
| `pd_clot/.../y_bifurcation_sources.py::YSources.__init__` 的实际所需部分 | `vessel.py::VesselSources` | 仅加载固定管腔查询所需网格、原 VTU、壁面/端口和小配置；移除历史资产注册路径解析 |
| `background_flow_sampler.BackgroundFlowSampler.locate` | `vessel.py::LumenSampler.locate` | VTK 定位、严格重心坐标、候选重查及拒绝外推原样保留；未迁入未调用的速度采样接口 |
| `y_mural_clot_geometry.locator/hit` | `vessel.py::locator/hit` | 仅提取实际使用的几何查询函数 |
| `native_y_wall_query.NativeYWallQuery` | `wall.py::NativeYWallQuery` | 最近点、真实相邻四面体定向、端口判断保留 |
| `native_y_wall_contact.CONTACT` 的当前有效参数 | `configs/model.toml` + `contact.CONTACT` | 当前壁面力在 surface 中按各点实际质量计算；未搬入旧 WallCallback |
| L1/L2/后处理重复的 connected-components 图构造 | `fragmentation.py::components` | 一份权威连通性函数；surface/analysis/visualization 显式调用 |
| `run_case.py` 的能量残差表达式 | `analysis.py::energy_residual` | 保留运算次序，在运行和测试中共用 |
| `analysis/y_clot_fracture_reference_style_v2/code/render.py::Scene/basis/lines/encoder` | `visualization.py` 同名实现 | 保留最近固定相机红色风格；输入和输出改为参数，源状态不可覆盖 |
| 当前 V4 生成的 6 份材料/接触/协调器 C++ | `native/src/` 同名文件 | 原字节复制；新构建不复用旧 `.o` |
| `V4/mapped_geometry.cpp` | `native/src/mapped_geometry.cpp` | 原字节复制 |

完整新函数/方法索引见 `archive/reports/NEW_SYMBOL_INDEX.json`。

## 实际执行依赖与历史辅助内容

**当前运行必需：** `clot_sim/`、三个物理输入 NPZ、必要 Y 资产、`configs/model.toml`、两个已编译原生库、其构建记录和正常安装的第三方包。

**只在构建时需要：** `native/src`、固定的 LAMMPS 核心源包、CMake/编译器。`native/work/` 中展开源码和对象文件可重建，不是权威维护位置。

**只在测试时需要：** `tests/fixtures/` 的三个完整原保存状态、旧实现独立输出、短步输出和一个无损压缩的已损伤原生检查点。它们不是正式轨迹全集，未改动其原始文件；压缩重启在临时测试目录解压并验证原始 SHA。

**参考输入的离线来源：** 旧 `prepare_reference.py` 曾使用旧 Y 规则图、连续形状和材料场生成当前不规则参考。旧 `prepare_boundary.py` 曾生成固定边界类别。`run_case.py` 不执行它们，本轮直接迁入其被冻结且校验过的输出。没有在新工程偷偷调用旧几何生成器；也没有声称已经迁入一套适用于任意新形状的采样/重新标定工具。

**删除的历史间接性：** 历史 native.py 的模块级 `External` 导入实际上没有被 V4 的 AST 类提取执行；旧 `Native.dynamics` 被 V4 覆盖；旧 `WallCallback` 没有被当前 CurrentSurface 路线调用；pd_clot 包初始化引入的一批旧 material/runner/export 模块仅是间接导入，不是当前材料动力学。没有将这些内容复制到新项目。

**展示/历史验证内容：** 历次后处理脚本、论文图缓存、旧 ELF 库、全部 3.6 ms 帧和中间原生检查点仍在原目录。新项目只保留所需科学报告和最小回归样本，未复制整个 archive 或运行工程。

## 独立性与工程变化

新运行不依赖旧目录中的 Python、C++、数据、字体映射或检查点。测试增加了文件打开审计钩子：访问原本地或服务器历史根目录会直接失败。原始路径仅存在于迁移/基线报告中，作为出处，不能作为运行时寻址。

变化限于模块组织、普通 import、相对路径、数据导出、显式编译、可选短步停止和验证报告。短步停止不会更改时间包络或推进算法。为了保留原始资产 SHA，正式保存状态、参考 NPZ、材料 NPZ、边界 NPZ 与原 Y 文件保持原字节；新原生文本输入允许采用说明当前工程的文件头，数值由原 NPZ 导出。原版 SHA 与新版构建 SHA 分开记录。

## 旧目录清理候选——尚未删除，须用户确认

这里的“可清理”只针对新项目运行已不再需要它们；不等于所有旧研究都不需要这些文件。

| 原位置 | 建议 | 前提 |
|---|---|---|
| 原 V4 `code/common.py/energy.py/native.py/run_case.py/core_surface.py/current_surface.py/supported_surface.py/mapped_helpers.py/shear_field_library.py` | 可从日常开发工作区退休 | 新项目验收；先确认原未提交内容已被 Git 或当前 SHA 清单覆盖，不逐个误删历史实验需要的文件 |
| 原 V4 `code/mapped_geometry.cpp` 及 `analysis/clot_irregular_lps_discretization_v4/code/` 下对应六份当前 C++ | 可不再维护旧副本 | 新 native 构建和回归已通过；原已接受库及其来源记录仍保留 |
| 原 `analysis/v3_softening_pilot_v1/native.py` | 新项目不需要；暂不单独删除 | 其他历史实验仍可能直接依赖，须按整个旧实验的保留策略处理 |
| 原 `pd_clot/lammps_fragmentation/` 和 `portable_paths.py` | 暂不删除目录 | 多条旧研究路线共享；只有本轮所需函数被提取 |
| 原 `__pycache__`、可重建的 mapped-geometry 缓存、重复临时渲染帧 | 优先清理候选 | 确认没有进程使用，确认不在正式交付/原始证据清单内 |
| 原临时 native 编译对象和编译暂存 | 候选，不自动删除 | 保留已接受原库、输入、构建审计及原始检查点；确认不再需要逐字节历史重建 |
| 各正式 verification/cases/checkpoints、完整轨迹、断键日志、账本、SHA 文件 | **保留** | 本轮只迁入少量回归样本，不能据此清理全部科研数据 |
| 其他历史任务、论文、未知用途文件 | **先保留** | `SOURCE_BASELINE.json` 已保留未提交清单；需按研究用途逐项判断 |

没有执行原目录删除、Git clean、重置或正式轨迹覆盖；也没有向 GitHub 推送本次新项目。
