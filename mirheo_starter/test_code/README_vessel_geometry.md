# 血管几何及标签迁移核查

本次仅核查血管几何和标签，尚未生成 SDF 或运行流体。
此输出是后续 Mirheo 适配使用的项目数据；未 import Mirheo、初始化 CUDA/MPI，
也不是求解器已执行的血管边界。

## 一键重新核查

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -m test_code.review_vessel_geometry --config py_scripts/vessel_geometry_import.yaml
```

程序打印正式数据包及 HTML 的绝对路径，不需要手动猜 run_id。
输入文件、配置字节、相关源码和必要依赖版本完全一致时，校验并加载已有结果，
不重新生成包或页面。不一致时分配新的 `migrate_<UTC时间>_<随机后缀>`；不覆盖旧结果。
损坏的已有数据包不会被复用，保留它并产生新结果。已有核查页被修改或损坏时会明确报错，
不会覆盖原页。缺少真实输入、LFS 指针、来源哈希不匹配时为 BLOCKED，不回退到示意模型。

正式导入入口（不生成页面）：

```bash
.venv/bin/python -m py_scripts.import_vessel_geometry --config py_scripts/vessel_geometry_import.yaml
```

相对配置文件路径和新迁移配置中的所有相对路径均以
`/home/lzy/projects/mirheo_starter` 项目根目录为基准，与启动时的工作目录无关。
使用 `-m` 时 Python 本身仍需能找到目标工程，所以示例先 `cd`。也支持从其他目录直接运行：

```bash
/home/lzy/projects/mirheo_starter/.venv/bin/python /home/lzy/projects/mirheo_starter/test_code/review_vessel_geometry.py
```

`--help` 只显示帮助，不读取数据包、不处理几何、不生成输出。

## 真实输入与依据

唯一指定包：

`/home/lzy/projects/bloodflow_starter/data/geometry/import_20260907T205342_653302Z_865c6f03`

先在旧项目 `py_scripts/lbm_grid_prepare.yaml`、`musubi_grid_prepare.yaml` 及数据包
`import_manifest.json` 中找到这个准确 run_id，再核对第一阶段正式
`py_scripts/vessel_geometry/export.py:load_package`。未根据最新时间或第一个 glob 选择模型。
固定来源清单 SHA-256 为：

`785c601845d1525f2e7116cfa1f9d891f6938c1c9abf778a86bc03b8acd3162e`

其唯一配套原运行是：

`/home/lzy/projects/ulm_3D_vascular/outputs/cfd_surface_prepare/vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611`

VTP 是该运行 `geometry/cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp`，
端口身份来自同运行 `boundaries/boundary_manifest.csv` 和 `qc/final_surface_qc.json`。
单位依据是同运行 `input/cfd_surface_prepare.yaml` 的 `geometry.input_unit: um`，
另有 `qc/meter_scale_qc.json`。配套米单位 STL 与各端口 STL 由已验收清单明确引用。
旧 Windows 端口路径原样保留，新迁移清单新增已确认本地路径及目标相对路径。

本次所有原输入和包文件 SHA-256 均与已验收记录对应。源包早期核查记录的 PENDING
保持原字节；用户本次声明及旧项目稍后 `stage1_acceptance: ACCEPTED_BY_USER`
另外记录为**来源模型历史人工验收**，不等同于本次迁移人工核查通过。

## 复用与目录

源码详细表及许可证说明见 `../py_scripts/vessel_geometry/SOURCES.md`，来源哈希见
`../py_scripts/vessel_geometry/code_reuse_manifest.json`。
复用了旧第一阶段的数据模型、VTP 读取、无序标签映射、单位检查、三角几何和拓扑诊断、
保存/加载/回读比较、40 项自动测试，以及完整 Plotly 离线核查页面。
没有复用 RBC PLY/轨迹查看器作为血管导入器，也没有修改其代码。

- `py_scripts/vessel_geometry/`：正式读取、检查、迁移、加载逻辑。
- `py_scripts/import_vessel_geometry.py` 和 `vessel_geometry_import.yaml`：正式 CLI 与真实配置。
- `data/geometry/<run_id>/`：正式数据包。
- `test_code/`：测试、显示与用户核查入口。
- `test_code/outputs/vessel_geometry/<run_id>/`：离线 HTML、JSON/中文报告及可选截图。
- `test_code/outputs/vessel_geometry/setup_*/`：本次环境盘点、依赖变更、实际测试日志及受保护文件审计。

数据包保留 `vessel_geometry.npz`、`boundary_patches.json`、`geometry_metadata.json`、
`source_surface.vtp`、`source_records/`、原 `import_*` 文件；这些原文件全部字节一致。
`source_companions/` 是配套 STL 副本，`historical_review/` 是对应历史核查和后续验收依据，
不复制旧流场或旧 outputs 全目录。
新增 `migration_intent.json` 表示迁移起始、`migration_config.yaml` 为配置原文副本、
`migration_checks.json` 为本次检查、`migration_manifest.json` 为完成状态和输入输出哈希。
`migration_intent.json` 保留在最终包里，只有同时存在有效 PASS 完成清单时才允许加载。
中断/失败包不能借用历史 `import_manifest.json` 的 PASS。
哈希清单不能包含自己的哈希，其余正式包文件均覆盖；核查报告另记录完成清单和 HTML 哈希。

## 独立加载与几何约定

```python
from pathlib import Path
from py_scripts.vessel_geometry.export import load_package

g = load_package(Path("data/geometry/实际打印的_run_id"))
points_m = g.points_m
triangles = g.triangles
labels = g.entity_ids
original_point_ids = g.original_point_ids
original_face_ids = g.original_face_ids
for boundary in g.patches:
    print(boundary.entity_id, boundary.port_id, boundary.original_role,
          boundary.boundary_origin, boundary.face_ids)
    faces = g.triangles[boundary.face_ids]  # 连接仍指向完整 points_m
```

加载器只访问给定目录内的相对文件；源项目离线、目录搬走时仍可加载。
上面的 `实际打印的_run_id` 仅演示 API；日常核查用一键命令即可，不要求手动填它。
原始编号继承已验收包定义，为原 VTP 点/面行编号；保留顺序、类型和全部未引用点。

本包有 **73,416 个原始点、67,262 个三角面**；33,633 个点被面引用，39,783 个点未引用。
33,673 是历史唯一坐标诊断数，不能用于替换原始点数。
`points_source` 保留原始 float32 微米坐标；`points_m` 保留已验收 float64 米数组。
本次不再次乘 1e-6，不转 Mirheo 模拟单位，不平移、旋转、归一化、修补或重网格化。
完整表面及所有封盖均保留；wall 和各端口以原始面编号集合分开，不把端口封盖设置为永久管壁。

本模型恰有下列边界；代码没有固定一入口三出口的假设：

| 名称 | entity_id | port_id 末尾（完整值保留在包和页面） | role | origin | 面数 | 面积 / µm² |
|---|---:|---|---|---|---:|---:|
| wall | 1 | null | WALL | WALL | 67071 | 2032.059490243 |
| inlet | 4 | cut_000 | ASSUMED_INLET | CUT_PORT | 56 | 7.819752112 |
| outlet_01 | 3 | cut_001 | ASSUMED_OUTLET | CUT_PORT | 44 | 4.416141896 |
| outlet_02 | 5 | cut_002 | ASSUMED_OUTLET | CUT_PORT | 42 | 4.315070430 |
| outlet_03 | 2 | terminal_000 | ASSUMED_OUTLET | TRUE_TERMINAL | 49 | 5.110621776 |

坐标、连接、编号、标签、分区、role/origin/port_id、原法向和来源身份的迁移比较是精确比较，
容差为 0。四个端口的重新计算面积、面积加权中心与单位法向误差也为 0。
跨 NumPy 版本复算管壁面积的差约为 4.14e-25 m²，中心差约为 4.07e-20 m；
原保存值保持不变。统计采用原验收容差：面积 `1e-18 m² + 1e-10 × |参考面积|`，
中心 `5e-11 m`。派生单位法向重核固定容差为 `1e-12`，无量纲；不因失败放宽。
原 QC 中心是 float32 三角质心的算术平均，复核同定义；不与面积加权中心混淆。

## 自动测试

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m unittest test_code.test_vessel_geometry test_code.test_vessel_geometry_migration -v
```

全部是独立临时目录中的合成测试，包括原 40 项测试及新增的包复制、缓存复用、配置/代码
改变分配新 run、损坏/缺失/LFS 阻断、源目录符号链接保护、独立加载、失败包拒绝、
米单位不重复转换、多 wall 标签、CPU 模块导入及帮助无处理副作用。
真实血管检查单独记录在正式包 `migration_checks.json`，合成测试不替代真实血管验收。
本次没有重新测试 RBC/Mirheo Hello World 功能；采用文件哈希保护这些已有成果。

复用当前目标 `.venv`。保留 NumPy 1.26.4、Plotly 6.5.2；仅安装 PyYAML 6.0.3、VTK 9.6.2
及 VTK 必要的传递依赖。所有实际新增版本、已有依赖未改变的对照见 setup 目录
`dependency_changes.json` 和 `dependency_install.log`。不要复制旧环境或安装旧 requirements。

## Windows 打开与人工观察

在核查命令末尾添加 `--open`，只把本地 HTML 交给 Windows 默认程序：

```bash
.venv/bin/python -m test_code.review_vessel_geometry --config py_scripts/vessel_geometry_import.yaml --open
```

也可把打印的 HTML 路径传给 `wslpath -w`，将得到的 `\\wsl.localhost\\...` 路径
粘贴到 Windows 资源管理器，再用 Edge/Chrome 打开。页面内嵌所需 JavaScript，离线可用；
不启动常驻服务器，不依赖 WSL 图形桌面、GPU 或 MPI。

页面以 µm 显示全部 67,262 个真实三角面，保持真实长宽比例，无降采样。
拖动旋转，滚轮缩放；勾选框隐藏管壁，滑条调整壁面透明度；下拉框突出一个端口，
图例可以隐藏分区。显示中心和继承的几何法向，法向注明闭合表面外向约定及自交未复查，
这些箭头不是流速。管壁没有强加全局唯一法向。

重点检查：

1. 整段形状、各开口位置、方向和旧验收模型是否一致，有无意外平移或形变。
2. `entity 4 → inlet`、`3 → outlet_01`、`5 → outlet_02`、`2 → outlet_03` 是否保持；不是按标签大小或 CSV 行序重新命名。
3. 管壁与封盖分色是否正确；隐藏 wall 后仍能单独看到四个端口。
4. 端口表完整 port_id、role、origin、面数、面积和面积加权中心是否符合预期。
5. 实际表面尺寸约为 `97.813232 × 92.801716 × 80.090073 µm`；坐标系未移到原点。
6. 警告和未检查项：ASSUMED 身份不是生理流向证明；wall 子集 199 条端口开边是正常分区，未补洞；复杂自交 NOT_CHECKED。

状态分别记录：来源历史人工验收 ACCEPTED_BY_USER；本次自动迁移按真实结果；
HTML 生成、浏览器交互按实际结果；**本次用户人工核查 PENDING**。
可选 `--browser-test` 在首次生成时调用当前机器已有的 Windows Chrome 与 Node，
软件渲染并实测隐藏、透明度、端口选择、鼠标旋转、缩放、恢复视角；截图和精确页面哈希
位于该 run 核查目录 `browser_final/`。不是把 HTML 生成当浏览器测试。
`browser_preview/` 是同一目标包的发布前测试，最终 `browser_final/` 还会再次执行最终页面。
普通重新核查只复用已记录页面；需要查看本次实际结果时看 `review_report.json`。
本次实际测试的是 Windows Chrome 152.0.7977.82，未单独实测 Edge。
Chrome 的临时 profile 位于 WSL UNC 路径时，stderr 有缓存、数据库与 Crashpad 锁告警；
日志保留在 `browser_final/browser_stderr.log`。实际页面没有未捕获 JavaScript 异常，
没有远程请求，全部规定交互通过。测试进程已关闭，临时 profile 已清理。
三维文字可能随视角被遮挡，可旋转查看，并以文字图例和完整端口表核对身份。

本次仅完成已验收血管几何及标签向 mirheo_starter 的迁移和核查。
未修改原血管，未生成 SDF、粒子或数值边界，未运行真实血流。
本次迁移的用户人工核查状态为 PENDING。
