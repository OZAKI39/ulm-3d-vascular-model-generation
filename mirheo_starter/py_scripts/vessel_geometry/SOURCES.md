# 本地源码与数据复用说明

本模块迁自本地 `/home/lzy/projects/bloodflow_starter` 已改进的第一阶段实现。
未从远程覆盖本地文件；不在运行时导入来源工程，也没有复制旧 Python 环境。
逐文件来源路径与迁移前源码 SHA-256 见 `code_reuse_manifest.json`。

| 能力 | 来源 | 迁入位置与修改 |
|---|---|---|
| GeometryError、BoundaryPatch、VesselGeometry | `py_scripts/vessel_geometry/model.py` | 同名文件，保持原文 |
| 无处理 VTP 读取、明确单位、CSV、原 Windows 路径解析、独占写入 | `py_scripts/vessel_geometry/io.py` | 同名文件；补两个只读工程保护、迁移配置路由、JSON 指数溢出拒绝 |
| 无序标签、可变端口数、多 wall 标签、面积加权中心、拓扑诊断、回读比较 | `py_scripts/vessel_geometry/validation.py` | 同名文件，保持原文 |
| schema 1 NPZ/JSON 保存和正式 `load_package` | `py_scripts/vessel_geometry/export.py` | 保留旧功能供迁入测试；补迁移清单校验、失败/中断拒绝、身份完整性检查；新配置转入 `migration.py` |
| 首阶段合成测试 | `test_code/test_vessel_geometry.py` | 同名文件，40 项原测试保持原文 |
| 全三角面 Plotly、端口中心/几何法向、交互控件、离线 JavaScript、HTML 核对 | `test_code/review_vessel_geometry.py` | 同名文件；保留显示主体，补结果复用、阶段状态、可选 Windows 打开和浏览器检查 |
| 入口结构 | `py_scripts/import_vessel_geometry.py` | 同名文件，调整阶段说明；使用正式导入路由 |

新增 `migration.py` 仅处理固定来源数据包、原文件副本、迁移清单、相关版本匹配和复核；
`test_vessel_geometry_migration.py` 补充迁移边界情况；`check_vessel_geometry_browser.cjs`
用现有 Windows Chrome/Node 检查离线页面，通过临时进程管道通信，结束后关闭。

检查过更早的 `/home/lzy/projects/ulm_3D_vascular/utils/cfd_flow/io.py`、`geometry.py`
及当前 `configs/cfd_flow.yaml`。旧 `geometry.py` 有固定一入口三出口、CSV 行序命名和
无条件三角化，因此采用 bloodflow_starter 的已改进实现，没有迁入这些旧限制，
也不运行旧重建或 CFD pipeline。

在所选 bloodflow_starter 项目根和源码中未发现项目级 LICENSE 或单独版权声明。
保留全部被复用文件的原有文件头、说明和来源，不赋予新的许可证或声称重新授权。
第三方 NumPy、VTK、PyYAML、Plotly 及其依赖保留安装包内各自的许可证。
更早来源项目中的 Ultraliser / LBPM 许可证不属于这些代码，未误套用。

正式迁移保留原 `import_manifest.json` 及它覆盖的全部文件。
新增 `migration_manifest.json` 记录目标相对路径映射；旧绝对路径及原 Windows 路径只是
历史元数据。独立加载器仅访问目标包内文件，无需来源目录或源项目环境。
