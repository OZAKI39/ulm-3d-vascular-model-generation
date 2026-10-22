# 科学结果与输入索引

这些内容具有不同用途，不能相互替代。路径均相对新工程根目录。

| 内容 | 路径 | 身份与本轮行为 |
|---|---|---|
| 原始 BraVa 血管 | `inputs/geometry/source_core_mm.stl` | 打印工程提供的四端口封闭内腔；源坐标 mm |
| 中心线、半径、端口与姿态 | `inputs/geometry/` | 原始 SWC、端口 STL、几何记录及 candidate 0/15 变换 |
| 权威 SI 网格 | `inputs/mesh/SV_MESH/` | 135852 节点、692313 四面体；网格与边界只保留一套权威输入 |
| 历史压力出口场 | `cases/balanced_pressure_final/` | 原 XML、PETSc 选项、日志、checkpoint、冻结场与原始 WSS 均保留 |
| 历史等流量校准场 | `cases/calibration_equal_split/` | 强制等流量，用途是校准压力；与最终压力出口场分开 |
| 旧正式群体 | `results/population/realization_0000_contactfix/` | 原 145 个文件，含全部状态块、三个压缩账本和身份清单 |
| 新固定协议 | `configs/population/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json` | 原协议逐字节保留；当前默认准备入口；尚无新正式群体 |
| 固定协议原验证 | `results/reports/fixed_window/` | 旧验证记录原样保留；本轮测试另在 `validation/` |
| 原 CFD 身份 | `results/reports/original_cfd/ACTIVE_FLOW.json` | 不会由新算例自动覆盖 |
| 速度急变与背景场诊断 | `diagnostics/` 下事件诊断目录 | 原始证据、图表和判断保留；当前计算函数在 `src/human_vessel_flow/event_diagnostics.py` |
| 参考直管与 Gmsh V2 | `diagnostics/reference_pipe/` | 包含 V2 参考网格、冻结场、误差分解数组与失败判定 |
| P2/P1 | `diagnostics/p2p1/` | 资源停止、曲面映射/压力回填、边界和原生补丁证据；不并入生产求解器 |
| 四视角微泡动画 | `results/visualization/population_four_views/` | 原 MP4、展示参数、QC、联系图逐字节保留 |
| 原流场可视化 | `results/visualization/` | 保留有价值成品，未搬运每一帧中间图片 |
| 本次总览 | `results/reviews/migration_final/OPEN_RESULTS.html` | 新生成的只读核查页，链接原四视角动画 |

## 最重要的 SHA256

| 文件 | SHA256 |
|---|---|
| 正式冻结流场 `cases/balanced_pressure_final/frozen_flow/steady_flow.vtu` | `73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc` |
| 校准冻结流场 `cases/calibration_equal_split/frozen_flow/steady_flow.vtu` | `f5d7a2ccf7cdcc08cbfebd725d0a17e23e5780df23e4ef6488a26ffa5983d4d2` |
| 权威网格 `inputs/mesh/SV_MESH/mesh-complete.mesh.vtu` | `e755a6ce04eca518fc02df5cd6139c98887e94617bdf9410d72f55b99f443c94` |
| 原粒径来源契约 | `5c0842d6e6c149e129519615257537fa53a95fb19214adda9105346a9185b429` |

完整复制前身份：[SOURCE_MANIFEST_BEFORE_COPY.json](provenance/SOURCE_MANIFEST_BEFORE_COPY.json)；补充诊断身份：[DIAGNOSTIC_SOURCES_BEFORE_COPY.json](provenance/DIAGNOSTIC_SOURCES_BEFORE_COPY.json)。最终逐文件对应关系见 [FINAL_FILE_MAPPING.json](provenance/FINAL_FILE_MAPPING.json)，最终不可变内容核查见 [FINAL_IMMUTABLE_DATA.json](../validation/FINAL_IMMUTABLE_DATA.json)。

旧日志里的服务器绝对路径和历史 Git commit 保留原文，属于溯源记录。当前读取代码使用 `configs/cfd/cases.json` 和 `configs/population/environment.json` 的新工程内部路径。历史 XML 没有改写成新实验；未来算例从模板生成新的 XML。
