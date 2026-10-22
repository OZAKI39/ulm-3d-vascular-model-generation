# BraVa landmark 辅助 MeVO 提取

本工具已实现 `inspect → annotate → extract → 可选 model/surface`，复用项目现有 SWC 校验、导出及 VascularMD 建模代码。严格 ROI 的坐标、半径和 TYPE 均保留原值。输出 SWC、原始编号映射、覆盖检查、JSON 清单及可交互检查的 VTP。

**当前真实 BG001 尚未标注。** 用户已说明没有现成 landmark；两个实际标注模板均保留空边界、`enabled: false` 和 `manually_verified: false`。human 入口已接入 MeVO 数据模式，在标注完成前返回 `NEEDS_MANUAL_REVIEW`。这表示缺少解剖输入，不是安装或渲染错误。

## 解剖范围与判断依据

本项目定义左右两侧 MCA 的 M2+M3、ACA 的 A2+A3、PCA 的 P2+P3，共六类 ROI。它们由人工或可靠影像/图谱支持的近端、远端边界决定；分支代次、分叉次数、半径、直径和坐标阈值均不作为解剖分类依据。原始 BG001 普通 SWC 的 TYPE 仅有 1、3；颜色编码版本使用另一套编号和类型值，不能把其中的 ID 直接移植到普通或优化文件。

联合 ROI 只要求起点和全部远端终点，不强制标注内部 M2/M3、A2/A3、P2/P3 转换点。项目中的 MCA 近端取主分叉/三分叉处，远端界定 M3 结束；ACA、PCA 按任务规定的 A2、P2 起点及 A3、P3 结束确定。BraVa 论文说明其重建省略了 ACom 连接，因此 ACA 近端是否可用需要额外核实。若所需近端边界在源数据中不可直接表示，应明确声明截断，不能把现有树根自动当作完整解剖起点。

## 使用当前 human 界面的 MeVO 数据

为保持左侧全局血管和现有显示效果，默认 human 仍读取原来的 `BG001.CNG_vmd_smooth.swc`，仅将 ROI 来源改为 landmark 限定子图。

- 当前可视化对应标注：[BG001.CNG_vmd_smooth.landmarks.yaml](../configs/mevo/BG001.CNG_vmd_smooth.landmarks.yaml)。
- 原始 BraVa 对应标注：[BG001.CNG.landmarks.yaml](../configs/mevo/BG001.CNG.landmarks.yaml)。
- 原空间 ROI 配置备份：[swc_roi_generate_human_spatial.yaml](../configs/swc_roi_generate_human_spatial.yaml)。

两个标注文件不能互换。节点编号及几何经过原生建模后发生变化，工具按具体 SWC 的 SHA-256 校验，不做按距离匹配或猜测转移。对原始 SWC 提取时，“原值”指原始 BraVa；对优化 SWC 提取时，“原值”指该优化文件，不能声称它仍逐点等于原始 BraVa。

从项目根目录执行以下命令。每次 inspect/extract/batch 的输出目录必须是新目录，已有结果不会被覆盖。

```bash
MEVO_PY=/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python
MEVO_SOURCE='vessel_model/T - Brava/swc_files_vmd/BG001.CNG_vmd_smooth.swc'
MEVO_LANDMARKS='configs/mevo/BG001.CNG_vmd_smooth.landmarks.yaml'

# 输出拓扑、原始编号映射和六类 ROI 的空模板。
"$MEVO_PY" tools/vascularmd_extract_mevo.py inspect "$MEVO_SOURCE" \
  --output-dir outputs/mevo/inspect_BG001

# 单独的轻量标注工具；不会修改现有双视窗渲染器。
"$MEVO_PY" tools/vascularmd_extract_mevo.py annotate "$MEVO_SOURCE" \
  --roi RMCA_M2M3 --output "$MEVO_LANDMARKS"

# 填好边界后精确提取，保留原几何。
"$MEVO_PY" tools/vascularmd_extract_mevo.py extract "$MEVO_SOURCE" \
  --landmarks "$MEVO_LANDMARKS" --output-dir outputs/mevo/extracted_BG001

# 另一次独立输出：严格 ROI + 原始上游上下文 + 官方参数化模型/表面。
"$MEVO_PY" tools/vascularmd_extract_mevo.py extract "$MEVO_SOURCE" \
  --landmarks "$MEVO_LANDMARKS" --output-dir outputs/mevo/modeled_BG001 \
  --proximal-context upstream-edge --model --surface

# 完成标注后，原入口直接将 MeVO 数据交给原双视窗。
"$MEVO_PY" s1-2_swc_roi_generate_human.py
```

上述原始数据路径可以替换成 `vessel_model/T - Brava/swc_files/BG001.CNG.swc`，同时必须使用其对应的原始标注文件。尚未标注时 extract 明确拒绝生成 ROI；不能直接运行空模板得到真实解剖分区。

标注器中，数字 1～6 选择 ROI，P/D/E 切换近端、远端、排除支路；右键点击中心线，将最近的原始 full-graph 节点 ID 加入或移出列表。B 声明近端解剖边界已确定，T 声明近端截断。V 在覆盖检查通过后由操作者显式切换人工确认状态，W 保存 YAML，并备份此前版本。修改边界会撤销该 ROI 的确认状态；没有任何操作会自动把 `manually_verified` 提升为 true。关闭窗口不会替代 W 保存。

没有 DISPLAY 时 annotate 保存或保留模板，返回 `NEEDS_MANUAL_REVIEW` 并提示使用 inspect 的 VTP 在 ParaView 中读取 `original_swc_id` 后编辑 YAML。它不会调用依赖图形桌面的窗口。点选功能用于定位节点，不能代替对影像/解剖依据的判断。

## Landmark schema

实际可编辑的完整 YAML 已在上述两个配置文件中生成，包含六个固定 ROI 名称。主要字段如下：

| 字段 | 约定 |
| --- | --- |
| `version` | 当前为 1 |
| `source.path / filename / sha256 / subject_id / units` | 路径相对标注文件；实际提取同时校验文件名和哈希；单位为 mm 或 um |
| `annotation_date / annotator` | 标注时间和操作者；允许增加影像/图谱依据等顶层 provenance，完整写入清单 |
| `rois.<name>.enabled` | 只输出实际启用的 ROI |
| `side / territory / segments` | 必须与固定名称一致，如 RMCA_M2M3 对应 R、MCA、[M2, M3] |
| `proximal_nodes` | 一个或多个原始/full SWC ID，不使用重新编号的 topo ID |
| `distal_nodes` | 所有远端切点；包含切点本身，停止访问其下游 |
| `exclude_subtree_roots` | 不包含这些节点及其下游子树 |
| `proximal_boundary_status` | `confirmed`、`truncated_to_available_data` 或 `unknown`；unknown 不允许提取 |
| `manually_verified` | 严格布尔值。false 可导出有完整边界的候选，但标为 UNVERIFIED |
| `internal_landmarks` | 可选名称到节点列表的映射，不是联合 ROI 提取的必要条件 |
| `notes` | 操作者说明、解剖依据及不确定性 |

默认禁止未经切点限定的路径延伸到自然末端。漏标时报告完整节点路径、对应 topology edges 和末端 ID，返回 `NEEDS_MANUAL_REVIEW`。显式 `--allow-natural-terminal` 可接受操作者确认的自然末端，并留下警告。`--allow-source-mismatch` 仅用于明确指定文件的 extract，输出强提示且有效解剖状态强制为 UNVERIFIED；批处理不以这一选项猜选标注。

多个近端形成不连通的 ROI 会被拒绝；可以选择共同近端分叉节点。反向边界、未被访问的切点/排除根、不存在的 ID、源图非法、非正半径等均会显式失败。不会自动反转血流方向，也不会删除所谓过大、过小或异常扩张的血管。

## 图语义与实现复用

当前固定 VascularMD 提交为 `770feb8fbfb591d6d43f966db57b1465010b9a13`。

- `full_graph` 节点保存四维 `coords=(x,y,z,r)`；适配层另外保存 `swc_type` 和 `original_swc_id`。普通完整图边的 `coords` 为空的 `(0,4)` 数组。
- `topo_graph` 压缩普通点，节点含 `coords/type/full_id`，边含中间点 `coords/full_id`。节点编号在上游重新生成；用节点和边的 `full_id` 建立显式双向定位记录。
- `model_graph` 进一步保存样条、分叉对象、切向和参考方向，`crsec_graph` 及网格属于后续建模阶段。提取不从这些网格切割。

检查/提取实际调用 `ArterialTree(..., filename=None)`、`set_full_graph()`、`get_full_graph()`、`get_topo_graph()` 和 `check_full_graph()`，不给上游自动重采样机会。当前适配层核对完整图的节点与连接，以及覆盖全部节点的 topo 映射。`check_full_graph()` 主要检查度数，DAG、父编号、根、连通性和数值合法性由项目公共 SWC 校验补足。

已核对 `set_topo_graph()`、`topo_to_full()`、`model_to_full()`：当前活动的转换实现会重建/重新编号节点，模型转换还会重新采样，所以精确提取时不调用它们。遍历使用 successors 的 DFS，复杂度为 O(V+E)，失败路径输出额外消耗与实际报告路径长度相关；没有为每对节点枚举路径。上游自身构建 topo 映射时仍使用它原有的路径查询，未修改这一第三方实现。

可选建模调用现有 `vascular_processing.pipeline.process_file`，由其调用官方 `model_network(radius_model=True, criterion="AIC")`、Spline/Model/Nfurcation、`compute_cross_sections()` 和 `mesh_surface()`。STL 仅为官方表面三角化后的序列化，不使用 tube 或圆柱拼接生成替代表面。

官方 Editor 使用 VPython，并把选择与拖拽/几何编辑耦合。本次仅复用已安装 PyVista 实现小型独立节点选择器，没有引入第二套血管编辑器，也没有改动 Editor 或其他上游源码。

## 输出与严格 ROI / 建模上下文

每个 ROI 独立输出 `*_roi.swc`、`*_roi_node_mapping.csv`。需要上游上下文时另存 `*_roi_modelable.swc` 及其映射；严格 ROI 永不被覆盖。映射包含原始父节点、重新编号后的父节点、TYPE、XYZR、近远端标志及 `is_anatomical_roi/is_context`。清单记录严格节点数、上下文节点 ID、长度及原生加载/建模兼容性。

共用分叉真正共用一个节点。MCA 近端三分叉可以是合法严格 SWC，但原生建模入口可能要求非分叉根；此时使用上游 topology edge 的原始路径作为上下文，不包含其旁支。不足的上游结构不会被补造。

总清单 `*_mevo_manifest.json` 记录状态、边界、覆盖计数、几何差值、半径/直径统计、人工确认和解剖完整性。`*_mevo_preview.vtp` 保存完整源树上的 ROI/边界/上下文标记；重叠 ROI 同时保留各 ROI 数组和位掩码，公共 `roi_id=-1` 明确表示多重归属。

显示适配使用原 `ROIRecord`、NPZ/CSV 存储和 `show_saved_run`。毫米只在数据进入显示时换成微米一次。ROI 包围盒是显示范围，绝不作为解剖裁剪条件；所有启用区域均显示，C 键的分组编号表示固定解剖 ROI 分组，未重新运行 K-means。未确认候选的名称带 `UNVERIFIED`。颜色、透明度、箭头、管面、相机及按键行为沿用现有渲染器，人脑左视窗 Y 轴方向保持。

## 批处理与状态

```bash
"$MEVO_PY" tools/vascularmd_extract_mevo.py batch 'vessel_model/T - Brava/swc_files' \
  --landmarks-dir configs/mevo --output-dir outputs/mevo/batch_raw --recursive
```

按源 SHA-256 匹配，随后校验文件名；不按排序位置。全禁用的空模板不参与匹配，同一哈希有多个活动标注时拒绝猜选。没有标注的文件返回 NEEDS_MANUAL_REVIEW，不影响其他已匹配文件的处理。

`PASS/PASS_WITH_WARNINGS` 退出码为 0；`NEEDS_MANUAL_REVIEW` 为 2；`FAILED` 为 1。模型失败时保留已经通过精确提取验证的严格 SWC，但整体状态为 FAILED，不伪造替代表面。原输出目录拒绝覆盖。原 SWC 始终不写入。

## 本次实际验证

真实检查目录为 [raw_BG001_inspect](../outputs/mevo_validation/raw_BG001_inspect) 和 [current_optimized_inspect](../outputs/mevo_validation/current_optimized_inspect)。原始 BG001：2 810 节点、198 topo 节点、197 分支；当前优化 BG001：2 816 节点、196 topo 节点、195 分支。真实 proximal/distal ID 尚未填写，因此没有声称提取成功的真实 MeVO 节点数或分支数，也不能认定其远端覆盖已完成。

合成数据的完整 CLI 验证见 [synthetic_model](../outputs/mevo_validation/synthetic_model)：严格 ROI 41 节点、2 分支，测试近端 21、远端 41/61，未限定路径为 0；另有 20 个上下文节点，官方模型、VTK/STL 均已实际生成。这些 ID 仅属于测试 fixture，不是任何 BraVa 受试者的解剖标注。

逐项测试与数据保护证据见 [验证目录](../outputs/mevo_validation)；中文研究报告见[实现与验证报告](../references/文本记录%20-%20BraVa%20MeVO%20解剖标注与精确子图提取.md)。
