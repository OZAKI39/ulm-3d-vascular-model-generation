# TopBrain 到 BraVa BG001 的最近邻语义迁移与 MeVO 区域建模

## 摘要

本研究旨在将 TopBrain 的大脑中动脉分段信息迁移到 BraVa BG001 的原始中心线，并在不改变原始血管几何和连接关系的条件下提取 MeVO 候选区域。根据用户明确选定的官方长度表对应关系，将 BG001 ColorCoded 文件的 TYPE 3 视为左侧 MCA、TYPE 4 视为右侧 MCA；该约定保留来源冲突说明。方法采用既有 Open3D 相似配准、同侧多供体最近邻投票，以及固定的支持范围拒绝规则。右侧有 5 个供体通过配准检查，得到 44 条 MeVO 分支，组成 5 个独立区域；其中 4 个区域完成原生 VascularMD 建模并生成平滑 SWC、VTK 和 STL。左侧未获得合格配准，全部保留为 UNKNOWN。两个显示模式的真实交互窗口测试通过，154 项回归及核心测试全部通过。结果支持所要求的单侧概念验证，但不构成人工解剖真值或制造适用性认证。

## 1 引言

TopBrain 与 BraVa 提供的信息互为补充。前者具有明确的血管分段标签，后者提供可直接用于几何处理的血管中心线、半径和树结构。本次工作的目标是借用前者的解剖语义，在后者中选择相应分支，而非使用配准后的 TopBrain 几何替换 BraVa 血管。

前期十侧 Pilot 的 MeVO 二分类比较支持最近邻方法：NN 的 MeVO F1 为 0.977445，Partial FGW 为 0.945324，NN 在十侧中胜出九侧。因此，本次保留已验证的相似配准和多供体 NN，退休 FGW 生产代码，不增加新的匹配方法。

## 2 资料与方法

### 2.1 数据版本及左右命名

迁移作用于官方标准化、未经 VascularMD 平滑的 `BG001.CNG.swc`。普通文件与 ColorCoded 文件均包含 2810 个节点，二者根和编号不同，因此使用已有的精度约束、一一节点对应及无向边核对建立版本映射。颜色编码仅用于确定 MCA 所属，不覆盖普通 SWC 的 TYPE、半径、坐标或父节点。

普通文件 SHA256 为 `5b69c6c406e4724ffc15a2a4c6b5958cbb169b822a05098d6ebedeecbbb3afbe`；ColorCoded 文件 SHA256 为 `733b255aeeb7795de45799f64e73b19a93ada45ca63518e9aa479c203658ef0e`。两者均与所取得的官方同版本文件逐字节相同。正式输入路径为 `vessel_model/T - Brava/swc_files/` 下的上述两个文件。

用户选用的命名依据是 BG001 官方具名长度表：TYPE 3 的同类边长总和四舍五入后为 1433.76 mm，对应 LEFT MCA；TYPE 4 为 1775.72 mm，对应 RIGHT MCA。这是长度对应推断，不是官方直接发布的 TYPE 字典。关联研究作者代码中相反的命名仍作为历史证据保留。本次约定只绑定所列 BG001 文件，不自动推广到其他受试者。

### 2.2 冻结的配准和语义支持

配准保留现有初始化、四个旋转候选、ICP 和质量阈值。完整同侧供体池参与配准，合格者按原规则排序，最多取五个。每个供体的 M2、M3 先合并为 MeVO，再进行等权投票。M1 与 MeVO 的概率为相应票数除以有效供体数；没有合格供体时，两者均记录为零，标签为 UNKNOWN，不解释为已估计的分类概率。

支持距离为最近语义样本距离除以该供体配准后的 inlier RMSE。UNKNOWN 规则只使用三个固定条件：至少三个有效供体，获胜类别票比不低于 0.60，供体归一化距离中位数不超过 Pilot 正确点分布的第 99 百分位，即 12.7949772843。阈值为无量纲数，不能解释为毫米距离。

十侧 Pilot 的原有硬预测、逐供体投票和距离已经冻结。本轮读取此前四组校准结果，没有重新运行 Pilot 配准、供体排序或参数搜索。新运行仅针对 BG001 的左右侧。

### 2.3 分支聚合与严格区域

正式区域以原有 VascularMD `topo_graph` 分支为单位。节点概率按照相邻中心线段长度的一半构成梯形权重；已知长度占比低于 0.60 时，该分支为 UNKNOWN。否则，只有平均 MeVO 或 M1 概率达到 0.60，才赋予对应标签。拓扑检查仅给出警告，不自动纠正标签。

严格 ROI 只保留 MeVO 分支的原有边。共享分叉端点属于保留分支的一部分，所以 ROI 节点数可以不同于点级 MeVO 标签数。UNKNOWN 间隔不补连，每个连通区域分别写入 SWC，允许重新编号及将新根的父编号设为 −1。节点映射和分支映射能够追溯到原始文件。坐标、半径、TYPE 和保留边均在 SWC 回读时逐项核验。

### 2.4 建模与显示

每个区域只在其直接上游分支已标为 M1 时加入一条真实 M1 context，形成独立的 modelable 输入。context 在分支及新增节点元数据中与解剖 ROI 分开；共享连接点仍属于严格 ROI。既不人为画直线，也不跨过 UNKNOWN 向上寻找入口。

建模调用现有 `vascular_processing.pipeline.process_file`，由上游 VascularMD 执行 `model_network(radius_model=True, criterion="AIC")`、截面计算和表面生成。本次保持严格拓扑政策，不接受原生分叉合并，也不降低入口校验。严格 SWC 始终单独保留；平滑 SWC 是派生建模结果，并不宣称其几何与原始输入完全相同。STL 从实际生成的 VTK 表面导出并回读检查。

显示层通过 s1-3 的 `topbrain-brava` 数据源读取预计算结果。左侧显示完整 MCA 并按 M1、MeVO、UNKNOWN 着色，右侧默认显示严格 ROI，也可显示成功建模的 VascularMD 表面。原有双视窗、黑色背景、坐标轴、相机、Trackball、框选和快捷键继续复用；左侧向上方向保持 Y 轴。毫米到微米仅在显示适配器中转换一次，窗口回调不执行配准或语义查询。

## 3 结果

### 3.1 四组支持阈值校准

表 1 为此前一次性校准的十侧平均指标；未知点仍计入 MeVO 总漏纳，未从评价分母删除。

| agreement | 距离分位数 | MeVO F1 | MeVO recall | M1 错纳为 MeVO | MeVO 总漏纳 | UNKNOWN 比例 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.60 | 95 | 0.950333 | 0.913994 | 0.088139 | 0.086006 | 0.053684 |
| 0.60 | 99 | 0.972496 | 0.954109 | 0.088139 | 0.045891 | 0.015527 |
| 0.80 | 95 | 0.933339 | 0.879609 | 0.038702 | 0.120391 | 0.119116 |
| 0.80 | 99 | 0.955848 | 0.919724 | 0.038702 | 0.080276 | 0.080959 |

没有组合同时满足降低 M1 错纳率和 F1 相对原 NN 下降不超过 0.01。故生产规则采用预先允许的 0.60 / p99 回退，状态为 `OOD_SUPPORT_GATE_ONLY` 和 `NOT_CLASSIFICATION_OPTIMIZED`。该结果用于识别缺少几何支持的点，不宣称改善域内分类，也不是独立测试集验证。

### 3.2 BG001 两侧迁移

| 指标 | LMCA：TYPE 3 | RMCA：TYPE 4 |
|---|---:|---:|
| 同侧缓存供体数 | 25 | 23 |
| 通过检查／实际投票的供体数 | 0／0 | 5／5 |
| 原始目标节点数 | 626 | 797 |
| 原始拓扑分支数 | 39 | 55 |
| M1 点数 | 0 | 90 |
| MeVO 点数 | 0 | 678 |
| UNKNOWN 点数 | 626 | 29 |
| M1 分支数 | 0 | 8 |
| MeVO 分支数 | 0 | 44 |
| UNKNOWN 分支数 | 39 | 3 |
| 严格 ROI 连通区域数 | 0 | 5 |
| 严格 ROI 节点／原始分支总数 | 0／0 | 687／44 |
| 高可信拓扑逆序警告数 | 0 | 0 |

RMCA 采用的供体为 MRA024、MRA013、MRA015、MRA022 和 MRA011。五者 fitness 的最小值／中位数／最大值为 0.563636／0.785455／0.861818；inlier RMSE 为 3.785730／4.101103／4.211019 mm；归一化 RMSE 为 0.047282／0.051221／0.052594。两侧实际配准池来自现有可用语义缓存，没有补造缺失供体。

LMCA 没有有效配准。其全部尝试的 fitness 范围为 0.912727–1.000000，中位数为 1.000000；inlier RMSE 范围为 0.931157–4.171934 mm，中位数为 2.873885 mm。较好的局部拟合分数不能代替完整 QC。所有尝试的根位置偏差均超过既有允许范围，因此未让这些供体投票，状态为 `REGISTRATION_UNSUPPORTED`。点和分支 CSV 仍完整输出，并全部标为 UNKNOWN。没有更改左右约定、配准阈值或初始化以取得成功。

RMCA 没有检测到 MeVO 下游重新出现 M1 的拓扑逆序，也没有检测到 MeVO→UNKNOWN→MeVO 的间隔模式。五个区域保持各自原有连接，不因数量多于一而自动合并。

### 3.3 严格 ROI 与 VascularMD

| RMCA 区域 | 严格节点数 | 原始分支数 | 真实 M1 context 分支 | 原生建模及 SWC/VTK/STL |
|---|---:|---:|---|---|
| part01 | 142 | 11 | 94_95 | 成功 |
| part02 | 44 | 1 | 95_107 | 成功 |
| part03 | 330 | 24 | 无 | 入口检查失败，严格 SWC 保留 |
| part04 | 130 | 6 | 93_137 | 成功 |
| part05 | 41 | 2 | 92_144 | 成功 |

part03 的严格根节点具有多个下游分支，且不具备符合本次规则的直接上游 M1 context。现有 VascularMD 适配器会将这类入口解释为 sink，故在原生建模前拒绝输入。没有人为延长、增加虚拟节点、移除分支或采用其他拟合方法掩盖这一失败。其 330 个节点、24 条分支仍完整保存在严格 ROI 和 modelable SWC 中；本次成功建模的四个区域不包括这一较大区域。

实际文件位于 `outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/`。五个严格输入为 `ROI/roi_part01.swc` 至 `ROI/roi_part05.swc`，modelable 输入对应 `ROI/roi_partNN_modelable.swc`。成功区域的平滑 SWC 和 VTK 位于 `VascularMD/roi_partNN/roi_partNN_modelable_vmd_smooth.swc`、`roi_partNN_modelable_vmd_surface.vtk`，STL 为同目录的 `roi_partNN_surface.stl`。各目录的 QC JSON 记录实际原生调用、失败阶段和几何比较；没有将示例表面或管状显示网格作为 VascularMD 输出。

### 3.4 界面、回归与资料保护

strict 和 surface 两种模式均在实际显示环境中创建了 1800×900 双视窗，使用 `vtkInteractorStyleTrackballCamera`，左侧 camera view-up 为 [0,1,0]。测试触发 ROI 切换、A/R/S/C 框选显示、F12 截图并自动关闭窗口，两种模式均通过。strict 模式可检查全部五个严格区域；surface 模式仅列出四个成功建模区域。默认 TopBrain 数据源仍保留，s1-2 及共享 renderer 未修改。

当前稳定回归与核心测试共 154 项通过，其中原稳定测试 146 项、核心测试 8 项，失败和错误均为零；61 条警告来自既有依赖弃用提示。另一次真实缓存复用核验将配准、NN 查询和建模入口替换为禁止调用断言，左右侧仍成功读取结果，表明重新打开已有生产结果不重复计算这些阶段。

Partial FGW 专用实现文件已删除，其评价搜索、安装逻辑和测试中的相应代码共退休 734 行；该数字为删除或替换的旧行数，不是项目净减行数。POT 的生产依赖声明已移除，活跃源码未检出 `import ot`。当前环境已安装的 POT 没有卸载，也没有安装新软件。保留的核心包括 `similarity_registration.py`、`semantic_ensemble.py`、`topbrain_semantic_points.py`、`nn_support_gate.py` 及 BraVa 精确分支、ROI 模块。历史最终比较结果与 NN、配准、供体缓存保留。

资料保护结果以 `continued_protection_verification.json` 为准，覆盖 BraVa 输入、TopBrain 原始影像及标签、上游 VascularMD、冻结配准和既有证据。s1-2 SHA256 保持 `3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729`。新的输入溯源文件还记录本次左右命名约定和全部参与语义缓存的 SHA256。

## 4 讨论

本次完成了至少一个 BG001 侧的真实“TopBrain 语义→BraVa 原始分支→严格 ROI→VascularMD→STL”链路。结果也显示，迁移成功、区域可提取与区域可建模是不同条件：左侧因配准失败而没有可提取区域；右侧虽有五个候选区域，其中一个仍因入口结构限制不能建模。分别报告这些状态，比把所有区域合并成一个成功状态更能说明方法的适用范围。

两类主要不确定性仍然存在。首先，左右命名采用用户选择的官方长度表推断，来源间的历史冲突未被新的直接证据消除。其次，BraVa 未提供本次区域的人工分段真值，TopBrain Pilot 的 F1 不能直接当作 BG001 的准确率。所有输出应称为 TopBrain 引导的解剖迁移候选区域。STL 也尚未经过制造或 CFD 所需的完整验证，状态保持 `MANUFACTURING_NOT_VALIDATED`。

本次未扩大至全部 BraVa 受试者，也未进行新的二十五病例验证。现有结果可用于检查候选边界、原始半径及模型几何；后续应在人工确认的具体分支上评估迁移可靠性，再决定是否扩大范围。

## 5 结论

在选定且明确记录左右命名约定的前提下，既有多供体 NN 方法完成了 BG001 RMCA 的 MeVO 候选分支提取。四个独立区域获得真实 VascularMD 平滑中心线及表面输出，未通过检查的左侧和右侧 part03 均保留明确失败状态。原始几何、受保护代码和现有显示风格得到保留，达到了本轮单侧真实链路的主要目标。

## 参考资料

1. [BraVa 官方 BG001 长度表](http://cng.gmu.edu/brava/query_subject.php?id=1)。
2. [BraVa 原始研究全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3971907/)。
3. [TopBrain 2025 v2 官方数据记录](https://zenodo.org/records/16878417)。
4. [VascularMD 上游项目](https://github.com/megdec/vascularmd)。本项目运行所用固定提交及源码哈希见各原生 QC 文件和资料保护审计。
5. 本地研究证据：`unknown_gate_summary.json`、`BG001/summary.json`、`BG001/input_provenance.json`、两份 `brava_*_ui_compatibility.json`、`continued_regression.xml`、`cache_reuse_verification.json`。上述文件均位于 `outputs/topbrain_brava_transfer/nn_production/` 之下。
