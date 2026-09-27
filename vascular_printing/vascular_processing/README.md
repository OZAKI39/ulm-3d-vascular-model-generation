# VascularMD 原生建模与 SWC 导出

BraVa 的解剖 landmark 辅助 MeVO 提取及现有 human 双视窗数据接入见 [MeVO 工作流](MEVO_README.md)。该步骤保持所选源 SWC 的 XYZ、半径及 TYPE；可选的后续建模继续使用本文的公共流程。当前真实 BG001 尚待人工标注，程序会明确报告 NEEDS_MANUAL_REVIEW。

本工具直接调用 [megdec/vascularmd](https://github.com/megdec/vascularmd)，固定于提交 `770feb8fbfb591d6d43f966db57b1465010b9a13`。核心拟合来自上游 `ArterialTree`、`Spline`、`Model` 与 `Nfurcation`，没有添加移动平均、Gaussian、LOWESS、自行编写的样条拟合、半径截断或生理约束。

默认 `--model-mode network` 使用整网 Nfurcation 模型，上游会合并距离很近的分叉。该模式默认严格保持原始拓扑；拓扑变化时拒绝导出 SWC，能够生成的表面仅供诊断。显式指定 `--accept-native-merges` 后，允许官方合并近邻分叉，但仍要求保留原始根和每一个原始末端，且所有连接均能追溯到原树上的路径。任何末端丢失或其他重连仍会失败。QC 分别记录原始拓扑是否保持、是否忠实保留原生模型拓扑及实际合并的原始节点。

当前 `s1-2_swc_roi_generate_human.py` 已按用户确认采用合并后的 BG001 SWC；配置为 `configs/swc_roi_generate_human.yaml`，输入目录为 `vessel_model/T - Brava/swc_files_vmd`，样本为 `BG001.CNG_vmd_smooth`，结果另存 `outputs/human_brava_vmd`。原始 SWC 与可视化实现保持不变。

## 运行

从项目根目录运行；当前使用 `/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。

接受原生分叉合并并重新生成完整 BraVa 样本（选择尚不存在输出文件的新目录）：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python tools/vascularmd_smooth_swc.py \
  "vessel_model/T - Brava/swc_files/BG001.CNG.swc" \
  --output-dir outputs/vascularmd/BG001_accepted --accept-native-merges --qc-plot
```

BG001 中官方算法将原始分叉 10、12 合并到 8 对应的分叉，99 个末端全部保留。不指定 `--accept-native-merges` 时，仍按原有严格策略返回状态码 1。只需要 SWC 时可加入 `--no-surface`；当前 human 默认派生数据采用此方式生成。

本次已成功验证的真实 BraVa ROI 输入示例：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python tools/vascularmd_smooth_swc.py \
  outputs/vascularmd_validation/brava_roi_inputs/BG001_anchor000647.swc \
  --output-dir outputs/vascularmd/roi647 --qc-plot
```

此输入是此前选定 ROI 内的原始 SWC 子树，不是完整脑血管文件。坐标、半径、类型、ID 保留原值，仅把裁剪后唯一根节点的父编号设为 `-1`；文件头记录来源和 SHA-256。成功后输出：

```text
outputs/vascularmd/roi647/
  BG001_anchor000647_vmd_smooth.swc
  BG001_anchor000647_vmd_surface.vtk
  BG001_anchor000647_vmd_qc.json
  BG001_anchor000647_vmd_radius.csv
  BG001_anchor000647_vmd_radius.png
  BG001_anchor000647_vmd_run.log
  vascularmd_batch_时间戳.json
```

输出文件已存在时拒绝覆盖，请使用新的输出目录。输入文件永远不覆盖；QC 校验处理前后的源文件 SHA-256。

目录递归批处理：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python tools/vascularmd_smooth_swc.py \
  "vessel_model/T - Brava/swc_files" \
  --output-dir outputs/vascularmd_batch --recursive --qc-plot
```

保留输入的相对子目录，跳过输出目录和已有的 `_vmd_smooth.swc`。逐文件输出日志，单文件失败不会阻止后续文件；不自动删除失败血管，不自动选取“看起来可用”的成分。完整数据集包含非法记录和零半径文件，可能出现输入检查失败。

## 参数

| 参数 | 默认值与含义 |
| --- | --- |
| `--model-mode network` / `--model-mode branches` | 默认整网 Nfurcation；branches 为官方逐分支 AIC，固定原始端点 XYZ 和半径，保留全部连接；必须配合 `--no-surface` 和 original-count，不能启用合并 |
| `--sample-mode original-count` | 默认；每个原始拓扑段输出点数与原段相同，值来自原生连续模型 |
| `--sample-mode spacing --spacing-mm 0.5` | 依据模型弧长均匀采样；步长不大于指定值；不对参数 `t` 等距采样 |
| `--auto-resample` / `--no-auto-resample` | 默认采用上游行为，即开启。先保存原始段映射，再调用原生 `automatic_resampling()` |
| `--surface` / `--no-surface` | 默认生成原生 VTK 表面 |
| `--circumferential-n 24` | 截面圆周节点数；要求至少 8，且为 4 的倍数 |
| `--longitudinal-density 0.2` | 上游 `d`，截面纵向间距相对半径的比例；**不是**每毫米截面数；越小通常越密 |
| `--qc-plot` | 生成原始跳变最明显的若干血管段诊断图 |
| `--plot-branches 6` | 诊断图中最多显示的原始段数 |
| `--blas-threads 1` | 限制小型线性代数求解的 BLAS 并行度，避免大量线程造成开销 |
| `--accept-native-merges` | 默认关闭；显式接受官方近邻分叉合并，保留原始根及全部末端，记录原始分叉编号与合并关系 |
| `--input-units mm` / `--input-units um` | 声明输入 XYZ 和半径已经共同使用的物理单位，不执行单位换算；默认 mm。um 输入采用 original-count 模式，不能配合 --spacing-mm |

允许合并时，`original-count` 指每条保留下来的原始路径的点数；被合并的短主干可能关联到多条输出路径，因此总节点数可以变化。原始总体 QC 始终在原始拓扑上计算一次，不重复计入这些共享主干。合并后 SWC 编号重新生成，human 流程会重新采样 ROI，不能直接沿用旧的锚点编号。

数值单位完全不变，第六列始终为**半径**。`--spacing-mm` 针对本次 BraVa 的毫米数值，既不把输入放大 1000 倍，也不把半径当作直径换算。上游若干原生密度参数按毫米尺度设计；对其他单位的数据，不能据此认为其物理尺度已经适配。

默认网格参数 `N=24、d=0.2` 来自所给论文第 6.3 节的人脑实验及 `Nfurcation` 当前默认设置。当前仓库 `main.py` 的演示采用 `N=48、d=0.25`，可以通过 CLI 显式设置。

## 小鼠数据适配

小鼠样本 `fMOST_0_5_6_0_0_6_0001_02_01` 的原始节点 5707 附近曾触发整网 Nfurcation 的 apex 搜索失败。经用户确认，新增显式的官方逐分支 AIC 模式，保持原始分叉连接和全部末端。先前整网失败目录中的中间文件仍不能当作成功结果。

当前 mouse 默认配置已切换到成功生成的 `vessel_model/T - A high-resolution dataset of mouse brain vasculature_vmd_branch_aic`，结果另存 `outputs/mouse_vmd_branch_aic`。242 条分析分支均已完成；图结构及 ROI 验收通过。全局验收保留原始 Mask 配准提示：拟合后 2 / 9 828 个节点未落在 Mask 前景中，整体状态为 WARNING；没有为消除提示而改写 Mask 或中心线。原始输入对照可运行 `python s1-1_swc_roi_generate_mouse.py configs/swc_roi_generate_mouse_raw.yaml`。半径变化的改善与剩余局限见[逐分支 AIC 实现报告](../references/文本记录%20-%20小鼠%20VascularMD%20逐分支%20AIC%20拟合与半径优化.md)。

小鼠原始 XYZ 为体素坐标，半径已经是微米，不能直接将原 SWC 作为等比例物理坐标交给 VascularMD。使用独立适配入口：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python tools/vascularmd_prepare_mouse.py \
  --config configs/swc_roi_generate_mouse_raw.yaml \
  --output-dir "vessel_model/T - A high-resolution dataset of mouse brain vasculature_vmd_branch_aic_new" \
  --model-mode branches --no-accept-native-merges
```

该命令复用既有规则，按物理中心线总长度选择当前分析分量（或采用配置中显式指定的分量）。先令 XYZ 乘以原体素间距，半径保持不变，再在微米坐标系运行原生 AIC 模型；输出 XYZ 除以原体素间距恢复到原坐标系，半径仍为微米。此处原生重采样密度 0.4～0.6 的单位为点/μm。默认 `max_distance=6` 是无量纲倍数，上游内部乘以平均半径后才成为空间距离阈值；不应将其解释为 6 μm 或 6 mm，也不宣称默认数值经过小鼠影像标定。

逐分支模式直接调用官方 `Spline.approximation(..., radius_model=True, criterion="AIC", akaike=False, max_distance=6)`，独立选择空间与半径的平滑参数。约束为 `[True, False, False, True]`：原始共享分叉、根及末端的位置和半径固定，不施加端点切向约束。沿模型弧长采样，保留每条原始分支的节点数、顺序和 TYPE；内部位置及半径来自拟合模型，不是原始数据插值。短分支沿用官方补足 6 个拟合点的规则。该模式不调用 Nfurcation、不合并分叉，也不生成或声称获得原生分叉表面；非正半径直接报错，不做截断。QC 逐分支记录两个平滑参数、控制点数量、拟合点数及原节点到输出编号的映射。

小鼠适配入口的 network 模式仍默认接受官方近邻分叉合并，可通过 `--no-accept-native-merges` 恢复严格模式；branches 模式始终关闭合并。只有选定分析分量接受建模，其余分量在派生 SWC 中保留原节点、类型、坐标、半径与连接，明确作为参考数据。新建模节点分配到原最大编号之外，避免与其他分量冲突。原始图像和掩膜采用符号链接复用，保持坐标配准及既有显示方式。任何建模或导出失败均不发布替代 SWC，不混入其他平滑算法。

输出保持原数据的目录结构和样本名称，并包含 `vascularmd_mouse_manifest.json`、物理单位输入、原生建模 QC、半径曲线和节点编号映射。拟合前将原始 SWC 及输入配置逐字节复制到派生目录的 `backups/`，校验 SHA-256；备份是独立文件，原始数据也继续保留。不能在原始数据目录内创建输出，已有派生目录也不能覆盖。原始配置保存在 `configs/swc_roi_generate_mouse_raw.yaml`，可继续用原 mouse 脚本传入该配置检查优化前数据。

适配层另增加两项执行保护：当原生 apex 搜索已经使用全部可用下游样本，仍返回“未找到分离位置”的端点结果时，终止不可推进的重复拟合并记录失败；原生 `distance()` 调用期间复用同一份离散点数组，调用结束立即恢复 getter，避免每次投影都重复转换整个数组。两者均不修改上游惩罚样条、AIC 或分叉数学，失败时也不生成替代几何。

## 源码核查与实际复用

| 上游对象/API | 本工具的使用方式 |
| --- | --- |
| `ArterialTree(patient_name, database_name, filename=None, automatic_resampling=True)` | 原生 SWC 读取；为记录原始映射，先暂缓默认重采样，再按用户选项调用同一个官方方法 |
| `get_full_graph()`、`get_topo_graph()`、`check_full_graph()` | 获取完整节点图、压缩拓扑图并验证原生要求 |
| `model_network(radius_model=True, criterion="AIC", akaike=False, max_distance=6, max_distance_radius=np.inf)` | network 模式每个文件调用一次；branches 模式不调用 |
| `Spline.approximation()` → `Model` | 上游先建空间模型，再以其参数化建立半径模型，分别优化两个平滑参数；branches 模式逐条原始拓扑分支直接调用 |
| `Model.quality("AIC")` | 上游 AIC 目标函数；由上游黄金分割搜索选择平滑参数，没有自行重写 |
| `Nfurcation` | 官方分叉形状、截面、apex、分离区域及轨迹建模 |
| `Spline.length()`、`length_to_time()`、`point(..., radius=True)` | 从普通血管样条按弧长读取连续模型 |
| `Nfurcation.get_tspl()`、`get_spl()`、`get_X()` | 分叉轨迹、形状与共用拓扑节点；相关半径表达限制见下文 |
| `compute_cross_sections(N, d, parallel=False)`、`mesh_surface()` | 原生表面生成，直接保存 VTK |

`criterion="AIC"` 用于优化平滑参数。`akaike` 是另一个用于选择控制点数量的选项，保持上游默认 `False`，不应混为一谈。`max_distance_radius` 虽出现在当前函数签名中，但函数体未使用，因此没有为其虚构 CLI 控制效果。

`full_graph` 节点保存 `(x,y,z,r)`；`topo_graph` 压缩度为 2 的普通节点，内部坐标及源 ID 存入边的 `coords/full_id`；`model_graph` 保存普通边样条，以及 `end/sep/bif` 节点和 `bifurcation` 对象。上游 SWC loader 丢弃 TYPE，因此包装层另读类型元数据，使用原生拓扑段关联恢复类型，不推断新的动脉分类。

普通血管段不采用跨分叉的均值滤波。`Nfurcation` 本身会联合使用母段与某个子段的样本拟合分叉形状，这是官方算法的组成部分，不能称为“所有段完全隔离拟合”。

### 为什么没有直接调用官方 SWC exporter

已检查 `model_to_full()` 和 `write_swc()`：前者使用样条默认离散点，不支持所需的原始段点数/弧长采样控制，且直接写入分叉轨迹的零半径占位值与仅含三维坐标的分叉中心；后者把根和其他节点的 TYPE 分别统一写为 1、3。分叉输出方向和端点去重也不能直接满足本任务。因此只新增导出适配层，不新增拟合算法。

### 分叉转换为 SWC 的限制

原生 `Nfurcation.get_tspl()` 是网格轨迹，第四维为**零占位**；不能将其视为血管半径。工具使用此原生轨迹的空间坐标，把位置投影到对应的原生形状样条，再通过其 `radius()` 读取半径。入口对应 shape 0，出口对应各自的 shape。

SWC 的一个共享节点只能存一个半径，无法完整表示非圆形的三维分叉。共享节点位置采用原生 `get_X()`，代表半径来自入口 shape 0 在该位置的投影值，**不平均母血管与子血管半径**。子段内部继续使用各自模型的半径。该约定会逐文件写入 QC 和 SWC 文件头说明，不能把转换后的 SWC 当作原生分叉表面的无损编码；精确检查分叉表面应使用官方 VTK。

类型继承按同一原始段的归一化弧长最近位置映射，内部采样点优先继承原段内部节点类型，分叉和端点使用原节点 TYPE。段内原本存在多种内部 TYPE 时记录 warning。

## QC、失败与可重复性

输入检查覆盖重复 ID、缺失父节点、多个根、环、自环、非有限数值和非正半径；不进行自动修补。SWC 父子方向按给定血流方向解释，文件本身不足以证明真实生理流向。

采用带原始末端标识的有根拓扑签名匹配原始段和模型段，忽略普通度为 2 节点；检查的不仅是分叉数量。共享分叉只生成一个节点，输出 ID 连续、父节点先于子节点。写出后重新读取并逐项核对连接、坐标、半径和 TYPE。

每条原始拓扑段分别统计 `abs(log(r[i+1]/r[i]))` 与对称相对半径变化，排除触及共享分叉节点的点对。输出样本数、median、P90、P95、P99、max，并报告半径范围、总中心线长度及节点/分支/分叉/末端数量。另提供每单位长度的对数半径变化统计，提醒 `spacing` 模式增加密度本身也可能降低相邻点跳变。

半径 CSV 使用 `series=raw/smooth` 的长表：两种曲线各用自己的弧长和样本索引，不虚构重采样后的逐点对应。`raw_radius` 与 `smooth_radius` 在对应的行填值；邻接变化是否计入总体 QC 由 `pair_evaluated` 指示。诊断图只显示原始半径与 VascularMD 模型采样值，标注分叉附近区域，没有第三条自行滤波的曲线。

QC 保存原始/预处理节点数、源 SHA-256、上游提交及源码摘要、依赖版本、全部参数、warning、错误阶段和可确定的 branch/node。若原生全网建模中途抛错，保留已建模边清单，报告人工复核；官方接口不保证失败后可以安全恢复全网，故不会伪造“其余段已完成”的整树 SWC。批处理继续处理下一个文件。

状态码 `0` 表示全部文件成功，`1` 表示至少有失败或部分输出，参数错误返回 `2`。`success` 只表示建模、导出和上述检查完成，不代表恢复了真实解剖半径，也不代表通过了 CFD 体网格质量验收。

没有额外调用表面的 Taubin/Laplacian 滤波。上游 `Nfurcation` 自身包含论文规定的网格松弛与回投影，这一原生步骤予以保留；其 PyVista 参数兼容修改没有改变迭代数或松弛系数。

## 验证与依赖

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -m pip install -r requirements-vascularmd.txt
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -m pytest -q \
  tests/test_vascularmd_smoothing.py tests/test_swc_export.py
```

附加依赖为 `geomdl==5.4.0` 和 `threadpoolctl==3.6.0`；其余使用项目已有的 NumPy、SciPy、NetworkX、PyVista、VTK、Matplotlib。本次切换 human 默认输入配置，保留原 ROI 可视化实现。

第三方来源、安装复现与最小兼容补丁见 [third_party/vascularmd_integration.md](../third_party/vascularmd_integration.md)。实测结果及人工复核事项见 [实现与验证报告](../references/文本记录%20-%20VascularMD%20原生建模与%20SWC%20导出验证.md)。
