# BraVa ColorCoded TYPE 的来源追踪与左右侧命名冲突

## 摘要

本次追踪旨在确定 BG001 的 ColorCoded SWC 中哪些 TYPE 对应左、右大脑中动脉，以避免将未经确认的左右侧定义带入 TopBrain 语义迁移。结果证实，本地 ColorCoded 文件和普通非平滑文件分别与 BraVa 官方同版本下载文件逐字节一致，版本错配不是当前问题的原因。关联研究的作者代码明确规定 TYPE 3 为右 MCA、TYPE 4 为左 MCA；然而，官方 BG001 网页标注的左右 MCA 长度与本地两类子树的长度反向对应。因此，目前能够确认 TYPE 3、4 属于 MCA 相关标注，但不能将其左右侧命名视为已经无歧义地确认。生产配置保持未赋值，状态记录为 `BRAVA_TYPE_LATERALITY_CONFLICT`。

## 资料与方法

首先，从 Wright 等人的 BraVa 原论文与 NITRC 登记页追踪至官方站点。HTTPS 入口出现证书域名不匹配；HTTP 入口正常返回页面，首页的跳转进一步导向 `home.php`。通过 Browse 页面进入 BG001 的独立记录，读取 LEFT/RIGHT MCA 表格以及普通、颜色编码两种查看器的 JNLP 文件。JNLP 中的真实参数给出了 SWC 下载地址，因此文件来源并非根据名称猜测。所有取回页面和文件均保存 SHA256。

其次，检查 Padmos 等人发表的脑血流耦合研究所公开的作者代码。在固定提交 `cadcae820299487a9b9ffc07c744e82065c22e8a` 中，`ProcessingBravaset.py` 将 SWC 第二列直接传入 `MajorVesselID`；`GeneralFunctions.py` 的名称字典明确将 3、4 分别命名为右、左 MCA。这是后续研究作者的明确实现约定，不等同于 BraVa 原作者发布的编码字典。

最后，将本地普通 SWC 与 ColorCoded SWC 建立版本对应。两者均包含 2810 个节点，但根和编号不同。对应必须同时满足：坐标处于 ColorCoded 四位有效数字的舍入误差内、半径一致到其存储精度、节点一一对应、全部无向边一致。随后只将 TYPE 元数据对应到普通文件，以普通文件的原始精度计算同 TYPE 节点之间的边长总和。没有更改任何坐标、半径、TYPE、父节点或左右侧标签，也没有进行配准或解剖分类。

## 结果

本地 `BG001_ColorCoded.CNG.swc` 与官方 `files/file_swc_color/` 文件逐字节一致，SHA256 为 `733b255aeeb7795de45799f64e73b19a93ada45ca63518e9aa479c203658ef0e`。本地 `BG001.CNG.swc` 与官方 `files/file_swc/` 标准化文件逐字节一致，SHA256 为 `5b69c6c406e4724ffc15a2a4c6b5958cbb169b822a05098d6ebedeecbbb3afbe`。

另一个官方地址 `files/file_swc/original/` 返回不同版本，其 Y 坐标与标准化版本符号相反；其他数值逐项一致。本次仅核对这一事实，没有对项目输入进行翻转。本地所谓“raw”更精确的名称应为“官方标准化、未经 VascularMD 平滑的 SWC”。

2810 个节点全部获得唯一对应，全部无向连接一致。两个版本最大坐标分量差为 0.05 mm，最大半径差为 0.0004 mm，均处于 ColorCoded 的十进制存储精度范围。这些差异仅用于版本对应，不写回原始几何。

| 项目 | TYPE 3 | TYPE 4 |
|---|---:|---:|
| ColorCoded 节点数 | 626 | 796 |
| 按普通 SWC 原始精度计算的同 TYPE 边长总和/mm | 1433.76408695 | 1775.71520751 |
| 四舍五入至官方表的两位小数/mm | 1433.76 | 1775.72 |
| BG001 官方表中相同长度所在列 | LEFT MCA | RIGHT MCA |
| 关联研究作者代码的明确名称 | RIGHT MCA | LEFT MCA |

官方页面的 LEFT/RIGHT 表头位于所保存 HTML 的第 523–524 行，相应长度位于第 544–545 行；后续作者代码的明确字典位于 `GeneralFunctions.py` 第 20–26 行。前者是由官方具名长度进行的对应推断，后者是明确的数值命名约定，两者的证据性质不同，不能将前者写成官方直接给出的 TYPE 字典。

辅助核查 BG0002 也呈相同的左右对应趋势，但该普通文件包含一个零半径节点，因此仅作只读来源审计，没有进入生产校验、迁移或建模。原论文图 1 提供动脉颜色示意而非 TYPE 数值表；不同渲染工具的颜色不能据此直接互换。官方查看器 JAR 仅作文件内容检查，未执行、安装，也未检出 MCA/cerebral 命名文字。

## 讨论与结论

来源追踪已将问题从“缺少映射文件”缩小到可复核的左右命名冲突。TYPE 3、4 的 MCA 身份有相互支持的证据，但目前尚不能判断冲突来自官方统计表列名、后续代码约定，还是历史坐标约定。单凭 X 坐标、图像左右位置或配准分数选择一种解释，会重新引入本轮明确禁止的解剖推断。

因此，左右 MCA 的生产配置保持空值；不自动选用其中一种约定，不运行 BG001 同侧 donor 配准，不生成带有已确认左右名称的 ROI，也不接入 s1-3。已经完成的 NN 冻结、四组支持阈值校准和精确子集代码予以保留。后续需要获得 BraVa 原作者的明确 TYPE 说明、带方向标记的原始影像对应，或对当前具体文件的人工确认；本次未向任何人发送消息。

## 后续采用约定

在本次追踪之后，用户明确决定采用官方 BG001 长度表的对应推断继续开发。生产约定因此设为 TYPE 3 对应 LMCA、TYPE 4 对应 RMCA，并绑定本文列出的两个 BG001 源文件 SHA256。该决定不改变上述来源冲突的事实，也不把推断提升为原作者的明确编码字典。后续真实迁移、ROI 提取和建模结果另行记录于《TopBrain 到 BraVa BG001 的 NN 语义迁移与 MeVO 建模》。本文前述“未赋值、未启动”描述的是来源追踪结束时的历史状态。

## 来源

1. Wright 等，原始 BraVa 研究：[论文全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3971907/)。
2. [NITRC BraVa 登记](https://www.nitrc.org/projects/breva/)。
3. [BraVa 官方 BG001 记录](http://cng.gmu.edu/brava/query_subject.php?id=1)。
4. [官方 ColorCoded 查看器参数](http://cng.gmu.edu/brava/files/file_jnlp_color/BG001_color.jnlp)。
5. Padmos 等，脑血流耦合研究：[论文及代码来源说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC7739918/)。
6. 作者代码：[TYPE 名称字典](https://github.com/rmpadmos/Coupling1DBF/blob/cadcae820299487a9b9ffc07c744e82065c22e8a/python/GeneralFunctions.py#L20-L27)，[SWC TYPE 读取](https://github.com/rmpadmos/Coupling1DBF/blob/cadcae820299487a9b9ffc07c744e82065c22e8a/python/ProcessingBravaset.py#L23-L39)。

机器可读证据、下载页面、文件哈希及节点映射保存在 `outputs/topbrain_brava_transfer/nn_production/type_mapping_trace/`。
