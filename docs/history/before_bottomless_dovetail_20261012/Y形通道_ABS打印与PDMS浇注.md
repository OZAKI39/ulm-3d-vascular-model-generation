# Y 形通道：ABS 打印与 PDMS 浇注

本流程将 `labelled_lumen_m.vtp` 中的 Y 形腔道做成 ABS 牺牲芯，配套一个有底开口盒和三个端口定位支架。打印设备为 **Bambu Lab P1S、0.4 mm 喷嘴、ABS**。打印完成后先去掉血管芯的打印支撑，再把干净的血管芯安装到盒中浇注 PDMS。三个端口从 PDMS 顶面竖直引出，中央通道靠近底部，便于后续从底部观察。

**2026-10-12 ABS 余量修订：血管芯 CAD 直径 1.60 mm；成品 PDMS 高度/净宽均要求至少 9.00 mm，按材料与误差预算计算 CAD 下限为 9.60 mm。Y 形腔体宽约 22.60 mm，浇注设计厚度为 9.60 mm。旧名义 9 mm 版本保留为历史。**

这是一条独立的 Y 通道制造流程。现有 BG001、candidate 0/15、s1-10/s1-11 的模型和可视化没有改动。

## 1. 输入是什么

本地工程根目录是 `/home/lzy/projects/vascular_printing`；服务器对应目录是 `/workspace/vascular_printing`。

输入保存在 `data/y_channel/source/`：

| 文件 | 用途 |
| --- | --- |
| `labelled_lumen_m.vtp` | 带壁面、入口、两个出口标签的封闭 Y 形表面 |
| `PLANES.json` | 三个原始端面的位置、法向和半径 |
| `SUMMARY.json`、`IDENTITY.json` | 原始生成参数、单位与来源记录 |
| `centerlines.vtp` | 原始中心线，供检查参考 |
| `provenance.json` | 本次复制来源及文件 SHA256 |

选用原 Y 通道项目中最后一次成功的镜像对称版本：

`/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/geometry/20261002T135639758501Z/`

它也是 Y-clot 验证流程使用的版本。源 VTP 的 SHA256 为：

```text
e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b
```

**源坐标为米，输出模型为毫米。** 代码根据 VTP 中的 `coordinate_scale_to_m` 和源文件校验值确认单位，进行 1000 倍单位换算。新制造版本保留分叉位置、14/12/12 mm 分支长度及原轴向布局，按用户要求把原来的 2.60/2.00 mm 管径重建为统一 1.60 mm；这不是全模型缩放，也不再声称保留原壁面。原始 VTP 不修改。这是一个理想化 Y 形通道，不是 BraVa 解剖血管。

## 2. 本次尺寸与构造

### ABS 尺寸余量的依据与计算

`minimum_pdms_height_mm=9` 和 `minimum_pdms_width_mm=9` 现在表示**最终 PDMS 的验收下限**，不是直接用于建模的尺寸。[Prusa 官方 ABS 说明](https://help.prusa3d.com/article/abs_2058)给出的通用冷却收缩参考约为 1%–2%，并建议通过试件测量校准。[Bambu ABS 技术数据表 V3](https://store.bblcdn.com/s7/default/23b4cf2b83d5470bb96d19970b5f3ae8/Bambu_ABS_Technical_Data_Sheet_V3.pdf)说明材料存在翘曲倾向，但未给出可直接用于本项目的固定成形收缩率。

本次采用 **2% 的收缩设计预算**，在宽度和高度上分别再预留 **0.30 mm 的总尺寸误差预算**，覆盖未校准时的制造、测量和液面控制误差。0.30 mm 是本项目选择的工程余量，不是厂家保证的误差，也不是每侧各加 0.30 mm。为与 0.20 mm 层高兼容，向上取整：

```text
CAD 最小尺寸 = 向上取整到 0.20 mm [(9.00 + 0.30) / (1 - 0.02)]
             = 9.60 mm
预算下限     = 9.60 × 0.98 - 0.30 = 9.108 mm
```

宽度与高度各自保存收缩参数，暂用同一保守预算，不声称 ABS 各方向的实测收缩相同。宽于 9.60 mm 的天然腔体不会缩小。只调整腔体边界、浇注高度及相关竖直杆和支架，不整体缩放模型；1.60 mm 芯径、1.50 mm 底层、2.00 mm ABS 壁/底板均保留其 CAD 设计值。

配置中的 `abs_dimensional_allowance` 保存材料、两方向收缩预算、误差余量和取整步长；共享实现为 `vascular_processing/abs_dimensional_allowance.py`。实际 STL 内壁射线与内底/盒口测量见 `geometry_qc.json` 的 `pdms_width_qc`、`pdms_height_qc`，其中 `allowance` 明确区分成品下限、CAD 下限、实际 CAD 尺寸及预算下限。名义 9 mm 的旧几何不能通过新检查。

3MF 中保持 `filament_shrink=100%`、XY 轮廓/孔补偿为 0、圆孔补偿关闭，避免对已有余量再次缩放，也避免改变血管芯直径。`print_qc.json` 会核对实际切片工程中的这些参数。不要再手动整体放大本版本工程。

**这些是未校准的设计预算，不是实物尺寸合格证明。** 尚未指定 PDMS 配方、固化温度与时间，代码没有虚构固化收缩率。浇注时可参考冷却后盒口下方约 2.4 mm 的位置，同时从实际内底面量取液厚，保留余量并避免浸入支架；不能按旧版 9.00 mm 直接控液面。成品宽度和高度仍以固化、冷却后实测均不少于 9.00 mm 为验收要求。首次试打记录实际内宽和内高，用于后续单独校准两方向预算。

Y 形原有腔体净宽约 **22.60 mm**，已大于带余量的 9.60 mm 下限，平面尺寸保持不变。盒体仅增高 0.60 mm；竖直杆上端、肩台和三个定位支架相应升高。

用户确认的 **1.5 mm 是通道外表面到 PDMS 底面的最小距离**。它与 2 mm 的 ABS 盒底厚度不同。

| 参数 | 本次值 |
| --- | --- |
| 入口、两出口、主体各分支管径 | 均为 1.60 mm |
| 原始轴向分支长度 | 入口 14 mm，两出口各 12 mm |
| PDMS 最薄底层 | 1.50 mm，主体保持靠近内底面 |
| 盒壁与底板厚度 | 均为 2.00 mm |
| 浇注腔平面内尺寸 | 37.59 × 22.60 mm |
| 盒体外尺寸 | 41.59 × 26.60 × 14.00 mm |
| 通道外包围盒到四侧内壁余量 | 各 3.00 mm |
| PDMS CAD 浇注厚度 | 9.60 mm，从内底面计算；不含 ABS 底板和液面上方余量 |
| 液面到盒口余量 | 2.40 mm |
| 目标液面下净容积 | 约 8.02 mL，不含混合、脱泡和转移损耗 |
| 入口 / 出口弯曲中心线半径 | 统一为 2.40 mm |

三条分支使用半径 0.80 mm 的圆截面 CAD 扫掠。端部先沿原端面法向延伸 0.6 mm，再通过中心线半径 2.4 mm 的四分之一圆弧转为竖直；直线和圆弧切向连续。分叉处按顺序施加 0.10 mm 和 0.05 mm 局部 CAD 圆角，消除交汇锐边；不使用全局平滑。分支汇合处自然形成非圆形连通区域，不能将其解释为独立圆管截面。42 处主体、弯段、竖直段截面及 3 个端面均从实际网格测量，直径约 1.5998–1.5999 mm，离散误差小于 0.003 mm。

盒外底面 z=0，盒内底面 z=2，主体中心线 z=4.3，CAD 目标液面 z=11.6，盒口 z=14.0，竖直杆顶面 z=19.0。因此 CAD 浇注厚度为 11.6−2=9.6 mm，不含盒底和液面上方余量。高于液面的定位肩台外径为 4 mm，它不参与形成 1.6 mm 通道。

三个竖直芯杆在液面上方带有定位肩台。独立 C 形支架承托肩台，并以两条小脚跨住盒壁，确定血管芯高度。支架最低处比目标液面高 0.4 mm；定位肩台也在液面上方。因此，这些定位结构不会被浇入有效 PDMS 区域。定位配合单侧间隙为 0.25 mm，实际配合仍需要试装。

这里“侧向余量 3 mm”是通道整体外包围盒与盒壁之间的尺寸约束。矩形盒在 Y 字凹侧会有额外空白，并不意味着每个壁面点到通道的距离都恰好为 3 mm。

## 3. 生成入口与修改参数

入口为 `generate_y_channel_casting.py`，调度与打印实现为 `vascular_processing/y_channel_casting.py`，等径几何实现为 `vascular_processing/y_uniform_channel.py`，配置为 `config/y_channel_p1s_abs.yaml`。

本地生成模型和检查图：

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh generate_y_channel_casting.py --stage geometry
```

服务器生成并使用已配置的 Linux Bambu Studio 完成真实切片：

```bash
cd /workspace/vascular_printing
./run_server_python.sh generate_y_channel_casting.py --stage all
```

只重新切片已有模型：

```bash
./run_server_python.sh generate_y_channel_casting.py --stage slice
```

`--stage slice` 会核对配置和 STL 校验值。如果改过底层厚度、盒体尺寸等参数，需要先重新执行 `geometry` 或 `all`，不能将旧 STL 当成新模型切片。

最常修改的参数：

| 配置项 | 含义 |
| --- | --- |
| `geometry_mode`、`diameter_mm` | 当前为 `uniform_diameter`、1.6；保留旧模式仅供历史复现 |
| `minimum_pdms_width_mm` | 成品 Y 方向净宽下限 9；先计算带 ABS 余量的 CAD 下限，原有约 22.60 不缩小 |
| `minimum_pdms_height_mm` | 成品高度下限 9；经材料预算生成 CAD 下限 9.6，净空不足时继续增高 |
| `pdms_floor_mm` | 通道到 PDMS 底面的最小距离，本次 1.5 |
| `base_mm`、`wall_mm` | ABS 盒底和盒壁厚度，本次均为 2 |
| `side_clearance_mm` | 通道整体包围盒到盒内壁的侧向余量 |
| `bend_radius_to_lumen_radius` | 弯管中心线半径 / 当前通道半径，本次为 3 |
| `pdms_above_elbow_mm` | 液面高于最高弯管末端中心的距离 |
| `freeboard_mm` | 液面到盒口的高度，必须给支架留出空间 |
| `retainer_fit_clearance_mm` | 支架与杆、盒壁的配合间隙 |
| `output` | 输出目录；保留旧版本时应使用新的目录名 |

其它血管不能只改文件名直接套用：该入口明确检查本 Y 模型的来源、四类标签、平面分叉与三个端面，以免误判单位和入口方向。

## 4. 文件在哪里、应该打开哪个

输出根目录：

```text
outputs/y_channel/D1p6_PDMS_min9_ABS_allowance_H9p6_Wmin9p6_P1S_04_floor1p5/
```

| 路径 | 内容与用途 |
| --- | --- |
| `Y_P1S_04_ABS_all_parts_project.3mf` | 可编辑 Bambu 工程，同一盘的血管芯、盒体、三个支架；保留 P1S/ABS 和分对象支撑设置 |
| `Y_P1S_04_ABS_all_parts_sliced.3mf` | 经真实 Bambu Studio 切片的打印盘，含实际打印支撑和刀路；导入后先检查预览 |
| `assembly/Y_lumen_core.stl` | 带三个竖直端口和定位肩台的血管芯，处于浇注装配坐标 |
| `assembly/Y_casting_box_base_2mm.stl` | 2 mm 盒底与四壁的一体式盒体 |
| `assembly/Y_retainer_INLET.stl` | 入口定位支架 |
| `assembly/Y_retainer_OUTLET_A.stl`、`Y_retainer_OUTLET_B.stl` | 两出口定位支架 |
| `assembly/Y_lumen_core.step`、`labelled_y_D1p6_lumen_mm.vtp` | 精确 CAD 与带端口标签、明确毫米单位的参考腔道 |
| `assembly/*_reference.stl` | 原始主体与不带肩台的通道参考模型，不是额外打印组件 |
| `print_parts/` | 已安排打印姿态与盘上位置的五个 STL |
| `print_layout.json` | 从装配坐标到打印坐标的变换及逆变换 |
| `QC/assembly.png`、`top.png`、`side.png` | 装配、俯视、侧视检查图；淡青平面表示目标液面 |
| `QC/core.png`、`exploded.png` | 血管芯单独显示和装配分解图 |
| `QC/native_support_plate.png` | 模型与实际打印支撑刀路叠加图 |
| `QC/native_print_support_toolpaths_mm.vtp` | 实际 G-code 支撑线段，可在 ParaView 中检查；不是可打印支撑 STL |
| `geometry_qc.json` | 实测截面、PDMS 高度、封闭性、碰撞和原始文件校验 |
| `print_qc.json`、`manufacturing_summary.json` | 实际切片检查与交付状态 |
| `profiles/profiles.json` | 使用的真实打印配置、来源校验值和本次覆盖项 |
| `bambu/command.json`、`stdout.log`、`stderr.log`、`slice_result.json` | 实际切片命令、日志与结果 |
| `bambu/plate_1.gcode` | Bambu 原生输出的 G-code，代码没有手工拼接刀路 |

**首次使用建议打开 `Y_P1S_04_ABS_all_parts_sliced.3mf`。** 五个对象已经分开排好；无需再次把装配坐标下的 STL 导入并重新排盘。参考 STL 不应与生产组件一起打印。

两类“支撑”不要混淆：黄色的三个定位支架是需要单独打印、用于浇注装配的实体；绿色刀路线是 Bambu 给血管芯生成的临时打印支撑，包含在已切片 3MF/G-code 中。临时支撑没有被伪装成独立实体 STL，也不能带入 PDMS 浇注。

## 5. 打印与装配顺序

1. 在 Bambu Studio 中检查实际机型为 **P1S**、喷嘴 **0.4 mm**、材料 **ABS**，使用纹理 PEI 板对应配置。层高 0.20 mm，轮廓线宽 0.42 mm。原生共享工艺名称中含 `@BBL X1C`，但机器配置及起止 G-code 已核验为 P1S。
2. 血管芯使用贴合式普通支撑和两层 raft；支撑顶隙 0.20 mm、侧隙 0.35 mm。盒体正向打印，三个定位支架翻面打印，均不启用自动支撑。外周有 5 mm brim。细小轮廓降速，血管尺寸没有为打印便利而放大。
3. 打印后先拆去 brim、raft 和血管芯的全部临时支撑。握持近端粗直段，分段剥离支撑，避免借支架或长管当杠杆扭动分叉。确认没有残余支撑连到将被 PDMS 包住的表面。
4. 在盒外将三个 C 形定位支架从侧面套到对应竖直杆上，使肩台坐在支架上。随后整体下放到盒口，让定位脚跨住对应盒壁。支架和血管芯作为一组下放，避免在盒内强行挤入支架。检查血管芯最低处离盒内底面 1.5 mm，三个端口方向与盒底垂直。
5. 浇注前确认支架不会晃动，肩台贴合承托面。这些 C 形支架负责承托和定位，没有封闭的上压锁扣；若浇注或脱泡会抬起血管芯，需要在液面外临时压住肩台，不能用贯穿 PDMS 区域的附加支柱固定。
6. CAD 液面相对内底面为 9.60 mm，设计净容积约 8.02 mL；实物按上文余量与控液说明核对。不要加至盒口，否则支架会嵌入 PDMS，并使保留的端口位置改变。可用深度尺或外部水平标尺辅助控制。
7. 固化后按实验室已有的 ABS 牺牲芯去除流程去掉 ABS；在有效 PDMS 之外拆解或切除定位件，不能强拔整个 Y 形芯。盒体也属于待移除的 ABS 部分。确认腔道通畅、无残留和泄漏后，再接管并进行观察。

## 6. 检查结论与边界

程序实际检查了单位、源文件 SHA256、单连通封闭实体、三端口方向、42 处截面与三个端面直径、带材料余量的 PDMS 宽高、组件体积相交、定位件与目标 PDMS 区域的分离，以及排盘间距。3MF 另做 ZIP CRC、真实机型/材料、组件几何和实际挤出刀路检查。

Bambu 的 raft 在本次刀路中把血管芯抬高了 0.6 mm，但 3MF 中存储的 STL 坐标没有同步抬高。检查代码从血管芯第一层实际模型挤出恢复该偏移，再对比血管表面。这个打印偏移不会改变浇注装配中的 1.5 mm 底层。

通过上述检查意味着数字模型和真实切片可以进入实物试打。**尚未实际验证支撑拆除、ABS 收缩、盒体密封、PDMS 固化/脱模、溶解残留及成像质量。** 当前交付状态会在 `manufacturing_summary.json` 中如实记录。

ABS 牺牲模板制备 PDMS 腔道有公开实验依据，见 [Saggiomo 与 Velders，Advanced Science，2015](https://doi.org/10.1002/advs.201500125)。本文采用这一制造思路，但 1.5 mm 底层、3 mm 侧余量和竖直引出是本项目的设计参数，并不是所有 PDMS 实验通用的“规范值”。

**显微镜观察发生在 ABS 底座去除之后。** 1.5 mm 只是最薄 PDMS 覆盖层，血管中心线位于 PDMS 底面上方2.30 mm；还需要考虑观察目标在管内的深度、物镜工作距离和底面平整度。不能据此保证高倍率物镜能对焦。[Nikon 的物镜说明](https://www.microscope.healthcare.nikon.com/literature/Brochures/2ce-mphh-2.pdf)明确区分工作距离与盖玻片校正厚度。若后续加玻璃衬片、改用玻璃键合底面或抛光处理，应把这些层厚纳入新的几何版本，不能直接加到现有盒底上却仍声称 PDMS 厚度为 1.5 mm。

## 7. 回归测试

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh -m pytest -q tests/test_abs_dimensional_allowance.py tests/test_y_channel_casting.py tests/test_y_uniform_channel.py
```

测试使用真实源 VTP，独立核对端口截面积、顶面法向、盒底射线交点、新模型截面、PDMS 高度与坐标变换；旧模式另保留原壁面回归检查；同时验证改变底层厚度后不能误切旧几何。服务器环境可将 `run_python.sh` 换成 `run_server_python.sh`。

当前完整交付包：`deliverables/Y_channel_D1p6_PDMS_min9_ABS_allowance_H9p6_Wmin9p6_P1S_04_floor1p5.zip`。旧文件保留在原输出目录；旧配置位于 `config/history/before_D1p6_H9_20261012/`，旧说明位于 `docs/history/before_D1p6_H9_20261012/`。

本轮宽度检查前的配置和说明保存在各自的 `history/before_minimum_pdms_width_9_20261012/`。本次增高后应使用全套新组件，不再与旧 STL 混用。ABS 余量修订前配置、说明位于各自 `history/before_abs_dimensional_allowance_20261012/`。
