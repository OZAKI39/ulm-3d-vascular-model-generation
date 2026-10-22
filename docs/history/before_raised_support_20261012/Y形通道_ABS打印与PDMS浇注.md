# Y 形通道：无底四壁盒、外置燕尾支架与横向锁键

2026-10-12 高度上限 15 mm 修订版，打印工艺为 **Bambu P1S / 0.4 mm / ABS**。

Y 形共有 **8 个生产组件**：血管芯、无底框架、三个燕尾支架、三枚锁键。三个端口均竖直向上。原始 `data/y_channel/source/labelled_lumen_m.vtp` 保持不变，SHA256 为 `e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b`；仍采用原分叉轴向的 14/12/12 mm 等径重建版本。

| 参数 | CAD 设计值 |
| --- | --- |
| 血管芯直径 | 1.60 mm |
| 最薄 PDMS 底层 | 1.50 mm，参考外部底板上表面 |
| 腔体平面净尺寸 | 37.59 × 22.60 mm |
| PDMS 设计浇注高度 | 9.60 mm；成品宽、高下限仍各为 9.00 mm |
| 四壁厚度 / 内部 ABS 底板 | 2.00 mm / 0 mm |
| 外伸安装脚厚度 | 2.00 mm，仅位于四壁外侧 |
| 框架外包围盒，含外置安装座 | 50.79 × 42.60 × 10.40 mm |
| CAD 净浇注容积 | 约 8.02 mL，不含转移等损耗 |

## 本轮 3dm 文件与培养皿高度

先用 Rhino 打开 **`Y_H15_casting_assembly.3dm`** 检查浇注装配；**`Y_H15_print_layout.3dm`** 是同一套组件的分件打印排版。两文件位于下方当前输出目录。每个组件是独立命名、独立图层、带颜色的封闭网格，单位为 **mm**，保存为 Rhino 7 格式。它们不是 NURBS 曲面；血管芯的解析 CAD 另存为 STEP。导出使用 [McNeel 官方 rhino3dm](https://mcneel.github.io/rhino3dm/python/api/File3dm.html)，并已通过库读回检查；没有宣称运行过 Rhino 桌面软件。

| 组件 | 装配顶面相对 Z0 | 独立零件高度 |
| --- | --- | --- |
| 四壁框架，含外置燕尾座 | 10.40 mm | 10.40 mm |
| 燕尾支架 | 14.80 mm | 12.80 mm |
| 横向锁键 | 7.10 mm | 3.20 mm |
| 血管芯，包含允许略高的竖直端口 | 15.40 mm | 13.90 mm |

非血管芯结构按 **15.00 mm 上限减去 0.20 mm 设计余量** 控制，未利用 ABS 收缩来降低名义高度；本轮实测 CAD 最高为 14.80 mm。血管芯端口高出 15 mm 限制 0.40 mm。此次降低的是液面上方空余高度，PDMS 设计浇注高度仍为 9.60 mm，支架下臂仍厚 2 mm、上叉仍厚 1.2 mm，侧壁仍厚 2 mm。燕尾有效啮合长度变为 8.40 mm，装配与锁定运动重新校验。

高度基准是培养皿内底面或独立浇注底板上表面。若在培养皿内另垫玻璃等底板，其厚度要额外加入装配高度，不能视为已包含在 14.80 mm 内。培养皿实际内径尚未提供，本轮验证高度，不宣称已完成培养皿实物放入或盖合测试；平面外形请按下方尺寸及 3dm 全部外伸件核对。

3dm 文件用于模型检查和后续几何处理，不含 Bambu 工艺、raft 或临时支撑刀路。实际打印仍打开重新生成的 P1S / 0.4 mm / ABS 已切片 3MF；无需把浇注装配态作为一个整体打印。

## 当前文件与入口

本地工程为 `/home/lzy/projects/vascular_printing`，服务器工程为 `/workspace/vascular_printing`。

```text
outputs/y_channel/D1p6_PDMSmin9_ABS_bottomless_dovetail_H15_3dm_P1S04/
```

优先打开 `Y_P1S_04_ABS_all_parts_sliced.3mf`；`Y_P1S_04_ABS_all_parts_project.3mf` 为可编辑工程。完整打印包为 `deliverables/Y_channel_D1p6_PDMSmin9_ABS_bottomless_dovetail_H15_3dm_P1S04.zip`。

| 文件 | 内容 |
| --- | --- |
| `assembly/Y_lumen_core.stl` | 血管芯及定位肩台 |
| `assembly/Y_bottomless_box_wall2mm.stl` | 无底四壁、外侧安装脚和燕尾座的一体框架 |
| `assembly/Y_dovetail_retainer_*.stl` | INLET、OUTLET_A、OUTLET_B 的肩台捕获支架 |
| `assembly/Y_dovetail_key_*.stl` | 对应横向锁键，分别打印一枚 |
| `assembly/Y_lumen_core.step` | 血管芯 CAD |
| `assembly/*reference.stl` | 非生产参考表面，不作为额外组件打印 |
| `print_parts/` | 已按打印姿态与盘位转换的独立 STL |
| `labelled_*_mm.vtp` | 带端口标签、单位明确为毫米的参考通道 |

```bash
# 本地重新生成几何与检查图
cd /home/lzy/projects/vascular_printing
./run_python.sh generate_y_channel_casting.py --stage geometry

# 服务器完整生成、真实 Bambu 切片与打包
cd /workspace/vascular_printing
./run_server_python.sh generate_y_channel_casting.py --stage all
```

配置为 `config/y_channel_p1s_abs.yaml`。STL 单位按 mm 导入；直接打开整盘 3MF 可保留对象排版及打印工艺。

## 结构为什么这样修改

盒体内部的 ABS 底板已完全移除，四壁内缘从 z=0 开始，向下、向上均敞开。外部的平整玻璃或其他独立底板上表面定义为 z=0。外置安装脚厚 2 mm，仅在燕尾座附近向外伸出 8 mm，不伸入浇注腔体。外部底板不是本次打印组件；浇注前需密封四壁下沿与该底板的接缝，不应直接将无底框架放空后浇注。

固定结构分为三部分：

| 结构 | 功能 |
| --- | --- |
| 与四壁一体打印的外置燕尾座 | 竖直导槽限制侧向移动，底部止挡确定支架高度；2 mm 腔体侧壁保持完整 |
| 独立燕尾支架 | 滑块从上方插入导槽；上部双层叉口从侧面接收血管芯定位肩台，同时限制下落和上浮 |
| 独立横向楔形锁键 | 穿过固定座和滑块的共同开口，限制支架向上退出；外侧有抓取头，可拆卸 |

原来的跨墙夹脚已经取消。燕尾座位于侧壁外侧；跨过盒口的支架臂位于目标液面上方。盒内没有用于固定支架的底板、安装脚或支柱。所有安装件与设计 PDMS 填充体积均无相交。

燕尾有效配合长度为 8.40 mm，后端宽 8 mm、喉口宽 5 mm、深 4.4 mm；两侧斜面的法向装配间隙为 0.20 mm。肩台的上方与径向余量各为 0.20 mm。横向锁键为带抓取头的渐薄键，设计开口为 2.4 × 2.4 mm，键体宽 2 mm。它通过几何阻挡限制支架上滑；实际松紧、承载和反复拆装性能仍须试打确认，不能把数字碰撞检查当作实物力学试验。

## 保留的尺寸和材料余量

通道 CAD 直径仍为 1.60 mm，主体最低点距外部底板上表面 1.50 mm。与上一无底版相比，液面以下血管主体、弯头、直径和端口轴向位置保持不变；只缩短液面上方的竖直端口，并降低定位肩台与支架。血管芯与支架必须使用本轮同一套文件，不能混用原 17 mm 高芯或 16.4 mm 高支架。

PDMS 成品净宽和高度仍要求均不少于 9.00 mm。沿用上一轮的 2% ABS 收缩设计预算、0.30 mm 总尺寸误差余量，按 0.20 mm 向上取整，得到 **9.60 mm 的 CAD 最小尺寸**：

```text
CAD 最小值 = ceil(((9.00 + 0.30) / 0.98) / 0.20) × 0.20 = 9.60 mm
预算下限   = 9.60 × 0.98 - 0.30 = 9.108 mm
```

[Prusa 官方 ABS 说明](https://help.prusa3d.com/article/abs_2058)给出的通用冷却收缩参考约为 1%–2%；[Bambu ABS 技术数据表 V3](https://store.bblcdn.com/s7/default/23b4cf2b83d5470bb96d19970b5f3ae8/Bambu_ABS_Technical_Data_Sheet_V3.pdf)未给出本项目可直接套用的固定收缩率。本预算尚未用本机实物校准，也没有虚构 PDMS 固化收缩率。最终宽高仍在固化、冷却后实测验收。

Z 坐标：外部底板上表面 0；通道主体中心线 2.30；目标液面 9.60；盒口 10.40；肩台下沿 12.40；肩台上沿 13.40；上叉口下表面 13.60；上叉口顶面 14.80；竖直端口顶面 15.40，单位均为 mm。最薄 PDMS 底层为 1.50 mm，独立玻璃等底板的厚度不计入其中。

## 打印和装配顺序

1. 打开包内 `_all_parts_sliced.3mf` 查看原生切片，机型为 **Bambu P1S、0.4 mm 喷嘴、ABS**。同盘所有组件已分开排列，层高 0.20 mm。只保留单个变体的一套组件进行装配。
2. 血管芯和肩台捕获支架启用原生普通贴合支撑；只有血管芯使用两层 raft。框架和锁键不启用支撑。先去除临时支撑、brim 和 raft，清理燕尾配合面、横向键孔和肩台槽。原生临时支撑位于切片刀路中，不是要装入浇注腔体的构件。
3. 在盒外，把各支架从对应竖直芯杆的外侧向内推入，使芯杆进入开口，肩台位于上下两层叉口之间。不要把细血管杆当作撬动手柄。
4. 将血管芯和支架作为一组，从上方沿各燕尾导槽垂直下放，直至滑块抵住底部止挡。不要逐个强扭已经入槽的支架。几何检查覆盖了侧向接收肩台和整组下放的路径。
5. 从各固定座外侧的键孔方向推入对应横向锁键，抓取头留在外侧。锁键就位后再检查芯高度及各支架的活动量。拆卸时先抽出全部锁键，再沿导槽上提整组；不要强拉血管芯。
6. 将无底框架放在独立平整底板上，确认下沿连续贴合并从外部密封接缝。不要让密封材料垫高芯高度基准或流入有效通道区域。浇注液面参考盒口下方约 0.80 mm，并从实际底板上表面复核液厚；保留 ABS 余量，不按旧名义 9.00 mm 直接控液。
7. 固化后按已有 ABS 牺牲芯去除流程处理。未进行实物密封、溶解残留或成像验证。显微观察仍需结合实际底板厚度和物镜工作距离。

3MF 保持 `filament_shrink=100%`、XY 补偿为 0，避免将现有腔体尺寸余量重复放大。共享工艺名可能含 `@BBL X1C`，实际机型配置和 G-code 均检查为 P1S。

## 检查图、日志与开发约定

`QC/assembly.png` 是装配图，`bottom.png` 从下方显示完全开放的内腔，`exploded.png` 表示支架向上退出、锁键侧向退出的顺序；`dovetail_detail.png` 为外置接口近景，`dovetail_section.png` 在 z=8 mm 截开接口，显示真实燕尾横截面。`native_support_plate.png` 展示实际模型与原生临时支撑刀路。

`rhino_qc.json` 记录两个真实 3dm 的毫米单位、对象/图层、顶点和面片读回、封闭性、高度及 SHA256；`geometry_qc.json` 的 `manufacturing_height_qc` 对每个组件分别记录装配顶面与自身高度。

`geometry_qc.json` 记录全内腔体积检查、底部射线、四壁连续性、管径、PDMS 宽高预算、所有装配碰撞，以及锁定和释放测试。`print_qc.json` 记录真实机型、喷嘴、材料、3MF 几何保持、各组件实际挤出刀路和补偿设置。`bambu/` 保存真实切片命令和日志；`delivery_receipt.json` 记录交付 ZIP 的 CRC 与 SHA256 核验结果。

共享结构实现为 `vascular_processing/bottomless_channel_fixture.py`，通过原有 Y/直管入口接入；`fixture_mode: bottomless_external_dovetail` 开启本模式，`base_mm: 0` 强制内部无底，`dovetail.external_foot_mm: 2` 是外伸脚厚度。独立底板不建模为打印件。`height_limits.fixture_max_mm: 15`、`cad_reserve_mm: 0.2` 为装配高度门槛，`freeboard_mm: 0.8` 为液面至盒口的空余高度，`rhino_export: true` 开启 3dm 导出。修改参数后需重新运行 `geometry` 或 `all`；`slice` 会拒绝配置与几何不一致的旧文件。

上一无底高支架版的配置与说明保存在各自 `history/before_height15_3dm_20261012/`；更早的有底版在 `history/before_bottomless_dovetail_20261012/`。旧模型和 ZIP 原位保留。默认入口生成本次低高度无底版本，不应混用旧支架、血管芯或切片。本轮测试与审计集中在 `logs/channel_height15_3dm_20261012/`。
